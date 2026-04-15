"""
Phase 4 ETL pipeline: extract → validate → enrich → transform → validate → load.

Supports SQLite source (publication) or MySQL (when tables mirror Starbase schema).
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import yaml

from starbase_validation.extractors import extract_from_sqlite
from starbase_validation.extractors.mysql_extractor import extract_from_mysql
from starbase_validation.loaders import DatabaseVersionManager, VersionInfo, load_to_sqlite, read_version_from_db
from starbase_validation.reporters import generate_report
from starbase_validation.transformers import enrich_records, normalize_records
from starbase_validation.validators import validate_all

logger = logging.getLogger(__name__)


def _default_config_path() -> Path:
    return Path(__file__).resolve().parent / "config" / "quality_rules.yaml"


def load_quality_rules_yaml(path: Optional[str] = None) -> Dict[str, Any]:
    """Load YAML and flatten quality_gates for validate_all()."""
    p = Path(path) if path else _default_config_path()
    with open(p, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    gates = raw.get("quality_gates", {})
    merged = {
        "taxonomy": gates.get("taxonomy", {}),
        "genome": gates.get("genome", {}),
        "classification": gates.get("classification", {}),
        "sequences": gates.get("sequences", {}),
        "blocking_issues": raw.get("blocking_issues", []),
        "warning_issues": raw.get("warning_issues", []),
    }
    return merged


class ValidationPipeline:
    """
    End-to-end validation and publication pipeline.

    Typical flow (SQLite source):
        pipeline = ValidationPipeline(quality_rules_path="...")
        result = pipeline.run(sqlite_path="starbase.sqlite", output_path="out.sqlite")
    """

    def __init__(
        self,
        mysql_config: Optional[Dict[str, Any]] = None,
        quality_rules: Optional[Dict[str, Any]] = None,
        quality_rules_path: Optional[str] = None,
    ):
        self.mysql_config = mysql_config
        if quality_rules is not None:
            self.quality_rules = quality_rules
        else:
            self.quality_rules = load_quality_rules_yaml(quality_rules_path)

    def extract_from_sqlite(self, sqlite_path: str) -> Dict[str, Any]:
        logger.info("Extracting data from SQLite: %s", sqlite_path)
        return extract_from_sqlite(sqlite_path)

    def extract_from_mysql(self) -> Dict[str, Any]:
        if not self.mysql_config:
            raise ValueError("mysql_config is required for extract_from_mysql()")
        logger.info("Extracting data from MySQL")
        return extract_from_mysql(self.mysql_config)

    def validate(self, records: Dict[str, Any]):
        blocking = set(self.quality_rules.get("blocking_issues", []))
        return validate_all(records, self.quality_rules, blocking)

    def validate_twice(self, records: Dict[str, Any]) -> Tuple[Any, Any]:
        """Run validation before and after enrich/transform (plan Phase 4)."""
        first = self.validate(records)
        enriched = enrich_records(records, {})
        normalized = normalize_records(enriched)
        second = self.validate(normalized)
        return first, second

    def run(
        self,
        sqlite_path: Optional[str] = None,
        output_path: Optional[str] = None,
        dry_run: bool = True,
        skip_final_validation: bool = False,
        version_bump: str = "minor",
        previous_published_db: Optional[str] = None,
        notes: str = "",
    ) -> Dict[str, Any]:
        """
        Run full pipeline.

        Args:
            sqlite_path: Source Starbase SQLite database
            output_path: Destination for published DB (required if not dry_run)
            dry_run: If True, do not write output file
            skip_final_validation: If True, skip second validation pass after transform
            version_bump: major | minor | patch for next version
            previous_published_db: Optional path to last published DB to read version from
            notes: Stored in database_version.notes

        Returns:
            Dict with keys: success, first_validation, second_validation, report,
            version, output_path, error (if any)
        """
        if not sqlite_path:
            raise ValueError("sqlite_path is required for SQLite-based ETL")

        records = self.extract_from_sqlite(sqlite_path)
        first = self.validate(records)

        if first.has_blocking_issues():
            report = generate_report(first)
            logger.error("Blocking issues: %s", len(first.blocking_issues))
            return {
                "success": False,
                "first_validation": first,
                "second_validation": None,
                "report": report,
                "version": None,
                "output_path": None,
                "error": f"{len(first.blocking_issues)} blocking issues",
            }

        enriched = enrich_records(records, {})
        transformed = normalize_records(enriched)

        if skip_final_validation:
            second = first
        else:
            second = self.validate(transformed)

        if not second.is_acceptable():
            report = generate_report(second)
            return {
                "success": False,
                "first_validation": first,
                "second_validation": second,
                "report": report,
                "version": None,
                "output_path": None,
                "error": "Final validation failed",
            }

        report = generate_report(second)
        joined_n = len(transformed.get("joined_ships", []))

        last_ver = None
        if previous_published_db:
            last_ver = read_version_from_db(previous_published_db)
        mgr = DatabaseVersionManager(last_version=last_ver, bump=version_bump)
        version = mgr.next_version()

        if dry_run:
            logger.info("[DRY RUN] Would create database version %s", version)
            return {
                "success": True,
                "first_validation": first,
                "second_validation": second,
                "report": report,
                "version": version,
                "output_path": None,
                "error": None,
            }

        if not output_path:
            raise ValueError("output_path is required when dry_run=False")

        info = VersionInfo(
            version=version,
            created_at=datetime.now(timezone.utc).isoformat(),
            validation_passed=True,
            total_joined_ships=joined_n,
            blocking_count=len(second.blocking_issues),
            warning_count=len(second.warning_issues),
            source="ValidationPipeline",
            notes=notes or "",
        )

        load_to_sqlite(
            sqlite_path,
            output_path,
            info,
            extra_metadata={"summary": second.summary},
        )
        logger.info("Published database version %s -> %s", version, output_path)

        return {
            "success": True,
            "first_validation": first,
            "second_validation": second,
            "report": report,
            "version": version,
            "output_path": output_path,
            "error": None,
        }


def run_cli() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="Starbase validation ETL pipeline")
    parser.add_argument("sqlite_path", help="Path to starbase.sqlite")
    parser.add_argument("--output", "-o", help="Output path for published DB")
    parser.add_argument("--apply", action="store_true", help="Write output (requires --output)")
    parser.add_argument("--config", "-c", help="Path to quality_rules.yaml")
    parser.add_argument("--previous", help="Path to previous published DB (for version bump)")
    parser.add_argument("--bump", default="minor", choices=("major", "minor", "patch"))
    parser.add_argument("--notes", default="", help="Notes stored in database_version")
    args = parser.parse_args()

    if args.apply and not args.output:
        parser.error("--apply requires --output /path/to/published.sqlite")

    # pylint: disable=import-outside-toplevel
    pipeline = ValidationPipeline(quality_rules_path=args.config)

    out = pipeline.run(
        sqlite_path=args.sqlite_path,
        output_path=args.output,
        dry_run=not args.apply,
        version_bump=args.bump,
        previous_published_db=args.previous,
        notes=args.notes,
    )

    print(out["report"])
    if out.get("error"):
        print(f"\nFAILED: {out['error']}", file=sys.stderr)
        return 1
    if args.apply:
        print(f"\nPublished: {out['version']} -> {out['output_path']}")
    else:
        print(f"\n[DRY RUN] Would publish version {out.get('version')}")
    return 0


if __name__ == "__main__":
    sys.exit(run_cli())
