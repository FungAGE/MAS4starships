"""
Main ETL pipeline orchestration.

Runs: Extract -> Validate -> (optional) Transform -> (optional) Load
"""

import argparse
import os
import sys
from pathlib import Path

# Ensure package is importable when run as script
_PKG_ROOT = Path(__file__).parent
_PKG_PARENT = _PKG_ROOT.parent
if str(_PKG_PARENT) not in sys.path:
    sys.path.insert(0, str(_PKG_PARENT))

import yaml

from starbase_validation.extractors import extract_from_sqlite
from starbase_validation.validators import validate_all
from starbase_validation.transformers import normalize_records
from starbase_validation.reporters import generate_report


def load_quality_rules(config_path: str = None) -> dict:
    """Load quality rules from YAML."""
    path = config_path or _PKG_ROOT / "config" / "quality_rules.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def run_pipeline(
    db_path: str,
    dry_run: bool = True,
    output_path: str = None,
    config_path: str = None,
) -> dict:
    """
    Run validation pipeline on SQLite database.

    Args:
        db_path: Path to starbase.sqlite
        dry_run: If True, only validate and report (no load)
        output_path: Optional path for output (ignored in dry_run)
        config_path: Optional path to quality_rules.yaml

    Returns:
        Dict with 'success', 'results', 'report'
    """
    if not os.path.exists(db_path):
        return {
            "success": False,
            "error": f"Database not found: {db_path}",
            "results": None,
            "report": None,
        }

    # 1. Extract
    records = extract_from_sqlite(db_path)

    # 2. Normalize
    records = normalize_records(records)

    # 3. Load rules and validate
    quality_rules = load_quality_rules(config_path)
    results = validate_all(records, quality_rules)

    # 4. Generate report
    report = generate_report(results)

    # 5. Check for blocking issues
    success = not results.has_blocking_issues()

    # 6. Load (if not dry_run and success)
    if success and not dry_run and output_path:
        # TODO: Implement in Phase 4
        pass

    return {
        "success": success,
        "results": results,
        "report": report,
        "blocking_count": len(results.blocking_issues),
        "warning_count": len(results.warning_issues),
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Starbase validation pipeline - validate SQLite database quality"
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        default="starbase.sqlite",
        help="Path to starbase SQLite database",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=True,
        help="Only validate, don't create output (default)",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply mode (create output if validation passes)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output path for validated database",
    )
    parser.add_argument(
        "--config",
        "-c",
        type=str,
        help="Path to quality_rules.yaml",
    )

    args = parser.parse_args()
    dry_run = not args.apply

    result = run_pipeline(
        db_path=args.db_path,
        dry_run=dry_run,
        output_path=args.output,
        config_path=args.config,
    )

    if result.get("error"):
        print(f"ERROR: {result['error']}", file=sys.stderr)
        sys.exit(1)

    print(result["report"])

    if not result["success"]:
        print(f"\nValidation FAILED: {result['blocking_count']} blocking issues found.")
        sys.exit(1)

    print("\nValidation PASSED.")
    if result["warning_count"]:
        print(f"  ({result['warning_count']} warnings - review recommended)")
    sys.exit(0)


if __name__ == "__main__":
    main()
