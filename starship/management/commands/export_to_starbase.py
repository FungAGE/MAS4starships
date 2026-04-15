"""
Export / publish Starbase SQLite using the quality-gate ValidationPipeline (Phase 5).

Reads the Starbase DB path from ``DATABASES['starbase']['NAME']``, runs
``starbase_validation.pipeline.ValidationPipeline``, optionally writes a versioned
SQLite file, and records :class:`~starship.models.ValidationRun` / ``ValidationIssue``
rows in MySQL (unless ``--no-tracking``).

Examples::

    # Validate only (default)
    python manage.py export_to_starbase

    # Publish to a file (requires passing validation)
    python manage.py export_to_starbase --apply --output /data/publish/starbase_v1.2.0.sqlite

    # Write issues CSV on failure
    python manage.py export_to_starbase --csv-export /tmp/issues.csv
"""

from __future__ import annotations

import os
import traceback

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from starship.models import ValidationRun


def _default_quality_rules_path() -> str:
    base = getattr(settings, "BASE_DIR", None)
    if base:
        p = os.path.join(base, "starbase_validation", "config", "quality_rules.yaml")
        if os.path.isfile(p):
            return p
    return os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
        "starbase_validation",
        "config",
        "quality_rules.yaml",
    )


class Command(BaseCommand):
    help = (
        "Run Starbase quality validation (ValidationPipeline) and optionally publish "
        "a versioned SQLite file. Uses DATABASES['starbase'] path."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--apply",
            action="store_true",
            help="Write output SQLite (requires --output). Default is dry-run.",
        )
        parser.add_argument(
            "--output",
            "-o",
            type=str,
            default=None,
            help="Destination path for published SQLite (e.g. starbase_v1.2.0.sqlite).",
        )
        parser.add_argument(
            "--config",
            "-c",
            type=str,
            default=None,
            help="Path to quality_rules.yaml (default: starbase_validation/config/quality_rules.yaml).",
        )
        parser.add_argument(
            "--previous",
            type=str,
            default=None,
            help="Path to previous published DB (for semver bump from database_version).",
        )
        parser.add_argument(
            "--bump",
            default="minor",
            choices=("major", "minor", "patch"),
            help="Semantic version bump when publishing.",
        )
        parser.add_argument(
            "--notes",
            default="",
            help="Notes stored in database_version metadata.",
        )
        parser.add_argument(
            "--sqlite-path",
            type=str,
            default=None,
            help="Override Starbase SQLite path (default: DATABASES['starbase']['NAME']).",
        )
        parser.add_argument(
            "--report",
            type=str,
            default=None,
            help="Write text validation report to this path.",
        )
        parser.add_argument(
            "--csv-export",
            type=str,
            default=None,
            help="On failure, write all issues to this CSV path.",
        )
        parser.add_argument(
            "--max-issues",
            type=int,
            default=2000,
            help="Max ValidationIssue rows to store in MySQL when tracking (default 2000).",
        )
        parser.add_argument(
            "--no-tracking",
            action="store_true",
            help="Do not create ValidationRun / ValidationIssue rows.",
        )
        parser.add_argument(
            "--skip-final-validation",
            action="store_true",
            help="Skip second validation pass after transform (same as pipeline).",
        )

    def handle(self, *args, **options):
        from starbase_validation.pipeline import ValidationPipeline
        from starbase_validation.reporters import (
            persist_issues_to_django_run,
            write_issues_csv,
        )

        sqlite_path = options["sqlite_path"]
        if not sqlite_path:
            sqlite_path = settings.DATABASES["starbase"]["NAME"]
        sqlite_path = os.path.abspath(sqlite_path)

        if not os.path.isfile(sqlite_path):
            raise CommandError(f"Starbase SQLite not found: {sqlite_path}")

        quality_rules_path = options["config"] or getattr(
            settings, "STARBASE_QUALITY_RULES_PATH", None
        ) or _default_quality_rules_path()
        if not os.path.isfile(quality_rules_path):
            self.stderr.write(
                self.style.WARNING(
                    f"Quality rules not found at {quality_rules_path}; using pipeline defaults."
                )
            )
            quality_rules_path = None

        dry_run = not options["apply"]
        output_path = options["output"]
        if options["apply"] and not output_path:
            raise CommandError("--apply requires --output /path/to/published.sqlite")

        if output_path:
            output_path = os.path.abspath(output_path)

        report_path = options["report"]
        if report_path:
            report_path = os.path.abspath(report_path)

        csv_export = options["csv_export"]
        if csv_export:
            csv_export = os.path.abspath(csv_export)

        no_tracking = options["no_tracking"]
        max_issues = options["max_issues"]

        run = None
        if not no_tracking:
            run = ValidationRun.objects.create(
                status="running",
                dry_run=dry_run,
                db_version="",
                report_path=report_path or "",
            )

        pipeline = ValidationPipeline(quality_rules_path=quality_rules_path)

        try:
            result = pipeline.run(
                sqlite_path=sqlite_path,
                output_path=output_path,
                dry_run=dry_run,
                skip_final_validation=options["skip_final_validation"],
                version_bump=options["bump"],
                previous_published_db=options["previous"],
                notes=options["notes"] or "",
            )
        except Exception as exc:
            if run:
                run.status = "failed"
                run.finished_at = timezone.now()
                run.error_message = f"{exc}\n{traceback.format_exc()}"
                run.save()
            self.stderr.write(self.style.ERROR(f"Pipeline error: {exc}"))
            raise

        report = result.get("report") or ""
        self.stdout.write(report)

        if report_path:
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(report)

        first = result.get("first_validation")
        second = result.get("second_validation")
        effective = second if second is not None else first

        if not result.get("success"):
            if csv_export and effective is not None:
                n = write_issues_csv(csv_export, effective)
                self.stdout.write(self.style.WARNING(f"Wrote {n} issues to {csv_export}"))

            if run:
                run.status = "failed"
                run.finished_at = timezone.now()
                run.total_issues = len(effective.issues) if effective else 0
                run.blocking_issues = len(effective.blocking_issues) if effective else 0
                run.error_message = result.get("error") or "Validation failed"
                run.save()
                if effective is not None:
                    n = persist_issues_to_django_run(
                        run, effective, max_issues=max_issues
                    )
                    self.stdout.write(
                        self.style.WARNING(f"Recorded {n} issues in MySQL (capped).")
                    )

            msg = result.get("error") or "Validation failed."
            self.stderr.write(self.style.ERROR(msg))
            raise CommandError(msg)

        # Success path
        version = result.get("version") or ""
        if run:
            run.status = "passed"
            run.finished_at = timezone.now()
            run.db_version = version
            run.total_issues = len(effective.issues) if effective else 0
            run.blocking_issues = len(effective.blocking_issues) if effective else 0
            run.error_message = ""
            run.save()

        if dry_run:
            self.stdout.write(
                self.style.SUCCESS(
                    f"[DRY RUN] Would publish version {version} (validation passed)."
                )
            )
        else:
            self.stdout.write(
                self.style.SUCCESS(
                    f"Published {version} -> {result.get('output_path')}"
                )
            )
