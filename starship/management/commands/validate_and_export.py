"""
Run the Starbase comprehensive validation / cleanup pipeline (SQLite starbase DB).

Wraps ``starbase_validation.cleanup.run_comprehensive_cleanup.main`` and records
a :class:`starship.models.ValidationRun` in the default (MySQL) database.
"""

import os
import traceback

from django.core.management.base import BaseCommand
from django.utils import timezone

from starship.models import ValidationRun


class Command(BaseCommand):
    help = (
        "Run Starbase SQLite validation/cleanup (comprehensive pipeline). "
        "Uses DATABASES['starbase'] path via starbase_validation config."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--apply",
            action="store_true",
            help="Run with fixes enabled (not dry-run for main workflow steps).",
        )
        parser.add_argument(
            "--apply-fixes",
            action="store_true",
            help="Apply database fixes after analysis (same as comprehensive script).",
        )
        parser.add_argument(
            "--report",
            type=str,
            default=None,
            help="Path to write the text report (optional).",
        )
        parser.add_argument(
            "--skip-accessions",
            action="store_true",
            help="Skip expensive accession cleanup step.",
        )
        parser.add_argument(
            "--no-record-issues",
            action="store_true",
            help="Do not write rows to cleanup_issues table in SQLite.",
        )
        parser.add_argument(
            "--fasta",
            type=str,
            default=None,
            help="Optional FASTA path for ship_id validation.",
        )
        parser.add_argument(
            "--db-version",
            type=str,
            default="",
            help="Label for this run (e.g. semantic version of the published DB).",
        )
        parser.add_argument(
            "--no-tracking",
            action="store_true",
            help="Do not create ValidationRun / ValidationIssue rows in MySQL.",
        )

    def handle(self, *args, **options):
        dry_run = not options["apply"]
        apply_fixes = options["apply_fixes"]
        report_path = options["report"]
        skip_accessions = options["skip_accessions"]
        no_record_issues = options["no_record_issues"]
        fasta_path = options["fasta"]
        db_version = options["db_version"] or ""
        no_tracking = options["no_tracking"]

        if report_path:
            report_path = os.path.abspath(report_path)

        run = None
        if not no_tracking:
            run = ValidationRun.objects.create(
                status="running",
                dry_run=dry_run,
                db_version=db_version,
                report_path=report_path or "",
            )

        try:
            from starbase_validation.cleanup.run_comprehensive_cleanup import main

            # main() uses dry_run=True by default; --apply sets dry_run False
            main(
                dry_run=dry_run,
                output_report=report_path,
                apply_fixes=apply_fixes,
                skip_accessions=skip_accessions,
                no_record_issues=no_record_issues,
                fasta_path=fasta_path,
            )

            if run:
                run.status = "passed"
                run.finished_at = timezone.now()
                run.total_issues, run.blocking_issues = 0, 0
                run.save()
            self.stdout.write(self.style.SUCCESS("Validation pipeline completed."))

        except Exception as exc:
            if run:
                run.status = "failed"
                run.finished_at = timezone.now()
                run.error_message = f"{exc}\n{traceback.format_exc()}"
                run.save()
            self.stderr.write(self.style.ERROR(f"Validation pipeline failed: {exc}"))
            raise
