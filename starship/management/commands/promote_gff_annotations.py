"""
Promote approved :class:`~starship.models.StagingGffFeature` rows into Starbase SQLite ``gff``.

Ensure Starbase schema includes ``gff.reviewed_run_id`` (see
``starbase_validation.cleanup.utils.migrate_gff_starbase_schema``).

Examples::

    python manage.py promote_gff_annotations --dry-run
    python manage.py promote_gff_annotations --run-id 1
    python manage.py promote_gff_annotations --accession-tags MYTAG
"""

from django.core.management.base import BaseCommand

from starship.gff_promotion import promote_approved_staging_to_starbase


class Command(BaseCommand):
    help = (
        "Insert approved staging GFF features into Starbase SQLite gff "
        "(publication layer)."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--run-id",
            type=int,
            action="append",
            dest="run_ids",
            help="GeneAnnotationRun id (repeatable).",
        )
        parser.add_argument(
            "--accession-tags",
            nargs="*",
            default=None,
            help="Only staging rows with these accession_tag values.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Count rows that would be promoted without writing.",
        )

    def handle(self, *args, **options):
        run_ids = options["run_ids"]
        accession_tags = options["accession_tags"]
        dry_run = options["dry_run"]

        n, log = promote_approved_staging_to_starbase(
            run_ids=run_ids,
            accession_tags=accession_tags,
            dry_run=dry_run,
        )
        for line in log:
            self.stdout.write(line)
        if dry_run:
            self.stdout.write(self.style.NOTICE(f"[dry-run] count={n}"))
        else:
            self.stdout.write(self.style.SUCCESS(f"Inserted {n} gff row(s)."))
