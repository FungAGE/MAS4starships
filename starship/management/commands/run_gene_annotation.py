"""
Enqueue or run MetaEuk gene annotation for Starbase ships that lack ``gff`` rows.

Examples::

    python manage.py run_gene_annotation --dry-run
    python manage.py run_gene_annotation --target-db /data/uniref50 --sync
    python manage.py run_gene_annotation --accession-tags TAG1 TAG2 --sync
"""

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from starship.metaeuk import ships_needing_annotation


class Command(BaseCommand):
    help = (
        "Find joined ships without Starbase gff rows and run MetaEuk "
        "(Celery async or --sync)."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--tool",
            default="metaeuk",
            choices=("metaeuk",),
            help="Annotation tool (only metaeuk is implemented).",
        )
        parser.add_argument(
            "--target-db",
            type=str,
            default=None,
            help="MetaEuk protein target DB path (default: settings.METAEUK_TARGET_DB).",
        )
        parser.add_argument(
            "--accession-tags",
            nargs="*",
            default=None,
            help="Restrict to these Starbase accession tags (optional).",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="List accession tags that would be processed; do not run tasks.",
        )
        parser.add_argument(
            "--sync",
            action="store_true",
            help="Run MetaEuk in-process instead of Celery (for debugging).",
        )

    def handle(self, *args, **options):
        if options["tool"] != "metaeuk":
            raise CommandError("Only --tool metaeuk is supported.")

        target_db = options["target_db"] or getattr(
            settings, "METAEUK_TARGET_DB", ""
        )
        if not target_db and not options["dry_run"]:
            raise CommandError(
                "Pass --target-db or set METAEUK_TARGET_DB / settings.METAEUK_TARGET_DB."
            )

        qs = None
        tags_filter = options["accession_tags"]
        if tags_filter:
            from starship.starbase_models import JoinedShips

            qs = JoinedShips.objects.filter(
                accession_id__accession_tag__in=tags_filter
            )

        need = ships_needing_annotation(joinedships_qs=qs)
        tags = []
        for js in need:
            if js.accession_id_id:
                tags.append(js.accession_id.accession_tag)

        if options["dry_run"]:
            self.stdout.write(
                self.style.NOTICE(
                    f"[dry-run] Would process {len(tags)} accession(s): {', '.join(tags) or '(none)'}"
                )
            )
            return

        if not tags:
            self.stdout.write(self.style.WARNING("No ships need annotation."))
            return

        if options["sync"]:
            from starship import metaeuk

            for tag in tags:
                self.stdout.write(f"Running MetaEuk for {tag!r}...")
                run, n = metaeuk.run_metaeuk_for_accession(tag, target_db=target_db)
                self.stdout.write(
                    self.style.SUCCESS(
                        f"Done run id={run.id}, staging_features={n}"
                    )
                )
            return

        from starship.tasks import run_metaeuk_annotation

        for tag in tags:
            run_metaeuk_annotation.delay(tag, target_db=target_db)
            self.stdout.write(f"Queued MetaEuk for {tag!r}")

        self.stdout.write(
            self.style.SUCCESS(f"Queued {len(tags)} MetaEuk task(s).")
        )
