"""
Promote approved :class:`~starship.models.StagingGffFeature` rows into Starbase SQLite ``gff``.
"""

from __future__ import annotations

import logging
from typing import List, Optional, Sequence, Tuple

from django.db import connections, transaction
from django.utils import timezone

logger = logging.getLogger(__name__)


def promote_approved_staging_to_starbase(
    *,
    run_ids: Optional[Sequence[int]] = None,
    accession_tags: Optional[Sequence[str]] = None,
    dry_run: bool = False,
) -> Tuple[int, List[str]]:
    """
    Insert approved, not-yet-promoted staging features into Starbase ``gff``.

    Returns:
        ``(inserted_count, log_messages)``
    """
    from starship.metaeuk import starbase_gff_has_reviewed_run_column
    from starship.models import StagingGffFeature
    from starship.starbase_models import Accessions, JoinedShips

    qs = StagingGffFeature.objects.filter(
        review_status="approved",
        promoted_at__isnull=True,
    ).select_related("run")
    if run_ids:
        qs = qs.filter(run_id__in=run_ids)
    if accession_tags:
        qs = qs.filter(accession_tag__in=accession_tags)

    log: List[str] = []
    has_reviewed_col = starbase_gff_has_reviewed_run_column()
    if not has_reviewed_col:
        log.append(
            "Warning: gff.reviewed_run_id missing in Starbase SQLite. "
            "Run: python -m starbase_validation.cleanup.utils.migrate_gff_starbase_schema"
        )

    inserted = 0
    for row in qs.iterator(chunk_size=200):
        acc = Accessions.objects.using("starbase").filter(
            accession_tag=row.accession_tag
        ).first()
        if not acc:
            log.append(f"Skip id={row.id}: no accessions row for {row.accession_tag!r}")
            continue
        js = (
            JoinedShips.objects.using("starbase")
            .filter(accession_id=acc.id)
            .first()
        )
        if not js or not js.ship_id_id:
            log.append(
                f"Skip id={row.id}: no joined_ships with ship for {row.accession_tag!r}"
            )
            continue

        run_pk = row.run_id
        if dry_run:
            inserted += 1
            continue

        with transaction.atomic(using="starbase"):
            with connections["starbase"].cursor() as cursor:
                cursor.execute("SELECT COALESCE(MAX(id), 0) FROM gff")
                next_id = cursor.fetchone()[0] + 1
                phase = row.phase
                score = row.score or ""
                if has_reviewed_col:
                    cursor.execute(
                        """
                        INSERT INTO gff (
                            id, accession_id, source, "type", start, end, phase,
                            strand, score, attributes, ship_id, reviewed_run_id
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        [
                            next_id,
                            acc.id,
                            row.source,
                            row.feature_type,
                            row.start,
                            row.end,
                            phase,
                            row.strand or None,
                            score,
                            row.attributes or "",
                            js.ship_id_id,
                            run_pk,
                        ],
                    )
                else:
                    cursor.execute(
                        """
                        INSERT INTO gff (
                            id, accession_id, source, "type", start, end, phase,
                            strand, score, attributes, ship_id
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        [
                            next_id,
                            acc.id,
                            row.source,
                            row.feature_type,
                            row.start,
                            row.end,
                            phase,
                            row.strand or None,
                            score,
                            row.attributes or "",
                            js.ship_id_id,
                        ],
                    )
        row.promoted_at = timezone.now()
        row.save(update_fields=["promoted_at"])
        inserted += 1

    if dry_run:
        log.append(f"[dry-run] Would promote {inserted} staging rows.")
    else:
        log.append(f"Promoted {inserted} staging rows into Starbase gff.")
    return inserted, log
