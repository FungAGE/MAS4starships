"""
MetaEuk integration: discover ships needing annotation, run ``metaeuk easy-predict``,
and load raw GFF3 into :class:`~starship.models.StagingGffFeature`.
"""

from __future__ import annotations

import glob
import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from django.conf import settings

logger = logging.getLogger(__name__)


def ships_needing_annotation(
    joinedships_qs=None,
) -> List[Any]:
    """
    Return ``JoinedShips`` rows that:
    - have a linked ``ship_id`` with sequence
    - have no ``gff`` rows for that ship in Starbase
    - have no blocking :class:`~starship.models.GeneAnnotationRun` (pending, running,
      or completed) for the accession tag

    Failed runs are not blocking — a new run may be scheduled.

    ``joinedships_qs`` defaults to all joined ships (Starbase DB).
    """
    # Late imports so ``manage.py check`` works without full DB
    from starship.models import GeneAnnotationRun
    from starship.starbase_models import Gff, JoinedShips

    qs = joinedships_qs or JoinedShips.objects.all()
    qs = qs.select_related("ship_id", "accession_id")
    need: List[Any] = []

    for js in qs:
        if not js.ship_id_id or not js.accession_id_id:
            continue
        seq = getattr(js.ship_id, "sequence", None)
        if not seq or not str(seq).strip():
            continue
        if Gff.objects.filter(ship_id=js.ship_id_id).exists():
            continue
        atag = js.accession_id.accession_tag
        if GeneAnnotationRun.objects.filter(
            accession_tag=atag,
            status__in=("pending", "running", "completed"),
        ).exists():
            continue
        need.append(js)
    return need


def parse_gff3_file(path: str) -> List[Dict[str, Any]]:
    """Parse a GFF3 file into dicts suitable for :class:`StagingGffFeature`."""
    rows: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            score = parts[5]
            if score == ".":
                score = ""
            phase_raw = parts[7]
            phase: Optional[int]
            if phase_raw in ("", "."):
                phase = None
            else:
                try:
                    phase = int(phase_raw)
                except ValueError:
                    phase = None
            rows.append(
                {
                    "source": parts[1],
                    "feature_type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": score,
                    "strand": parts[6] if parts[6] != "." else "",
                    "phase": phase,
                    "attributes": parts[8],
                }
            )
    return rows


def store_staging_features(
    run,
    gff3_path: str,
    *,
    accession_tag: str,
    batch_size: int = 500,
) -> int:
    """
    Parse ``gff3_path`` and bulk-create :class:`~starship.models.StagingGffFeature`
    rows for ``run``.

    Returns:
        Number of rows created.
    """
    from starship.models import StagingGffFeature

    records = parse_gff3_file(gff3_path)
    if not records:
        return 0
    objs = []
    for r in records:
        objs.append(
            StagingGffFeature(
                run=run,
                accession_tag=accession_tag,
                source=r["source"],
                feature_type=r["feature_type"],
                start=r["start"],
                end=r["end"],
                strand=r["strand"][:1] if r["strand"] else "",
                phase=r["phase"],
                score=r["score"][:64] if r["score"] else "",
                attributes=r["attributes"] or "",
                review_status="pending",
            )
        )
    StagingGffFeature.objects.bulk_create(objs, batch_size=batch_size)
    return len(objs)


def _metaeuk_version(binary: str) -> str:
    try:
        p = subprocess.run(
            [binary, "-h"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        out = (p.stdout or "") + (p.stderr or "")
        if out:
            m = re.search(r"Version\s*:\s*(\S+)", out)
            if m:
                return m.group(1)
    except OSError:
        pass
    return ""


def run_metaeuk(
    query_fasta: str,
    target_db: str,
    work_dir: str,
    out_prefix: str,
    *,
    threads: Optional[int] = None,
    binary: Optional[str] = None,
) -> Tuple[str, str]:
    """
    Run ``metaeuk easy-predict`` (MetaEuk 4+).

    Args:
        query_fasta: Nucleotide FASTA path (one contig per ship is typical).
        target_db: MMseqs2/MetaEuk protein target database prefix.
        work_dir: Writable directory for temp and outputs.
        out_prefix: Output basename (no path) for MetaEuk products.

    Returns:
        ``(gff3_path, faa_path)`` — paths to primary GFF3 and FAA if found (else "").
    """
    binary = binary or getattr(settings, "METAEUK_BINARY", None) or "metaeuk"
    threads = threads if threads is not None else int(
        getattr(settings, "METAEUK_THREADS", 4)
    )
    os.makedirs(work_dir, exist_ok=True)
    prefix_path = os.path.join(work_dir, out_prefix)

    cmd = [
        binary,
        "easy-predict",
        os.path.abspath(query_fasta),
        os.path.abspath(target_db),
        os.path.abspath(work_dir),
        os.path.basename(prefix_path),
        "--threads",
        str(threads),
    ]
    logger.info("Running MetaEuk: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=work_dir)

    gff_candidates = sorted(glob.glob(os.path.join(work_dir, "**", "*.gff"), recursive=True))
    faa_candidates = sorted(
        glob.glob(os.path.join(work_dir, "**", "*.faa"), recursive=True)
    ) + sorted(glob.glob(os.path.join(work_dir, "**", "*.fas"), recursive=True))

    gff_path = gff_candidates[0] if gff_candidates else ""
    faa_path = faa_candidates[0] if faa_candidates else ""
    if not gff_path:
        # Some builds write TSV/GFF next to prefix
        for pat in (f"{prefix_path}*.gff", f"{prefix_path}*.gff3"):
            found = glob.glob(pat)
            if found:
                gff_path = found[0]
                break
    return gff_path, faa_path


def run_metaeuk_for_accession(
    accession_tag: str,
    target_db: Optional[str] = None,
    *,
    work_root: Optional[str] = None,
    celery_task_id: str = "",
) -> Tuple[object, int]:
    """
    Load sequence from Starbase, run MetaEuk, populate staging features.

    Creates and updates a :class:`~starship.models.GeneAnnotationRun`.

    Returns:
        ``(run, n_staging_features)``
    """
    from django.utils import timezone

    from starship.models import GeneAnnotationRun, StagingGffFeature
    from starship.starbase_models import Accessions, JoinedShips

    target_db = target_db or getattr(settings, "METAEUK_TARGET_DB", "") or ""
    if not target_db:
        raise ValueError(
            "MetaEuk target DB path is required (argument target_db or settings.METAEUK_TARGET_DB)."
        )

    acc = Accessions.objects.using("starbase").get(accession_tag=accession_tag)
    js = (
        JoinedShips.objects.using("starbase")
        .filter(accession_id=acc.id)
        .select_related("ship_id")
        .first()
    )
    if not js or not js.ship_id_id:
        raise ValueError(f"No joined_ships row with accession_tag={accession_tag!r}")
    seq = js.ship_id.sequence
    if not seq:
        raise ValueError(f"Ship has no sequence for accession_tag={accession_tag!r}")

    tool_version = _metaeuk_version(
        getattr(settings, "METAEUK_BINARY", None) or "metaeuk"
    )
    run = GeneAnnotationRun.objects.create(
        accession_tag=accession_tag,
        tool="metaeuk",
        tool_version=tool_version[:32],
        target_db=target_db[:512],
        status="running",
        celery_task_id=(celery_task_id or "")[:255],
    )
    tmp_root = work_root or tempfile.mkdtemp(prefix="metaeuk_")
    try:
        fasta = os.path.join(tmp_root, "query.fna")
        with open(fasta, "w", encoding="utf-8") as out:
            out.write(f">{accession_tag}\n{seq}\n")

        out_prefix = f"mas_{run.id}_{accession_tag}".replace(" ", "_")[:200]
        gff_path, faa_path = run_metaeuk(
            fasta,
            target_db,
            tmp_root,
            out_prefix,
        )
        if not gff_path or not os.path.isfile(gff_path):
            raise RuntimeError("MetaEuk finished but no .gff output was found.")

        persistent = os.path.join(
            settings.MEDIA_ROOT, "gene_annotation_runs", str(run.id)
        )
        os.makedirs(persistent, exist_ok=True)
        dest_gff = os.path.join(persistent, "metaeuk_output.gff")
        shutil.copy2(gff_path, dest_gff)
        run.gff3_path = dest_gff[:512]
        dest_faa = ""
        if faa_path and os.path.isfile(faa_path):
            dest_faa = os.path.join(persistent, "metaeuk_output.faa")
            shutil.copy2(faa_path, dest_faa)
            run.faa_path = dest_faa[:512]
        else:
            run.faa_path = ""
        run.save(update_fields=["gff3_path", "faa_path"])

        StagingGffFeature.objects.filter(run=run).delete()
        n = store_staging_features(run, dest_gff, accession_tag=accession_tag)
        run.status = "completed"
        run.completed_at = timezone.now()
        run.error_message = ""
        run.save(update_fields=["status", "completed_at", "error_message"])
        return run, n
    except Exception as exc:
        run.status = "failed"
        run.completed_at = timezone.now()
        run.error_message = str(exc)[:10000]
        run.save(update_fields=["status", "completed_at", "error_message"])
        raise
    finally:
        if work_root is None and os.path.isdir(tmp_root):
            try:
                shutil.rmtree(tmp_root)
            except OSError:
                logger.warning("Could not remove temp dir %s", tmp_root)


def starbase_gff_has_reviewed_run_column() -> bool:
    """Return True if Starbase SQLite ``gff`` has ``reviewed_run_id`` column."""
    from django.db import connections

    conn = connections["starbase"]
    with conn.cursor() as cursor:
        cursor.execute("PRAGMA table_info(gff)")
        cols = {row[1] for row in cursor.fetchall()}
    return "reviewed_run_id" in cols
