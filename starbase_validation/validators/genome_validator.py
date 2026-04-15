"""Genome validation — from quality_tags.check_genomic_context."""

from datetime import datetime
from typing import Any, Dict, Optional, Set, Union

from .results import ValidationIssue, ValidationResults


def _g(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_genome_record(genome: Union[Dict, Any]) -> list:
    if not genome:
        return ["missing_genome_context"]
    issues = []
    if not _g(genome, "assembly_accession"):
        issues.append("missing_assembly")
    if not _g(genome, "biosample"):
        issues.append("missing_biosample")
    tax = _g(genome, "taxonomy_id") or _g(genome, "tax_id")
    if tax is None:
        issues.append("broken_genome_taxonomy")
    ad = _g(genome, "acquisition_date")
    if ad:
        try:
            if datetime.now().year - int(ad) > 5:
                issues.append("old_assembly")
        except (ValueError, TypeError):
            pass
    return issues


def validate_genome(
    records: Dict,
    rules: Dict,
    blocking_issues: Optional[Set[str]] = None,
) -> ValidationResults:
    blocking = blocking_issues or {"broken_genome_taxonomy"}
    results = ValidationResults()
    genomes = records.get("genomes", [])
    by_id = {
        _g(g, "id") or _g(g, "genome_id"): g
        for g in genomes
        if (_g(g, "id") or _g(g, "genome_id")) is not None
    }
    joined = records.get("joined_ships", [])
    # One pass per genome id (avoid duplicate rows when many joined_ships share a genome)
    seen_genome_checks = set()

    for js in joined:
        gid = _g(js, "genome_id") or _g(js, "genome")
        if gid is None:
            continue
        st = _g(js, "starshipID") or _g(js, "starship_id")
        jid = _g(js, "id")
        g = by_id.get(gid)
        if g is None:
            results.add_issue(
                ValidationIssue(
                    "orphaned_genome_ref",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    f"genome {gid} missing",
                ),
                blocking,
            )
            continue
        if gid in seen_genome_checks:
            continue
        seen_genome_checks.add(gid)
        for it in check_genome_record(g):
            results.add_issue(
                ValidationIssue(
                    it,
                    "blocking" if it == "broken_genome_taxonomy" else "warning",
                    "genomes",
                    gid,
                    st,
                    it,
                ),
                blocking,
            )

    ref = {_g(js, "genome_id") or _g(js, "genome") for js in joined}
    ref.discard(None)
    for g in genomes:
        gid = _g(g, "id") or _g(g, "genome_id")
        if gid is not None and gid not in ref:
            results.add_issue(
                ValidationIssue(
                    "orphaned_genomes",
                    "warning",
                    "genomes",
                    gid,
                    None,
                    "Not referenced by joined_ships",
                ),
                blocking,
            )
    return results
