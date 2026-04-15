"""Foreign key and relationship integrity."""

from typing import Any, Dict, Optional, Set

from .results import ValidationIssue, ValidationResults


def _g(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def validate_relationships(
    records: Dict,
    blocking_issues: Optional[Set[str]] = None,
) -> ValidationResults:
    blocking = blocking_issues or {
        "orphaned_foreign_keys",
        "duplicate_accessions",
        "broken_genome_taxonomy",
    }
    results = ValidationResults()
    ships = {
        _g(s, "id") or _g(s, "ship_id"): s
        for s in records.get("ships", [])
        if (_g(s, "id") or _g(s, "ship_id")) is not None
    }
    acc = {
        _g(a, "id") or _g(a, "accession_id"): a
        for a in records.get("accessions", [])
        if (_g(a, "id") or _g(a, "accession_id")) is not None
    }
    caps = {
        _g(c, "id") or _g(c, "captain_id"): c
        for c in records.get("captains", [])
        if (_g(c, "id") or _g(c, "captain_id")) is not None
    }
    genomes = {
        _g(g, "id") or _g(g, "genome_id"): g
        for g in records.get("genomes", [])
        if (_g(g, "id") or _g(g, "genome_id")) is not None
    }
    tax = {
        _g(t, "id") or _g(t, "taxonomy_id"): t
        for t in records.get("taxonomy", [])
        if (_g(t, "id") or _g(t, "taxonomy_id")) is not None
    }

    for js in records.get("joined_ships", []):
        jid = _g(js, "id")
        st = _g(js, "starshipID") or _g(js, "starship_id")
        sid = _g(js, "ship_id")
        if sid is not None and sid not in ships:
            results.add_issue(
                ValidationIssue(
                    "orphaned_foreign_keys",
                    "blocking",
                    "joined_ships",
                    jid,
                    st,
                    f"ship_id {sid} invalid",
                ),
                blocking,
            )
        cid = _g(js, "captain_id")
        if cid is not None and cid not in caps:
            results.add_issue(
                ValidationIssue(
                    "orphaned_foreign_keys",
                    "blocking",
                    "joined_ships",
                    jid,
                    st,
                    f"captain_id {cid} invalid",
                ),
                blocking,
            )
        gid = _g(js, "genome_id")
        if gid is not None and gid not in genomes:
            results.add_issue(
                ValidationIssue(
                    "orphaned_foreign_keys",
                    "blocking",
                    "joined_ships",
                    jid,
                    st,
                    f"genome_id {gid} invalid",
                ),
                blocking,
            )

    for g in records.get("genomes", []):
        tid = _g(g, "taxonomy_id") or _g(g, "tax_id")
        gid = _g(g, "id")
        if tid is not None and tid not in tax:
            results.add_issue(
                ValidationIssue(
                    "broken_genome_taxonomy",
                    "blocking",
                    "genomes",
                    gid,
                    None,
                    f"taxonomy_id {tid} invalid",
                ),
                blocking,
            )

    for f in records.get("starship_features", []):
        sid = _g(f, "ship_id")
        if sid is not None and sid not in ships:
            results.add_issue(
                ValidationIssue(
                    "orphaned_features",
                    "warning",
                    "starship_features",
                    _g(f, "id"),
                    None,
                    f"ship_id {sid} invalid",
                ),
                blocking,
            )

    for gf in records.get("gff", []):
        sid = _g(gf, "ship_id")
        if sid is not None and sid not in ships:
            results.add_issue(
                ValidationIssue(
                    "orphaned_gff",
                    "warning",
                    "gff",
                    _g(gf, "id"),
                    None,
                    f"ship_id {sid} invalid",
                ),
                blocking,
            )

    tags = {}
    for a in records.get("accessions", []):
        t = _g(a, "accession_tag")
        if t:
            tags[t] = tags.get(t, 0) + 1
    for tag, cnt in tags.items():
        if cnt > 1:
            results.add_issue(
                ValidationIssue(
                    "duplicate_accessions",
                    "blocking",
                    "accessions",
                    None,
                    None,
                    f"duplicate tag {tag}",
                ),
                blocking,
            )

    return results
