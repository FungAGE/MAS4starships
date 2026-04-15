"""Taxonomy validation — extracted from cleanup/utils/quality_tags.py."""

from typing import Any, Dict, List, Optional, Set, Union

from .results import ValidationIssue, ValidationResults

UNCERTAIN = ("sp.", "cf.", "aff.", "uncertain", "unknown")


def _g(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_taxonomy_record(taxonomy: Union[Dict, Any], rules: Dict) -> List[str]:
    if not taxonomy:
        return ["missing_taxonomy"]
    issues = []
    genus, species = _g(taxonomy, "genus"), _g(taxonomy, "species")
    if not genus or not species:
        issues.append("incomplete_taxonomy")
    strain = _g(taxonomy, "strain")
    if not strain and rules.get("require_strain"):
        issues.append("missing_strain")
    if species:
        s = str(species).lower()
        for ind in UNCERTAIN:
            if ind in s:
                issues.append("taxonomy_uncertain")
                break
    return issues


def validate_taxonomy(
    records: Dict,
    rules: Dict,
    blocking_issues: Optional[Set[str]] = None,
) -> ValidationResults:
    blocking = blocking_issues or set()
    results = ValidationResults()
    taxonomy_list = records.get("taxonomy", [])
    joined_ships = records.get("joined_ships", [])
    by_id = {}
    for t in taxonomy_list:
        tid = _g(t, "id") or _g(t, "taxonomy_id")
        if tid is not None:
            by_id[tid] = t
    ref_tax = set()
    ref_tax.update(_g(js, "tax_id") or _g(js, "taxonomy_id") for js in joined_ships)
    ref_tax.update(_g(g, "taxonomy_id") or _g(g, "tax_id") for g in records.get("genomes", []))
    ref_tax.discard(None)

    for js in joined_ships:
        tax_id = _g(js, "tax_id") or _g(js, "taxonomy_id")
        sid = _g(js, "starshipID") or _g(js, "starship_id")
        jid = _g(js, "id")
        if tax_id is None:
            results.add_issue(
                ValidationIssue(
                    "missing_taxonomy",
                    "warning",
                    "joined_ships",
                    jid,
                    sid,
                    "No taxonomy reference",
                ),
                blocking,
            )
            continue
        tax = by_id.get(tax_id)
        if tax is None:
            results.add_issue(
                ValidationIssue(
                    "orphaned_taxonomy_ref",
                    "warning",
                    "joined_ships",
                    jid,
                    sid,
                    f"taxonomy id {tax_id} missing",
                ),
                blocking,
            )
            continue
        for it in check_taxonomy_record(tax, rules):
            results.add_issue(
                ValidationIssue(
                    it,
                    "warning",
                    "taxonomy",
                    tax_id,
                    sid,
                    it,
                ),
                blocking,
            )

    for t in taxonomy_list:
        tid = _g(t, "id") or _g(t, "taxonomy_id")
        if tid is not None and tid not in ref_tax:
            results.add_issue(
                ValidationIssue(
                    "orphaned_taxonomy",
                    "warning",
                    "taxonomy",
                    tid,
                    None,
                    "Not referenced by joined_ships or genomes",
                ),
                blocking,
            )
    return results
