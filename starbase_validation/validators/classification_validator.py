"""Classification (family/navis/haplotype) validation."""

from typing import Any, Dict, Optional, Set

from .results import ValidationIssue, ValidationResults


def _g(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def validate_classification(
    records: Dict,
    rules: Dict,
    blocking_issues: Optional[Set[str]] = None,
) -> ValidationResults:
    blocking = blocking_issues or {"missing_family"}
    results = ValidationResults()
    navis = {_g(n, "id"): n for n in records.get("navis", []) if _g(n, "id") is not None}
    hap = {_g(h, "id"): h for h in records.get("haplotype", []) if _g(h, "id") is not None}

    for js in records.get("joined_ships", []):
        fam = _g(js, "ship_family_id") or _g(js, "family_id")
        nid = _g(js, "ship_navis_id") or _g(js, "navis_id")
        hid = _g(js, "ship_haplotype_id") or _g(js, "haplotype_id")
        st = _g(js, "starshipID") or _g(js, "starship_id")
        jid = _g(js, "id")

        if fam is None:
            results.add_issue(
                ValidationIssue(
                    "missing_family",
                    "blocking",
                    "joined_ships",
                    jid,
                    st,
                    "No family classification",
                ),
                blocking,
            )
        if fam is not None and nid is None:
            results.add_issue(
                ValidationIssue(
                    "missing_navis",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    "No navis",
                ),
                blocking,
            )
        elif nid and nid in navis and _g(navis[nid], "previous_navis_name"):
            results.add_issue(
                ValidationIssue(
                    "classification_outdated",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    "Outdated navis name",
                ),
                blocking,
            )
        if fam is not None and hid is None:
            results.add_issue(
                ValidationIssue(
                    "missing_haplotype",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    "No haplotype",
                ),
                blocking,
            )
        elif hid and hid in hap and _g(hap[hid], "previous_haplotype_name"):
            results.add_issue(
                ValidationIssue(
                    "classification_outdated",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    "Outdated haplotype name",
                ),
                blocking,
            )
    return results
