"""Sequence and boundary validation."""

import re
from typing import Any, Dict, List, Optional, Set, Union

from .results import ValidationIssue, ValidationResults


def _g(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_sequence_quality(
    sequence: str,
    md5: Optional[str],
    rev_comp_md5: Optional[str],
    md5_counts: Dict[str, int],
    rev_comp_counts: Dict[str, int],
) -> List[str]:
    issues = []
    if not sequence:
        return ["missing_sequence"]
    seq = str(sequence).upper()
    if re.search(r"[^ATCGN\-]", seq):
        issues.append("invalid_sequence")
    slen = len(seq.replace("-", "").replace("N", ""))
    if slen < 5000:
        issues.append("short_sequence")
    elif slen > 300000:
        issues.append("long_sequence")
    if md5 and md5_counts.get(md5, 0) > 1:
        issues.append("duplicate_sequence")
    if rev_comp_md5 and rev_comp_counts.get(rev_comp_md5, 0) > 1:
        issues.append("duplicate_rev_comp")
    return issues


def check_boundary_features(features: Union[Dict, Any]) -> List[str]:
    if not features:
        return ["missing_boundaries"]
    issues = []
    eb = _g(features, "elementBegin") or _g(features, "element_begin")
    ee = _g(features, "elementEnd") or _g(features, "element_end")
    if not eb or not ee:
        issues.append("missing_boundaries")
    if not (_g(features, "upDR") or _g(features, "up_dr")) or not (
        _g(features, "downDR") or _g(features, "down_dr")
    ):
        issues.append("missing_direct_repeats")
    if not (_g(features, "upTIR") or _g(features, "up_tir")) or not (
        _g(features, "downTIR") or _g(features, "down_tir")
    ):
        issues.append("missing_tir")
    es = _g(features, "emptySiteID") or _g(features, "empty_site_id")
    esq = _g(features, "emptySeq") or _g(features, "empty_seq")
    if not es or not esq:
        issues.append("missing_empty_site")
    return issues


def validate_sequence(
    records: Dict,
    rules: Dict,
    blocking_issues: Optional[Set[str]] = None,
) -> ValidationResults:
    blocking = blocking_issues or set()
    results = ValidationResults()
    ships = records.get("ships", [])
    joined = records.get("joined_ships", [])
    by_ship = {}
    for f in records.get("starship_features", []):
        sid = _g(f, "ship_id")
        if sid is not None:
            by_ship[sid] = f
    md5_counts, rc_counts = {}, {}
    for s in ships:
        m = _g(s, "md5")
        if m:
            md5_counts[m] = md5_counts.get(m, 0) + 1
        r = _g(s, "rev_comp_md5")
        if r:
            rc_counts[r] = rc_counts.get(r, 0) + 1

    for js in joined:
        ship_id = _g(js, "ship_id")
        st = _g(js, "starshipID") or _g(js, "starship_id")
        jid = _g(js, "id")
        if ship_id is None:
            results.add_issue(
                ValidationIssue(
                    "missing_sequence",
                    "warning",
                    "joined_ships",
                    jid,
                    st,
                    "No ship linked",
                ),
                blocking,
            )
            continue
        ship = next(
            (x for x in ships if (_g(x, "id") or _g(x, "ship_id")) == ship_id),
            None,
        )
        if not ship:
            results.add_issue(
                ValidationIssue(
                    "orphaned_ship_ref",
                    "blocking",
                    "joined_ships",
                    jid,
                    st,
                    f"ship_id {ship_id} missing",
                ),
                blocking,
            )
            continue
        for it in check_sequence_quality(
            _g(ship, "sequence"),
            _g(ship, "md5"),
            _g(ship, "rev_comp_md5"),
            md5_counts,
            rc_counts,
        ):
            results.add_issue(
                ValidationIssue(it, "warning", "ships", ship_id, st, it),
                blocking,
            )
        feat = by_ship.get(ship_id)
        for it in check_boundary_features(feat) if feat else ["missing_boundaries"]:
            results.add_issue(
                ValidationIssue(
                    it,
                    "warning",
                    "starship_features",
                    _g(feat, "id") if feat else None,
                    st,
                    it,
                ),
                blocking,
            )
    return results
