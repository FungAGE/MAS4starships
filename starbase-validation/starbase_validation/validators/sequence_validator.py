"""
Sequence validator - checks sequence quality and boundaries.

Extracted from cleanup/utils/quality_tags.py check_sequence_quality(),
check_boundary_features(), and check_annotation_quality().
"""

import re
from typing import Dict, Any, List, Union

from .results import ValidationResults, ValidationIssue


def _get_val(obj: Any, key: str, default: Any = None) -> Any:
    """Get value from dict or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_sequence_quality(sequence: str, md5: str = None, rev_comp_md5: str = None,
                          md5_counts: Dict[str, int] = None,
                          rev_comp_counts: Dict[str, int] = None,
                          ship_id: int = None) -> List[str]:
    """
    Check sequence-related quality issues. Pure function.

    Args:
        sequence: DNA sequence
        md5: MD5 hash of sequence
        rev_comp_md5: MD5 of reverse complement
        md5_counts: Optional dict md5 -> count (for duplicate detection)
        rev_comp_counts: Optional dict rev_comp_md5 -> count
        ship_id: For duplicate check exclusion

    Returns:
        List of issue type strings
    """
    issues = []

    if not sequence:
        return ['missing_sequence']

    seq_upper = str(sequence).upper()
    if re.search(r'[^ATCGN\-]', seq_upper):
        issues.append('invalid_sequence')

    seq_length = len(seq_upper.replace('-', '').replace('N', ''))
    if seq_length < 5000:
        issues.append('short_sequence')
    elif seq_length > 300000:
        issues.append('long_sequence')

    if md5 and md5_counts and md5_counts.get(md5, 0) > 1:
        issues.append('duplicate_sequence')
    if rev_comp_md5 and rev_comp_counts and rev_comp_counts.get(rev_comp_md5, 0) > 1:
        issues.append('duplicate_rev_comp')

    return issues


def check_boundary_features(features: Union[Dict, Any]) -> List[str]:
    """
    Check boundary and feature-related quality issues. Pure function.
    """
    issues = []

    if not features:
        return ['missing_boundaries']

    elem_begin = _get_val(features, 'elementBegin') or _get_val(features, 'element_begin')
    elem_end = _get_val(features, 'elementEnd') or _get_val(features, 'element_end')
    up_dr = _get_val(features, 'upDR') or _get_val(features, 'up_dr')
    down_dr = _get_val(features, 'downDR') or _get_val(features, 'down_dr')
    up_tir = _get_val(features, 'upTIR') or _get_val(features, 'up_tir')
    down_tir = _get_val(features, 'downTIR') or _get_val(features, 'down_tir')
    dr_edit = _get_val(features, 'DRedit') or _get_val(features, 'dr_edit')
    tir_edit = _get_val(features, 'TIRedit') or _get_val(features, 'tir_edit')
    empty_site = _get_val(features, 'emptySiteID') or _get_val(features, 'empty_site_id')
    empty_seq = _get_val(features, 'emptySeq') or _get_val(features, 'empty_seq')

    if not elem_begin or not elem_end:
        issues.append('missing_boundaries')
    if not up_dr or not down_dr:
        issues.append('missing_direct_repeats')
    if not up_tir or not down_tir:
        issues.append('missing_tir')
    if dr_edit or tir_edit:
        issues.append('boundary_edited')
        if dr_edit and 'inconsistent' in str(dr_edit).lower():
            issues.append('boundary_inconsistent')
        if tir_edit and 'inconsistent' in str(tir_edit).lower():
            issues.append('boundary_inconsistent')
    if not empty_site or not empty_seq:
        issues.append('missing_empty_site')

    return issues


def validate_sequence(
    records: Dict,
    rules: Dict,
    blocking_issues: set = None,
) -> ValidationResults:
    """
    Validate sequence and boundary data.

    Args:
        records: Dict with 'ships', 'joined_ships', 'starship_features', 'gff'
        rules: Quality rules
        blocking_issues: Set of blocking issue types

    Returns:
        ValidationResults
    """
    blocking = blocking_issues or set()
    results = ValidationResults()

    ships = records.get('ships', [])
    joined_ships = records.get('joined_ships', [])
    features_by_ship = {}
    for f in records.get('starship_features', []):
        sid = _get_val(f, 'ship_id')
        if sid is not None:
            features_by_ship[sid] = f

    # Build MD5 counts for duplicate detection
    md5_counts = {}
    rev_comp_counts = {}
    for s in ships:
        m = _get_val(s, 'md5')
        if m:
            md5_counts[m] = md5_counts.get(m, 0) + 1
        r = _get_val(s, 'rev_comp_md5')
        if r:
            rev_comp_counts[r] = rev_comp_counts.get(r, 0) + 1

    for js in joined_ships:
        ship_id = _get_val(js, 'ship_id')
        starship_id = _get_val(js, 'starshipID') or _get_val(js, 'starship_id')
        joined_id = _get_val(js, 'id') or _get_val(js, 'joined_ship_id')

        if ship_id is None:
            results.add_issue(
                ValidationIssue(
                    issue_type='missing_sequence',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details='No ship/sequence linked',
                ),
                blocking,
            )
            continue

        ship = next((s for s in ships if (_get_val(s, 'id') or _get_val(s, 'ship_id')) == ship_id), None)
        if not ship:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_ship_ref',
                    category='blocking',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details=f'Ship id {ship_id} not found',
                ),
                blocking,
            )
            continue

        seq = _get_val(ship, 'sequence')
        seq_issues = check_sequence_quality(
            seq,
            _get_val(ship, 'md5'),
            _get_val(ship, 'rev_comp_md5'),
            md5_counts,
            rev_comp_counts,
        )
        for issue_type in seq_issues:
            results.add_issue(
                ValidationIssue(
                    issue_type=issue_type,
                    category='warning',
                    table_name='ships',
                    record_id=ship_id,
                    starship_id=starship_id,
                    details=f'Sequence issue: {issue_type}',
                ),
                blocking,
            )

        features = features_by_ship.get(ship_id)
        boundary_issues = check_boundary_features(features) if features else ['missing_boundaries']
        for issue_type in boundary_issues:
            results.add_issue(
                ValidationIssue(
                    issue_type=issue_type,
                    category='warning',
                    table_name='starship_features',
                    record_id=_get_val(features, 'id') if features else None,
                    starship_id=starship_id,
                    details=f'Boundary issue: {issue_type}',
                ),
                blocking,
            )

    return results
