"""
Classification validator - checks family/navis/haplotype completeness.

Extracted from cleanup/utils/quality_tags.py check_classification_hierarchy()
and quality_issues_inventory.md - missing_family is BLOCKING.
"""

from typing import Dict, Any, Union

from .results import ValidationResults, ValidationIssue


def _get_val(obj: Any, key: str, default: Any = None) -> Any:
    """Get value from dict or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def validate_classification(
    records: Dict,
    rules: Dict,
    blocking_issues: set = None,
) -> ValidationResults:
    """
    Validate classification (family, navis, haplotype) across joined_ships.

    Args:
        records: Dict with 'joined_ships', optionally 'navis', 'haplotype'
        rules: Quality rules (max_missing_family: 0 for blocking)
        blocking_issues: Set of issue types that block publication

    Returns:
        ValidationResults
    """
    blocking = blocking_issues or {'missing_family'}
    results = ValidationResults()

    joined_ships = records.get('joined_ships', [])
    navis_by_id = {_get_val(n, 'id'): n for n in records.get('navis', []) if _get_val(n, 'id') is not None}
    haplotype_by_id = {_get_val(h, 'id'): h for h in records.get('haplotype', []) if _get_val(h, 'id') is not None}

    for js in joined_ships:
        ship_family_id = _get_val(js, 'ship_family_id') or _get_val(js, 'family_id')
        navis_id = _get_val(js, 'ship_navis_id') or _get_val(js, 'navis_id')
        haplotype_id = _get_val(js, 'ship_haplotype_id') or _get_val(js, 'haplotype_id')
        starship_id = _get_val(js, 'starshipID') or _get_val(js, 'starship_id')
        joined_id = _get_val(js, 'id') or _get_val(js, 'joined_ship_id')

        # Missing family - BLOCKING
        if ship_family_id is None:
            results.add_issue(
                ValidationIssue(
                    issue_type='missing_family',
                    category='blocking',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details='Ship has no family classification',
                ),
                blocking,
            )

        # Missing navis - warning
        if navis_id is None and ship_family_id is not None:
            results.add_issue(
                ValidationIssue(
                    issue_type='missing_navis',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details='Ship has no navis classification',
                ),
                blocking,
            )
        elif navis_id and navis_id in navis_by_id:
            navis = navis_by_id[navis_id]
            if _get_val(navis, 'previous_navis_name'):
                results.add_issue(
                    ValidationIssue(
                        issue_type='classification_outdated',
                        category='warning',
                        table_name='joined_ships',
                        record_id=joined_id,
                        starship_id=starship_id,
                        details='Using outdated navis name',
                    ),
                    blocking,
                )

        # Missing haplotype - warning
        if haplotype_id is None and ship_family_id is not None:
            results.add_issue(
                ValidationIssue(
                    issue_type='missing_haplotype',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details='Ship has no haplotype classification',
                ),
                blocking,
            )
        elif haplotype_id and haplotype_id in haplotype_by_id:
            haplotype = haplotype_by_id[haplotype_id]
            if _get_val(haplotype, 'previous_haplotype_name'):
                results.add_issue(
                    ValidationIssue(
                        issue_type='classification_outdated',
                        category='warning',
                        table_name='joined_ships',
                        record_id=joined_id,
                        starship_id=starship_id,
                        details='Using outdated haplotype name',
                    ),
                    blocking,
                )

    return results
