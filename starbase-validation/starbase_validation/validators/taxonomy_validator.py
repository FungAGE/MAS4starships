"""
Taxonomy validator - checks taxonomy completeness and consistency.

Extracted from cleanup/utils/quality_tags.py check_taxonomic_quality()
and quality_issues_inventory.md.
"""

import re
from typing import Dict, List, Any, Union

from .results import ValidationResults, ValidationIssue


# Uncertainty indicators in species names
UNCERTAIN_INDICATORS = ['sp.', 'cf.', 'aff.', 'uncertain', 'unknown']


def _get_val(obj: Any, key: str, default: Any = None) -> Any:
    """Get value from dict or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_taxonomy_record(
    taxonomy: Union[Dict, Any],
    taxonomy_id: int,
    rules: Dict,
) -> List[str]:
    """
    Check a single taxonomy record for quality issues.
    Pure function - no database access.

    Returns list of issue type strings.
    """
    issues = []

    if not taxonomy:
        return ['missing_taxonomy']

    genus = _get_val(taxonomy, 'genus')
    species = _get_val(taxonomy, 'species')
    strain = _get_val(taxonomy, 'strain')

    # Incomplete taxonomy
    if not genus or not species:
        issues.append('incomplete_taxonomy')

    # Missing strain (warning)
    if not strain and rules.get('require_strain', False):
        issues.append('missing_strain')

    # Uncertain taxonomy
    if species:
        species_str = str(species).lower()
        for indicator in UNCERTAIN_INDICATORS:
            if indicator in species_str:
                issues.append('taxonomy_uncertain')
                break

    return issues


def validate_taxonomy(
    records: Dict,
    rules: Dict,
    blocking_issues: set = None,
) -> ValidationResults:
    """
    Validate taxonomy data across all records.

    Args:
        records: Dict with 'taxonomy' list and 'joined_ships' (for tax_id refs)
        rules: Quality rules (max_missing_species, max_missing_strain, require_genus)
        blocking_issues: Set of issue types that block publication

    Returns:
        ValidationResults
    """
    blocking = blocking_issues or set()
    results = ValidationResults()

    taxonomy_list = records.get('taxonomy', [])
    joined_ships = records.get('joined_ships', [])

    # Build taxonomy lookup
    taxonomy_by_id = {}
    for t in taxonomy_list:
        tid = _get_val(t, 'id') or _get_val(t, 'taxonomy_id')
        if tid is not None:
            taxonomy_by_id[tid] = t

    # Check taxonomy referenced by joined_ships
    missing_taxonomy = 0
    incomplete_count = 0
    missing_strain_count = 0

    for js in joined_ships:
        tax_id = _get_val(js, 'tax_id') or _get_val(js, 'taxonomy_id')
        starship_id = _get_val(js, 'starshipID') or _get_val(js, 'starship_id')
        joined_id = _get_val(js, 'id') or _get_val(js, 'joined_ship_id')

        if tax_id is None:
            results.add_issue(
                ValidationIssue(
                    issue_type='missing_taxonomy',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details='Joined ship has no taxonomy reference',
                ),
                blocking,
            )
            missing_taxonomy += 1
            continue

        taxonomy = taxonomy_by_id.get(tax_id)
        if taxonomy is None:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_taxonomy_ref',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details=f'Taxonomy id {tax_id} not found',
                ),
                blocking,
            )
            continue

        issue_types = check_taxonomy_record(taxonomy, tax_id, rules)
        for issue_type in issue_types:
            results.add_issue(
                ValidationIssue(
                    issue_type=issue_type,
                    category='warning',
                    table_name='taxonomy',
                    record_id=tax_id,
                    starship_id=starship_id,
                    details=f'Taxonomy issue: {issue_type}',
                ),
                blocking,
            )
            if issue_type == 'incomplete_taxonomy':
                incomplete_count += 1
            elif issue_type == 'missing_strain':
                missing_strain_count += 1

    # Check for orphaned taxonomy (not referenced by joined_ships or genomes)
    genome_tax_ids = {_get_val(g, 'taxonomy_id') or _get_val(g, 'tax_id') 
                      for g in records.get('genomes', [])}
    ref_tax_ids = {_get_val(js, 'tax_id') or _get_val(js, 'taxonomy_id') 
                   for js in joined_ships}
    ref_tax_ids.update(genome_tax_ids)

    for t in taxonomy_list:
        tid = _get_val(t, 'id') or _get_val(t, 'taxonomy_id')
        if tid is not None and tid not in ref_tax_ids:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_taxonomy',
                    category='warning',
                    table_name='taxonomy',
                    record_id=tid,
                    details='Taxonomy not referenced by joined_ships or genomes',
                ),
                blocking,
            )

    total = len(joined_ships)
    if total > 0:
        results.summary['missing_taxonomy'] = results.summary.get('missing_taxonomy', 0) + missing_taxonomy
        results.summary['incomplete_taxonomy'] = results.summary.get('incomplete_taxonomy', 0) + incomplete_count
        results.summary['missing_strain'] = results.summary.get('missing_strain', 0) + missing_strain_count

    return results
