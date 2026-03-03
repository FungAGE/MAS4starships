"""
Genome validator - checks genome metadata completeness.

Extracted from cleanup/utils/quality_tags.py check_genomic_context()
and quality_issues_inventory.md.
"""

from datetime import datetime
from typing import Dict, Any, Union

from .results import ValidationResults, ValidationIssue


def _get_val(obj: Any, key: str, default: Any = None) -> Any:
    """Get value from dict or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def check_genome_record(genome: Union[Dict, Any], genome_id: int) -> list:
    """
    Check a single genome record for quality issues.
    Pure function - no database access.

    Returns list of issue type strings.
    """
    issues = []

    if not genome:
        return ['missing_genome_context']

    assembly = _get_val(genome, 'assembly_accession')
    biosample = _get_val(genome, 'biosample')
    taxonomy_id = _get_val(genome, 'taxonomy_id') or _get_val(genome, 'tax_id')
    acquisition_date = _get_val(genome, 'acquisition_date')

    if not assembly:
        issues.append('missing_assembly')
    if not biosample:
        issues.append('missing_biosample')
    if taxonomy_id is None:
        issues.append('broken_genome_taxonomy')  # BLOCKING

    # Old assembly (>5 years)
    if acquisition_date:
        try:
            year = int(acquisition_date)
            if datetime.now().year - year > 5:
                issues.append('old_assembly')
        except (ValueError, TypeError):
            pass

    return issues


def validate_genome(
    records: Dict,
    rules: Dict,
    blocking_issues: set = None,
) -> ValidationResults:
    """
    Validate genome data across all records.

    Args:
        records: Dict with 'genomes' list and 'joined_ships' (for genome_id refs)
        rules: Quality rules (require_taxonomy_link, max_missing_assembly)
        blocking_issues: Set of issue types that block publication

    Returns:
        ValidationResults
    """
    blocking = blocking_issues or {'broken_genome_taxonomy'}
    results = ValidationResults()

    genomes = records.get('genomes', [])
    genome_by_id = {_get_val(g, 'id') or _get_val(g, 'genome_id'): g for g in genomes if (_get_val(g, 'id') or _get_val(g, 'genome_id')) is not None}
    joined_ships = records.get('joined_ships', [])

    for js in joined_ships:
        genome_id = _get_val(js, 'genome_id') or _get_val(js, 'genome')
        starship_id = _get_val(js, 'starshipID') or _get_val(js, 'starship_id')
        joined_id = _get_val(js, 'id') or _get_val(js, 'joined_ship_id')

        if genome_id is None:
            continue  # No genome ref - might be OK

        genome = genome_by_id.get(genome_id)
        if genome is None:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_genome_ref',
                    category='warning',
                    table_name='joined_ships',
                    record_id=joined_id,
                    starship_id=starship_id,
                    details=f'Genome id {genome_id} not found',
                ),
                blocking,
            )
            continue

        issue_types = check_genome_record(genome, genome_id)
        for issue_type in issue_types:
            results.add_issue(
                ValidationIssue(
                    issue_type=issue_type,
                    category='blocking' if issue_type == 'broken_genome_taxonomy' else 'warning',
                    table_name='genomes',
                    record_id=genome_id,
                    starship_id=starship_id,
                    details=f'Genome issue: {issue_type}',
                ),
                blocking,
            )

    # Orphaned genomes (not referenced by joined_ships)
    ref_genome_ids = {_get_val(js, 'genome_id') or _get_val(js, 'genome') for js in joined_ships}
    for g in genomes:
        gid = _get_val(g, 'id') or _get_val(g, 'genome_id')
        if gid is not None and gid not in ref_genome_ids:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_genomes',
                    category='warning',
                    table_name='genomes',
                    record_id=gid,
                    details='Genome not referenced by joined_ships',
                ),
                blocking,
            )

    return results
