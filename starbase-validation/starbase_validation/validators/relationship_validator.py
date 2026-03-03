"""
Relationship validator - checks foreign key integrity.

Extracted from cleanup/utils/database_cleanup.py and quality_issues_inventory.md.
Checks: orphaned_foreign_keys, ship_id_mismatches, duplicate_accessions.
"""

from typing import Dict, Any

from .results import ValidationResults, ValidationIssue


def _get_val(obj: Any, key: str, default: Any = None) -> Any:
    """Get value from dict or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def validate_relationships(
    records: Dict,
    blocking_issues: set = None,
) -> ValidationResults:
    """
    Validate foreign key relationships across tables.

    Args:
        records: Dict with joined_ships, ships, accessions, captains,
                 starship_features, gff, genomes, taxonomy
        blocking_issues: Set of blocking issue types

    Returns:
        ValidationResults
    """
    blocking = blocking_issues or {
        'orphaned_foreign_keys',
        'duplicate_accessions',
        'ship_id_mismatches',
    }
    results = ValidationResults()

    joined_ships = records.get('joined_ships', [])
    ships = {_get_val(s, 'id') or _get_val(s, 'ship_id'): s for s in records.get('ships', []) if (_get_val(s, 'id') or _get_val(s, 'ship_id')) is not None}
    accessions = {_get_val(a, 'id') or _get_val(a, 'accession_id'): a for a in records.get('accessions', []) if (_get_val(a, 'id') or _get_val(a, 'accession_id')) is not None}
    captains = {_get_val(c, 'id') or _get_val(c, 'captain_id'): c for c in records.get('captains', []) if (_get_val(c, 'id') or _get_val(c, 'captain_id')) is not None}
    genomes = {_get_val(g, 'id') or _get_val(g, 'genome_id'): g for g in records.get('genomes', []) if (_get_val(g, 'id') or _get_val(g, 'genome_id')) is not None}
    taxonomy = {_get_val(t, 'id') or _get_val(t, 'taxonomy_id'): t for t in records.get('taxonomy', []) if (_get_val(t, 'id') or _get_val(t, 'taxonomy_id')) is not None}

    # Check joined_ships -> ship_id
    broken_ship_links = 0
    for js in joined_ships:
        ship_id = _get_val(js, 'ship_id')
        if ship_id is not None and ship_id not in ships:
            broken_ship_links += 1
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_foreign_keys',
                    category='blocking',
                    table_name='joined_ships',
                    record_id=_get_val(js, 'id'),
                    starship_id=_get_val(js, 'starshipID'),
                    details=f'joined_ships.ship_id={ship_id} references non-existent ship',
                ),
                blocking,
            )

    # Check joined_ships -> captain_id
    for js in joined_ships:
        captain_id = _get_val(js, 'captain_id')
        if captain_id is not None and captain_id not in captains:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_foreign_keys',
                    category='blocking',
                    table_name='joined_ships',
                    record_id=_get_val(js, 'id'),
                    starship_id=_get_val(js, 'starshipID'),
                    details=f'joined_ships.captain_id={captain_id} references non-existent captain',
                ),
                blocking,
            )

    # Check joined_ships -> genome_id
    for js in joined_ships:
        genome_id = _get_val(js, 'genome_id')
        if genome_id is not None and genome_id not in genomes:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_foreign_keys',
                    category='blocking',
                    table_name='joined_ships',
                    record_id=_get_val(js, 'id'),
                    starship_id=_get_val(js, 'starshipID'),
                    details=f'joined_ships.genome_id={genome_id} references non-existent genome',
                ),
                blocking,
            )

    # Check genomes -> taxonomy_id
    for g in records.get('genomes', []):
        tax_id = _get_val(g, 'taxonomy_id') or _get_val(g, 'tax_id')
        if tax_id is not None and tax_id not in taxonomy:
            results.add_issue(
                ValidationIssue(
                    issue_type='broken_genome_taxonomy',
                    category='blocking',
                    table_name='genomes',
                    record_id=_get_val(g, 'id'),
                    details=f'genome.taxonomy_id={tax_id} references non-existent taxonomy',
                ),
                blocking,
            )

    # Orphaned starship_features (ship_id not in ships)
    for f in records.get('starship_features', []):
        ship_id = _get_val(f, 'ship_id')
        if ship_id is not None and ship_id not in ships:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_features',
                    category='warning',
                    table_name='starship_features',
                    record_id=_get_val(f, 'id'),
                    details=f'starship_features.ship_id={ship_id} references non-existent ship',
                ),
                blocking,
            )

    # Orphaned GFF (ship_id not in ships)
    for g in records.get('gff', []):
        ship_id = _get_val(g, 'ship_id')
        if ship_id is not None and ship_id not in ships:
            results.add_issue(
                ValidationIssue(
                    issue_type='orphaned_gff',
                    category='warning',
                    table_name='gff',
                    record_id=_get_val(g, 'id'),
                    details=f'gff.ship_id={ship_id} references non-existent ship',
                ),
                blocking,
            )

    # Duplicate accessions
    accession_tags = {}
    for a in records.get('accessions', []):
        tag = _get_val(a, 'accession_tag')
        if tag:
            accession_tags[tag] = accession_tags.get(tag, 0) + 1
    for tag, count in accession_tags.items():
        if count > 1:
            results.add_issue(
                ValidationIssue(
                    issue_type='duplicate_accessions',
                    category='blocking',
                    table_name='accessions',
                    details=f'Duplicate accession_tag: {tag} ({count} times)',
                ),
                blocking,
            )

    return results
