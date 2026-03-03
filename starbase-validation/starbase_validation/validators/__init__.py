"""
Reusable validators - database-agnostic quality checks.

Validators work on record dicts, not database sessions.
Extracted from cleanup/utils/quality_tags.py and database_cleanup.py.
"""

from typing import Dict, Optional

from .results import ValidationResults, ValidationIssue
from .taxonomy_validator import validate_taxonomy
from .genome_validator import validate_genome
from .classification_validator import validate_classification
from .sequence_validator import validate_sequence
from .relationship_validator import validate_relationships


def validate_all(
    records: Dict,
    quality_rules: Dict,
    blocking_issues: Optional[set] = None,
) -> ValidationResults:
    """
    Run all validators on extracted records.

    Args:
        records: Dict with keys 'joined_ships', 'ships', 'taxonomy', 'genomes',
                 'starship_features', 'gff', 'accessions', 'captains', etc.
        quality_rules: Config from quality_rules.yaml
        blocking_issues: Set of issue types that block publication

    Returns:
        ValidationResults with all issues found
    """
    blocking = blocking_issues or set(
        quality_rules.get("blocking_issues", [])
    )
    results = ValidationResults()

    # Run validators
    results.add(validate_taxonomy(records, quality_rules.get("taxonomy", {}), blocking))
    results.add(validate_genome(records, quality_rules.get("genome", {}), blocking))
    results.add(validate_classification(records, quality_rules.get("classification", {}), blocking))
    results.add(validate_sequence(records, quality_rules.get("sequences", {}), blocking))
    results.add(validate_relationships(records, blocking))

    return results


__all__ = [
    "ValidationResults",
    "ValidationIssue",
    "validate_all",
    "validate_taxonomy",
    "validate_genome",
    "validate_classification",
    "validate_sequence",
    "validate_relationships",
]
