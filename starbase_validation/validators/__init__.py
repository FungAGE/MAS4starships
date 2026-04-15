"""Validators — database-agnostic checks on record dicts."""

from typing import Dict, Optional

from .results import ValidationError, ValidationIssue, ValidationResults
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
    blocking = blocking_issues or set(quality_rules.get("blocking_issues", []))
    results = ValidationResults()
    results.add(validate_taxonomy(records, quality_rules.get("taxonomy", {}), blocking))
    results.add(validate_genome(records, quality_rules.get("genome", {}), blocking))
    results.add(
        validate_classification(records, quality_rules.get("classification", {}), blocking)
    )
    results.add(validate_sequence(records, quality_rules.get("sequences", {}), blocking))
    results.add(validate_relationships(records, blocking))
    return results


__all__ = [
    "ValidationResults",
    "ValidationIssue",
    "ValidationError",
    "validate_all",
    "validate_taxonomy",
    "validate_genome",
    "validate_classification",
    "validate_sequence",
    "validate_relationships",
]
