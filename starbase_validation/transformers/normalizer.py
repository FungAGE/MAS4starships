"""Normalize records (placeholder for enrichment pipeline)."""

from typing import Any, Dict


def normalize_records(records: Dict[str, Any]) -> Dict[str, Any]:
    """Pass-through; extend for column renames and null handling."""
    return records


def enrich_records(records: Dict[str, Any], _enrichment_config: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Enrich metadata from external sources (NCBI, etc.) — Phase 4 stub.

    Returns records unchanged until enrichment is implemented.
    """
    return records
