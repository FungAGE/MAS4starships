"""
SQLite extractor - extract records from SQLite starbase database.

Returns record dicts suitable for validators (database-agnostic).
"""

import sqlite3
from typing import Dict, List, Any, Optional


def _row_to_dict(cursor: sqlite3.Cursor, row: tuple) -> Dict[str, Any]:
    """Convert row to dict using column names."""
    return dict(zip([c[0] for c in cursor.description], row))


def extract_from_sqlite(db_path: str) -> Dict[str, List[Dict]]:
    """
    Extract all relevant tables from SQLite starbase database.

    Args:
        db_path: Path to starbase.sqlite file

    Returns:
        Dict with keys: joined_ships, ships, accessions, taxonomy, genomes,
        captains, family_names, navis, haplotype, starship_features, gff
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    tables = [
        'joined_ships', 'ships', 'accessions', 'taxonomy', 'genomes',
        'captains', 'family_names', 'navis_names', 'haplotype_names',
        'starship_features', 'gff'
    ]

    # Map table names to record keys (some have different names)
    key_map = {
        'navis_names': 'navis',
        'haplotype_names': 'haplotype',
    }

    records = {}
    for table in tables:
        try:
            cursor.execute(f"SELECT * FROM {table}")
            rows = cursor.fetchall()
            columns = [desc[0] for desc in cursor.description]
            key = key_map.get(table, table)
            records[key] = [dict(zip(columns, row)) for row in rows]
        except sqlite3.OperationalError:
            records[key_map.get(table, table)] = []

    conn.close()
    return records
