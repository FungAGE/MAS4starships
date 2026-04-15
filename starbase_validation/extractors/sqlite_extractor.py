"""Extract record dicts from Starbase SQLite."""

import sqlite3
from typing import Any, Dict, List


def extract_from_sqlite(db_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Load all starbase tables into lists of row dicts.

    Table names match typical Starbase schema (see starship/starbase_models.py).
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    tables = [
        "joined_ships",
        "ships",
        "accessions",
        "taxonomy",
        "genomes",
        "captains",
        "family_names",
        "navis_names",
        "haplotype_names",
        "starship_features",
        "gff",
    ]
    key_map = {
        "navis_names": "navis",
        "haplotype_names": "haplotype",
    }

    records: Dict[str, List[Dict[str, Any]]] = {}
    for table in tables:
        key = key_map.get(table, table)
        try:
            cursor.execute(f"SELECT * FROM {table}")
            rows = cursor.fetchall()
            columns = [d[0] for d in cursor.description]
            records[key] = [dict(zip(columns, row)) for row in rows]
        except sqlite3.OperationalError:
            records[key] = []

    conn.close()
    return records
