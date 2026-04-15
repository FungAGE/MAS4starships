"""
Extract Starbase-shaped tables from MySQL (optional).

Use when curated data lives in MySQL before SQLite publication.
Requires: pymysql
"""

from typing import Any, Dict, List, Optional

# Same logical tables as SQLite extractor
_TABLES = [
    ("joined_ships", "joined_ships"),
    ("ships", "ships"),
    ("accessions", "accessions"),
    ("taxonomy", "taxonomy"),
    ("genomes", "genomes"),
    ("captains", "captains"),
    ("family_names", "family_names"),
    ("navis_names", "navis"),
    ("haplotype_names", "haplotype"),
    ("starship_features", "starship_features"),
    ("gff", "gff"),
]


def extract_from_mysql(mysql_config: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Connect to MySQL and read tables into the same record dict shape as SQLite.

    mysql_config keys: user, password, host, port (optional, default 3306), database

    Raises:
        ImportError: if pymysql is not installed
        Exception: on connection/query errors
    """
    try:
        import pymysql
        from pymysql.cursors import DictCursor
    except ImportError as e:
        raise ImportError(
            "extract_from_mysql requires pymysql. Install with: pip install pymysql"
        ) from e

    port = int(mysql_config.get("port", 3306))
    conn = pymysql.connect(
        host=mysql_config["host"],
        user=mysql_config["user"],
        password=mysql_config["password"],
        database=mysql_config["database"],
        port=port,
        cursorclass=DictCursor,
    )
    records: Dict[str, List[Dict[str, Any]]] = {}
    try:
        with conn.cursor() as cur:
            for sql_name, key in _TABLES:
                try:
                    cur.execute(f"SELECT * FROM `{sql_name}`")
                    rows = cur.fetchall()
                    records[key] = [dict(r) for r in rows]
                except Exception:
                    records[key] = []
    finally:
        conn.close()
    return records


def mysql_config_from_django(settings_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Build mysql_config from Django DATABASES['default'] style dict."""
    return {
        "user": settings_dict.get("USER") or settings_dict.get("user"),
        "password": settings_dict.get("PASSWORD") or settings_dict.get("password"),
        "host": settings_dict.get("HOST") or settings_dict.get("host", "127.0.0.1"),
        "port": settings_dict.get("PORT") or settings_dict.get("port", 3306),
        "database": settings_dict.get("NAME") or settings_dict.get("database"),
    }
