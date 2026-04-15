"""
Publish validated Starbase data to a versioned SQLite file.

Default strategy: copy the validated source database (byte-for-byte schema + data),
then append `database_version` metadata. This avoids re-inserting every row in Python.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Dict, Optional

from .version_manager import VersionInfo, insert_version_row
import sqlite3


def load_to_sqlite(
    source_db_path: str,
    output_path: str,
    version_info: VersionInfo,
    extra_metadata: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Copy source SQLite to output_path and stamp version metadata.

    Args:
        source_db_path: Validated starbase.sqlite (or temp copy)
        output_path: Destination path (e.g. starbase_v1.2.0.sqlite)
        version_info: Version and validation summary
        extra_metadata: Optional JSON-serializable dict stored in extra_json

    Returns:
        Version string written to database_version
    """
    src = Path(source_db_path)
    dst = Path(output_path)
    if not src.is_file():
        raise FileNotFoundError(f"Source database not found: {source_db_path}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    conn = sqlite3.connect(str(dst))
    try:
        insert_version_row(conn, version_info, extra_metadata)
        conn.commit()
    finally:
        conn.close()

    return version_info.version


def write_records_to_sqlite(
    _records: Dict[str, Any],
    _output_path: str,
    _version_info: VersionInfo,
) -> str:
    """
    Future: materialize SQLite from in-memory record dicts (full ETL).

    Not implemented — use file copy + stamp for Phase 4.
    """
    raise NotImplementedError(
        "write_records_to_sqlite is reserved for a future row-by-row loader; "
        "use load_to_sqlite(source_db_path, ...) with a validated SQLite file."
    )
