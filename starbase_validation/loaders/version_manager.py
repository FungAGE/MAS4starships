"""
Semantic versioning for published Starbase SQLite databases.

Adds `database_version` metadata table to the published file.
"""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


@dataclass
class VersionInfo:
    version: str
    created_at: str
    validation_passed: bool
    total_joined_ships: int
    blocking_count: int
    warning_count: int
    source: str
    notes: str = ""


def parse_semver(s: str) -> Tuple[int, int, int]:
    m = re.match(r"^(\d+)\.(\d+)\.(\d+)$", s.strip())
    if not m:
        return (0, 1, 0)
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def bump_version(
    current: Optional[str],
    bump: str = "minor",
) -> str:
    """
    Increment semantic version.

    bump: 'major' | 'minor' | 'patch'
    """
    if not current:
        return "1.0.0"
    major, minor, patch = parse_semver(current)
    if bump == "major":
        return f"{major + 1}.0.0"
    if bump == "minor":
        return f"{major}.{minor + 1}.0"
    return f"{major}.{minor}.{patch + 1}"


def read_version_from_db(db_path: str) -> Optional[str]:
    """Read version string from database_version table if present."""
    if not Path(db_path).is_file():
        return None
    conn = sqlite3.connect(db_path)
    try:
        row = conn.execute(
            "SELECT version FROM database_version ORDER BY created_at DESC LIMIT 1"
        ).fetchone()
        return row[0] if row else None
    except sqlite3.OperationalError:
        return None
    finally:
        conn.close()


def ensure_version_table(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS database_version (
            version TEXT NOT NULL,
            created_at TEXT NOT NULL,
            validation_passed INTEGER NOT NULL,
            total_records INTEGER,
            blocking_count INTEGER,
            warning_count INTEGER,
            quality_score REAL,
            source TEXT,
            notes TEXT,
            extra_json TEXT
        )
        """
    )


def insert_version_row(conn: sqlite3.Connection, info: VersionInfo, extra: Optional[Dict] = None) -> None:
    ensure_version_table(conn)
    extra_json = json.dumps(extra) if extra else None
    total = info.total_joined_ships
    qs = None
    if total > 0 and info.warning_count >= 0:
        qs = max(0.0, 1.0 - (info.blocking_count + info.warning_count * 0.01) / max(total, 1))

    conn.execute(
        """
        INSERT INTO database_version (
            version, created_at, validation_passed, total_records,
            blocking_count, warning_count, quality_score, source, notes, extra_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            info.version,
            info.created_at,
            1 if info.validation_passed else 0,
            total,
            info.blocking_count,
            info.warning_count,
            qs,
            info.source,
            info.notes,
            extra_json,
        ),
    )


class DatabaseVersionManager:
    """Coordinates version strings and metadata for published databases."""

    def __init__(self, last_version: Optional[str] = None, bump: str = "minor"):
        self.last_version = last_version
        self.bump = bump

    def next_version(self) -> str:
        return bump_version(self.last_version, self.bump)

    def stamp_database(
        self,
        db_path: str,
        info: VersionInfo,
        extra: Optional[Dict] = None,
    ) -> str:
        """Append version row to db_path. Returns version string."""
        conn = sqlite3.connect(db_path)
        try:
            insert_version_row(conn, info, extra)
            conn.commit()
        finally:
            conn.close()
        return info.version
