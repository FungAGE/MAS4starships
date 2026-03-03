"""Extractors - read data from source databases."""

from .sqlite_extractor import extract_from_sqlite

__all__ = ["extract_from_sqlite"]
