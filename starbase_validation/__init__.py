"""
Starbase validation and cleanup utilities for MAS.

This package provides database hygiene, validation, and ETL helpers that operate
on the Starbase SQLite database (and related workflows) without depending on a
separate Starbase application checkout.

Use ``from starbase_validation.pipeline import ValidationPipeline`` to avoid
import cycles when running ``python -m starbase_validation.pipeline``.
"""

__all__ = ["cleanup"]
