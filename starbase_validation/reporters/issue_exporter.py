"""
Export validation issues to CSV and persist to Django (Phase 5).

Used by ``export_to_starbase`` when validation fails or for audit trails.
"""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING, Any, List, Optional, Sequence, Union

if TYPE_CHECKING:
    from starbase_validation.validators.results import ValidationIssue, ValidationResults


def _format_details(issue: "ValidationIssue") -> str:
    parts = []
    if getattr(issue, "starship_id", None):
        parts.append(f"starshipID={issue.starship_id}")
    if issue.details:
        parts.append(issue.details)
    return "; ".join(parts) if parts else ""


def validation_results_to_rows(results: "ValidationResults") -> List[dict]:
    """Flatten ValidationResults into row dicts for CSV or inspection."""
    rows: List[dict] = []
    for issue in results.issues:
        rows.append(
            {
                "issue_type": issue.issue_type,
                "category": issue.category,
                "severity": "blocking" if issue.category == "blocking" else "warning",
                "table_name": issue.table_name,
                "record_id": issue.record_id,
                "starship_id": getattr(issue, "starship_id", None) or "",
                "details": _format_details(issue),
            }
        )
    return rows


def write_issues_csv(
    path: str,
    results: "ValidationResults",
    issues: Optional[Sequence["ValidationIssue"]] = None,
) -> int:
    """
    Write issues to CSV. If ``issues`` is None, writes all ``results.issues``.

    Returns:
        Number of rows written.
    """
    rows = issues if issues is not None else results.issues
    fieldnames = [
        "issue_type",
        "category",
        "severity",
        "table_name",
        "record_id",
        "starship_id",
        "details",
    ]
    count = 0
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for issue in rows:
            w.writerow(
                {
                    "issue_type": issue.issue_type,
                    "category": issue.category,
                    "severity": "blocking" if issue.category == "blocking" else "warning",
                    "table_name": issue.table_name,
                    "record_id": issue.record_id,
                    "starship_id": getattr(issue, "starship_id", None) or "",
                    "details": _format_details(issue),
                }
            )
            count += 1
    return count


def persist_issues_to_django_run(
    run: Any,
    results: "ValidationResults",
    *,
    max_issues: int = 2000,
    prefer_blocking: bool = True,
) -> int:
    """
    Bulk-create :class:`starship.models.ValidationIssue` rows for a ValidationRun.

    Caps total rows at ``max_issues``. If ``prefer_blocking``, store blocking issues
    first, then warnings until the cap.

    Returns:
        Number of rows created.

    Requires Django; import ``ValidationIssue`` inside the function to avoid
    importing Django when running the CLI outside Django.
    """
    from starship.models import ValidationIssue as DjangoValidationIssue

    blocking = list(results.blocking_issues)
    warnings = list(results.warning_issues)
    to_store: List = []
    if prefer_blocking:
        to_store.extend(blocking)
        remaining = max_issues - len(to_store)
        if remaining > 0:
            to_store.extend(warnings[:remaining])
    else:
        combined = blocking + warnings
        to_store = combined[:max_issues]

    objs = []
    for issue in to_store:
        sev = "blocking" if issue.category == "blocking" else "warning"
        objs.append(
            DjangoValidationIssue(
                run=run,
                issue_type=issue.issue_type[:200],
                table_name=(issue.table_name or "")[:200],
                record_id=issue.record_id,
                severity=sev,
                details=_format_details(issue)[:10000],
            )
        )
    DjangoValidationIssue.objects.bulk_create(objs, batch_size=500)
    return len(objs)
