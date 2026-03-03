"""
Validation report generator.
"""

from typing import Dict, Any

from ..validators.results import ValidationResults


def generate_report(results: ValidationResults) -> str:
    """
    Generate human-readable validation report.

    Args:
        results: ValidationResults from validate_all()

    Returns:
        Formatted report string
    """
    lines = [
        "=" * 60,
        "STARBASE VALIDATION REPORT",
        "=" * 60,
        "",
        f"Total issues: {len(results.issues)}",
        f"Blocking issues: {len(results.blocking_issues)}",
        f"Warning issues: {len(results.warning_issues)}",
        "",
    ]

    if results.blocking_issues:
        lines.append("BLOCKING ISSUES (must fix before pushing to starbase):")
        lines.append("-" * 40)
        for issue in results.blocking_issues[:20]:
            lines.append(f"  [{issue.issue_type}] {issue.table_name} id={issue.record_id}: {issue.details}")
        if len(results.blocking_issues) > 20:
            lines.append(f"  ... and {len(results.blocking_issues) - 20} more")
        lines.append("")

    if results.warning_issues:
        lines.append("WARNING ISSUES:")
        lines.append("-" * 40)
        for issue_type, count in sorted(results.summary.items(), key=lambda x: -x[1])[:15]:
            if issue_type not in {i.issue_type for i in results.blocking_issues}:
                lines.append(f"  {issue_type}: {count}")
        lines.append("")

    lines.append("=" * 60)
    return "\n".join(lines)
