"""Human-readable validation reports."""

from ..validators.results import ValidationResults


def generate_report(results: ValidationResults) -> str:
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
        lines.append("BLOCKING ISSUES:")
        lines.append("-" * 40)
        for issue in results.blocking_issues[:25]:
            lines.append(
                f"  [{issue.issue_type}] {issue.table_name} id={issue.record_id}: {issue.details}"
            )
        if len(results.blocking_issues) > 25:
            lines.append(f"  ... and {len(results.blocking_issues) - 25} more")
        lines.append("")
    if results.summary:
        lines.append("ISSUE COUNTS (top):")
        for k, v in sorted(results.summary.items(), key=lambda x: -x[1])[:20]:
            lines.append(f"  {k}: {v}")
    lines.append("=" * 60)
    return "\n".join(lines)
