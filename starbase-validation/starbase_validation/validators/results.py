"""
Validation results container - aggregates issues from all validators.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set


@dataclass
class ValidationIssue:
    """Single validation issue."""
    issue_type: str
    category: str  # 'blocking' or 'warning'
    table_name: str
    record_id: Optional[int] = None
    starship_id: Optional[str] = None
    details: str = ""
    severity: str = "warning"  # 'critical', 'major', 'minor'


@dataclass
class ValidationResults:
    """Aggregated results from all validators."""
    issues: List[ValidationIssue] = field(default_factory=list)
    blocking_issues: List[ValidationIssue] = field(default_factory=list)
    warning_issues: List[ValidationIssue] = field(default_factory=list)
    summary: Dict[str, int] = field(default_factory=dict)

    def add(self, result: "ValidationResults") -> None:
        """Merge another ValidationResults into this one."""
        self.issues.extend(result.issues)
        self.blocking_issues.extend(result.blocking_issues)
        self.warning_issues.extend(result.warning_issues)
        for key, count in result.summary.items():
            self.summary[key] = self.summary.get(key, 0) + count

    def add_issue(self, issue: ValidationIssue, blocking_issues: Set[str]) -> None:
        """Add an issue, categorizing as blocking or warning."""
        self.issues.append(issue)
        if issue.issue_type in blocking_issues:
            issue.category = "blocking"
            self.blocking_issues.append(issue)
        else:
            issue.category = "warning"
            self.warning_issues.append(issue)
        self.summary[issue.issue_type] = self.summary.get(issue.issue_type, 0) + 1

    def has_blocking_issues(self) -> bool:
        """Return True if any blocking issues were found."""
        return len(self.blocking_issues) > 0

    def is_acceptable(self, allow_warnings: bool = True) -> bool:
        """Return True if validation passed (no blocking issues)."""
        return not self.has_blocking_issues()

    def get_issue_types(self) -> Set[str]:
        """Return set of all issue types found."""
        return {i.issue_type for i in self.issues}
