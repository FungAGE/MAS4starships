"""Validation results aggregation."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


@dataclass
class ValidationIssue:
    issue_type: str
    category: str
    table_name: str
    record_id: Optional[int] = None
    starship_id: Optional[str] = None
    details: str = ""
    severity: str = "warning"


@dataclass
class ValidationResults:
    issues: List[ValidationIssue] = field(default_factory=list)
    blocking_issues: List[ValidationIssue] = field(default_factory=list)
    warning_issues: List[ValidationIssue] = field(default_factory=list)
    summary: Dict[str, int] = field(default_factory=dict)

    def add(self, result: "ValidationResults") -> None:
        self.issues.extend(result.issues)
        self.blocking_issues.extend(result.blocking_issues)
        self.warning_issues.extend(result.warning_issues)
        for key, count in result.summary.items():
            self.summary[key] = self.summary.get(key, 0) + count

    def add_issue(self, issue: ValidationIssue, blocking_issues: Set[str]) -> None:
        self.issues.append(issue)
        if issue.issue_type in blocking_issues:
            issue.category = "blocking"
            self.blocking_issues.append(issue)
        else:
            issue.category = "warning"
            self.warning_issues.append(issue)
        self.summary[issue.issue_type] = self.summary.get(issue.issue_type, 0) + 1

    def has_blocking_issues(self) -> bool:
        return len(self.blocking_issues) > 0

    def is_acceptable(self) -> bool:
        return not self.has_blocking_issues()


class ValidationError(Exception):
    """Raised when validation fails with blocking issues."""

    def __init__(self, message: str, results: Optional[ValidationResults] = None):
        super().__init__(message)
        self.results = results
