# Validators

Database-agnostic validation logic extracted from `cleanup/utils/`.

## Source Mapping

| Validator | Extracted From | Functions |
|-----------|----------------|-----------|
| taxonomy_validator | quality_tags.py | check_taxonomic_quality() |
| genome_validator | quality_tags.py | check_genomic_context() |
| classification_validator | quality_tags.py | check_classification_hierarchy() |
| sequence_validator | quality_tags.py | check_sequence_quality(), check_boundary_features() |
| relationship_validator | database_cleanup.py | check_*_relationships(), identify_duplicate_* |

## Design Principles

- **Pure functions**: Validators work on record dicts, no database sessions
- **Reusable**: Same logic for SQLite validation and MySQL export validation
- **Configurable**: Blocking vs warning issues from quality_rules.yaml
