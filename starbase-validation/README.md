# Starbase Validation

Standalone validation/ETL app for Starbase database quality. Sits between Django curation (MySQL) and the Starbase web app (SQLite), ensuring data quality before publication.

## Architecture

```
Django (MySQL) → Validation App (Quality Gate) → Versioned SQLite
                        ↓ (if issues found)
                   Report to Curators
```

## Installation

```bash
cd starbase-validation
pip install -e .
```

## Usage

### Validate SQLite database

```bash
# Validate starbase.sqlite (dry-run, default)
starbase-validate /path/to/starbase.sqlite

# Or as Python module
python -m starbase_validation.pipeline /path/to/starbase.sqlite

# With custom config
starbase-validate starbase.sqlite --config quality_rules.yaml
```

### Programmatic use

```python
from starbase_validation.extractors import extract_from_sqlite
from starbase_validation.validators import validate_all
from starbase_validation.reporters import generate_report
import yaml

# Load data
records = extract_from_sqlite("starbase.sqlite")
with open("config/quality_rules.yaml") as f:
    rules = yaml.safe_load(f)

# Validate
results = validate_all(records, rules)
print(generate_report(results))

if results.has_blocking_issues():
    print(f"Blocking: {len(results.blocking_issues)}")
```

## Structure

```
starbase-validation/
├── starbase_validation/
│   ├── config/
│   │   ├── quality_rules.yaml      # Quality gates
│   │   └── enrichment_sources.yaml
│   ├── validators/                 # Database-agnostic validators
│   │   ├── taxonomy_validator.py
│   │   ├── genome_validator.py
│   │   ├── classification_validator.py
│   │   ├── sequence_validator.py
│   │   └── relationship_validator.py
│   ├── extractors/
│   │   └── sqlite_extractor.py
│   ├── transformers/
│   ├── loaders/
│   ├── reporters/
│   └── pipeline.py
├── cleanup/                        # Legacy cleanup scripts
│   └── utils/
│       ├── quality_tags.py        # Source of validation logic
│       └── ...
└── setup.py
```

## Validators (extracted from cleanup)

| Validator | Source | Checks |
|-----------|--------|--------|
| taxonomy_validator | quality_tags.check_taxonomic_quality | incomplete_taxonomy, missing_strain, taxonomy_uncertain |
| genome_validator | quality_tags.check_genomic_context | missing_assembly, missing_biosample, broken_genome_taxonomy |
| classification_validator | quality_tags.check_classification_hierarchy | missing_family (BLOCKING), missing_navis, missing_haplotype |
| sequence_validator | quality_tags.check_sequence_quality, check_boundary_features | invalid_sequence, short/long, duplicates, boundaries |
| relationship_validator | database_cleanup | orphaned_foreign_keys, duplicate_accessions |

## Quality Rules

Edit `config/quality_rules.yaml` to adjust:
- **blocking_issues**: Issue types that prevent publication
- **warning_issues**: Reported but don't block
- **quality_gates**: Thresholds per category

## Next Steps (from plan)

- Phase 4: Implement ETL pipeline with load step
- Phase 5: Django integration (export_to_starbase command)
- Phase 6: Simplify cleanup utilities to diagnostics only
