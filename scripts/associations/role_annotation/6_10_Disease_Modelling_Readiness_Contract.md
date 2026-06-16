# Disease Modelling Readiness Contract

This document records the current input/output contract for
`scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R`.
It is a refactor guardrail: C11 helper extraction must preserve the schemas,
row counts, output paths, and evidence-interpretation boundaries described here.

## Scope

The readiness build is a planning and collaborator handoff surface. It combines
WHO analysis-unit rules with role-annotation QA, species roster, direct vector
evidence, country-evidence summaries, and SDM availability inventories.

It must not:

- assign new host/vector biological roles.
- infer vector evidence from host-vector observations alone.
- merge country-evidence layers into a single truth set.
- change upstream inclusion, review, or modelling-scope decisions.

## Required Inputs

The build stops if either required table is missing:

- `who_master_plus_analysis_units_path()`
  - current target:
    `pathogen_association_data/evidence/who_diseases/master_expansion/master_plus_who_analysis_units.csv`
- `who_diseases_transmission_rules_path("master_plus_who_transmission_rules_manual_reviewed_v2.csv")`
  - current target:
    `pathogen_association_data/manual/who_diseases/transmission_rules/master_plus_who_transmission_rules_manual_reviewed_v2.csv`

## Optional Inputs

Missing optional inputs are treated as empty evidence layers, not as failures:

- `who_master_disease_analysis_units_path()`
- `pathogen_association_data/evidence/role_annotation/qa/disease_evidence_readiness.csv`
- `pathogen_association_data/evidence/role_annotation/qa/vector_evidence_readiness_by_disease.csv`
- `pathogen_association_data/evidence/role_annotation/species_host_vector_roster.csv`
- `pathogen_association_data/evidence/who_don_v2/final/who_don_modelling_ready.csv`
  with legacy fallback under `pathogen_association_data/WHO/disease_outbreak_news_v2/final/`
- `pathogen_association_data/evidence/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv`
  with fallback to the standard GenBank summary.
- `sdms/outputs/catalog/accessible_sdm_species.csv`
- `sdms/outputs/projections/projection_manifest.csv`
- `sdms/outputs/comparisons/comparison_manifest.csv`
- `pathogen_association_data/staged/virion/outputs/who_pathogens_virion_taxid.csv`
- `pathogen_association_data/staged/clover/outputs/who_bacteria_clover_taxid.csv`

## Outputs

The script writes these tracked readiness surfaces under
`pathogen_association_data/readiness/`:

- `disease_modelling_readiness.csv`
- `disease_modelling_pilot.csv`
- `disease_modelling_readiness_full.csv`
- `disease_modelling_pilot_package/manifest.csv`
- `disease_modelling_pilot_package/disease_modelling_pilot.csv`
- `disease_modelling_pilot_package/pilot_hosts.csv`
- `disease_modelling_pilot_package/pilot_vectors.csv`
- `disease_modelling_pilot_package/pilot_countries.csv`
- `disease_modelling_pilot_package/pilot_sdm_species.csv`
- `disease_modelling_pilot_package/pilot_evidence_summary.csv`
- `disease_modelling_pilot_package/README.md`
- `disease_modelling_pilot_package.rds`
- `disease_modelling_pilot_package.xlsx`, when `writexl` is available.

## Current Invariant Baseline

Current schema and row-count baseline after the C11 helper split and the
master-plus compatibility migration for role-roster inputs:

| Output | Rows | Columns |
| --- | ---: | ---: |
| `disease_modelling_readiness.csv` | 88 | 30 |
| `disease_modelling_readiness_full.csv` | 90 | 135 |
| `disease_modelling_pilot.csv` | 31 | 34 |
| `disease_modelling_pilot_package/manifest.csv` | 6 | 6 |
| `disease_modelling_pilot_package/pilot_hosts.csv` | 2235 | 18 |
| `disease_modelling_pilot_package/pilot_vectors.csv` | 448 | 26 |
| `disease_modelling_pilot_package/pilot_countries.csv` | 1188 | 16 |
| `disease_modelling_pilot_package/pilot_sdm_species.csv` | 2683 | 17 |
| `disease_modelling_pilot_package/pilot_evidence_summary.csv` | 31 | 73 |

Host-role fields are part of the current baseline for `pilot_hosts.csv` and
`pilot_sdm_species.csv`:

- `host_role_assignment`
- `host_role_confidence`
- `host_role_needs_manual_review`
- `host_role_assignment_status`

Vector-role fields are part of the current baseline for `pilot_vectors.csv` and
`pilot_sdm_species.csv`:

- `vector_role_assignment`
- `vector_role_confidence`
- `vector_role_needs_manual_review`
- `vector_role_assignment_status`

## Package Manifest

The package manifest must contain exactly six rows, one for each package data
table:

- `disease_modelling_pilot`
- `pilot_hosts`
- `pilot_vectors`
- `pilot_countries`
- `pilot_sdm_species`
- `pilot_evidence_summary`

Manifest columns must remain:

- `generated_at_utc`
- `table_name`
- `file_name`
- `rows`
- `columns`
- `source_description`

`generated_at_utc` is intentionally volatile on each rerun. Refactors may
change only that timestamp unless upstream inputs have changed.

## Validation Commands

Run from the repository root:

```sh
Rscript scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R
```

Parse-check affected R files:

```sh
Rscript -e 'files <- c("scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R", "scripts/associations/role_annotation/disease_modelling_readiness_helpers.R"); bad <- files[!vapply(files, function(f) is.expression(parse(f)), logical(1))]; print(bad)'
```

Review whitespace and patch hygiene:

```sh
git --no-pager diff --check
```

After rerun, compare row counts, column names, column order, and deterministic
CSV content against the pre-refactor baseline. Ignore only `generated_at_utc` in
`disease_modelling_pilot_package/manifest.csv`.
