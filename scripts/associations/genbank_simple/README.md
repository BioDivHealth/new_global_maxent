# GenBank Simple Workflow

This folder contains the active GenBank pathogen-country enrichment workflow for
WHO disease modelling readiness. It replaces the older broad/adaptive GenBank
workflow with a manifest-driven process that keeps retrieval targets reviewable
and writes standardized disease-country evidence for downstream readiness
tables.

## Target Surfaces

The workflow supports two target surfaces:

- `genbank_simple_manifest.csv`: the original 19-target point-data-backed WHO
  zoonotic manifest, retained as a reference/control surface.
- `genbank_simple_readiness_manifest.csv`: the expanded readiness manifest
  built from
  `pathogen_association_data/WHO/role_annotation/qa/disease_modelling_readiness.csv`.

Readiness mode is the current main path. It starts from non-held readiness rows,
joins the full readiness audit table for query/provenance fields, and builds one
retrieval target per unique species-level query label.

## Guardrails

- Coronavirus rows remain deferred unless narrower species/strain retrieval
  targets are reviewed later.
- Broad influenza labels are not queried unless a concrete subtype is available
  in the source label, such as H5N1 or H7N9.
- Salmonella is skipped for readiness GenBank retrieval because the record
  volume is too broad for the current modelling-use case.
- Targets above `GENBANK_SIMPLE_MAX_RECORDS_FOUND` are deferred rather than
  partially downloaded as if complete.

## Scripts

Run scripts from the repository root.

1. `01_build_manifest.R`
   Builds the original 19-target GenBank-simple manifest.

2. `01b_build_readiness_manifest.R`
   Builds the expanded readiness manifest and row-level manifest QA table.

3. `02_run_genbank_full_retrieval.R`
   Retrieves NCBI nuccore records with deterministic pagination and per-target
   checkpoints.

4. `03_summarize_country_metadata.R`
   Binds checkpoint outputs and writes pathogen-country and disease-country
   summaries. With `GENBANK_SIMPLE_SUMMARY_KIND=readiness_combined`, it combines
   the original and readiness runs.

5. `04_quality_checks.R`
   Writes search-log, target-level, and summary QA tables.

6. `05_standardize_countries.R`
   Standardizes country names and writes standardization QA.

7. `06_map_disease_countries.R`
   Writes map-control CSVs and generated disease-country PNG maps.

## Commit Policy

Commit lightweight, reviewable readiness artifacts:

- top-level manifests, query overrides, and standardized disease-country
  summaries;
- `qa/` control and QA tables;
- `intermediate/` aggregate summaries;
- compact `maps_readiness/*.csv` map-control files.

Do not commit bulky or transient retrieval outputs:

- `pathogen_runs/`
- `pathogen_runs_readiness/`
- record-level `*_country_records*.csv`
- generated PNG maps under `maps*/disease_country_records/`
- local retrieval logs.

## Downstream Handoff

Downstream modelling-readiness scripts should prefer
`pathogen_association_data/WHO/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv`
when it exists.
