# GenBank Simple Workflow

This folder contains the active GenBank pathogen-country enrichment workflow for
WHO disease modelling readiness. It replaces the older broad/adaptive GenBank
workflow with a manifest-driven process that keeps retrieval targets reviewable
and writes standardized disease-country evidence for downstream readiness
tables.

## Target Surface

The active target surface is `genbank_simple_readiness_manifest.csv`, the
expanded readiness manifest built from
`pathogen_association_data/readiness/disease_modelling_readiness.csv` and written
under `pathogen_association_data/staged/genbank_simple/manifests/`.

Readiness mode starts from non-held readiness rows,
joins the full readiness audit table for query/provenance fields, and builds one
retrieval target per unique species-level query label.

`genbank_simple_manifest.csv` is the frozen original 19-target point-data-backed
WHO zoonotic manifest. It is retained as legacy provenance/control data, not as
the active pipeline scope. Set `GENBANK_SIMPLE_USE_LEGACY_19_MANIFEST=TRUE`
only when `01b_build_readiness_manifest.R` needs temporary old-manifest
comparison fields.

The old `pathogen_runs/` checkpoint files are also frozen cache evidence. They
are still combined into readiness summaries because rerunning GenBank retrieval
is slow and external-state-sensitive. Do not delete or ignore those checkpoint
files unless a full explicit retrieval refresh is planned.

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
   Legacy builder for the frozen original 19-target GenBank-simple manifest.

2. `01b_build_readiness_manifest.R`
   Builds the active expanded readiness manifest and row-level manifest QA
   table.

3. `02_run_genbank_full_retrieval.R`
   Retrieves NCBI nuccore records from the readiness manifest by default with
   deterministic pagination and per-target checkpoints. Set
   `GENBANK_SIMPLE_MANIFEST_KIND=standard` only for explicit legacy 19-target
   runs.

4. `03_summarize_country_metadata.R`
   Binds checkpoint outputs and writes readiness pathogen-country and
   disease-country summaries by default. It combines frozen 19-target
   checkpoints with readiness checkpoints as cached evidence, while keeping the
   readiness manifest as the target surface.

5. `04_quality_checks.R`
   Writes readiness search-log, target-level, and summary QA tables by default.

6. `05_standardize_countries.R`
   Standardizes readiness country names and writes standardization QA by
   default.

7. `06_map_disease_countries.R`
   Writes readiness map-control CSVs and generated disease-country PNG maps by
   default.

## Commit Policy

Commit lightweight, reviewable readiness artifacts:

- generated manifests under `pathogen_association_data/staged/genbank_simple/manifests/`;
- manual query overrides under `pathogen_association_data/manual/genbank_simple/`;
- standardized readiness disease-country summaries under
  `pathogen_association_data/evidence/genbank_simple/`;
- QA tables under `pathogen_association_data/evidence/genbank_simple/qa/`;
- aggregate readiness summaries under
  `pathogen_association_data/staged/genbank_simple/intermediate/`;
- compact readiness map-control CSVs under
  `pathogen_association_data/staged/genbank_simple/maps/readiness/`.

Do not commit bulky or transient retrieval outputs:

- `staged/genbank_simple/local_runs/pathogen_runs/`
- `staged/genbank_simple/local_runs/pathogen_runs_readiness/`
- record-level `*_country_records*.csv`
- generated PNG maps under `staged/genbank_simple/maps/*/disease_country_records/`
- standard-mode map-control CSVs under `staged/genbank_simple/maps/standard/`
- local retrieval logs.

## Downstream Handoff

Downstream modelling-readiness scripts should prefer
`pathogen_association_data/evidence/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv`
when it exists.
