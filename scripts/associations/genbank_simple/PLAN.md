# GenBank Simple Workflow Plan

## Purpose

Build a simpler GenBank pathogen-country metadata workflow that only queries the
current WHO zoonotic analysis universe and avoids broad pathogen searches that do
not match the refined network targets.

This workflow should replace the current broad/adaptive GenBank behavior for the
next clean run, while leaving `scripts/associations/genbank/` available as a
reference archive.

## Core Scope

The simple workflow should only consider pathogens represented in both the
current WHO disease/pathogen list and the canonical zoonotic network:

- `pathogen_association_data/WHO/who_diseases/who_pathogens_diseases_zoonotic.csv`
- `pathogen_association_data/WHO/networks/combined_who_network_canonical_zoonotic.csv`

Within that universe, keep only rows where point-data support is present in at
least one source:

```r
in_gibb_etal == TRUE | in_empres_i == TRUE
```

These columns indicate whether we already have point occurrence data for the
disease in the Gibb et al. or EMPRES-i datasets. GenBank should complement that
point-data-backed disease set rather than expanding to every possible WHO
pathogen.

## Explicit Exclusions

Exclude coronavirus targets for now because the current WHO/network labels are
too broad for meaningful GenBank retrieval and need a separate refined
species/strain decision first.

Exclude rows matching any of the following:

- `Subgenus Sarbecovirus`
- `subgenus Merbecovirus`
- SARS
- MERS
- SARS-like
- MERS-like
- COVID
- SARS-CoV
- SARS-CoV-2
- MERS-CoV

This means no SARS/COVID special disease-level query profile and no broad
Sarbecovirus/Merbecovirus retrieval in the simple workflow.

## Influenza Policy

Do not query broad influenza records.

Keep only the influenza analysis units that are present in the refined WHO
network:

- `Alphainfluenzavirus influenzae (H5N1)`
- `Alphainfluenzavirus influenzae (H7N9)`

Do not run broad fallback queries such as:

- `Alphainfluenzavirus influenzae`
- `Influenza A virus`
- `H5Nx`
- `H7Nx`
- other subtype expansions not present in the filtered network

## Retrieval Policy

For included pathogens, retrieve all matched records through deterministic
pagination and checkpointing.

Do not use:

- random sampling
- interval sampling
- adaptive plateau stopping
- first-hit-only country summaries

The workflow can still use batching, rate limiting, retries, and per-pathogen
checkpoints so long as the intended final state is full retrieval for every
included query.

## Proposed Script Layout

Keep this folder small and sequential:

1. `01_build_manifest.R`
   - read the two allowed WHO/network input files
   - normalize pathogen and disease columns
   - keep only rows with `in_gibb_etal == TRUE | in_empres_i == TRUE`
   - exclude broad coronavirus targets
   - exclude broad influenza targets
   - write a reviewable manifest

2. `02_run_genbank_full_retrieval.R`
   - read the approved manifest
   - query NCBI with deterministic pagination
   - write per-pathogen checkpoint files
   - avoid sampling and adaptive stopping

3. `03_summarize_country_metadata.R`
   - bind checkpoint outputs
   - parse and standardize country fields
   - summarize pathogen-country and disease-country coverage
   - preserve links back to canonical WHO network labels

4. `04_quality_checks.R`
   - report excluded targets
   - report pathogens with zero records
   - report pathogens with records but no usable country metadata
   - flag unusually broad queries before interpretation

## Suggested Outputs

Suggested outputs should live under:

```text
pathogen_association_data/WHO/genbank_simple/
```

Proposed files:

- `genbank_simple_manifest.csv`
- `excluded_targets.csv`
- `pathogen_runs/search_logs/*.csv`
- `pathogen_runs/country_records/*.csv`
- `genbank_country_records.csv`
- `genbank_pathogen_country_summary.csv`
- `genbank_disease_country_summary.csv`
- `genbank_simple_qa_summary.csv`

## Design Notes

The old GenBank scripts contain useful parsing and checkpointing ideas, but they
also contain logic that is no longer aligned with the current goal: broad
influenza expansion, SARS/COVID disease profiles, adaptive stopping, interval
sampling, and diagnostic rerun branches.

For this reason, the simple workflow should reuse only small helper patterns
where they remain useful, rather than continuing to add exceptions to the older
scripts.
