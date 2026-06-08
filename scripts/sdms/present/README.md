# Present-Day SDM Fitting

This folder contains the present-day SDM workflow for regenerating calibration
models and building vector SDMs. The current bulk workflow is disease-aware at
the manifest stage, but fits each vector species once so the same SDM can be
reused across diseases.

## Folder Layout

```text
scripts/sdms/present/
  utils.R
  manifests/
    01_build_vector_sdm_target_manifests.R
  occurrences/
    01_prepare_one_gbif_species.R
    01_prepare_gbif_occurrences.R
    02_extract_local_vector_occurrences.R
    03_prepare_species_occurrences_batch.R
    03_prepare_chikungunya_occurrences_batch.R
    04_submit_gbif_download_requests.R
    05_fetch_gbif_download_requests.R
    06_combine_vector_occurrences.R
    07_plot_combined_vector_occurrence_maps.R
    08_plot_occurrence_period_maps.R
    09_audit_gbif_synonyms.R
  models/
    01_run_present_model.R
    02a_run_species_models_batch.R
    02b_run_manifest_model_batch.R
  calibration/
    01_prepare_host_regeneration_manifest.R
    02_compare_occurrences_to_existing_model.R
    03_compare_rousettus_models.R
```

`manifests/` builds SDM target manifests. `occurrences/` downloads, imports,
extracts, and combines occurrence records. `models/` runs or dry-runs present-day
AutoMaxent models. `calibration/` contains diagnostics used while comparing our
regenerated host models against Gonzalo's saved host SDMs.

## Vector SDM Push

The all-disease vector workflow writes to:

```text
sdms/runs/vector_sdm_push/
```

Build target manifests first:

```r
manifest_config <- list(
  recommended_next_action = "find_or_build_vector_sdm",
  roles = "vector",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  overwrite = TRUE
)

source("scripts/sdms/present/manifests/01_build_vector_sdm_target_manifests.R")
```

This creates:

```text
sdms/runs/vector_sdm_push/disease_vector_sdm_targets.csv
sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv
sdms/runs/vector_sdm_push/manifest_build_summary.csv
```

`disease_vector_sdm_targets.csv` keeps disease-vector context. Use it later for
disease proxy stacking. `vector_species_sdm_targets.csv` is the operational
species manifest used for occurrence downloads and model fitting.

## Occurrence Records

For one species, use `occurrences/01_prepare_one_gbif_species.R`. It supports:

- `direct-gbif`: month-by-month `rgbif::occ_search()` download for smaller jobs;
- `spatial-spp`: SDM_Pipeline synonym-expanded download, requiring
  `IUCN_REDLIST_KEY` or `IUCN_API_KEY`;
- `gbif-download`: asynchronous GBIF download API for record-rich species.

For vector batches, use `gbif-download` from 1970 through the current year.
GBIF credentials are read from environment variables named `GBIF_USER`,
`GBIF_PASSWORD`, and `GBIF_EMAIL`, or from repo-ignored `.env` entries named
`gbif_username`, `gbif_password`, and `gbif_email`.

Raw asynchronous GBIF downloads are not filtered by `occurrenceStatus`, so raw
GBIF source files can retain explicit `ABSENT` rows for audit. GBIF inputs used
for modelling are filtered later: the preparation/combination scripts remove
explicit non-present rows plus rows with `individualCount == 0` before
deduplication and cleaning.

Extract local VectorMap and MapVEu records before combining sources:

```r
source("scripts/sdms/present/occurrences/02_extract_local_vector_occurrences.R")
```

The extractor writes exact species matches only:

```text
sdms/runs/vector_sdm_push/occurrences/<Species_safe>/vectormap/raw/
sdms/runs/vector_sdm_push/occurrences/<Species_safe>/mapveu/raw/
sdms/runs/vector_sdm_push/local_vector_occurrence_sources_manifest.csv
```

Current bulk vector occurrence files may be stored outside the repo on the
LaCie drive:

```text
/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/occurrences/
```

To audit whether GBIF downloads likely cover older names and synonyms, run the
diagnostic synonym audit. It reads existing GBIF files and optional GBIF
taxonomy/count metadata, but it does not submit downloads or change occurrence
inputs:

```r
batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  occurrence_root = "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/occurrences",
  request_manifest_path = "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/gbif_download_requests.csv",
  audit_run_root = "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/gbif_synonym_audit_runs",
  roles = "vector",
  query_gbif_api = FALSE,
  query_gbif_counts = FALSE,
  gbif_api_timeout_seconds = 20,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/09_audit_gbif_synonyms.R")
```

## Two-Phase GBIF Downloads

Submit GBIF download requests first:

```r
batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  request_manifest_path = "sdms/runs/vector_sdm_push/gbif_download_requests.csv",
  occurrence_root = "sdms/runs/vector_sdm_push/occurrences",
  request_run_root = "sdms/runs/vector_sdm_push/gbif_download_request_runs",
  roles = "vector",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  max_new_submissions = 3,
  refresh_existing_status = TRUE,
  resubmit_existing = FALSE,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/04_submit_gbif_download_requests.R")
```

The submit script refreshes saved GBIF statuses before submitting and only fills
available GBIF download slots. Re-run it after earlier downloads finish.

Fetch and clean completed downloads later:

```r
batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  request_manifest_path = "sdms/runs/vector_sdm_push/gbif_download_requests.csv",
  occurrence_root = "sdms/runs/vector_sdm_push/occurrences",
  fetch_run_root = "sdms/runs/vector_sdm_push/gbif_download_fetch_runs",
  roles = "vector",
  fetch_statuses = "SUCCEEDED",
  redownload_occurrences = FALSE,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/05_fetch_gbif_download_requests.R")
```

## Combined Vector Occurrences

After GBIF records have been fetched and local records extracted, combine all
available sources into a `combined` occurrence method:

```r
batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  local_source_manifest_path = "sdms/runs/vector_sdm_push/local_vector_occurrence_sources_manifest.csv",
  occurrence_root = "sdms/runs/vector_sdm_push/occurrences",
  combined_run_root = "sdms/runs/vector_sdm_push/combined_vector_occurrence_runs",
  roles = "vector",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  coordinate_round_digits = 5,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/06_combine_vector_occurrences.R")
```

The preferred modelling input is `combined`. If local records are absent, this
can still be GBIF-only after combination. Local-only modelling should be treated
as a fallback when GBIF has too few usable records and local sources have enough
clean records. The combined run summary records how many GBIF rows were removed
by the presence filter before source merging.

To map the model-facing occurrence inputs by source provenance, plot the cleaned
combined records:

```r
batch_config <- list(
  roles = "vector",
  occurrence_method = "combined",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/07_plot_combined_vector_occurrence_maps.R")
```

These maps use the deduplicated/cleaned combined layer, so they show the points
the model workflow will see rather than raw pre-combination source records.
Shared coordinate groups are labelled as mixed-source points.

Map scripts write to stable, human-readable run folders by default so rerunning
the same diagnostic replaces the previous map set. Chikungunya map diagnostics
live under `sdms/runs/chikungunya/maps/`, not under `calibration/`. Set
`timestamped_run_dir = TRUE` only when you intentionally want to keep every
exploratory map run.

## Model Runs

Use `models/01_run_present_model.R` for one species,
`models/02a_run_species_models_batch.R` as the RStudio-friendly batch config,
and `models/02b_run_manifest_model_batch.R` as the generic batch implementation.

The batch script does not prepare occurrences. It expects cleaned occurrence
files under the configured occurrence method folder. By default it is a
status/preflight run only.

For vector reruns intended to align with Gonzalo's saved host SDM catalogue,
use the host-catalog-style settings: `start_year = 2000`, `end_year = 2026`,
`n_background = "dynamic"`, `beta_values = "4,8,12"`, `random_features = TRUE`,
`n_models = 25`, `n_selected_models = 10`, and `use_boyce = 0.5`. The default
`test_percent = "dynamic"` uses `20%` test data below 60 model records and `30%`
otherwise. Set `n_background` to an integer, such as `8000`, for a
fixed-background sensitivity.

```r
batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  occurrence_root = "sdms/runs/vector_sdm_push/occurrences",
  model_output_root = "sdms/runs/vector_sdm_push/models",
  model_batch_run_root = "sdms/runs/vector_sdm_push/model_batch_runs",
  roles = "vector",
  occurrence_method = "combined",
  fit_models = FALSE,
  dry_run_models = TRUE,
  start_year = 2000,
  end_year = 2026
)

source("scripts/sdms/present/models/02a_run_species_models_batch.R")
```

Model outputs from the bulk vector push stay under
`sdms/runs/vector_sdm_push/models/` until they have been reviewed. Do not promote
them into the shared `sdms/models/` catalogue during the first batch run.

## Chikungunya Compatibility

Existing Chikungunya outputs under `sdms/runs/chikungunya/` are not moved or
deleted. The Chikungunya occurrence/model batch defaults now read and write
occurrences from the shared vector workspace:

```text
sdms/runs/vector_sdm_push/occurrences/
```

This avoids duplicating Chikungunya vector records in
`sdms/runs/chikungunya/calibration/occurrences/`. The old model batch filename
has been replaced by the generic `02a`/`02b` pair:

```text
occurrences/03_prepare_chikungunya_occurrences_batch.R
models/02a_run_species_models_batch.R
models/02b_run_manifest_model_batch.R
```

Those scripts are path-configurable. Their run summaries and regenerated model
outputs still default to the existing Chikungunya calibration workspace, but
their occurrence root defaults to the shared vector-push occurrence folder.

The older `occurrences/01_prepare_gbif_occurrences.R` filename is kept as a
compatibility wrapper around `01_prepare_one_gbif_species.R`.
