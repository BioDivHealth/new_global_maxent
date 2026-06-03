# Present-Day SDM Fitting

This folder contains the present-day SDM workflow for Chikungunya host/vector
calibration and vector model generation.

## Folder Layout

```text
scripts/sdms/present/
  utils.R
  occurrences/
    01_prepare_gbif_occurrences.R
    02_extract_local_vector_occurrences.R
    03_prepare_chikungunya_occurrences_batch.R
    04_submit_gbif_download_requests.R
    05_fetch_gbif_download_requests.R
    06_combine_vector_occurrences.R
  models/
    01_run_present_model.R
    02_run_chikungunya_models_batch.R
  calibration/
    01_prepare_host_regeneration_manifest.R
    02_compare_occurrences_to_existing_model.R
    03_compare_rousettus_models.R
```

`occurrences/` scripts create occurrence inputs. `models/` scripts run or
dry-run present-day AutoMaxent models. `calibration/` scripts are comparison and
diagnostic scripts used while checking our regenerated models against Gonzalo's
saved host SDMs.

## Occurrence Preparation

Use `occurrences/01_prepare_gbif_occurrences.R` for one species. It supports
three GBIF pathways:

- `direct-gbif`: month-by-month `rgbif::occ_search()` download for smaller jobs;
- `spatial-spp`: SDM_Pipeline synonym-expanded download, requiring
  `IUCN_REDLIST_KEY` or `IUCN_API_KEY`;
- `gbif-download`: asynchronous GBIF download API for record-rich species.

For Chikungunya vectors, the working default is `gbif-download` with records from
1970 through the current calendar year.

```sh
Rscript scripts/sdms/present/occurrences/01_prepare_gbif_occurrences.R \
  --species "Aedes albopictus" \
  --method gbif-download \
  --manifest sdms/runs/chikungunya/sdm_target_manifest.csv \
  --start-year 1970 \
  --end-year 2026 \
  --redownload
```

The script reads GBIF credentials from environment variables named `GBIF_USER`,
`GBIF_PASSWORD`, and `GBIF_EMAIL`, or from repo-ignored `.env` entries named
`gbif_username`, `gbif_password`, and `gbif_email`.

Outputs are separated by species and method:

```text
sdms/runs/chikungunya/calibration/occurrences/<Species_safe>/<method>/
```

Before cleaning, synonym-expanded downloads are deduplicated by GBIF key where
that key is available. The occurrence summary records raw, year-filtered,
deduplicated, cleaned, and unique-coordinate counts.

## Two-Phase GBIF Download Runs

Use the two-phase workflow for Chikungunya vectors and other record-rich
species. It avoids waiting for every GBIF download inside one long serial run.

Submit requests first:

```r
batch_config <- list(
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

This writes a durable request ledger:

```text
sdms/runs/chikungunya/calibration/gbif_download_requests.csv
```

By default, the submit script first seeds that ledger from any existing
per-species `gbif-download/raw/raw_download_manifest.csv` files, so previously
downloaded species are recorded before new GBIF requests are submitted.
It also defaults to `max_new_submissions = 3`, matching GBIF's simultaneous
download limit for the account. Before submitting, it refreshes saved GBIF
download statuses and only uses the free slots. Re-run the submit script after
earlier jobs have finished to submit the next batch.

Later, after GBIF has finished preparing the downloads, fetch and clean ready
requests:

```r
batch_config <- list(
  roles = "vector",
  fetch_statuses = "SUCCEEDED",
  redownload_occurrences = FALSE,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/05_fetch_gbif_download_requests.R")
```

`05_fetch_gbif_download_requests.R` refreshes GBIF status metadata, skips jobs
that are not ready, and calls `01_prepare_gbif_occurrences.R` with the saved
`gbif_download_key` for ready jobs. The one-species script then imports the ZIP,
deduplicates, cleans, writes occurrence summaries, and updates
`sdm_target_manifest.csv` when requested.

Submit/fetch run summaries are written under:

```text
sdms/runs/chikungunya/calibration/gbif_download_request_runs/
sdms/runs/chikungunya/calibration/gbif_download_fetch_runs/
```

## One-Pass Occurrence Runs

Use `occurrences/03_prepare_chikungunya_occurrences_batch.R` only for small
one-pass jobs where it is acceptable for each species to submit, wait, import,
and clean before the next species starts. For `gbif-download` vector batches,
prefer the two-phase workflow above.

For RStudio use, either edit the top `batch_config` block in the script, or
define `batch_config` in the console immediately before sourcing:

```r
batch_config <- list(
  roles = "vector",
  occurrence_method = "direct-gbif",
  prepare_occurrences = FALSE,
  redownload_occurrences = FALSE,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y"))
)

source("scripts/sdms/present/occurrences/03_prepare_chikungunya_occurrences_batch.R")
```

The one-pass batch script retries occurrence preparation using
`occurrence_download_attempts` and `occurrence_retry_sleep_seconds` from its
internal defaults. Logs and summaries are written under:

```text
sdms/runs/chikungunya/calibration/occurrence_batch_runs/
```

## Local Vector Occurrence Sources

Use `occurrences/02_extract_local_vector_occurrences.R` to copy exact species
matches from local VectorMap and MapVEu raw tables into the same occurrence
workspace. These are source-specific raw folders only; they are not yet the
combined cleaned occurrence input for model fitting.

```text
sdms/runs/chikungunya/calibration/occurrences/<Species_safe>/vectormap/raw/
sdms/runs/chikungunya/calibration/occurrences/<Species_safe>/mapveu/raw/
```

## Combined Vector Occurrences

After GBIF records have been fetched and local VectorMap/MapVEu records have
been extracted, use `occurrences/06_combine_vector_occurrences.R` to write a
new `combined` occurrence method per species. The script reads the source
folders only; it does not modify `gbif-download`, `vectormap`, `mapveu`,
`direct-gbif`, or `spatial-spp` inputs.

```r
batch_config <- list(
  roles = "vector",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  coordinate_round_digits = 5,
  dry_run = FALSE
)

source("scripts/sdms/present/occurrences/06_combine_vector_occurrences.R")
```

Combined outputs are written under:

```text
sdms/runs/chikungunya/calibration/occurrences/<Species_safe>/combined/
```

Each species gets standardized, coordinate-deduplicated, and cleaned CSVs plus
an occurrence-preparation summary. Run-level summaries are written under:

```text
sdms/runs/chikungunya/calibration/combined_vector_occurrence_runs/
```

## Model Runs

Use `models/01_run_present_model.R` for one species and
`models/02_run_chikungunya_models_batch.R` for manifest-driven batches.

The model batch script does not prepare occurrences. It expects cleaned
occurrence files to already exist under the configured occurrence method folder.
By default it is a status/preflight run only. Set either `dry_run_models = TRUE`
or `fit_models = TRUE` in the top `batch_config` block, or define
`batch_config` in the console immediately before sourcing:

```r
batch_config <- list(
  roles = "vector",
  occurrence_method = "combined",
  fit_models = TRUE,
  dry_run_models = FALSE,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y"))
)

source("scripts/sdms/present/models/02_run_chikungunya_models_batch.R")
```

Current model defaults use `predictor_mode = "bio-elev"` and
`candidate_set = "iucn_complete_all"`. `range_filter = "auto"` applies IUCN
ranges when available and falls back to the equivalent no-range candidate set
when a species has no matching range polygon.

Logs and summaries are written under:

```text
sdms/runs/chikungunya/calibration/model_batch_runs/
```
