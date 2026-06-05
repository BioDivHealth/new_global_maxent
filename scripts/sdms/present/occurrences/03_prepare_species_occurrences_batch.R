#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 03_prepare_species_occurrences_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: Generic one-pass occurrence batch entrypoint for present-day SDM
#          target manifests.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

generic_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
  occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  batch_run_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrence_batch_runs"),
  roles = "vector",
  occurrence_method = "direct-gbif",
  prepare_occurrences = FALSE,
  redownload_occurrences = FALSE,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y"))
)

batch_config <- if (exists("batch_config", inherits = FALSE)) {
  utils::modifyList(generic_batch_config, batch_config)
} else {
  generic_batch_config
}

source(file.path(
  here::here(),
  "scripts",
  "sdms",
  "present",
  "occurrences",
  "03_prepare_chikungunya_occurrences_batch.R"
))
