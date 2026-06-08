#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02_run_species_models_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: Generic manifest-driven present-day SDM batch entrypoint.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

batch_config <- list(
  target_manifest_path = "sdms/runs/vector_sdm_push/vector_species_sdm_targets.csv",
  occurrence_root = "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/occurrences",
  model_output_root = "/Volumes/LaCie/new_global_maxent/sdms/models_artur/vector_sdm_push",
  model_batch_run_root = "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push/model_batch_runs",
  
  roles = "vector",
  occurrence_method = "combined",
  
  species_filter = paste(c(
    "Opifex fuscus",
    "Aedes togoi",
    "Eretmapodites chrysogaster",
    "Aedes procax",
    "Verrallina funerea",
    "Aedes africanus",
    "Aedes vittatus",
    "Coquillettidia linealis",
    "Culex sitiens",
    "Culex annulirostris",
    "Aedes vigilax",
    "Aedes notoscriptus",
    "Aedes triseriatus",
    "Aedes vexans",
    "Aedes aegypti",
    "Aedes albopictus"
  ), collapse = ","),
  
  fit_models = TRUE,
  dry_run_models = FALSE,
  
  start_year = 1985,
  end_year = 2026,
  maxent_threads = 4
)



source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

generic_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
  occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  model_output_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "models"),
  model_batch_run_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "model_batch_runs"),
  roles = "vector",
  occurrence_method = "combined",
  fit_models = FALSE,
  dry_run_models = FALSE,
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
  "models",
  "02_run_chikungunya_models_batch.R"
))
