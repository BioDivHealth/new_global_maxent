#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02a_run_species_models_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: User-facing config wrapper for present-day species SDM batches.
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
    #"Aedes togoi",
    "Eretmapodites chrysogaster",
    #"Aedes procax",
    #"Verrallina funerea",
    "Aedes africanus"
    #"Aedes vittatus",
    #"Coquillettidia linealis",
    #"Culex sitiens",
    #"Culex annulirostris",
    #"Aedes vigilax",
    #"Aedes notoscriptus",
    #"Aedes triseriatus",
    #"Aedes vexans",
    #"Aedes aegypti",
    #"Aedes albopictus"
  ), collapse = ","),

  fit_models = TRUE,
  dry_run_models = FALSE,
  skip_existing_models = FALSE,

  start_year = 2000,
  end_year = 2026,

  candidate_set = "iucn_complete_all",
  predictor_mode = "bio-elev",
  range_filter = "auto",
  range_mode = "strict",
  range_buffer = 4,
  min_obs = 20,
  n_background = "dynamic",
  test_percent = "dynamic",
  beta_values = "4,8,12",
  #random_features = FALSE, 
  random_features = TRUE,
  n_models = 25,
  n_selected_models = 25,
  #n_selected_models = 10,
  #use_boyce = 0.5,
  use_boyce = -1,

  maxent_threads = 2,
  java_memory_gb = 8
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
  start_year = 2000,
  end_year = 2026,
  candidate_set = "iucn_complete_all",
  predictor_mode = "bio-elev",
  range_filter = "auto",
  range_mode = "strict",
  range_buffer = 4,
  min_obs = 20,
  n_background = "dynamic",
  test_percent = "dynamic",
  beta_values = "4,8,12",
  #random_features = TRUE,
  random_features = FALSE,
  n_models = 25,
  #n_selected_models = 10,
  n_selected_models = 25,
  #use_boyce = 0.5,
  use_boyce = -1,
  maxent_threads = 2,
  java_memory_gb = 8
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
  "02b_run_manifest_model_batch.R"
))
