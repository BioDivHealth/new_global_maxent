#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02b_run_manifest_model_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: Generic manifest-driven present-day SDM batch implementation.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

defaults <- list(
  target_manifest_path = file.path(here::here(), "sdms", "runs", "chikungunya", "sdm_target_manifest.csv"),
  occurrence_root = file.path(here::here(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  model_output_root = file.path(here::here(), "sdms", "runs", "chikungunya", "calibration", "regenerated_models"),
  model_batch_run_root = file.path(here::here(), "sdms", "runs", "chikungunya", "calibration", "model_batch_runs"),
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  occurrence_method = "gbif-download",
  fit_models = FALSE,
  dry_run_models = FALSE,
  skip_existing_models = TRUE,
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
  random_features = TRUE,
  n_models = 25,
  n_selected_models = 10,
  use_boyce = 0.5,
  maxent_threads = 2,
  java_memory_gb = 8,
  predictor_stack_path = file.path(here::here(), "sdms", "cache", "Resample_rast.tif"),
  iucn_range_path = file.path(here::here(), "sdms", "cache", "MAMMALS_TERRESTRIAL_ONLY", "MAMMALS_TERRESTRIAL_ONLY.shp"),
  automaxent_root = Sys.getenv("AUTOMAXENT_ROOT", unset = "/Users/arturtrebski/Coding_Projects/AutoMaxent")
)

if (!exists("batch_config", inherits = FALSE)) {
  batch_config <- list()
}

cfg <- utils::modifyList(defaults, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
get_cfg <- function(key, field = gsub("-", "_", key)) get_arg(args, key, cfg[[field]])

cfg$target_manifest_path <- get_cfg("target-manifest-path")
cfg$occurrence_root <- get_cfg("occurrence-root")
cfg$model_output_root <- get_cfg("model-output-root")
cfg$model_batch_run_root <- get_cfg("model-batch-run-root")
cfg$roles <- split_arg(get_cfg("roles"))
cfg$species_filter <- split_arg(get_cfg("species-filter"))
cfg$include_not_needed <- as_logical_arg(get_cfg("include-not-needed"))
cfg$include_already_available <- as_logical_arg(get_cfg("include-already-available"))
cfg$max_species <- as.numeric(get_cfg("max-species"))
cfg$occurrence_method <- get_cfg("occurrence-method")
cfg$fit_models <- as_logical_arg(get_cfg("fit-models")) || has_flag(args, "fit-models")
cfg$dry_run_models <- as_logical_arg(get_cfg("dry-run-models")) || has_flag(args, "dry-run-models")
cfg$skip_existing_models <- as_logical_arg(get_cfg("skip-existing-models"))
cfg$start_year <- as.integer(get_cfg("start-year"))
cfg$end_year <- as.integer(get_cfg("end-year"))
cfg$candidate_set <- get_cfg("candidate-set")
cfg$predictor_mode <- get_cfg("predictor-mode")
cfg$range_filter <- get_cfg("range-filter")
cfg$range_mode <- get_cfg("range-mode")
cfg$range_buffer <- as.numeric(get_cfg("range-buffer"))
cfg$min_obs <- as.integer(get_cfg("min-obs"))
cfg$n_background <- get_cfg("n-background")
cfg$test_percent <- get_cfg("test-percent")
cfg$beta_values <- get_cfg("beta-values")
cfg$random_features <- as_logical_arg(get_cfg("random-features")) || has_flag(args, "random-features")
if (has_flag(args, "feature-grid")) {
  cfg$random_features <- FALSE
}
cfg$n_models <- as.integer(get_cfg("n-models"))
cfg$n_selected_models <- as.integer(get_cfg("n-selected-models"))
cfg$use_boyce <- get_cfg("use-boyce")
cfg$maxent_threads <- as.integer(get_cfg("threads", "maxent_threads"))
cfg$java_memory_gb <- as.numeric(get_cfg("java-memory-gb"))
cfg$predictor_stack_path <- get_cfg("predictor-stack-path")
cfg$iucn_range_path <- get_cfg("iucn-range-path")
cfg$automaxent_root <- get_cfg("automaxent-root")

if (!file.exists(cfg$target_manifest_path)) {
  stop("Missing SDM target manifest: ", cfg$target_manifest_path, call. = FALSE)
}
if (cfg$fit_models && cfg$dry_run_models) {
  stop("Use only one of `fit_models = TRUE` or `dry_run_models = TRUE`.", call. = FALSE)
}

model_grid_tag <- if (cfg$random_features) {
  paste0("random", cfg$n_models, "_beta", gsub(",", "-", cfg$beta_values, fixed = TRUE))
} else {
  paste0("feature_grid_beta", gsub(",", "-", cfg$beta_values, fixed = TRUE))
}
boyce_tag <- if (cfg$use_boyce %in% c("", "NA", "NULL", "FALSE")) "boyceNA" else paste0("boyce", cfg$use_boyce)
background_tag <- if (cfg$n_background == "dynamic") "bk[0-9]+" else paste0("bk", cfg$n_background)
run_config_tag <- paste(
  background_tag,
  model_grid_tag,
  paste0("select", cfg$n_selected_models),
  boyce_tag,
  paste0("threads", cfg$maxent_threads),
  sep = "__"
)

find_existing_model <- function(species) {
  species_safe <- safe_species_name(species)
  output_dir <- file.path(cfg$model_output_root, species_safe)
  if (!dir.exists(output_dir)) {
    return(NA_character_)
  }

  pattern <- paste0(
    "^", species_safe, "__", cfg$occurrence_method, "__.*__",
    cfg$predictor_mode, "__", cfg$start_year, "_", cfg$end_year, "__",
    run_config_tag, "__model[.]rds$"
  )
  hits <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) NA_character_ else normalizePath(sort(hits)[[1]], winslash = "/", mustWork = TRUE)
}

run_model_script <- function(model_args, log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(model_args), stdout = log_path, stderr = log_path)
  if (is.null(status)) 0L else as.integer(status)
}

targets <- select_sdm_targets(
  target_manifest = read.csv(cfg$target_manifest_path, check.names = FALSE, stringsAsFactors = FALSE),
  roles = cfg$roles,
  species_filter = cfg$species_filter,
  include_not_needed = cfg$include_not_needed,
  include_already_available = cfg$include_already_available,
  max_species = cfg$max_species
)

cat("Selected target species:", nrow(targets), "\n")
if (nrow(targets) == 0) {
  warning("No target species selected. Check batch filters.", call. = FALSE)
}

batch_dir <- ensure_dir(file.path(
  cfg$model_batch_run_root,
  paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
))
log_dir <- ensure_dir(file.path(batch_dir, "logs"))
summary_path <- file.path(batch_dir, "model_batch_summary.csv")
model_script <- file.path(repo_root(), "scripts", "sdms", "present", "models", "01_run_present_model.R")
rows <- vector("list", nrow(targets))

for (i in seq_len(nrow(targets))) {
  target <- targets[i, , drop = FALSE]
  species <- target$species_name_canonical[[1]]
  species_safe <- safe_species_name(species)
  occurrence_path <- file.path(
    cfg$occurrence_root,
    species_safe,
    cfg$occurrence_method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )
  existing_model_path <- find_existing_model(species)
  occurrence_status <- if (file.exists(occurrence_path)) "ready" else "missing"
  model_status <- "not_requested"
  model_exit_status <- NA_integer_
  model_log <- NA_character_
  model_path <- existing_model_path

  if (cfg$skip_existing_models && !is.na(existing_model_path) && file.exists(existing_model_path)) {
    model_status <- "skipped_existing_model"
  } else if (!file.exists(occurrence_path)) {
    model_status <- "skipped_missing_occurrences"
  } else if (cfg$fit_models || cfg$dry_run_models) {
    model_log <- file.path(log_dir, paste0(species_safe, "__model.log"))
    model_args <- c(
      model_script,
      "--species", species,
      "--method", cfg$occurrence_method,
      "--occurrences", occurrence_path,
      "--candidate-set", cfg$candidate_set,
      "--predictor-mode", cfg$predictor_mode,
      "--range-filter", cfg$range_filter,
      "--range-mode", cfg$range_mode,
      "--range-buffer", as.character(cfg$range_buffer),
      "--start-year", as.character(cfg$start_year),
      "--end-year", as.character(cfg$end_year),
      "--min-obs", as.character(cfg$min_obs),
      "--n-background", as.character(cfg$n_background),
      "--test-percent", as.character(cfg$test_percent),
      "--beta-values", cfg$beta_values,
      "--n-models", as.character(cfg$n_models),
      "--n-selected-models", as.character(cfg$n_selected_models),
      "--use-boyce", as.character(cfg$use_boyce),
      "--threads", as.character(cfg$maxent_threads),
      "--java-memory-gb", as.character(cfg$java_memory_gb),
      "--predictor-stack", cfg$predictor_stack_path,
      "--iucn-range-path", cfg$iucn_range_path,
      "--automaxent-root", cfg$automaxent_root,
      "--output-root", cfg$model_output_root,
      if (cfg$random_features) "--random-features" else "--feature-grid",
      if (cfg$fit_models) "--run" else "--dry-run"
    )

    model_exit_status <- run_model_script(model_args, model_log)
    model_status <- if (model_exit_status == 0) {
      if (cfg$fit_models) "completed_or_prepared_by_model_script" else "dry_run_completed"
    } else {
      "failed"
    }
    model_path <- find_existing_model(species)
  }

  rows[[i]] <- data.frame(
    species_name = species,
    manifest_species_name = target$species_name[[1]],
    species_role = target$species_role[[1]],
    sdm_needed_for_disease = target$sdm_needed_for_disease[[1]],
    run_priority = target$run_priority[[1]],
    sdm_available = target$sdm_available[[1]],
    manifest_run_status = target$run_status[[1]],
    occurrence_method = cfg$occurrence_method,
    occurrence_path = occurrence_path,
    occurrence_status = occurrence_status,
    model_status = model_status,
    model_exit_status = model_exit_status,
    model_log = model_log,
    model_path = model_path,
    candidate_set = cfg$candidate_set,
    predictor_mode = cfg$predictor_mode,
    start_year = cfg$start_year,
    end_year = cfg$end_year,
    n_background = cfg$n_background,
    beta_values = cfg$beta_values,
    random_features = cfg$random_features,
    n_models = cfg$n_models,
    n_selected_models = cfg$n_selected_models,
    use_boyce = cfg$use_boyce,
    maxent_threads = cfg$maxent_threads,
    java_memory_gb = cfg$java_memory_gb,
    prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  write.csv(do.call(rbind, rows[seq_len(i)]), summary_path, row.names = FALSE, na = "")
  cat("[", i, "/", nrow(targets), "] ", species, ": ", occurrence_status, ", ", model_status, "\n", sep = "")
}

summary <- if (length(rows) == 0) {
  data.frame()
} else {
  do.call(rbind, rows)
}
write.csv(summary, summary_path, row.names = FALSE, na = "")
cat("Wrote model batch summary:", summary_path, "\n")
cat("Batch directory:", batch_dir, "\n")
