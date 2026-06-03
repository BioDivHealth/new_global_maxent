#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02_run_chikungunya_models_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run present-day SDMs across the Chikungunya SDM target manifest.
#
# Default behavior is status/preflight only. Edit the RStudio config block before
# sourcing to dry-run or fit models, or pass command-line arguments from Rscript.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# RStudio config: edit this block before sourcing the script ----
# -----------------------------------------------------------------------------|

if (!exists("batch_config", inherits = FALSE)) {
  batch_config <- list(
    roles = "vector",
    occurrence_method = "gbif-download",
    fit_models = FALSE,
    dry_run_models = FALSE,
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y"))
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_batch_config <- list(
  target_manifest_path = file.path(here::here(), "sdms", "runs", "chikungunya", "sdm_target_manifest.csv"),
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  occurrence_method = "gbif-download",
  fit_models = FALSE,
  dry_run_models = FALSE,
  skip_existing_models = TRUE,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  candidate_set = "iucn_complete_all",
  predictor_mode = "bio-elev",
  range_filter = "auto",
  range_mode = "strict",
  range_buffer = 4,
  min_obs = 20,
  n_background = 8000,
  test_percent = 30,
  beta_values = "4,8",
  n_selected_models = 10,
  maxent_threads = 4,
  predictor_stack_path = file.path(here::here(), "sdms", "cache", "Resample_rast.tif"),
  iucn_range_path = file.path(
    here::here(),
    "sdms",
    "cache",
    "MAMMALS_TERRESTRIAL_ONLY",
    "MAMMALS_TERRESTRIAL_ONLY.shp"
  )
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Config parsing helpers ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
}

expected_run_config_tag <- function(n_background, beta_values, n_selected_models, maxent_threads) {
  paste(
    paste0("bk", n_background),
    paste0("feature_grid_beta", gsub(",", "-", beta_values, fixed = TRUE)),
    paste0("select", n_selected_models),
    paste0("threads", maxent_threads),
    sep = "__"
  )
}

existing_model_for_config <- function(species_name, method, predictor_mode, start_year, end_year, run_config_tag) {
  species_safe <- safe_species_name(species_name)
  output_dir <- file.path(
    repo_root(),
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "regenerated_models",
    species_safe
  )
  if (!dir.exists(output_dir)) {
    return(NA_character_)
  }

  pattern <- paste0(
    "^",
    species_safe,
    "__",
    method,
    "__.*__",
    predictor_mode,
    "__",
    start_year,
    "_",
    end_year,
    "__",
    run_config_tag,
    "__model[.]rds$"
  )
  hits <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) {
    return(NA_character_)
  }

  normalizePath(sort(hits)[[1]], winslash = "/", mustWork = TRUE)
}

run_rscript <- function(script, script_args, log_path) {
  command <- file.path(R.home("bin"), "Rscript")
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(command, shQuote(c(script, script_args)), stdout = log_path, stderr = log_path)
  if (is.null(status)) {
    status <- 0L
  }

  as.integer(status)
}

# -----------------------------------------------------------------------------|
# Resolve batch settings ----
# -----------------------------------------------------------------------------|

target_manifest_path <- config_arg("target-manifest-path")
roles <- split_arg(config_arg("roles"))
species_filter <- split_arg(config_arg("species-filter"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
max_species <- as.numeric(config_arg("max-species"))
occurrence_method <- config_arg("occurrence-method")
fit_models <- as_logical_arg(config_arg("fit-models")) || has_flag(args, "fit-models")
dry_run_models <- as_logical_arg(config_arg("dry-run-models")) || has_flag(args, "dry-run-models")
skip_existing_models <- as_logical_arg(config_arg("skip-existing-models"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
candidate_set <- config_arg("candidate-set")
predictor_mode <- config_arg("predictor-mode")
range_filter <- config_arg("range-filter")
range_mode <- config_arg("range-mode")
range_buffer <- as.numeric(config_arg("range-buffer"))
min_obs <- as.integer(config_arg("min-obs"))
n_background <- as.integer(config_arg("n-background"))
test_percent <- as.integer(config_arg("test-percent"))
beta_values <- config_arg("beta-values")
n_selected_models <- as.integer(config_arg("n-selected-models"))
maxent_threads <- as.integer(config_arg("threads", "maxent_threads"))
predictor_stack_path <- config_arg("predictor-stack-path")
iucn_range_path <- config_arg("iucn-range-path")

if (!file.exists(target_manifest_path)) {
  stop("Missing Chikungunya SDM target manifest: ", target_manifest_path, call. = FALSE)
}

if (fit_models && dry_run_models) {
  stop("Use only one of `fit_models = TRUE` or `dry_run_models = TRUE`.", call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Select target species ----
# -----------------------------------------------------------------------------|

target_manifest <- read.csv(target_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
targets <- select_sdm_targets(
  target_manifest = target_manifest,
  roles = roles,
  species_filter = species_filter,
  include_not_needed = include_not_needed,
  include_already_available = include_already_available,
  max_species = max_species
)

cat("Selected target species:", nrow(targets), "\n")
if (nrow(targets) == 0) {
  warning(
    "No target species selected. Check roles/species_filter/include_not_needed/include_already_available settings.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------|
# Batch output paths ----
# -----------------------------------------------------------------------------|

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
batch_dir <- ensure_dir(file.path(
  repo_root(),
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "model_batch_runs",
  timestamp
))
log_dir <- ensure_dir(file.path(batch_dir, "logs"))
summary_path <- file.path(batch_dir, "model_batch_summary.csv")
run_config_tag <- expected_run_config_tag(n_background, beta_values, n_selected_models, maxent_threads)
script_model <- file.path(repo_root(), "scripts", "sdms", "present", "models", "01_run_present_model.R")

rows <- vector("list", nrow(targets))

# -----------------------------------------------------------------------------|
# Run optional model fitting ----
# -----------------------------------------------------------------------------|

for (i in seq_len(nrow(targets))) {
  target <- targets[i, , drop = FALSE]
  species <- target$species_name_canonical[[1]]
  species_safe <- safe_species_name(species)
  occurrence_path <- file.path(
    repo_root(),
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "occurrences",
    species_safe,
    occurrence_method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )
  existing_model_path <- existing_model_for_config(
    species,
    occurrence_method,
    predictor_mode,
    start_year,
    end_year,
    run_config_tag
  )

  occurrence_status <- if (file.exists(occurrence_path)) "ready" else "missing"
  model_status <- "not_requested"
  model_log <- NA_character_
  model_exit_status <- NA_integer_
  model_path <- existing_model_path

  if (skip_existing_models && !is.na(existing_model_path) && file.exists(existing_model_path)) {
    model_status <- "skipped_existing_model"
  } else if (!file.exists(occurrence_path)) {
    model_status <- "skipped_missing_occurrences"
  } else if (fit_models || dry_run_models) {
    model_log <- file.path(log_dir, paste0(species_safe, "__model.log"))
    model_args <- c(
      "--species", species,
      "--method", occurrence_method,
      "--candidate-set", candidate_set,
      "--predictor-mode", predictor_mode,
      "--range-filter", range_filter,
      "--range-mode", range_mode,
      "--range-buffer", as.character(range_buffer),
      "--start-year", as.character(start_year),
      "--end-year", as.character(end_year),
      "--min-obs", as.character(min_obs),
      "--n-background", as.character(n_background),
      "--test-percent", as.character(test_percent),
      "--beta-values", beta_values,
      "--n-selected-models", as.character(n_selected_models),
      "--threads", as.character(maxent_threads),
      "--predictor-stack", predictor_stack_path,
      "--iucn-range-path", iucn_range_path
    )
    if (fit_models) {
      model_args <- c(model_args, "--run")
    } else {
      model_args <- c(model_args, "--dry-run")
    }

    model_exit_status <- run_rscript(script_model, model_args, model_log)
    model_status <- if (model_exit_status == 0) {
      if (fit_models) "completed_or_prepared_by_model_script" else "dry_run_completed"
    } else {
      "failed"
    }
    model_path <- existing_model_for_config(
      species,
      occurrence_method,
      predictor_mode,
      start_year,
      end_year,
      run_config_tag
    )
  }

  rows[[i]] <- data.frame(
    species_name = species,
    manifest_species_name = target$species_name[[1]],
    species_role = target$species_role[[1]],
    sdm_needed_for_disease = target$sdm_needed_for_disease[[1]],
    run_priority = target$run_priority[[1]],
    sdm_available = target$sdm_available[[1]],
    manifest_run_status = target$run_status[[1]],
    occurrence_method = occurrence_method,
    occurrence_path = occurrence_path,
    occurrence_status = occurrence_status,
    model_status = model_status,
    model_exit_status = model_exit_status,
    model_log = model_log,
    model_path = model_path,
    candidate_set = candidate_set,
    predictor_mode = predictor_mode,
    start_year = start_year,
    end_year = end_year,
    n_background = n_background,
    n_selected_models = n_selected_models,
    maxent_threads = maxent_threads,
    prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  write.csv(do.call(rbind, rows[seq_len(i)]), summary_path, row.names = FALSE, na = "")
  cat("[", i, "/", nrow(targets), "] ", species, ": ", occurrence_status, ", ", model_status, "\n", sep = "")
}

# -----------------------------------------------------------------------------|
# Final batch summary ----
# -----------------------------------------------------------------------------|

summary <- if (length(rows) == 0) {
  data.frame()
} else {
  do.call(rbind, rows)
}
write.csv(summary, summary_path, row.names = FALSE, na = "")

cat("Wrote model batch summary:", summary_path, "\n")
cat("Batch directory:", batch_dir, "\n")
