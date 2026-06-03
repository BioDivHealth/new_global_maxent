#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_prepare_host_regeneration_manifest.R ----
# -----------------------------------------------------------------------------|
# Purpose: Build a small host-SDM calibration manifest for Chikungunya.
#
# The first calibration target is a host species with an existing Gonzalo SDM.
# This lets us compare a regenerated present-day SDM against the received model
# before attempting new vector SDMs.
#
# Output: sdms/runs/chikungunya/calibration/host_regeneration_manifest.csv
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package `terra` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# Arguments and paths ----
# -----------------------------------------------------------------------------|

species_filter <- commandArgs(trailingOnly = TRUE)
if (length(species_filter) == 0) {
  species_filter <- "Scotophilus kuhlii"
}

target_manifest_path <- file.path(
  repo_root(),
  "sdms",
  "runs",
  "chikungunya",
  "sdm_target_manifest.csv"
)
predictor_stack_path <- file.path(repo_root(), "sdms", "cache", "Resample_rast.tif")
calibration_dir <- ensure_dir(file.path(repo_root(), "sdms", "runs", "chikungunya", "calibration"))
output_path <- file.path(calibration_dir, "host_regeneration_manifest.csv")

# -----------------------------------------------------------------------------|
# Validate required inputs ----
# -----------------------------------------------------------------------------|

if (!file.exists(target_manifest_path)) {
  stop("Missing target manifest: ", target_manifest_path, call. = FALSE)
}

if (!file.exists(predictor_stack_path)) {
  stop("Missing predictor stack: ", predictor_stack_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Load target species ----
# -----------------------------------------------------------------------------|

target_manifest <- read.csv(
  target_manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

targets <- target_manifest[
  target_manifest$species_role == "host" &
    target_manifest$species_name %in% species_filter,
  ,
  drop = FALSE
]

if (nrow(targets) == 0) {
  stop(
    "No host rows found for requested species: ",
    paste(species_filter, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------|
# Summarise predictor stack ----
# -----------------------------------------------------------------------------|

predictors <- terra::rast(predictor_stack_path)
predictor_layer_names <- paste(names(predictors), collapse = "; ")
predictor_n_layers <- terra::nlyr(predictors)
predictor_crs <- terra::crs(predictors)
predictor_extent <- paste(as.vector(terra::ext(predictors)), collapse = "; ")
predictor_resolution <- paste(terra::res(predictors), collapse = "; ")
predictor_size_bytes <- file.info(predictor_stack_path)$size

# -----------------------------------------------------------------------------|
# Extract received-model calibration metadata ----
# -----------------------------------------------------------------------------|

rows <- lapply(seq_len(nrow(targets)), function(i) {
  row <- targets[i, , drop = FALSE]
  species <- row$species_name[[1]]
  existing_model_path <- file.path(
    repo_root(),
    "sdms",
    "models",
    species,
    paste0(species, ".rds")
  )

  if (!file.exists(existing_model_path)) {
    stop("Missing existing SDM for calibration species: ", existing_model_path, call. = FALSE)
  }

  # Read the received model only to record its settings; this script does not
  # alter or refit the model object.
  model_obj <- readRDS(existing_model_path)
  params <- model_obj$params
  aicc <- model_obj$AICc

  top_model <- if (is.data.frame(aicc) && nrow(aicc) > 0 && "mod" %in% names(aicc)) {
    aicc$mod[[1]]
  } else {
    NA_character_
  }

  data.frame(
    analysis_unit_id = row$analysis_unit_id,
    readiness_disease_name = row$readiness_disease_name,
    species_name = species,
    existing_model_path = existing_model_path,
    predictor_stack_path = predictor_stack_path,
    predictor_n_layers = predictor_n_layers,
    predictor_layer_names = predictor_layer_names,
    predictor_crs = predictor_crs,
    predictor_extent = predictor_extent,
    predictor_resolution = predictor_resolution,
    predictor_size_bytes = predictor_size_bytes,
    existing_model_n_presence = if (is.data.frame(params) && "n_presence" %in% names(params)) {
      coalesce_scalar(params$n_presence)
    } else {
      NA
    },
    existing_model_n_background = if (is.data.frame(params) && "n_background" %in% names(params)) {
      coalesce_scalar(params$n_background)
    } else {
      NA
    },
    existing_model_background_type = if (is.data.frame(params) && "bk_type" %in% names(params)) {
      collapse_unique(params$bk_type)
    } else {
      NA_character_
    },
    existing_model_selected_variables = collapse_unique(model_obj$variables),
    existing_model_count = if (!is.null(model_obj$mods)) length(model_obj$mods) else NA_integer_,
    existing_top_model = top_model,
    occurrence_input_path = NA_character_,
    occurrence_status = "missing_original_presence_points",
    planned_output_dir = file.path(
      repo_root(),
      "sdms",
      "runs",
      "chikungunya",
      "calibration",
      "regenerated_models",
      species
    ),
    calibration_status = "waiting_for_occurrence_points",
    notes = "Use the original presence coordinates if available; a fresh GBIF pull is a fallback and may not reproduce Gonzalo exactly.",
    stringsAsFactors = FALSE
  )
})

# -----------------------------------------------------------------------------|
# Write calibration manifest ----
# -----------------------------------------------------------------------------|

manifest <- do.call(rbind, rows)
write.csv(manifest, output_path, row.names = FALSE, na = "")

cat("Wrote calibration manifest:", output_path, "\n")
cat("Species:", paste(manifest$species_name, collapse = ", "), "\n")
cat("Predictor stack:", predictor_stack_path, "\n")
cat("Predictor layers:", predictor_n_layers, "\n")
