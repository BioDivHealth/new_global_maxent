#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02_compare_occurrences_to_existing_model.R ----
# -----------------------------------------------------------------------------|
# Purpose: Compare cleaned candidate occurrence records with an existing SDM.
#
# The existing AutoMaxent model does not keep original GBIF rows or coordinates,
# so the closest reproducible comparison is:
#   - occurrence-counts under candidate filters,
#   - model study-area bbox overlap,
#   - IUCN range overlap,
#   - selected-environment signature overlap against fitted MaxEnt presences.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package `terra` is required.", call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# Arguments and paths ----
# -----------------------------------------------------------------------------|

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
species <- get_arg(args, "species", "Artibeus jamaicensis")
methods <- strsplit(get_arg(args, "methods", "direct-gbif,spatial-spp"), ",", fixed = TRUE)[[1]]
methods <- trimws(methods)

species_safe <- safe_species_name(species)
repo <- repo_root()

occurrence_base <- file.path(
  repo,
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "occurrences",
  species_safe
)
diagnostics_dir <- ensure_dir(file.path(occurrence_base, "diagnostics"))
range_diagnostics_base <- ensure_dir(file.path(
  repo,
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "range_diagnostics",
  species_safe
))

model_path <- file.path(repo, "sdms", "models", species, paste0(species, ".rds"))
predictor_stack_path <- file.path(repo, "sdms", "cache", "Resample_rast.tif")
iucn_range_path <- get_arg(
  args,
  "iucn-range-path",
  file.path(repo, "sdms", "cache", "MAMMALS_TERRESTRIAL_ONLY", "MAMMALS_TERRESTRIAL_ONLY.shp")
)

if (!file.exists(model_path)) {
  stop("Missing existing SDM: ", model_path, call. = FALSE)
}
if (!file.exists(predictor_stack_path)) {
  stop("Missing predictor stack: ", predictor_stack_path, call. = FALSE)
}
if (!file.exists(iucn_range_path)) {
  stop("Missing IUCN range shapefile: ", iucn_range_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Load existing model and predictor context ----
# -----------------------------------------------------------------------------|

model <- readRDS(model_path)
predictors <- terra::rast(predictor_stack_path)
selected_variables <- model$variables
missing_variables <- setdiff(selected_variables, names(predictors))
if (length(missing_variables) > 0) {
  stop("Predictor stack is missing selected model variables: ", paste(missing_variables, collapse = ", "), call. = FALSE)
}

model_presence <- as.data.frame(model$mods[[1]]@presence)
model_n_presence <- nrow(model_presence)
model_bbox <- sf::st_bbox(model$study.area)
model_signature <- apply(round(model_presence[, selected_variables, drop = FALSE], 6), 1, paste, collapse = "|")
model_unique_signature <- unique(model_signature)

# -----------------------------------------------------------------------------|
# Load IUCN range used for spatial candidate filters ----
# -----------------------------------------------------------------------------|

iucn_ranges <- sf::st_read(iucn_range_path, quiet = TRUE)
species_range_all <- iucn_ranges[iucn_ranges$sci_name == species, ]
if (nrow(species_range_all) == 0) {
  stop("No IUCN range rows found for species: ", species, call. = FALSE)
}

species_range_strict <- species_range_all[
  species_range_all$presence == 1 &
    species_range_all$origin == 1 &
    species_range_all$seasonal == 1,
]
if (nrow(species_range_strict) == 0) {
  species_range_strict <- species_range_all
}
species_range_strict <- sf::st_make_valid(species_range_strict)
species_range_union <- sf::st_union(species_range_strict)
range_bbox <- sf::st_bbox(species_range_union)

bbox_summary <- data.frame(
  species_name = species,
  model_n_presence = model_n_presence,
  model_xmin = unname(model_bbox["xmin"]),
  model_ymin = unname(model_bbox["ymin"]),
  model_xmax = unname(model_bbox["xmax"]),
  model_ymax = unname(model_bbox["ymax"]),
  iucn_range_rows = nrow(species_range_all),
  iucn_strict_rows = nrow(species_range_strict),
  iucn_xmin = unname(range_bbox["xmin"]),
  iucn_ymin = unname(range_bbox["ymin"]),
  iucn_xmax = unname(range_bbox["xmax"]),
  iucn_ymax = unname(range_bbox["ymax"]),
  stringsAsFactors = FALSE
)
write.csv(
  bbox_summary,
  file.path(diagnostics_dir, paste0(species_safe, "_model_iucn_bbox_summary.csv")),
  row.names = FALSE,
  na = ""
)

date_to_flag <- function(x, cutoff) {
  # GBIF eventDate strings can include times, so compare on the YYYY-MM-DD part.
  dates <- suppressWarnings(as.Date(substr(x, 1, 10)))
  !is.na(dates) & dates <= as.Date(cutoff)
}

# -----------------------------------------------------------------------------|
# Build candidate-set diagnostics for each occurrence method ----
# -----------------------------------------------------------------------------|

summarise_candidate_sets <- function(method) {
  cleaned_path <- file.path(
    occurrence_base,
    method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )
  if (!file.exists(cleaned_path)) {
    stop("Missing cleaned occurrence file for method `", method, "`: ", cleaned_path, call. = FALSE)
  }

  cleaned <- read.csv(cleaned_path, check.names = FALSE, stringsAsFactors = FALSE)
  coordinate_key <- paste(cleaned$decimalLongitude, cleaned$decimalLatitude, sep = "|")
  cleaned <- cleaned[!duplicated(coordinate_key), , drop = FALSE]

  occurrence_points <- terra::vect(
    cleaned,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  )
  all_predictor_values <- terra::extract(predictors, occurrence_points, ID = FALSE)
  selected_predictor_values <- all_predictor_values[, selected_variables, drop = FALSE]

  # Keep both complete-all and complete-selected flags because Gonzalo's saved
  # model may only require a subset of the predictor stack.
  cleaned$complete_all_predictors <- complete.cases(all_predictor_values)
  cleaned$complete_selected_predictors <- complete.cases(selected_predictor_values)
  cleaned$inside_model_bbox <- cleaned$decimalLongitude >= model_bbox["xmin"] &
    cleaned$decimalLongitude <= model_bbox["xmax"] &
    cleaned$decimalLatitude >= model_bbox["ymin"] &
    cleaned$decimalLatitude <= model_bbox["ymax"]
  cleaned$basis_human_observation <- if ("basisOfRecord" %in% names(cleaned)) {
    cleaned$basisOfRecord == "HUMAN_OBSERVATION"
  } else {
    NA
  }
  cleaned$github_issues_ok <- if ("issues" %in% names(cleaned)) {
    is.na(cleaned$issues) | cleaned$issues != "cdc,cdround,gass84,muluriiv"
  } else {
    NA
  }

  years <- suppressWarnings(as.integer(cleaned$year))
  cleaned$year_to_2024 <- !is.na(years) & years <= 2024
  cleaned$year_to_2025 <- !is.na(years) & years <= 2025
  cleaned$date_to_2025_01_31 <- date_to_flag(cleaned$eventDate, "2025-01-31")

  occurrence_sf <- sf::st_as_sf(
    cleaned,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )
  cleaned$inside_iucn_strict_range <- lengths(sf::st_intersects(occurrence_sf, species_range_union)) > 0

  # The model object does not retain GBIF rows, but the fitted MaxEnt presences
  # do retain environmental values. Rounded signatures give a reproducible
  # approximate overlap check between candidate points and saved presences.
  selected_signature <- rep(NA_character_, nrow(cleaned))
  complete_selected <- cleaned$complete_selected_predictors
  selected_signature[complete_selected] <- apply(
    round(selected_predictor_values[complete_selected, selected_variables, drop = FALSE], 6),
    1,
    paste,
    collapse = "|"
  )
  cleaned$selected_env_signature <- selected_signature
  cleaned$matches_gonzalo_presence_signature <- !is.na(selected_signature) &
    selected_signature %in% model_unique_signature

  candidate_sets <- list(
    cleaned_unique = rep(TRUE, nrow(cleaned)),
    complete_all_predictors = cleaned$complete_all_predictors,
    complete_selected_predictors = cleaned$complete_selected_predictors,
    human_observation = cleaned$basis_human_observation %in% TRUE,
    github_human_and_issues = cleaned$basis_human_observation %in% TRUE &
      cleaned$github_issues_ok %in% TRUE,
    year_to_2024 = cleaned$year_to_2024,
    year_to_2025 = cleaned$year_to_2025,
    date_to_2025_01_31 = cleaned$date_to_2025_01_31,
    model_bbox = cleaned$inside_model_bbox,
    iucn_strict_range = cleaned$inside_iucn_strict_range,
    iucn_model_bbox = cleaned$inside_iucn_strict_range & cleaned$inside_model_bbox,
    iucn_complete_all = cleaned$inside_iucn_strict_range & cleaned$complete_all_predictors,
    iucn_complete_selected = cleaned$inside_iucn_strict_range & cleaned$complete_selected_predictors,
    model_bbox_complete_all = cleaned$inside_model_bbox & cleaned$complete_all_predictors,
    model_bbox_complete_selected = cleaned$inside_model_bbox & cleaned$complete_selected_predictors,
    iucn_model_bbox_complete_all = cleaned$inside_iucn_strict_range &
      cleaned$inside_model_bbox &
      cleaned$complete_all_predictors,
    iucn_model_bbox_complete_selected = cleaned$inside_iucn_strict_range &
      cleaned$inside_model_bbox &
      cleaned$complete_selected_predictors,
    iucn_human_observation = cleaned$inside_iucn_strict_range &
      cleaned$basis_human_observation %in% TRUE,
    iucn_human_complete_all = cleaned$inside_iucn_strict_range &
      cleaned$basis_human_observation %in% TRUE &
      cleaned$complete_all_predictors,
    iucn_human_complete_selected = cleaned$inside_iucn_strict_range &
      cleaned$basis_human_observation %in% TRUE &
      cleaned$complete_selected_predictors,
    iucn_year_to_2024 = cleaned$inside_iucn_strict_range & cleaned$year_to_2024,
    iucn_year_to_2024_complete_all = cleaned$inside_iucn_strict_range &
      cleaned$year_to_2024 &
      cleaned$complete_all_predictors,
    iucn_year_to_2024_complete_selected = cleaned$inside_iucn_strict_range &
      cleaned$year_to_2024 &
      cleaned$complete_selected_predictors
  )

  summary <- dplyr::bind_rows(lapply(names(candidate_sets), function(candidate_name) {
    index <- candidate_sets[[candidate_name]] %in% TRUE
    candidate_signatures <- unique(cleaned$selected_env_signature[index & !is.na(cleaned$selected_env_signature)])
    data.frame(
      method = method,
      candidate_set = candidate_name,
      candidate_unique_coords = sum(index),
      diff_vs_gonzalo_n = sum(index) - model_n_presence,
      model_rows_represented_by_env_signature = sum(model_signature %in% candidate_signatures),
      model_unique_signatures_represented = sum(model_unique_signature %in% candidate_signatures),
      candidate_rows_with_model_env_signature = sum(
        cleaned$selected_env_signature[index] %in% model_unique_signature,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::arrange(abs(diff_vs_gonzalo_n), dplyr::desc(model_rows_represented_by_env_signature))

  candidate_bbox_summary <- dplyr::bind_rows(lapply(names(candidate_sets), function(candidate_name) {
    index <- candidate_sets[[candidate_name]] %in% TRUE
    if (!any(index)) {
      return(data.frame(
        method = method,
        candidate_set = candidate_name,
        rows = 0,
        xmin = NA_real_,
        ymin = NA_real_,
        xmax = NA_real_,
        ymax = NA_real_
      ))
    }
    data.frame(
      method = method,
      candidate_set = candidate_name,
      rows = sum(index),
      xmin = min(cleaned$decimalLongitude[index], na.rm = TRUE),
      ymin = min(cleaned$decimalLatitude[index], na.rm = TRUE),
      xmax = max(cleaned$decimalLongitude[index], na.rm = TRUE),
      ymax = max(cleaned$decimalLatitude[index], na.rm = TRUE)
    )
  }))

  method_range_dir <- ensure_dir(file.path(range_diagnostics_base, method))
  write.csv(
    cleaned,
    file.path(diagnostics_dir, paste0(species_safe, "_", method, "_candidate_occurrences_with_overlap_flags.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    summary,
    file.path(diagnostics_dir, paste0(species_safe, "_", method, "_model_presence_environment_overlap_summary.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    candidate_bbox_summary,
    file.path(diagnostics_dir, paste0(species_safe, "_", method, "_candidate_bbox_summary.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    cleaned,
    file.path(method_range_dir, "candidate_occurrences_with_iucn_range_flags.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    summary,
    file.path(method_range_dir, "iucn_range_count_summary.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    candidate_bbox_summary,
    file.path(method_range_dir, "iucn_range_bbox_summary.csv"),
    row.names = FALSE,
    na = ""
  )

  summary
}

# -----------------------------------------------------------------------------|
# Write combined diagnostics ----
# -----------------------------------------------------------------------------|

all_summaries <- dplyr::bind_rows(lapply(methods, summarise_candidate_sets))
all_summary_path <- file.path(diagnostics_dir, paste0(species_safe, "_all_methods_model_overlap_summary.csv"))
write.csv(all_summaries, all_summary_path, row.names = FALSE, na = "")

cat("Model presences:", model_n_presence, "\n")
cat("Model bbox:", paste(round(unname(model_bbox), 5), collapse = ", "), "\n")
cat("IUCN strict bbox:", paste(round(unname(range_bbox), 5), collapse = ", "), "\n")
cat("Wrote:", all_summary_path, "\n")
print(dplyr::group_by(all_summaries, method) |> dplyr::slice_head(n = 12) |> dplyr::ungroup(), n = Inf)
