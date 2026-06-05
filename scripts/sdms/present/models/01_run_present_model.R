#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_run_present_model.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run or dry-run one present-day SDM with the AutoMaxent settings
#          inferred from Gonzalo's existing host models.
#
# Default behavior is controlled by the RStudio-friendly config block below.
# Command-line arguments can still override these values when running via Rscript.
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
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# RStudio config ----
# -----------------------------------------------------------------------------|
# Edit this block, then source the script from the repository root.

sdm_config <- list(
  species = "Rousettus aegyptiacus",
  method = "spatial-spp",
  candidate_set = "iucn_complete_all",
  start_year = 2015,
  end_year = 2025,
  run_model = TRUE,
  range_filter = "auto",
  range_mode = "strict",
  range_buffer = 4,
  predictor_mode = "bio-elev",
  select_var = "NUMERICAL",
  min_obs = 20,
  n_background = 8000,
  test_percent = 30,
  beta_values = "4,8",
  n_models = 10,
  n_selected_models = 10,
  seed = 185,
  random_features = FALSE,
  maxent_threads = 4,
  predictor_stack_path = file.path(here::here(), "sdms", "cache", "Resample_rast.tif"),
  iucn_range_path = file.path(
    here::here(),
    "sdms",
    "cache",
    "MAMMALS_TERRESTRIAL_ONLY",
    "MAMMALS_TERRESTRIAL_ONLY.shp"
  ),
  output_root = file.path(
    here::here(),
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "regenerated_models"
  ),
  automaxent_root = "/Users/arturtrebski/Coding_Projects/AutoMaxent"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Resolve config and command-line overrides ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, sdm_config[[config_key]])
}

species <- config_arg("species")
method <- config_arg("method")
candidate_set <- config_arg("candidate-set")
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
run_model <- as.logical(config_arg("run-model")) || has_flag(args, "run")
if (has_flag(args, "dry-run")) {
  run_model <- FALSE
}
run_started_at <- NA_character_
run_finished_at <- NA_character_
run_elapsed_seconds <- NA_real_

range_filter <- config_arg("range-filter")
range_mode <- config_arg("range-mode")
range_buffer <- as.numeric(config_arg("range-buffer"))
predictor_mode <- config_arg("predictor-mode")
select_var <- config_arg("select-var")
if (identical(select_var, "FALSE")) {
  select_var <- FALSE
}

min_obs <- as.integer(config_arg("min-obs"))
n_background <- as.integer(config_arg("n-background"))
test_percent <- as.integer(config_arg("test-percent"))
beta_values <- as.numeric(strsplit(config_arg("beta-values"), ",", fixed = TRUE)[[1]])
n_models <- as.integer(config_arg("n-models"))
n_selected_models <- as.integer(config_arg("n-selected-models"))
seed <- as.integer(config_arg("seed"))
random_features <- as.logical(config_arg("random-features")) || has_flag(args, "random-features")
if (has_flag(args, "feature-grid")) {
  random_features <- FALSE
}
maxent_threads <- as.integer(config_arg("threads", "maxent_threads"))
candidate_model_grid_size <- if (random_features) {
  n_models
} else {
  sum(vapply(1:4, function(k) choose(5, k), numeric(1))) * length(beta_values)
}

repo <- repo_root()
species_safe <- safe_species_name(species)

occurrence_path <- get_arg(
  args,
  "occurrences",
  file.path(
    repo,
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "occurrences",
    species_safe,
    method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )
)

predictor_stack_path <- get_arg(
  args,
  "predictor-stack",
  sdm_config$predictor_stack_path
)

iucn_range_path <- get_arg(
  args,
  "iucn-range-path",
  sdm_config$iucn_range_path
)

existing_model_path <- get_arg(
  args,
  "existing-model",
  file.path(repo, "sdms", "models", species, paste0(species, ".rds"))
)

automaxent_root <- get_arg(args, "automaxent-root", sdm_config$automaxent_root)
output_root <- get_arg(args, "output-root", sdm_config$output_root)
output_dir <- ensure_dir(file.path(output_root, species_safe))
model_work_dir <- ensure_dir(file.path(output_dir, "maxent_work"))

# -----------------------------------------------------------------------------|
# Validate requested run settings ----
# -----------------------------------------------------------------------------|

if (!candidate_set %in% c(
  "cleaned_unique",
  "complete_all_predictors",
  "complete_selected_predictors",
  "iucn_strict_range",
  "iucn_complete_all",
  "iucn_complete_selected"
)) {
  stop(
    "`--candidate-set` must be one of cleaned_unique, complete_all_predictors, ",
    "complete_selected_predictors, iucn_strict_range, iucn_complete_all, ",
    "iucn_complete_selected.",
    call. = FALSE
  )
}

if (!range_filter %in% c("auto", "apply", "none")) {
  stop("`--range-filter` must be one of: auto, apply, none", call. = FALSE)
}

if (!range_mode %in% c("strict", "all")) {
  stop("`--range-mode` must be one of: strict, all", call. = FALSE)
}

if (!predictor_mode %in% c("bio-elev", "all", "existing-selected")) {
  stop("`--predictor-mode` must be one of: bio-elev, all, existing-selected", call. = FALSE)
}

if (!file.exists(occurrence_path)) {
  stop("Missing cleaned occurrence file: ", occurrence_path, call. = FALSE)
}

if (!file.exists(predictor_stack_path)) {
  stop("Missing predictor stack: ", predictor_stack_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# AutoMaxent compatibility helpers ----
# -----------------------------------------------------------------------------|

load_automaxent <- function(root) {
  function_dir <- file.path(root, "Functions")
  needed <- c(
    "auto_MaxEnt_complementary.R",
    "BackgroundPOINTS.R",
    "Environmental_weigthing_random_points.R",
    "Time_matchine.R",
    "auto_MaxEnt.R"
  )
  paths <- file.path(function_dir, needed)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing AutoMaxent function files: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  invisible(lapply(paths, source))
}

date_to_year <- function(x) {
  # GBIF eventDate strings can include times, so compare on the YYYY-MM-DD part.
  event_date <- suppressWarnings(as.Date(substr(x, 1, 10)))
  as.integer(format(event_date, "%Y"))
}

pick_species_range <- function(path, species_name, mode) {
  # In auto mode, missing range data falls back to occurrence-only candidate
  # filters; in apply mode it is treated as a hard error.
  if (!file.exists(path)) {
    if (range_filter == "apply") {
      stop("Missing IUCN range shapefile: ", path, call. = FALSE)
    }
    return(NULL)
  }

  ranges <- sf::st_read(path, quiet = TRUE)
  if (!"sci_name" %in% names(ranges)) {
    stop("IUCN range file is missing `sci_name`: ", path, call. = FALSE)
  }

  species_range <- ranges[ranges$sci_name == species_name, ]
  if (nrow(species_range) == 0) {
    if (range_filter == "apply") {
      stop("No IUCN range rows found for species: ", species_name, call. = FALSE)
    }
    return(NULL)
  }

  if (mode == "strict" && all(c("presence", "origin", "seasonal") %in% names(species_range))) {
    strict_range <- species_range[
      species_range$presence == 1 &
        species_range$origin == 1 &
        species_range$seasonal == 1,
    ]
    if (nrow(strict_range) > 0) {
      species_range <- strict_range
    }
  }

  # AutoMaxent expects a single study-area geometry in the same CRS as points.
  species_range <- sf::st_transform(sf::st_make_valid(species_range), 4326)
  sf::st_as_sf(sf::st_union(species_range))
}

summarise_bbox <- function(x) {
  if (is.null(x)) {
    return(c(xmin = NA_real_, ymin = NA_real_, xmax = NA_real_, ymax = NA_real_))
  }
  unname(sf::st_bbox(x))[c(1, 2, 3, 4)] |>
    stats::setNames(c("xmin", "ymin", "xmax", "ymax"))
}

# -----------------------------------------------------------------------------|
# Load and pre-filter occurrences ----
# -----------------------------------------------------------------------------|

occurrences <- read.csv(occurrence_path, check.names = FALSE, stringsAsFactors = FALSE)
required_xy <- c("decimalLongitude", "decimalLatitude")
missing_xy <- setdiff(required_xy, names(occurrences))
if (length(missing_xy) > 0) {
  stop("Occurrence file is missing coordinate columns: ", paste(missing_xy, collapse = ", "), call. = FALSE)
}

coord_key <- paste(occurrences$decimalLongitude, occurrences$decimalLatitude, sep = "|")
occurrences <- occurrences[!duplicated(coord_key), , drop = FALSE]

years <- if ("year" %in% names(occurrences)) {
  suppressWarnings(as.integer(occurrences$year))
} else if ("eventDate" %in% names(occurrences)) {
  date_to_year(occurrences$eventDate)
} else {
  rep(NA_integer_, nrow(occurrences))
}
occurrences$sdm_year <- years
occurrences <- occurrences[
  !is.na(occurrences$sdm_year) &
    occurrences$sdm_year >= start_year &
    occurrences$sdm_year <= end_year,
  ,
  drop = FALSE
]

if (nrow(occurrences) < min_obs) {
  stop(
    "Fewer than `--min-obs` records remain after date filtering for ",
    species,
    ": ",
    nrow(occurrences),
    call. = FALSE
  )
}

predictors <- terra::rast(predictor_stack_path)
existing_model <- if (file.exists(existing_model_path)) readRDS(existing_model_path) else NULL
existing_selected_variables <- if (!is.null(existing_model)) existing_model$variables else character()

# -----------------------------------------------------------------------------|
# Select predictor layers passed to AutoMaxent ----
# -----------------------------------------------------------------------------|

if (predictor_mode == "bio-elev") {
  bio_elev_variables <- names(predictors)[
    grepl("wc2[.]1_2[.]5m_bio_", names(predictors)) |
      grepl("elev", names(predictors), ignore.case = TRUE)
  ]
  if (length(bio_elev_variables) == 0) {
    stop("`--predictor-mode bio-elev` found no bioclim/elevation layers in predictor stack.", call. = FALSE)
  }
  predictors_for_model <- predictors[[bio_elev_variables]]
} else if (predictor_mode == "existing-selected") {
  if (length(existing_selected_variables) == 0) {
    stop("`--predictor-mode existing-selected` requires an existing model with variables.", call. = FALSE)
  }
  missing_variables <- setdiff(existing_selected_variables, names(predictors))
  if (length(missing_variables) > 0) {
    stop("Predictor stack is missing existing selected variables: ", paste(missing_variables, collapse = ", "), call. = FALSE)
  }
  predictors_for_model <- predictors[[existing_selected_variables]]
} else {
  predictors_for_model <- predictors
}

# -----------------------------------------------------------------------------|
# Apply predictor completeness and optional range filters ----
# -----------------------------------------------------------------------------|

occurrence_points <- terra::vect(
  occurrences,
  geom = c("decimalLongitude", "decimalLatitude"),
  crs = "EPSG:4326"
)
predictor_values <- terra::extract(predictors_for_model, occurrence_points, ID = FALSE)
occurrences$complete_model_predictors <- complete.cases(predictor_values)

species_range <- NULL
range_status <- "not_requested"
if (range_filter %in% c("auto", "apply")) {
  species_range <- pick_species_range(iucn_range_path, species, range_mode)
  range_status <- if (is.null(species_range)) "not_available" else "available"
}

if (!is.null(species_range)) {
  occurrence_sf_for_filter <- sf::st_as_sf(
    occurrences,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )
  occurrences$inside_iucn_range <- lengths(sf::st_intersects(occurrence_sf_for_filter, species_range)) > 0
} else {
  occurrences$inside_iucn_range <- NA
}

# -----------------------------------------------------------------------------|
# Select candidate records for model fitting ----
# -----------------------------------------------------------------------------|

effective_candidate_set <- candidate_set
if (is.null(species_range) && candidate_set %in% c("iucn_strict_range", "iucn_complete_all", "iucn_complete_selected")) {
  effective_candidate_set <- switch(
    candidate_set,
    iucn_strict_range = "cleaned_unique",
    iucn_complete_all = "complete_all_predictors",
    iucn_complete_selected = "complete_selected_predictors"
  )
}

candidate_index <- switch(
  effective_candidate_set,
  cleaned_unique = rep(TRUE, nrow(occurrences)),
  complete_all_predictors = occurrences$complete_model_predictors,
  complete_selected_predictors = occurrences$complete_model_predictors,
  iucn_strict_range = occurrences$inside_iucn_range %in% TRUE,
  iucn_complete_all = occurrences$inside_iucn_range %in% TRUE & occurrences$complete_model_predictors,
  iucn_complete_selected = occurrences$inside_iucn_range %in% TRUE & occurrences$complete_model_predictors
)

model_occurrences <- occurrences[candidate_index %in% TRUE, , drop = FALSE]
if (nrow(model_occurrences) < min_obs) {
  stop(
    "Fewer than `--min-obs` records remain in candidate set `",
    effective_candidate_set,
    "` for ",
    species,
    ": ",
    nrow(model_occurrences),
    call. = FALSE
  )
}

model_occurrences$species <- species

# -----------------------------------------------------------------------------|
# Prepare deterministic output names ----
# -----------------------------------------------------------------------------|

model_grid_tag <- if (random_features) {
  paste0("random", n_models)
} else {
  paste0("feature_grid", "_beta", paste(beta_values, collapse = "-"))
}
run_config_tag <- paste(
  paste0("bk", n_background),
  model_grid_tag,
  paste0("select", n_selected_models),
  paste0("threads", maxent_threads),
  sep = "__"
)

run_prefix <- paste(
  species_safe,
  method,
  effective_candidate_set,
  predictor_mode,
  paste0(start_year, "_", end_year),
  run_config_tag,
  sep = "__"
)
used_occurrence_path <- file.path(output_dir, paste0(run_prefix, "__occurrences_used.csv"))
summary_path <- file.path(output_dir, paste0(run_prefix, "__run_summary.csv"))
model_path <- file.path(output_dir, paste0(run_prefix, "__model.rds"))

write.csv(model_occurrences, used_occurrence_path, row.names = FALSE, na = "")

# -----------------------------------------------------------------------------|
# Write run summary before optional fitting ----
# -----------------------------------------------------------------------------|

range_bbox <- summarise_bbox(species_range)
model_occurrence_bbox <- c(
  xmin = min(model_occurrences$decimalLongitude, na.rm = TRUE),
  ymin = min(model_occurrences$decimalLatitude, na.rm = TRUE),
  xmax = max(model_occurrences$decimalLongitude, na.rm = TRUE),
  ymax = max(model_occurrences$decimalLatitude, na.rm = TRUE)
)

summary <- data.frame(
  species_name = species,
  method = method,
  candidate_set = candidate_set,
  effective_candidate_set = effective_candidate_set,
  start_year = start_year,
  end_year = end_year,
  run_model = run_model,
  occurrence_path = occurrence_path,
  used_occurrence_path = used_occurrence_path,
  predictor_stack_path = predictor_stack_path,
  predictor_mode = predictor_mode,
  predictor_n_layers = terra::nlyr(predictors_for_model),
  predictor_layer_names = paste(names(predictors_for_model), collapse = "; "),
  select_var = if (identical(select_var, FALSE)) "FALSE" else select_var,
  iucn_range_path = iucn_range_path,
  range_filter = range_filter,
  range_mode = range_mode,
  range_status = range_status,
  range_buffer = range_buffer,
  input_unique_coordinates_after_date_filter = nrow(occurrences),
  iucn_range_coordinates = if (!is.null(species_range)) {
    sum(occurrences$inside_iucn_range %in% TRUE)
  } else {
    NA_integer_
  },
  complete_model_predictor_coordinates = sum(occurrences$complete_model_predictors),
  candidate_occurrence_rows = nrow(model_occurrences),
  candidate_complete_model_predictor_coordinates = sum(model_occurrences$complete_model_predictors),
  candidate_incomplete_model_predictor_coordinates = sum(!model_occurrences$complete_model_predictors),
  min_obs = min_obs,
  n_background = n_background,
  background_type = "BwData",
  test_percent = test_percent,
  beta_values = paste(beta_values, collapse = "; "),
  random_features = random_features,
  mod_select = TRUE,
  requested_n_models = n_models,
  candidate_model_grid_size = candidate_model_grid_size,
  n_selected_models = n_selected_models,
  seed = seed,
  maxent_threads = maxent_threads,
  model_occurrence_xmin = unname(model_occurrence_bbox["xmin"]),
  model_occurrence_ymin = unname(model_occurrence_bbox["ymin"]),
  model_occurrence_xmax = unname(model_occurrence_bbox["xmax"]),
  model_occurrence_ymax = unname(model_occurrence_bbox["ymax"]),
  iucn_xmin = unname(range_bbox["xmin"]),
  iucn_ymin = unname(range_bbox["ymin"]),
  iucn_xmax = unname(range_bbox["xmax"]),
  iucn_ymax = unname(range_bbox["ymax"]),
  existing_model_path = if (file.exists(existing_model_path)) existing_model_path else NA_character_,
  existing_model_n_presence = if (!is.null(existing_model)) unique(existing_model$params$n_presence)[[1]] else NA_integer_,
  existing_model_selected_variables = if (length(existing_selected_variables) > 0) {
    paste(existing_selected_variables, collapse = "; ")
  } else {
    NA_character_
  },
  output_model_path = if (run_model) model_path else NA_character_,
  run_status = if (run_model) "prepared" else "dry_run",
  error_message = NA_character_,
  run_started_at = run_started_at,
  run_finished_at = run_finished_at,
  run_elapsed_seconds = run_elapsed_seconds,
  run_elapsed_minutes = if (is.na(run_elapsed_seconds)) NA_real_ else run_elapsed_seconds / 60,
  prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  stringsAsFactors = FALSE
)
write.csv(summary, summary_path, row.names = FALSE, na = "")

cat("Prepared present SDM run for:", species, "\n")
cat("Candidate occurrences:", nrow(model_occurrences), "\n")
cat("Predictor mode:", predictor_mode, "\n")
cat("Range status:", range_status, "\n")
cat("Wrote occurrence input:", used_occurrence_path, "\n")
cat("Wrote run summary:", summary_path, "\n")

if (!run_model) {
  cat("Dry run only. Pass `--run` to fit and save the AutoMaxent model.\n")
} else {
  ##############################################################################
  # Fit AutoMaxent model and record elapsed time
  ##############################################################################

  load_automaxent(automaxent_root)

  # AutoMaxent calls install.packages() for missing dependencies. Override it so
  # the run fails clearly instead of trying to install packages mid-model.
  install.packages <- function(pkgs, ...) {
    pkgs <- pkgs[!is.na(pkgs) & nzchar(pkgs)]
    if (length(pkgs) == 0) {
      return(invisible(TRUE))
    }

    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
      stop(
        "Missing required R packages for AutoMaxent: ",
        paste(missing, collapse = ", "),
        ". Install them before running this script.",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }
  
  detectCores <- function(logical = TRUE) {
    maxent_threads + 2L
  }
  
  # Older AutoMaxent helper code can hand terra::rast() a list of Raster objects.
  # Convert that case explicitly while leaving normal terra calls unchanged.
  rast <- function(x, ...) {
    if (is.list(x) && length(x) > 0 && all(vapply(x, inherits, logical(1), "Raster"))) {
      return(terra::rast(raster::stack(x)))
    }
    
    terra::rast(x, ...)
  }
  
  run_start_time <- Sys.time()
  run_started_at <- format(run_start_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  
  model_result <- tryCatch(
    Auto_maxent(
      presence_dat = model_occurrences,
      predictors = predictors_for_model,
      coords.p = c("decimalLongitude", "decimalLatitude"),
      min_obs = min_obs,
      rm.dp = TRUE,
      name.mod = species,
      sp_range = species_range,
      crs.r = "EPSG:4326",
      buff_lim = if (is.null(species_range)) 0 else range_buffer,
      n_bk = n_background,
      type_bk = "BwData",
      Test_n = test_percent,
      time_macth = FALSE,
      select_var = select_var,
      random_features = random_features,
      seed.r = seed,
      beta.val = beta_values,
      n.m = n_models,
      Mod.route = model_work_dir,
      mod.select = TRUE,
      n.mods = n_selected_models,
      return.all = TRUE
    ),
    error = function(e) e
  )
  
  run_finish_time <- Sys.time()
  run_finished_at <- format(run_finish_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  run_elapsed_seconds <- as.numeric(difftime(run_finish_time, run_start_time, units = "secs"))
  
  if (inherits(model_result, "error")) {
    summary$run_status <- "failed"
    summary$error_message <- conditionMessage(model_result)
    summary$run_started_at <- run_started_at
    summary$run_finished_at <- run_finished_at
    summary$run_elapsed_seconds <- run_elapsed_seconds
    summary$run_elapsed_minutes <- run_elapsed_seconds / 60
    write.csv(summary, summary_path, row.names = FALSE, na = "")
    stop(conditionMessage(model_result), call. = FALSE)
  }
  
  model <- model_result
  saveRDS(model, model_path)
  
  summary$run_status <- "completed"
  summary$error_message <- NA_character_
  summary$run_started_at <- run_started_at
  summary$run_finished_at <- run_finished_at
  summary$run_elapsed_seconds <- run_elapsed_seconds
  summary$run_elapsed_minutes <- run_elapsed_seconds / 60
  summary$output_model_path <- model_path
  write.csv(summary, summary_path, row.names = FALSE, na = "")
  
  cat("Wrote model:", model_path, "\n")
  cat("Elapsed minutes:", round(run_elapsed_seconds / 60, 2), "\n")
}
