#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_run_present_model.R ----
# -----------------------------------------------------------------------------|
# Purpose: Prepare and optionally fit one present-day AutoMaxent SDM.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) stop("Package `here` is required.", call. = FALSE)
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package `terra` is required.", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package `sf` is required.", call. = FALSE)
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

sdm_config <- list(
  species = "Rousettus aegyptiacus",
  method = "spatial-spp",
  candidate_set = "iucn_complete_all",
  start_year = 2000,
  end_year = 2026,
  run_model = TRUE,
  range_filter = "auto",
  range_mode = "strict",
  range_buffer = 4,
  predictor_mode = "bio-elev",
  select_var = "NUMERICAL",
  min_obs = 20,
  n_background = "dynamic",
  test_percent = "dynamic",
  beta_values = "4,8,12",
  n_models = 25,
  n_selected_models = 10,
  use_boyce = 0.5,
  seed = 185,
  random_features = TRUE,
  maxent_threads = 2,
  java_memory_gb = 8,
  predictor_stack_path = file.path(here::here(), "sdms", "cache", "Resample_rast.tif"),
  iucn_range_path = file.path(here::here(), "sdms", "cache", "MAMMALS_TERRESTRIAL_ONLY", "MAMMALS_TERRESTRIAL_ONLY.shp"),
  output_root = file.path(here::here(), "sdms", "runs", "chikungunya", "calibration", "regenerated_models"),
  automaxent_root = "/Users/arturtrebski/Coding_Projects/AutoMaxent"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
cfg <- sdm_config
get_cfg <- function(key, field = gsub("-", "_", key)) get_arg(args, key, cfg[[field]])

cfg$species <- get_cfg("species")
cfg$method <- get_cfg("method")
cfg$candidate_set <- get_cfg("candidate-set")
cfg$start_year <- as.integer(get_cfg("start-year"))
cfg$end_year <- as.integer(get_cfg("end-year"))
cfg$run_model <- as_logical_arg(get_cfg("run-model")) || has_flag(args, "run")
if (has_flag(args, "dry-run")) cfg$run_model <- FALSE
cfg$range_filter <- get_cfg("range-filter")
cfg$range_mode <- get_cfg("range-mode")
cfg$range_buffer <- as.numeric(get_cfg("range-buffer"))
cfg$predictor_mode <- get_cfg("predictor-mode")
cfg$select_var <- get_cfg("select-var")
if (identical(cfg$select_var, "FALSE")) cfg$select_var <- FALSE
cfg$min_obs <- as.integer(get_cfg("min-obs"))
cfg$n_background <- get_cfg("n-background")
cfg$test_percent <- get_cfg("test-percent")
cfg$beta_values <- as.numeric(strsplit(get_cfg("beta-values"), ",", fixed = TRUE)[[1]])
cfg$n_models <- as.integer(get_cfg("n-models"))
cfg$n_selected_models <- as.integer(get_cfg("n-selected-models"))
cfg$use_boyce <- get_cfg("use-boyce")
if (cfg$use_boyce %in% c("", "NA", "NULL", "FALSE")) cfg$use_boyce <- NULL else cfg$use_boyce <- as.numeric(cfg$use_boyce)
cfg$seed <- as.integer(get_cfg("seed"))
cfg$random_features <- as_logical_arg(get_cfg("random-features")) || has_flag(args, "random-features")
if (has_flag(args, "feature-grid")) cfg$random_features <- FALSE
cfg$maxent_threads <- as.integer(get_cfg("threads", "maxent_threads"))
cfg$java_memory_gb <- as.numeric(get_cfg("java-memory-gb"))
cfg$predictor_stack_path <- get_arg(args, "predictor-stack", cfg$predictor_stack_path)
cfg$iucn_range_path <- get_arg(args, "iucn-range-path", cfg$iucn_range_path)
cfg$automaxent_root <- get_arg(args, "automaxent-root", cfg$automaxent_root)
cfg$output_root <- get_arg(args, "output-root", cfg$output_root)

repo <- repo_root()
species_safe <- safe_species_name(cfg$species)
occurrence_path <- get_arg(
  args,
  "occurrences",
  file.path(repo, "sdms", "runs", "chikungunya", "calibration", "occurrences", species_safe, cfg$method, "cleaned", paste0(species_safe, "_cleaned.csv"))
)
existing_model_path <- get_arg(args, "existing-model", file.path(repo, "sdms", "models", cfg$species, paste0(cfg$species, ".rds")))

if (!cfg$candidate_set %in% c("cleaned_unique", "complete_all_predictors", "complete_selected_predictors", "iucn_strict_range", "iucn_complete_all", "iucn_complete_selected")) {
  stop("Unsupported candidate_set: ", cfg$candidate_set, call. = FALSE)
}
if (!cfg$range_filter %in% c("auto", "apply", "none")) stop("range_filter must be auto, apply, or none.", call. = FALSE)
if (!cfg$range_mode %in% c("strict", "all")) stop("range_mode must be strict or all.", call. = FALSE)
if (!cfg$predictor_mode %in% c("bio-elev", "all", "existing-selected")) stop("predictor_mode must be bio-elev, all, or existing-selected.", call. = FALSE)
if (!is.null(cfg$use_boyce) && (is.na(cfg$use_boyce) || cfg$use_boyce < -1 || cfg$use_boyce > 1)) stop("use_boyce must be between -1 and 1.", call. = FALSE)
if (cfg$n_background != "dynamic" && is.na(suppressWarnings(as.integer(cfg$n_background)))) stop("n_background must be dynamic or an integer.", call. = FALSE)
if (!cfg$test_percent %in% c("dynamic", "auto") && is.na(suppressWarnings(as.integer(cfg$test_percent)))) stop("test_percent must be dynamic or an integer.", call. = FALSE)
if (!file.exists(occurrence_path)) stop("Missing cleaned occurrence file: ", occurrence_path, call. = FALSE)
if (!file.exists(cfg$predictor_stack_path)) stop("Missing predictor stack: ", cfg$predictor_stack_path, call. = FALSE)

load_automaxent <- function(root) {
  paths <- file.path(
    root,
    "Functions",
    c("auto_MaxEnt_complementary.R", "BackgroundPOINTS.R", "Environmental_weigthing_random_points.R", "Time_matchine.R", "auto_MaxEnt.R")
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) stop("Missing AutoMaxent function files: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(lapply(paths, source))
}

record_year <- function(x) {
  if ("year" %in% names(x)) return(suppressWarnings(as.integer(x$year)))
  if ("eventDate" %in% names(x)) return(as.integer(format(suppressWarnings(as.Date(substr(x$eventDate, 1, 10))), "%Y")))
  rep(NA_integer_, nrow(x))
}

dynamic_background <- function(setting, n_presence) {
  if (setting != "dynamic") return(as.integer(setting))
  if (n_presence > 10000) return(min(250000L, as.integer(1.5 * n_presence)))
  10000L
}

dynamic_test_percent <- function(setting, n_presence) {
  if (setting %in% c("dynamic", "auto")) {
    return(if (n_presence < 60) 20L else 30L)
  }

  as.integer(setting)
}

boyce_tag <- function(x) {
  if (is.null(x)) "boyceNA" else paste0("boyce", format(x, scientific = FALSE, trim = TRUE))
}

species_range <- function(path, species, mode, required) {
  if (!file.exists(path)) {
    if (required) stop("Missing IUCN range shapefile: ", path, call. = FALSE)
    return(NULL)
  }

  ranges <- sf::st_read(path, quiet = TRUE)
  if (!"sci_name" %in% names(ranges)) stop("IUCN range file is missing `sci_name`.", call. = FALSE)
  ranges <- ranges[ranges$sci_name == species, ]
  if (nrow(ranges) == 0) {
    if (required) stop("No IUCN range rows found for species: ", species, call. = FALSE)
    return(NULL)
  }

  if (mode == "strict" && all(c("presence", "origin", "seasonal") %in% names(ranges))) {
    strict <- ranges[ranges$presence == 1 & ranges$origin == 1 & ranges$seasonal == 1, ]
    if (nrow(strict) > 0) ranges <- strict
  }

  sf::st_as_sf(sf::st_union(sf::st_transform(sf::st_make_valid(ranges), 4326)))
}

pick_predictors <- function(predictors, mode, existing_model_path) {
  if (mode == "all") return(predictors)
  if (mode == "bio-elev") {
    keep <- grepl("wc2[.]1_2[.]5m_bio_", names(predictors)) | grepl("elev", names(predictors), ignore.case = TRUE)
    if (!any(keep)) stop("No bioclim/elevation layers found in predictor stack.", call. = FALSE)
    return(predictors[[keep]])
  }

  if (!file.exists(existing_model_path)) stop("existing-selected requires an existing model: ", existing_model_path, call. = FALSE)
  variables <- readRDS(existing_model_path)$variables
  missing <- setdiff(variables, names(predictors))
  if (length(missing) > 0) stop("Predictor stack is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  predictors[[variables]]
}

candidate_rows <- function(records, candidate_set, range) {
  if (is.null(range) && candidate_set %in% c("iucn_strict_range", "iucn_complete_all", "iucn_complete_selected")) {
    candidate_set <- switch(
      candidate_set,
      iucn_strict_range = "cleaned_unique",
      iucn_complete_all = "complete_all_predictors",
      iucn_complete_selected = "complete_selected_predictors"
    )
  }

  keep <- rep(TRUE, nrow(records))
  if (grepl("complete", candidate_set)) keep <- keep & records$complete_model_predictors
  if (grepl("iucn", candidate_set)) keep <- keep & records$inside_iucn_range %in% TRUE
  list(rows = records[keep, , drop = FALSE], effective_candidate_set = candidate_set)
}

occurrences <- read.csv(occurrence_path, check.names = FALSE, stringsAsFactors = FALSE)
missing_xy <- setdiff(c("decimalLongitude", "decimalLatitude"), names(occurrences))
if (length(missing_xy) > 0) stop("Occurrence file is missing: ", paste(missing_xy, collapse = ", "), call. = FALSE)

occurrences <- occurrences[!duplicated(paste(occurrences$decimalLongitude, occurrences$decimalLatitude, sep = "|")), , drop = FALSE]
occurrences$sdm_year <- record_year(occurrences)
occurrences <- occurrences[
  !is.na(occurrences$sdm_year) & occurrences$sdm_year >= cfg$start_year & occurrences$sdm_year <= cfg$end_year,
  ,
  drop = FALSE
]
if (nrow(occurrences) < cfg$min_obs) stop("Fewer than min_obs records remain after date filtering: ", nrow(occurrences), call. = FALSE)

predictors <- pick_predictors(terra::rast(cfg$predictor_stack_path), cfg$predictor_mode, existing_model_path)
points <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
occurrences$complete_model_predictors <- complete.cases(terra::extract(predictors, points, ID = FALSE))

range <- NULL
if (cfg$range_filter %in% c("auto", "apply")) {
  range <- species_range(cfg$iucn_range_path, cfg$species, cfg$range_mode, required = cfg$range_filter == "apply")
}
occurrences$inside_iucn_range <- if (is.null(range)) {
  NA
} else {
  sf_points <- sf::st_as_sf(occurrences, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
  lengths(sf::st_intersects(sf_points, range)) > 0
}

picked <- candidate_rows(occurrences, cfg$candidate_set, range)
model_occurrences <- picked$rows
if (nrow(model_occurrences) < cfg$min_obs) {
  stop("Fewer than min_obs records remain in candidate set `", picked$effective_candidate_set, "`: ", nrow(model_occurrences), call. = FALSE)
}
model_occurrences$species <- cfg$species

n_background <- dynamic_background(cfg$n_background, nrow(model_occurrences))
test_percent <- dynamic_test_percent(cfg$test_percent, nrow(model_occurrences))
model_grid_tag <- if (cfg$random_features) {
  paste0("random", cfg$n_models, "_beta", paste(cfg$beta_values, collapse = "-"))
} else {
  paste0("feature_grid_beta", paste(cfg$beta_values, collapse = "-"))
}
run_config_tag <- paste(
  paste0("bk", n_background),
  model_grid_tag,
  paste0("select", cfg$n_selected_models),
  boyce_tag(cfg$use_boyce),
  paste0("threads", cfg$maxent_threads),
  sep = "__"
)
run_prefix <- paste(
  species_safe,
  cfg$method,
  picked$effective_candidate_set,
  cfg$predictor_mode,
  paste0(cfg$start_year, "_", cfg$end_year),
  run_config_tag,
  sep = "__"
)

output_dir <- ensure_dir(file.path(cfg$output_root, species_safe))
model_work_dir <- ensure_dir(file.path(output_dir, "maxent_work"))
used_occurrence_path <- file.path(output_dir, paste0(run_prefix, "__occurrences_used.csv"))
summary_path <- file.path(output_dir, paste0(run_prefix, "__run_summary.csv"))
model_path <- file.path(output_dir, paste0(run_prefix, "__model.rds"))
write.csv(model_occurrences, used_occurrence_path, row.names = FALSE, na = "")

summary <- data.frame(
  species_name = cfg$species,
  method = cfg$method,
  candidate_set = cfg$candidate_set,
  effective_candidate_set = picked$effective_candidate_set,
  start_year = cfg$start_year,
  end_year = cfg$end_year,
  occurrence_path = occurrence_path,
  used_occurrence_path = used_occurrence_path,
  predictor_stack_path = cfg$predictor_stack_path,
  predictor_mode = cfg$predictor_mode,
  predictor_layer_names = paste(names(predictors), collapse = "; "),
  range_status = if (is.null(range)) "not_available" else "available",
  input_unique_coordinates_after_date_filter = nrow(occurrences),
  complete_model_predictor_coordinates = sum(occurrences$complete_model_predictors),
  candidate_occurrence_rows = nrow(model_occurrences),
  min_obs = cfg$min_obs,
  requested_n_background = cfg$n_background,
  n_background = n_background,
  background_type = "BwData",
  requested_test_percent = cfg$test_percent,
  test_percent = test_percent,
  beta_values = paste(cfg$beta_values, collapse = "; "),
  random_features = cfg$random_features,
  requested_n_models = cfg$n_models,
  n_selected_models = cfg$n_selected_models,
  use_boyce = if (is.null(cfg$use_boyce)) NA_real_ else cfg$use_boyce,
  seed = cfg$seed,
  maxent_threads = cfg$maxent_threads,
  java_memory_gb = cfg$java_memory_gb,
  output_model_path = if (cfg$run_model) model_path else NA_character_,
  run_status = if (cfg$run_model) "prepared" else "dry_run",
  error_message = NA_character_,
  run_started_at = NA_character_,
  run_finished_at = NA_character_,
  run_elapsed_seconds = NA_real_,
  run_elapsed_minutes = NA_real_,
  prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  stringsAsFactors = FALSE
)
write.csv(summary, summary_path, row.names = FALSE, na = "")

cat("Prepared present SDM run for:", cfg$species, "\n")
cat("Candidate occurrences:", nrow(model_occurrences), "\n")
cat("Wrote run summary:", summary_path, "\n")
if (!cfg$run_model) {
  cat("Dry run only. Pass `--run` to fit and save the AutoMaxent model.\n")
  quit(save = "no")
}

if (!is.na(cfg$java_memory_gb) && cfg$java_memory_gb > 0 && !"rJava" %in% loadedNamespaces()) {
  options(java.parameters = paste0("-Xmx", cfg$java_memory_gb, "g"))
}
load_automaxent(cfg$automaxent_root)

install.packages <- function(pkgs, ...) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
detectCores <- function(logical = TRUE) cfg$maxent_threads + 2L
rast <- function(x, ...) {
  if (is.list(x) && length(x) > 0 && all(vapply(x, inherits, logical(1), "Raster"))) return(terra::rast(raster::stack(x)))
  terra::rast(x, ...)
}

started_at <- Sys.time()
model_result <- tryCatch(
  Auto_maxent(
    presence_dat = model_occurrences,
    predictors = predictors,
    coords.p = c("decimalLongitude", "decimalLatitude"),
    min_obs = cfg$min_obs,
    rm.dp = TRUE,
    name.mod = cfg$species,
    sp_range = range,
    crs.r = "EPSG:4326",
    buff_lim = if (is.null(range)) 0 else cfg$range_buffer,
    n_bk = n_background,
    type_bk = "BwData",
    Test_n = test_percent,
    time_macth = FALSE,
    select_var = cfg$select_var,
    random_features = cfg$random_features,
    seed.r = cfg$seed,
    beta.val = cfg$beta_values,
    n.m = cfg$n_models,
    Mod.route = model_work_dir,
    mod.select = TRUE,
    n.mods = cfg$n_selected_models,
    use.boyce = cfg$use_boyce,
    return.all = TRUE
  ),
  error = function(e) e
)
finished_at <- Sys.time()
elapsed <- as.numeric(difftime(finished_at, started_at, units = "secs"))

summary$run_started_at <- format(started_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
summary$run_finished_at <- format(finished_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
summary$run_elapsed_seconds <- elapsed
summary$run_elapsed_minutes <- elapsed / 60

if (inherits(model_result, "error")) {
  summary$run_status <- "failed"
  summary$error_message <- conditionMessage(model_result)
  write.csv(summary, summary_path, row.names = FALSE, na = "")
  stop(conditionMessage(model_result), call. = FALSE)
}

saveRDS(model_result, model_path)
summary$run_status <- "completed"
summary$output_model_path <- model_path
write.csv(summary, summary_path, row.names = FALSE, na = "")
cat("Wrote model:", model_path, "\n")
cat("Elapsed minutes:", round(elapsed / 60, 2), "\n")
