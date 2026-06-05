#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 06_combine_vector_occurrences.R ----
# -----------------------------------------------------------------------------|
# Purpose: Combine GBIF, VectorMap, and MapVEu records into one cleaned
#          occurrence input per vector species.
#
# This script writes a new `combined` occurrence method and does not modify the
# source-specific GBIF, VectorMap, or MapVEu occurrence folders.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
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
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y")),
    species_filter = character(),
    max_species = Inf,
    min_points = 20,
    coordinate_round_digits = 5,
    dry_run = FALSE
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
  local_source_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "local_vector_occurrence_sources_manifest.csv"),
  occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  combined_run_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "combined_vector_occurrence_runs"),
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  min_points = 20,
  coordinate_round_digits = 5,
  sdm_pipeline_root = "/Users/arturtrebski/Coding_Projects/SDM_Pipeline",
  dry_run = FALSE
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Config helpers ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
}

target_manifest_path <- config_arg("target-manifest-path")
local_source_manifest_path <- config_arg("local-source-manifest-path")
occurrence_root <- config_arg("occurrence-root")
combined_run_root <- config_arg("combined-run-root")
roles <- split_arg(config_arg("roles"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
species_filter <- split_arg(config_arg("species-filter"))
max_species <- as.numeric(config_arg("max-species"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
min_points <- as.integer(config_arg("min-points"))
coordinate_round_digits <- as.integer(config_arg("coordinate-round-digits"))
sdm_pipeline_root <- config_arg("sdm-pipeline-root")
dry_run <- as_logical_arg(config_arg("dry-run")) || has_flag(args, "dry-run")

if (!file.exists(target_manifest_path)) {
  stop("Missing SDM target manifest: ", target_manifest_path, call. = FALSE)
}

if (!file.exists(local_source_manifest_path)) {
  stop("Missing local vector occurrence manifest: ", local_source_manifest_path, call. = FALSE)
}

cleaning_functions <- file.path(sdm_pipeline_root, "Functions", "CleanGBIF_Points.R")
if (!file.exists(cleaning_functions)) {
  stop("Missing SDM_Pipeline cleaning functions: ", cleaning_functions, call. = FALSE)
}
source(cleaning_functions)

# -----------------------------------------------------------------------------|
# Small helpers ----
# -----------------------------------------------------------------------------|

utc_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

update_if_present <- function(data, idx, col, value) {
  if (col %in% names(data)) {
    data[[col]][idx] <- value
  }
  data
}

first_existing_column <- function(data, candidates) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

column_or_na <- function(data, col) {
  if (is.na(col) || !col %in% names(data)) {
    return(rep(NA_character_, nrow(data)))
  }

  as.character(data[[col]])
}

numeric_column_or_na <- function(data, col) {
  if (is.na(col) || !col %in% names(data)) {
    return(rep(NA_real_, nrow(data)))
  }

  suppressWarnings(as.numeric(data[[col]]))
}

filter_gbif_present_records <- function(data) {
  rows_before <- nrow(data)
  keep <- rep(TRUE, nrow(data))
  occurrence_status_removed_rows <- 0L
  individual_count_zero_removed_rows <- 0L

  if ("occurrenceStatus" %in% names(data)) {
    status <- toupper(trimws(as.character(data$occurrenceStatus)))
    non_present <- !is.na(status) & nzchar(status) & status != "PRESENT"
    occurrence_status_removed_rows <- sum(non_present)
    keep <- keep & !non_present
  }

  if ("individualCount" %in% names(data)) {
    individual_count <- suppressWarnings(as.numeric(data$individualCount))
    zero_count <- !is.na(individual_count) & individual_count == 0
    individual_count_zero_removed_rows <- sum(zero_count & keep)
    keep <- keep & !zero_count
  }

  list(
    data = data[keep, , drop = FALSE],
    rows_before = rows_before,
    rows_after = sum(keep),
    occurrence_status_removed_rows = occurrence_status_removed_rows,
    individual_count_zero_removed_rows = individual_count_zero_removed_rows,
    rows_removed = rows_before - sum(keep)
  )
}

missing_year_value <- function(x) {
  value <- trimws(as.character(x))
  is.na(value) |
    !nzchar(value) |
    tolower(value) %in% c("na", "n/a", "no data", "null", "none", "unknown")
}

extract_plausible_years <- function(x) {
  value <- trimws(as.character(x))
  if (length(value) != 1 || missing_year_value(value)) {
    return(integer())
  }

  matches <- gregexpr("\\b(18|19|20)[0-9]{2}\\b", value, perl = TRUE)
  years <- regmatches(value, matches)[[1]]
  if (length(years) == 0 || identical(years, character(0))) {
    return(integer())
  }

  years <- suppressWarnings(as.integer(years))
  unique(years[!is.na(years) & years >= 1800 & years <= end_year])
}

parse_occurrence_year_value <- function(x) {
  if (missing_year_value(x)) {
    return(list(year = NA_integer_, status = "missing"))
  }

  years <- extract_plausible_years(x)
  if (length(years) == 0) {
    return(list(year = NA_integer_, status = "out_of_supported_range"))
  }

  if (length(years) > 1) {
    return(list(year = min(years), status = "ambiguous_multiple_years"))
  }

  list(year = years[[1]], status = "parsed")
}

parse_occurrence_years <- function(data, candidates) {
  out <- data.frame(
    year = rep(NA_integer_, nrow(data)),
    year_source_column = rep(NA_character_, nrow(data)),
    year_parse_status = rep("missing", nrow(data)),
    stringsAsFactors = FALSE
  )

  candidates <- candidates[candidates %in% names(data)]
  if (length(candidates) == 0 || nrow(data) == 0) {
    return(out)
  }

  for (col in candidates) {
    needs_year <- is.na(out$year)
    if (!any(needs_year)) {
      break
    }

    parsed <- lapply(data[[col]][needs_year], parse_occurrence_year_value)
    parsed_year <- vapply(parsed, `[[`, integer(1), "year")
    parsed_status <- vapply(parsed, `[[`, character(1), "status")
    target_rows <- which(needs_year)

    parsed_rows <- parsed_status %in% c("parsed", "ambiguous_multiple_years")
    if (any(parsed_rows)) {
      rows <- target_rows[parsed_rows]
      out$year[rows] <- parsed_year[parsed_rows]
      out$year_source_column[rows] <- col
      out$year_parse_status[rows] <- parsed_status[parsed_rows]
    }

    unresolved_rows <- target_rows[!parsed_rows & out$year_parse_status[target_rows] == "missing"]
    if (length(unresolved_rows) > 0) {
      unresolved_status <- parsed_status[!parsed_rows]
      out_of_range <- unresolved_status == "out_of_supported_range"
      if (any(out_of_range)) {
        rows <- unresolved_rows[out_of_range]
        out$year_source_column[rows] <- col
        out$year_parse_status[rows] <- "out_of_supported_range"
      }
    }
  }

  out
}

derive_year <- function(data, source_method) {
  candidates <- if (source_method == "vectormap") {
    c(
      "EarliestYearCollected",
      "LatestYearCollected",
      "EarliestDateCollected",
      "LatestDateCollected",
      "VerbatimCollectingDate",
      "eventDate"
    )
  } else if (source_method == "mapveu") {
    c(
      "year",
      "eventDate",
      "specimen collection start date [EUPATH_0043256]",
      "specimen collection end date [EUPATH_0043257]",
      "specimen collection date(s) (raw) [OBI_0001619]"
    )
  } else {
    c("year", "eventDate")
  }

  parse_occurrence_years(data, candidates)
}

derive_event_date <- function(data, source_method) {
  candidates <- if (source_method == "vectormap") {
    c("eventDate", "VerbatimCollectingDate", "EarliestDateCollected", "LatestDateCollected")
  } else if (source_method == "mapveu") {
    c(
      "eventDate",
      "specimen collection start date [EUPATH_0043256]",
      "specimen collection end date [EUPATH_0043257]",
      "specimen collection date(s) (raw) [OBI_0001619]"
    )
  } else {
    c("eventDate")
  }

  column_or_na(data, first_existing_column(data, candidates))
}

fallback_record_ids <- function(data, source_method, source_record_col) {
  record_ids <- column_or_na(data, source_record_col)
  missing <- is.na(record_ids) | !nzchar(record_ids)
  record_ids[missing] <- paste0(source_method, "_row_", which(missing))
  record_ids
}

collapse_compact <- function(x, max_values = 50) {
  x <- unique(as.character(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0) {
    return(NA_character_)
  }

  shown <- head(sort(x), max_values)
  suffix <- if (length(x) > max_values) {
    paste0(" ... (+", length(x) - max_values, " more)")
  } else {
    ""
  }
  paste0(paste(shown, collapse = "; "), suffix)
}

min_int_or_na <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_integer_)
  }

  min(x)
}

max_int_or_na <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_integer_)
  }

  max(x)
}

prepare_points_with_fallback <- function(points_sp, range_sp = NULL, xy_cols = c("decimalLongitude", "decimalLatitude")) {
  tryCatch(
    {
      list(
        points = Prepare_points(points_sp = points_sp, range_sp = range_sp, xy.c = xy_cols),
        cleaning_status = "sdm_pipeline_prepare_points"
      )
    },
    error = function(err) {
      message <- conditionMessage(err)
      known_sea_error <- grepl("unused argument .*full_url", message) ||
        grepl("ne_file_name.*full_url", message)

      if (!known_sea_error || !is.null(range_sp)) {
        stop(err)
      }

      if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
        stop("Package `CoordinateCleaner` is required for the fallback cleaner.", call. = FALSE)
      }

      warning(
        "Prepare_points() failed in the CoordinateCleaner sea test because of a ",
        "CoordinateCleaner/rnaturalearth API mismatch. Retrying without the `seas` test.",
        call. = FALSE
      )

      keep <- CoordinateCleaner::clean_coordinates(
        points_sp,
        lon = xy_cols[[1]],
        lat = xy_cols[[2]],
        tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "zeros")
      )$.summary

      list(
        points = points_sp[keep, , drop = FALSE],
        cleaning_status = "fallback_without_seas_due_rnaturalearth_api_mismatch"
      )
    }
  )
}

gbif_source_path <- function(species) {
  species_safe <- safe_species_name(species)
  candidates <- c(
    file.path(occurrence_root, species_safe, "gbif-download", "deduplicated", paste0(species_safe, "_gbif_key_deduplicated.csv")),
    file.path(occurrence_root, species_safe, "gbif-download", "cleaned", paste0(species_safe, "_cleaned.csv")),
    file.path(occurrence_root, species_safe, "gbif-download", "raw", paste0(species, ".csv"))
  )

  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    return(NA_character_)
  }

  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

resolve_local_source_output_path <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }

  if (file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  old_occurrence_root <- file.path(
    repo_root(),
    "sdms",
    "runs",
    "vector_sdm_push",
    "occurrences"
  )
  old_occurrence_root <- normalizePath(
    old_occurrence_root,
    winslash = "/",
    mustWork = FALSE
  )
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (!startsWith(normalized_path, paste0(old_occurrence_root, "/"))) {
    return(NA_character_)
  }

  relative_path <- substring(
    normalized_path,
    nchar(old_occurrence_root) + 2L
  )
  remapped_path <- file.path(occurrence_root, relative_path)

  if (!file.exists(remapped_path)) {
    return(NA_character_)
  }

  normalizePath(remapped_path, winslash = "/", mustWork = TRUE)
}

local_source_path <- function(local_manifest, species, source_method) {
  if (nrow(local_manifest) == 0) {
    return(NA_character_)
  }

  species_key <- canonical_species_name(species)
  local_species_key <- canonical_species_name(local_manifest$species_name)
  hits <- local_manifest[
    local_species_key == species_key &
      local_manifest$source_method == source_method &
      !is.na(local_manifest$output_path) &
      nzchar(local_manifest$output_path),
    ,
    drop = FALSE
  ]
  if (nrow(hits) == 0) {
    return(NA_character_)
  }

  resolved_paths <- vapply(
    hits$output_path,
    resolve_local_source_output_path,
    character(1)
  )
  resolved_paths <- resolved_paths[!is.na(resolved_paths) & nzchar(resolved_paths)]

  if (length(resolved_paths) == 0) {
    return(NA_character_)
  }

  resolved_paths[[1]]
}

standardize_source <- function(data, species, source_method, source_dataset, source_path) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }

  lon <- numeric_column_or_na(data, "decimalLongitude")
  lat <- numeric_column_or_na(data, "decimalLatitude")
  if (all(is.na(lon)) || all(is.na(lat))) {
    return(data.frame())
  }

  source_record_col <- switch(
    source_method,
    "gbif-download" = first_existing_column(data, c("gbifID", "key")),
    vectormap = first_existing_column(data, c("OBJECTID")),
    mapveu = first_existing_column(data, c("Sample_ID", "Collection_ID")),
    NA_character_
  )

  years <- derive_year(data, source_method)
  event_dates <- derive_event_date(data, source_method)
  source_record_ids <- fallback_record_ids(data, source_method, source_record_col)

  out <- data.frame(
    species_name = species,
    species = species,
    source_method = source_method,
    source_dataset = source_dataset,
    source_record_id = source_record_ids,
    decimalLongitude = lon,
    decimalLatitude = lat,
    year = years$year,
    year_source_column = years$year_source_column,
    year_parse_status = years$year_parse_status,
    eventDate = event_dates,
    source_path = source_path,
    source_row_number = seq_len(nrow(data)),
    stringsAsFactors = FALSE
  )

  out$source_record_column <- if (is.na(source_record_col)) NA_character_ else source_record_col
  out$coordinate_status <- ifelse(
    is.na(out$decimalLongitude) | is.na(out$decimalLatitude),
    "missing_coordinates",
    "has_coordinates"
  )
  out
}

update_target_manifest_from_combined <- function(target_manifest, summary_rows, target_manifest_path) {
  if (dry_run || nrow(summary_rows) == 0 || !file.exists(target_manifest_path)) {
    return(invisible(target_manifest))
  }

  target_manifest$species_name_canonical <- canonical_species_name(target_manifest$species_name)
  summary_rows$species_name_canonical <- canonical_species_name(summary_rows$species_name)

  for (i in seq_len(nrow(summary_rows))) {
    row <- summary_rows[i, , drop = FALSE]
    idx <- target_manifest$species_name_canonical == row$species_name_canonical[[1]]
    if (!any(idx)) {
      next
    }

    target_manifest <- update_if_present(target_manifest, idx, "occurrence_status", row$combine_status[[1]])
    target_manifest <- update_if_present(target_manifest, idx, "occurrence_rows_raw", row$raw_source_rows[[1]])
    target_manifest <- update_if_present(target_manifest, idx, "occurrence_rows_clean", row$cleaned_rows[[1]])
    target_manifest <- update_if_present(target_manifest, idx, "passes_min_points", row$passes_min_points[[1]])
    target_manifest <- update_if_present(target_manifest, idx, "combined_occurrence_path", row$cleaned_path[[1]])

    if ("run_status" %in% names(target_manifest) && row$combine_status[[1]] == "cleaned") {
      target_manifest$run_status[idx] <- if (isTRUE(row$passes_min_points[[1]])) {
        "occurrences_ready"
      } else {
        "occurrences_below_min_points"
      }
    }
  }

  write.csv(target_manifest, target_manifest_path, row.names = FALSE, na = "")
  invisible(target_manifest)
}

write_empty_summary <- function(path) {
  data.table::fwrite(data.frame(), path, na = "")
}

process_species <- function(species, local_manifest) {
  species_safe <- safe_species_name(species)
  combined_dir <- file.path(occurrence_root, species_safe, "combined")
  standardized_dir <- file.path(combined_dir, "standardized")
  deduplicated_dir <- file.path(combined_dir, "deduplicated")
  cleaned_dir <- file.path(combined_dir, "cleaned")
  standardized_path <- file.path(standardized_dir, paste0(species_safe, "_combined_standardized.csv"))
  deduplicated_path <- file.path(deduplicated_dir, paste0(species_safe, "_combined_coordinate_deduplicated.csv"))
  cleaned_path <- file.path(cleaned_dir, paste0(species_safe, "_cleaned.csv"))
  summary_path <- file.path(combined_dir, "occurrence_preparation_summary.csv")

  gbif_path <- gbif_source_path(species)
  vectormap_path <- local_source_path(local_manifest, species, "vectormap")
  mapveu_path <- local_source_path(local_manifest, species, "mapveu")

  sources <- list(
    list(method = "gbif-download", dataset = "GBIF occurrence download", path = gbif_path),
    list(method = "vectormap", dataset = "VectorMap MosquitoMap2", path = vectormap_path),
    list(method = "mapveu", dataset = "MapVEu sample + collection site + collection metadata", path = mapveu_path)
  )

  gbif_presence_stats <- data.frame(
    gbif_download_source_rows = 0L,
    gbif_download_rows_after_presence_filter = 0L,
    gbif_occurrence_status_removed_rows = 0L,
    gbif_individual_count_zero_removed_rows = 0L,
    gbif_presence_filter_removed_rows = 0L,
    stringsAsFactors = FALSE
  )

  standardized_sources <- lapply(sources, function(source) {
    if (is.na(source$path) || !file.exists(source$path)) {
      return(data.frame())
    }
    data <- read_csv_if_exists(source$path)
    if (source$method == "gbif-download") {
      gbif_filter <- filter_gbif_present_records(data)
      data <- gbif_filter$data
      gbif_presence_stats$gbif_download_source_rows <<- gbif_filter$rows_before
      gbif_presence_stats$gbif_download_rows_after_presence_filter <<- gbif_filter$rows_after
      gbif_presence_stats$gbif_occurrence_status_removed_rows <<- gbif_filter$occurrence_status_removed_rows
      gbif_presence_stats$gbif_individual_count_zero_removed_rows <<- gbif_filter$individual_count_zero_removed_rows
      gbif_presence_stats$gbif_presence_filter_removed_rows <<- gbif_filter$rows_removed
    }
    standardize_source(data, species, source$method, source$dataset, source$path)
  })
  standardized <- data.table::rbindlist(standardized_sources, use.names = TRUE, fill = TRUE, idcol = FALSE)

  raw_source_rows <- nrow(standardized)
  rows_by_source <- table(factor(standardized$source_method, levels = c("gbif-download", "vectormap", "mapveu")))

  if (raw_source_rows == 0) {
    summary <- data.frame(
      species_name = species,
      requested_start_year = start_year,
      requested_end_year = end_year,
      min_points = min_points,
      combine_status = "no_source_records",
      source_scope = "none",
      gbif_source_path = gbif_path,
      vectormap_source_path = vectormap_path,
      mapveu_source_path = mapveu_path,
      gbif_download_source_rows = gbif_presence_stats$gbif_download_source_rows,
      gbif_download_rows_after_presence_filter = gbif_presence_stats$gbif_download_rows_after_presence_filter,
      gbif_occurrence_status_removed_rows = gbif_presence_stats$gbif_occurrence_status_removed_rows,
      gbif_individual_count_zero_removed_rows = gbif_presence_stats$gbif_individual_count_zero_removed_rows,
      gbif_presence_filter_removed_rows = gbif_presence_stats$gbif_presence_filter_removed_rows,
      raw_source_rows = 0,
      standardized_rows = 0,
      year_missing_rows = 0,
      year_parsed_rows = 0,
      year_ambiguous_rows = 0,
      year_out_of_supported_range_rows = 0,
      year_filtered_rows = 0,
      within_source_deduplicated_rows = 0,
      coordinate_deduplicated_rows = 0,
      cleaned_rows = 0,
      cleaned_unique_coordinate_rows = 0,
      passes_min_points = FALSE,
      cleaning_status = NA_character_,
      standardized_path = NA_character_,
      deduplicated_path = NA_character_,
      cleaned_path = NA_character_,
      prepared_at = utc_now(),
      stringsAsFactors = FALSE
    )
    if (!dry_run) {
      ensure_dir(combined_dir)
      data.table::fwrite(summary, summary_path, na = "")
    }
    return(summary)
  }

  standardized <- standardized[standardized$coordinate_status == "has_coordinates", , drop = FALSE]
  year_missing_rows <- sum(is.na(standardized$year))
  year_parsed_rows <- sum(standardized$year_parse_status == "parsed", na.rm = TRUE)
  year_ambiguous_rows <- sum(standardized$year_parse_status == "ambiguous_multiple_years", na.rm = TRUE)
  year_out_of_supported_range_rows <- sum(standardized$year_parse_status == "out_of_supported_range", na.rm = TRUE)
  year_known <- !is.na(standardized$year)
  year_keep <- !year_known | (standardized$year >= start_year & standardized$year <= end_year)
  year_filtered_rows <- sum(!year_keep)
  standardized <- standardized[year_keep, , drop = FALSE]
  standardized_rows <- nrow(standardized)

  if (standardized_rows == 0) {
    combine_status <- "no_records_after_year_filter"
    within_deduplicated <- standardized
    coordinate_deduplicated <- standardized
    cleaned <- standardized
    cleaning_status <- NA_character_
  } else {
    within_key <- paste(standardized$source_method, standardized$source_record_id, sep = "|")
    within_deduplicated <- standardized[!duplicated(within_key), , drop = FALSE]

    within_deduplicated$rounded_longitude <- round(within_deduplicated$decimalLongitude, coordinate_round_digits)
    within_deduplicated$rounded_latitude <- round(within_deduplicated$decimalLatitude, coordinate_round_digits)
    within_deduplicated$coordinate_group_key <- paste(
      within_deduplicated$species_name,
      within_deduplicated$rounded_longitude,
      within_deduplicated$rounded_latitude,
      sep = "|"
    )
    within_deduplicated$source_priority <- match(
      within_deduplicated$source_method,
      c("gbif-download", "vectormap", "mapveu")
    )

    dt <- data.table::as.data.table(within_deduplicated)
    group_summary <- dt[
      ,
      .(
        source_methods = collapse_compact(source_method),
        source_record_ids = collapse_compact(source_record_id),
        source_row_count = .N,
        gbif_download_row_count = sum(source_method == "gbif-download"),
        vectormap_row_count = sum(source_method == "vectormap"),
        mapveu_row_count = sum(source_method == "mapveu"),
        year_min = min_int_or_na(year),
        year_max = max_int_or_na(year)
      ),
      by = coordinate_group_key
    ]

    data.table::setorder(dt, coordinate_group_key, source_priority, source_method, source_record_id)
    representatives <- dt[!duplicated(coordinate_group_key)]
    coordinate_deduplicated <- merge(
      as.data.frame(representatives),
      as.data.frame(group_summary),
      by = "coordinate_group_key",
      all.x = TRUE,
      sort = FALSE
    )
    coordinate_deduplicated$source_priority <- NULL

    cleaning_result <- prepare_points_with_fallback(
      points_sp = coordinate_deduplicated,
      range_sp = NULL,
      xy_cols = c("decimalLongitude", "decimalLatitude")
    )
    cleaned <- cleaning_result$points
    cleaning_status <- cleaning_result$cleaning_status
    combine_status <- if (is.null(cleaned) || nrow(cleaned) == 0) {
      "cleaning_failed_no_usable_records"
    } else {
      cleaned_unique <- nrow(unique(cleaned[, c("decimalLongitude", "decimalLatitude"), drop = FALSE]))
      if (cleaned_unique >= min_points) "cleaned" else "cleaned_below_min_points"
    }
  }

  if (nrow(cleaned) > 0) {
    cleaned_unique_coordinate_rows <- nrow(unique(cleaned[, c("decimalLongitude", "decimalLatitude"), drop = FALSE]))
    passes_min_points <- cleaned_unique_coordinate_rows >= min_points
  } else {
    cleaned_unique_coordinate_rows <- 0L
    passes_min_points <- FALSE
  }

  has_gbif <- sum(standardized$source_method == "gbif-download", na.rm = TRUE) > 0
  has_local <- any(standardized$source_method %in% c("vectormap", "mapveu"))
  source_scope <- if (has_gbif && has_local) {
    "gbif_and_local"
  } else if (has_gbif) {
    "gbif_only"
  } else if (has_local) {
    "local_only"
  } else {
    "none"
  }

  if (!dry_run) {
    ensure_dir(standardized_dir)
    ensure_dir(deduplicated_dir)
    ensure_dir(cleaned_dir)
    data.table::fwrite(standardized, standardized_path, na = "")
    data.table::fwrite(coordinate_deduplicated, deduplicated_path, na = "")
    if (nrow(cleaned) > 0) {
      data.table::fwrite(cleaned, cleaned_path, na = "")
    }
  }

  summary <- data.frame(
    species_name = species,
    requested_start_year = start_year,
    requested_end_year = end_year,
    min_points = min_points,
    coordinate_round_digits = coordinate_round_digits,
    combine_status = combine_status,
    source_scope = source_scope,
    gbif_source_path = gbif_path,
    vectormap_source_path = vectormap_path,
    mapveu_source_path = mapveu_path,
    gbif_download_source_rows = gbif_presence_stats$gbif_download_source_rows,
    gbif_download_rows_after_presence_filter = gbif_presence_stats$gbif_download_rows_after_presence_filter,
    gbif_occurrence_status_removed_rows = gbif_presence_stats$gbif_occurrence_status_removed_rows,
    gbif_individual_count_zero_removed_rows = gbif_presence_stats$gbif_individual_count_zero_removed_rows,
    gbif_presence_filter_removed_rows = gbif_presence_stats$gbif_presence_filter_removed_rows,
    gbif_download_raw_rows = as.integer(rows_by_source[["gbif-download"]]),
    vectormap_raw_rows = as.integer(rows_by_source[["vectormap"]]),
    mapveu_raw_rows = as.integer(rows_by_source[["mapveu"]]),
    raw_source_rows = raw_source_rows,
    standardized_rows = standardized_rows,
    year_missing_rows = year_missing_rows,
    year_parsed_rows = year_parsed_rows,
    year_ambiguous_rows = year_ambiguous_rows,
    year_out_of_supported_range_rows = year_out_of_supported_range_rows,
    year_filtered_rows = year_filtered_rows,
    within_source_deduplicated_rows = nrow(within_deduplicated),
    coordinate_deduplicated_rows = nrow(coordinate_deduplicated),
    cleaned_rows = nrow(cleaned),
    cleaned_unique_coordinate_rows = cleaned_unique_coordinate_rows,
    passes_min_points = passes_min_points,
    cleaning_status = cleaning_status,
    standardized_path = if (!dry_run && file.exists(standardized_path)) standardized_path else NA_character_,
    deduplicated_path = if (!dry_run && file.exists(deduplicated_path)) deduplicated_path else NA_character_,
    cleaned_path = if (!dry_run && file.exists(cleaned_path)) cleaned_path else NA_character_,
    prepared_at = utc_now(),
    stringsAsFactors = FALSE
  )

  if (!dry_run) {
    data.table::fwrite(summary, summary_path, na = "")
  }
  summary
}

# -----------------------------------------------------------------------------|
# Select targets and run combination ----
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

local_manifest <- read.csv(local_source_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
run_dir <- file.path(
  combined_run_root,
  timestamp
)
run_summary_path <- file.path(run_dir, "combined_vector_occurrence_summary.csv")

cat("Selected target species:", nrow(targets), "\n")
if (dry_run) {
  cat("Dry run: combined occurrence files and summaries will not be written.\n")
}

if (nrow(targets) == 0) {
  warning(
    "No target species selected. Check roles/species_filter/include_not_needed/include_already_available settings.",
    call. = FALSE
  )
  summary <- data.frame()
  if (!dry_run) {
    ensure_dir(run_dir)
    data.table::fwrite(summary, run_summary_path, na = "")
    cat("Wrote combined occurrence run summary:", run_summary_path, "\n")
  }
} else {
  rows <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    species <- targets$species_name_canonical[[i]]
    result <- tryCatch(
      process_species(species, local_manifest),
      error = function(err) {
        data.frame(
          species_name = species,
          requested_start_year = start_year,
          requested_end_year = end_year,
          min_points = min_points,
          coordinate_round_digits = coordinate_round_digits,
          combine_status = "failed_error",
          source_scope = NA_character_,
          gbif_download_source_rows = NA_integer_,
          gbif_download_rows_after_presence_filter = NA_integer_,
          gbif_occurrence_status_removed_rows = NA_integer_,
          gbif_individual_count_zero_removed_rows = NA_integer_,
          gbif_presence_filter_removed_rows = NA_integer_,
          raw_source_rows = NA_integer_,
          standardized_rows = NA_integer_,
          year_missing_rows = NA_integer_,
          year_parsed_rows = NA_integer_,
          year_ambiguous_rows = NA_integer_,
          year_out_of_supported_range_rows = NA_integer_,
          year_filtered_rows = NA_integer_,
          within_source_deduplicated_rows = NA_integer_,
          coordinate_deduplicated_rows = NA_integer_,
          cleaned_rows = NA_integer_,
          cleaned_unique_coordinate_rows = NA_integer_,
          passes_min_points = FALSE,
          cleaning_status = NA_character_,
          standardized_path = NA_character_,
          deduplicated_path = NA_character_,
          cleaned_path = NA_character_,
          prepared_at = utc_now(),
          notes = conditionMessage(err),
          stringsAsFactors = FALSE
        )
      }
    )

    rows[[i]] <- result
    if (!dry_run) {
      ensure_dir(run_dir)
      data.table::fwrite(data.table::rbindlist(rows[seq_len(i)], use.names = TRUE, fill = TRUE), run_summary_path, na = "")
    }

    cat(
      "[",
      i,
      "/",
      nrow(targets),
      "] ",
      species,
      ": ",
      result$combine_status[[1]],
      ", cleaned unique coords = ",
      coalesce_scalar(result$cleaned_unique_coordinate_rows, default = "NA"),
      "\n",
      sep = ""
    )
  }

  summary <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!dry_run) {
    ensure_dir(run_dir)
    data.table::fwrite(summary, run_summary_path, na = "")
    target_manifest <- update_target_manifest_from_combined(target_manifest, summary, target_manifest_path)
  }

  if (!dry_run) {
    cat("Wrote combined occurrence run summary:", run_summary_path, "\n")
  }
}
