#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_prepare_one_gbif_species.R ----
# -----------------------------------------------------------------------------|
# Purpose: Download and clean occurrence records for host/vector SDM calibration
#          using the available SDM_Pipeline functions.
#
# This mostly reuses the existing pipeline wrappers:
#   - Spatial_Information_sp.R: Download_gbif() / Spatial_spp()
#   - CleanGBIF_Points.R: Prepare_points()
# For bulk GBIF vector occurrence downloads, `gbif-download` uses rgbif's
# asynchronous occurrence download API instead of month-by-month occ_search().
#
# Default target:
#   sdms/runs/chikungunya/calibration/host_regeneration_manifest.csv
# -----------------------------------------------------------------------------|
suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})
source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# Arguments and method selection ----
# -----------------------------------------------------------------------------|

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

species <- get_arg(args, "species", "Rousettus aegyptiacus")
method <- get_arg(args, "method", "spatial-spp")
start_year <- as.integer(get_arg(args, "start-year", "1970"))
end_year_arg <- get_arg(args, "end-year", as.character(as.integer(format(Sys.Date(), "%Y"))))
end_year <- if (is.na(end_year_arg) || !nzchar(end_year_arg)) NA_integer_ else as.integer(end_year_arg)
time_limit_arg <- get_arg(args, "time-limit", NA_character_)
time_limit <- if (is.na(time_limit_arg) || !nzchar(time_limit_arg)) NA_real_ else as.numeric(time_limit_arg)
min_points <- as.integer(get_arg(args, "min-points", "20"))
gbif_download_key_arg <- get_arg(args, "gbif-download-key", NA_character_)
redownload <- has_flag(args, "redownload")
update_manifest <- !has_flag(args, "no-manifest-update")

if (!method %in% c("direct-gbif", "spatial-spp", "gbif-download")) {
  stop("`--method` must be one of: direct-gbif, spatial-spp, gbif-download", call. = FALSE)
}

manifest_path <- get_arg(
  args,
  "manifest",
  file.path(repo_root(), "sdms", "runs", "chikungunya", "calibration", "host_regeneration_manifest.csv")
)
occurrence_root <- get_arg(
  args,
  "occurrence-root",
  file.path(repo_root(), "sdms", "runs", "chikungunya", "calibration", "occurrences")
)

# -----------------------------------------------------------------------------|
# External SDM_Pipeline functions ----
# -----------------------------------------------------------------------------|

pipeline_root <- get_arg(args, "sdm-pipeline-root", "/Users/arturtrebski/Coding_Projects/SDM_Pipeline")
spatial_functions <- file.path(pipeline_root, "Functions", "Spatial_Information_sp.R")
cleaning_functions <- file.path(pipeline_root, "Functions", "CleanGBIF_Points.R")

if (!file.exists(spatial_functions)) {
  stop("Missing SDM_Pipeline spatial functions: ", spatial_functions, call. = FALSE)
}

if (!file.exists(cleaning_functions)) {
  stop("Missing SDM_Pipeline cleaning functions: ", cleaning_functions, call. = FALSE)
}

if (update_manifest && !file.exists(manifest_path)) {
  stop("Missing calibration manifest: ", manifest_path, call. = FALSE)
}

source(spatial_functions)
source(cleaning_functions)

# -----------------------------------------------------------------------------|
# Shared small helpers ----
# -----------------------------------------------------------------------------|

canonical_species_key <- function(x) {
  x <- trimws(tolower(as.character(x)))
  gsub("[[:space:]]+", " ", x)
}

update_if_present <- function(data, idx, col, value) {
  if (col %in% names(data)) {
    data[[col]][idx] <- value
  }
  data
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

      # The fallback is deliberately narrow: it only bypasses the known
      # CoordinateCleaner/rnaturalearth sea-test API mismatch.
      warning(
        "Prepare_points() failed in the CoordinateCleaner sea test because of a ",
        "CoordinateCleaner/rnaturalearth API mismatch. Retrying without the `seas` test.",
        call. = FALSE
      )

      if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
        stop("Package `CoordinateCleaner` is required for the fallback cleaner.", call. = FALSE)
      }

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

# -----------------------------------------------------------------------------|
# GBIF bulk-download helpers ----
# -----------------------------------------------------------------------------|

filter_occurrence_years <- function(data, start_year, end_year) {
  if (!"year" %in% names(data)) {
    return(list(data = data, status = "not_applied_no_year_column", rows_before = nrow(data)))
  }

  years <- suppressWarnings(as.integer(data$year))
  keep <- is.na(years) | years >= start_year
  if (!is.na(end_year)) {
    keep <- keep & (is.na(years) | years <= end_year)
  }

  list(
    data = data[keep, , drop = FALSE],
    status = "applied",
    rows_before = nrow(data)
  )
}

filter_gbif_present_records <- function(data) {
  rows_before <- nrow(data)
  occurrence_status_removed_rows <- 0L
  individual_count_zero_removed_rows <- 0L

  keep <- rep(TRUE, nrow(data))
  status <- rep(NA_character_, nrow(data))
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
    occurrence_status_removed_rows = occurrence_status_removed_rows,
    individual_count_zero_removed_rows = individual_count_zero_removed_rows,
    rows_removed = rows_before - sum(keep)
  )
}

read_dotenv_value <- function(keys, env_path = file.path(repo_root(), ".env")) {
  # Allows ignored local credentials to be used without requiring a shell export.
  if (!file.exists(env_path)) {
    return("")
  }

  lines <- readLines(env_path, warn = FALSE)
  for (key in keys) {
    pattern <- paste0("^[[:space:]]*", key, "[[:space:]]*=")
    hit <- grep(pattern, lines, value = TRUE)
    if (length(hit) == 0) {
      next
    }

    value <- sub(pattern, "", hit[[1]])
    value <- sub("[[:space:]]+#.*$", "", value)
    value <- trimws(value)
    value <- sub("^\"", "", value)
    value <- sub("\"$", "", value)
    value <- sub("^'", "", value)
    value <- sub("'$", "", value)
    if (nzchar(value)) {
      return(value)
    }
  }

  ""
}

gbif_credential <- function(keys) {
  for (key in keys) {
    value <- Sys.getenv(key, unset = "")
    if (nzchar(value)) {
      return(value)
    }
  }

  read_dotenv_value(keys)
}

with_gbif_retries <- function(label, expr, attempts = 10, sleep_seconds = 30) {
  # Retry only transient GBIF status/get calls; the download request itself is
  # still submitted once per script run unless a key is explicitly reused.
  last_error <- NULL
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(expr(), error = function(err) {
      last_error <<- err
      NULL
    })
    if (!is.null(result)) {
      return(result)
    }

    if (attempt < attempts) {
      message(label, " failed on attempt ", attempt, "/", attempts, ": ", conditionMessage(last_error))
      Sys.sleep(sleep_seconds)
    }
  }

  stop(last_error)
}

download_gbif_bulk <- function(species, start_year, end_year, raw_path, raw_dir, existing_download_key = NA_character_) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package `rgbif` is required for `--method gbif-download`.", call. = FALSE)
  }

  gbif_user <- gbif_credential(c("GBIF_USER", "GBIF_USERNAME", "gbif_user", "gbif_username"))
  gbif_password <- gbif_credential(c("GBIF_PASSWORD", "GBIF_PWD", "gbif_password", "gbif_pwd"))
  gbif_email <- gbif_credential(c("GBIF_EMAIL", "gbif_email"))
  if (!nzchar(gbif_email) && grepl("@", gbif_user, fixed = TRUE)) {
    gbif_email <- gbif_user
  }

  if (!nzchar(gbif_user) || !nzchar(gbif_password) || !nzchar(gbif_email)) {
    stop(
      "`--method gbif-download` requires GBIF_USER/GBIF_USERNAME/gbif_username, ",
      "GBIF_PASSWORD/gbif_password, and GBIF_EMAIL/gbif_email. ",
      "If gbif_username is an email address it will be reused as the email.",
      call. = FALSE
    )
  }

  if (!is.na(existing_download_key) && nzchar(existing_download_key)) {
    download_key <- existing_download_key
    message("Reusing GBIF occurrence download key for ", species, ": ", download_key)
  } else {
    backbone <- rgbif::name_backbone(name = species, rank = "species")
    taxon_key <- backbone$usageKey
    if (is.null(taxon_key) || is.na(taxon_key)) {
      stop("GBIF backbone lookup did not return a usageKey for ", species, call. = FALSE)
    }

    predicates <- list(
      rgbif::pred("taxonKey", taxon_key),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      rgbif::pred_gte("year", start_year)
    )
    if (!is.na(end_year)) {
      predicates <- c(predicates, list(rgbif::pred_lte("year", end_year)))
    }

    message("Submitting GBIF occurrence download for ", species, " (taxonKey ", taxon_key, ")")
    download_key <- do.call(
      rgbif::occ_download,
      c(
        predicates,
        list(
          format = "SIMPLE_CSV",
          user = gbif_user,
          pwd = gbif_password,
          email = gbif_email
        )
      )
    )
  }

  message("Waiting for GBIF download: ", download_key)
  with_gbif_retries(
    label = paste0("GBIF download wait for ", download_key),
    expr = function() rgbif::occ_download_wait(download_key, status_ping = 20, quiet = FALSE),
    attempts = 12,
    sleep_seconds = 30
  )

  download_archive <- with_gbif_retries(
    label = paste0("GBIF download get for ", download_key),
    expr = function() {
      rgbif::occ_download_get(
        download_key,
        path = raw_dir,
        overwrite = TRUE
      )
    },
    attempts = 6,
    sleep_seconds = 30
  )
  raw <- rgbif::occ_download_import(download_archive)
  raw <- as.data.frame(raw)

  if (!all(c("decimalLongitude", "decimalLatitude") %in% names(raw))) {
    stop(
      "GBIF download did not include decimalLongitude/decimalLatitude columns for ",
      species,
      call. = FALSE
    )
  }

  write.csv(raw, raw_path, row.names = FALSE, na = "")
  download_key
}

gbif_backbone_names <- function(species) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    return(species)
  }

  backbone <- tryCatch(rgbif::name_backbone(name = species, rank = "species"), error = function(err) NULL)
  if (is.null(backbone) || is.null(backbone$usageKey) || is.na(backbone$usageKey)) {
    return(species)
  }

  synonyms <- tryCatch(
    rgbif::name_usage(key = backbone$usageKey, data = "synonyms"),
    error = function(err) data.frame()
  )

  names <- unique(c(
    species,
    backbone$canonicalName,
    backbone$species,
    synonyms$canonicalName,
    synonyms$species
  ))
  names <- names[!is.na(names) & nzchar(names)]
  names <- names[!grepl(" ", names) | lengths(strsplit(names, " ")) >= 2]
  unique(names)
}

# -----------------------------------------------------------------------------|
# Spatial_spp synonym fallback ----
# -----------------------------------------------------------------------------|

pipeline_iucn_names_without_itis <- function(species, iucn_key) {
  taxonomy <- tryCatch(
    retrieve_syns(
      spp_name = species,
      IUCN_api = iucn_key,
      Gbif = FALSE,
      ITIS = FALSE
    ),
    error = function(err) NULL
  )

  if (is.null(taxonomy)) {
    return(list(names = species, taxonomy = data.frame(Or_name = species)))
  }

  list(names = taxonomy$Spp_syn, taxonomy = taxonomy$TaxDat)
}

download_spatial_spp_with_fallback <- function(species, raw_dir, taxonomy_dir, start_year, end_year, iucn_key, time_limit) {
  tryCatch(
    Spatial_spp(
      sci_sp = species,
      p.route = raw_dir,
      t.route = taxonomy_dir,
      start_date = start_year,
      end_date = if (is.na(end_year)) NULL else end_year,
      IUCN_api = iucn_key,
      time_limit = if (is.na(time_limit)) NULL else time_limit
    ),
    error = function(err) {
      message <- conditionMessage(err)
      if (!grepl("ITIS Search is likely down", message, fixed = TRUE)) {
        stop(err)
      }

      warning(
        "Spatial_spp() failed during the ITIS synonym lookup. Retrying with ITIS disabled ",
        "and GBIF/IUCN synonym matching enabled.",
        call. = FALSE
      )

      # Preserve the pipeline-style synonym expansion when ITIS is unavailable,
      # then pass the expanded name list through Download_gbif().
      iucn_taxonomy <- pipeline_iucn_names_without_itis(species, iucn_key)
      gbif_names <- gbif_backbone_names(species)
      synonym_names <- unique(c(species, iucn_taxonomy$names, gbif_names))
      synonym_names <- synonym_names[!is.na(synonym_names) & nzchar(synonym_names)]

      taxonomy <- iucn_taxonomy$taxonomy
      taxonomy$GBIF_backbone_names_used <- paste(gbif_names, collapse = ";")
      taxonomy$Synonym_names_used <- paste(synonym_names, collapse = ";")

      ensure_dir(taxonomy_dir)
      write.csv(taxonomy, file.path(taxonomy_dir, paste0(species, ".csv")), row.names = FALSE)

      Download_gbif(
        sp_list = synonym_names,
        initial_date = start_year,
        end_date = if (is.na(end_year)) NULL else end_year,
        exit_route = raw_dir,
        n_records = 150000,
        time_limit = if (is.na(time_limit)) NULL else time_limit
      )

      taxonomy
    }
  )
}

# -----------------------------------------------------------------------------|
# Occurrence output paths ----
# -----------------------------------------------------------------------------|

species_safe <- safe_species_name(species)
species_occurrence_dir <- ensure_dir(file.path(occurrence_root, species_safe))
occurrence_dir <- ensure_dir(file.path(species_occurrence_dir, method))
raw_dir <- ensure_dir(file.path(occurrence_dir, "raw"))
deduplicated_dir <- ensure_dir(file.path(occurrence_dir, "deduplicated"))
clean_dir <- ensure_dir(file.path(occurrence_dir, "cleaned"))

raw_path <- file.path(raw_dir, paste0(species, ".csv"))
raw_manifest_path <- file.path(raw_dir, "raw_download_manifest.csv")
deduplicated_path <- file.path(deduplicated_dir, paste0(species_safe, "_gbif_key_deduplicated.csv"))
cleaned_path <- file.path(clean_dir, paste0(species_safe, "_cleaned.csv"))
summary_path <- file.path(occurrence_dir, "occurrence_preparation_summary.csv")
raw_download_status <- if (!file.exists(raw_path) || redownload) "downloaded" else "reused_existing_raw"
raw_download_method <- NA_character_
raw_download_start_year <- NA_integer_
raw_download_end_year <- NA_integer_
raw_downloaded_at <- NA_character_
raw_download_key <- NA_character_

# -----------------------------------------------------------------------------|
# Download or reuse raw occurrences ----
# -----------------------------------------------------------------------------|

if (!file.exists(raw_path) || redownload) {
  message("Downloading occurrence records for ", species, " using method: ", method)

  if (method == "direct-gbif") {
    Download_gbif(
      sp_list = species,
      initial_date = start_year,
      end_date = if (is.na(end_year)) NULL else end_year,
      exit_route = raw_dir,
      n_records = 150000,
      time_limit = if (is.na(time_limit)) NULL else time_limit
    )
  } else if (method == "gbif-download") {
    raw_download_key <- download_gbif_bulk(
      species = species,
      start_year = start_year,
      end_year = end_year,
      raw_path = raw_path,
      raw_dir = raw_dir,
      existing_download_key = gbif_download_key_arg
    )
  } else {
    iucn_key <- Sys.getenv("IUCN_REDLIST_KEY", unset = Sys.getenv("IUCN_API_KEY", unset = ""))
    if (!nzchar(iucn_key)) {
      stop(
        "`--method spatial-spp` requires IUCN_REDLIST_KEY or IUCN_API_KEY in the environment.",
        call. = FALSE
      )
    }

    download_spatial_spp_with_fallback(
      species = species,
      raw_dir = raw_dir,
      taxonomy_dir = file.path(occurrence_dir, "taxonomy"),
      start_year = start_year,
      end_year = end_year,
      iucn_key = iucn_key,
      time_limit = time_limit
    )
  }
}

if (!file.exists(raw_path)) {
  stop("Expected raw occurrence file was not found after download: ", raw_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Filter, deduplicate, and clean occurrence records ----
# -----------------------------------------------------------------------------|

raw <- read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
# The extra year filter keeps reused or resumed raw files aligned with the
# requested window, even if the download function returned broader records.
year_filter <- filter_occurrence_years(raw, start_year, end_year)
raw <- year_filter$data
raw_rows_before_year_filter <- year_filter$rows_before
raw_year_filter_status <- year_filter$status
presence_filter <- filter_gbif_present_records(raw)
raw <- presence_filter$data
raw_rows_before_presence_filter <- presence_filter$rows_before
gbif_occurrence_status_removed_rows <- presence_filter$occurrence_status_removed_rows
gbif_individual_count_zero_removed_rows <- presence_filter$individual_count_zero_removed_rows
gbif_presence_filter_removed_rows <- presence_filter$rows_removed
gbif_key_col <- "key"
raw_unique_gbif_keys <- NA_integer_
raw_duplicate_gbif_key_rows <- NA_integer_
deduplication_key <- NA_character_
deduplication_status <- "not_applied_no_gbif_key_column"
occurrence_input <- raw

if (gbif_key_col %in% names(raw)) {
  gbif_keys <- as.character(raw[[gbif_key_col]])
  usable_gbif_keys <- !is.na(gbif_keys) & nzchar(gbif_keys)
  duplicate_gbif_key_rows <- duplicated(gbif_keys) & usable_gbif_keys

  raw_unique_gbif_keys <- length(unique(gbif_keys[usable_gbif_keys]))
  raw_duplicate_gbif_key_rows <- sum(duplicate_gbif_key_rows)
  deduplication_key <- gbif_key_col
  deduplication_status <- "applied"
  occurrence_input <- raw[!duplicate_gbif_key_rows, , drop = FALSE]
}

write.csv(occurrence_input, deduplicated_path, row.names = FALSE, na = "")

# -----------------------------------------------------------------------------|
# Record raw-download provenance ----
# -----------------------------------------------------------------------------|

if (raw_download_status == "downloaded") {
  raw_download_method <- method
  raw_download_start_year <- start_year
  raw_download_end_year <- if (is.na(end_year)) NA_integer_ else end_year
  raw_downloaded_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

  raw_manifest <- data.frame(
    species_name = species,
    method = raw_download_method,
    start_year = raw_download_start_year,
    end_year = raw_download_end_year,
    raw_path = raw_path,
    raw_rows_before_year_filter = raw_rows_before_year_filter,
    raw_year_filter_status = raw_year_filter_status,
    raw_rows_before_presence_filter = raw_rows_before_presence_filter,
    gbif_occurrence_status_removed_rows = gbif_occurrence_status_removed_rows,
    gbif_individual_count_zero_removed_rows = gbif_individual_count_zero_removed_rows,
    gbif_presence_filter_removed_rows = gbif_presence_filter_removed_rows,
    raw_rows = nrow(raw),
    raw_unique_gbif_keys = raw_unique_gbif_keys,
    raw_duplicate_gbif_key_rows = raw_duplicate_gbif_key_rows,
    gbif_download_key = raw_download_key,
    downloaded_at = raw_downloaded_at,
    spatial_function = spatial_functions,
    stringsAsFactors = FALSE
  )
  write.csv(raw_manifest, raw_manifest_path, row.names = FALSE, na = "")
} else if (file.exists(raw_manifest_path)) {
  raw_manifest <- read.csv(raw_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  raw_download_method <- coalesce_scalar(raw_manifest$method)
  raw_download_start_year <- suppressWarnings(as.integer(coalesce_scalar(raw_manifest$start_year)))
  raw_download_end_year <- suppressWarnings(as.integer(coalesce_scalar(raw_manifest$end_year)))
  raw_downloaded_at <- coalesce_scalar(raw_manifest$downloaded_at)
  raw_download_key <- coalesce_scalar(raw_manifest$gbif_download_key)

  if (!is.na(raw_download_method) && raw_download_method != method) {
    warning(
      "Reusing an existing raw occurrence file downloaded with method `",
      raw_download_method,
      "` while requested method is `",
      method,
      "`. Pass `--redownload` to fetch a new raw file.",
      call. = FALSE
    )
  }
} else {
  warning(
    "Reusing an existing raw occurrence file without a raw download manifest: ",
    raw_path,
    ". Pass `--redownload` to refresh it and record provenance.",
    call. = FALSE
  )
}

required_xy <- c("decimalLongitude", "decimalLatitude")
missing_xy <- setdiff(required_xy, names(occurrence_input))
if (length(missing_xy) > 0) {
  stop("Occurrence input is missing coordinate columns: ", paste(missing_xy, collapse = ", "), call. = FALSE)
}

cleaning_result <- prepare_points_with_fallback(points_sp = occurrence_input, range_sp = NULL, xy_cols = required_xy)
cleaned <- cleaning_result$points
cleaning_status <- cleaning_result$cleaning_status

if (is.null(cleaned) || nrow(cleaned) == 0) {
  stop("Cleaning produced no usable occurrence records for ", species, call. = FALSE)
}

write.csv(cleaned, cleaned_path, row.names = FALSE, na = "")

# -----------------------------------------------------------------------------|
# Write summary and update manifest ----
# -----------------------------------------------------------------------------|

unique_coord_count <- nrow(unique(cleaned[, required_xy, drop = FALSE]))
passes_min_points <- unique_coord_count >= min_points
summary <- data.frame(
  species_name = species,
  requested_method = method,
  requested_start_year = start_year,
  requested_end_year = if (is.na(end_year)) NA_integer_ else end_year,
  min_points = min_points,
  raw_download_status = raw_download_status,
  raw_download_method = raw_download_method,
  raw_download_start_year = raw_download_start_year,
  raw_download_end_year = raw_download_end_year,
  raw_downloaded_at = raw_downloaded_at,
  raw_download_key = raw_download_key,
  raw_path = raw_path,
  raw_manifest_path = raw_manifest_path,
  raw_rows_before_year_filter = raw_rows_before_year_filter,
  raw_year_filter_status = raw_year_filter_status,
  raw_rows_before_presence_filter = raw_rows_before_presence_filter,
  gbif_occurrence_status_removed_rows = gbif_occurrence_status_removed_rows,
  gbif_individual_count_zero_removed_rows = gbif_individual_count_zero_removed_rows,
  gbif_presence_filter_removed_rows = gbif_presence_filter_removed_rows,
  occurrence_input_path = deduplicated_path,
  deduplicated_path = deduplicated_path,
  cleaned_path = cleaned_path,
  raw_rows = nrow(raw),
  raw_unique_gbif_keys = raw_unique_gbif_keys,
  raw_duplicate_gbif_key_rows = raw_duplicate_gbif_key_rows,
  deduplication_status = deduplication_status,
  deduplication_key = deduplication_key,
  occurrence_input_rows = nrow(occurrence_input),
  cleaned_rows = nrow(cleaned),
  cleaned_unique_coordinate_rows = unique_coord_count,
  passes_min_points = passes_min_points,
  cleaning_status = cleaning_status,
  cleaning_function = cleaning_functions,
  spatial_function = spatial_functions,
  prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  stringsAsFactors = FALSE
)
write.csv(summary, summary_path, row.names = FALSE, na = "")

if (update_manifest) {
  manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  idx <- manifest$species_name == species
  if (!any(idx) && "species_name" %in% names(manifest)) {
    idx <- canonical_species_key(manifest$species_name) == canonical_species_key(species)
  }
  if (!any(idx)) {
    warning("Species not found in calibration manifest, so manifest was not updated: ", species)
  } else {
    manifest <- update_if_present(manifest, idx, "occurrence_input_path", cleaned_path)
    manifest <- update_if_present(manifest, idx, "occurrence_status", "cleaned_with_sdm_pipeline")
    manifest <- update_if_present(manifest, idx, "calibration_status", "ready_for_model_regeneration")
    manifest <- update_if_present(manifest, idx, "occurrence_source", method)
    manifest <- update_if_present(manifest, idx, "occurrence_rows_raw", nrow(raw))
    manifest <- update_if_present(manifest, idx, "occurrence_rows_clean", nrow(cleaned))
    manifest <- update_if_present(manifest, idx, "passes_min_points", passes_min_points)
    if ("run_status" %in% names(manifest)) {
      selected_rows <- which(idx)
      current_status <- manifest$run_status[selected_rows]
      rows_to_update <- selected_rows[is.na(current_status) | current_status != "already_available"]
      manifest$run_status[rows_to_update] <- "occurrences_ready"
    }
    write.csv(manifest, manifest_path, row.names = FALSE, na = "")
  }
}

message("Wrote cleaned occurrences: ", cleaned_path)
message("Wrote occurrence summary: ", summary_path)
