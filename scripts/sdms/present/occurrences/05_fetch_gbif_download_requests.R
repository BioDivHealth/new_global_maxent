#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 05_fetch_gbif_download_requests.R ----
# -----------------------------------------------------------------------------|
# Purpose: Check submitted GBIF download keys, fetch ready ZIPs, and clean the
#          resulting occurrence records.
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
    fetch_statuses = "SUCCEEDED",
    redownload_occurrences = FALSE,
    dry_run = FALSE
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "chikungunya", "sdm_target_manifest.csv"),
  request_manifest_path = file.path(repo_root(), "sdms", "runs", "chikungunya", "calibration", "gbif_download_requests.csv"),
  roles = "vector",
  species_filter = character(),
  max_species = Inf,
  include_imported = FALSE,
  fetch_statuses = "SUCCEEDED",
  refresh_status = TRUE,
  update_target_manifest = TRUE,
  redownload_occurrences = FALSE,
  min_points = 20,
  dry_run = FALSE
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Helpers ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
}

request_columns <- c(
  "species_name",
  "species_name_canonical",
  "species_role",
  "sdm_needed_for_disease",
  "run_priority",
  "occurrence_method",
  "start_year",
  "end_year",
  "taxon_key",
  "gbif_matched_name",
  "gbif_download_key",
  "request_status",
  "gbif_status",
  "gbif_doi",
  "gbif_download_link",
  "gbif_total_records",
  "gbif_size",
  "gbif_created",
  "gbif_modified",
  "submitted_at",
  "status_checked_at",
  "import_status",
  "cleaned_path",
  "occurrence_summary_path",
  "notes"
)

ensure_manifest_columns <- function(data) {
  for (col in request_columns) {
    if (!col %in% names(data)) {
      data[[col]] <- NA
    }
  }
  data[, request_columns, drop = FALSE]
}

write_request_manifest <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(ensure_manifest_columns(data), path, row.names = FALSE, na = "")
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

classify_import_failure <- function(log_path) {
  if (is.na(log_path) || !file.exists(log_path)) {
    return("failed")
  }

  log_text <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  if (grepl("Cleaning produced no usable occurrence records", log_text, fixed = TRUE)) {
    return("cleaning_failed_no_usable_records")
  }
  if (grepl("Occurrence input is missing coordinate columns", log_text, fixed = TRUE)) {
    return("cleaning_failed_missing_coordinate_columns")
  }

  "failed"
}

update_status_columns <- function(requests, row_idx, status_row) {
  for (col in intersect(names(status_row), names(requests))) {
    requests[[col]][row_idx] <- status_row[[col]][[1]]
  }
  requests
}

cleaned_occurrence_path <- function(species, method = "gbif-download") {
  species_safe <- safe_species_name(species)
  file.path(
    repo_root(),
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
}

occurrence_summary_path <- function(species, method = "gbif-download") {
  file.path(
    repo_root(),
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "occurrences",
    safe_species_name(species),
    method,
    "occurrence_preparation_summary.csv"
  )
}

# -----------------------------------------------------------------------------|
# Resolve config ----
# -----------------------------------------------------------------------------|

target_manifest_path <- config_arg("target-manifest-path")
request_manifest_path <- config_arg("request-manifest-path")
roles <- split_arg(config_arg("roles"))
species_filter <- canonical_species_name(split_arg(config_arg("species-filter")))
max_species <- as.numeric(config_arg("max-species"))
include_imported <- as_logical_arg(config_arg("include-imported"))
fetch_statuses <- split_arg(config_arg("fetch-statuses"))
refresh_status <- as_logical_arg(config_arg("refresh-status"))
update_target_manifest <- as_logical_arg(config_arg("update-target-manifest"))
redownload_occurrences <- as_logical_arg(config_arg("redownload-occurrences")) || has_flag(args, "redownload-occurrences")
min_points <- as.integer(config_arg("min-points"))
dry_run <- as_logical_arg(config_arg("dry-run")) || has_flag(args, "dry-run")

if (!file.exists(request_manifest_path)) {
  stop("Missing GBIF request manifest: ", request_manifest_path, call. = FALSE)
}

if (update_target_manifest && !file.exists(target_manifest_path)) {
  stop("Missing Chikungunya SDM target manifest: ", target_manifest_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Select request rows ----
# -----------------------------------------------------------------------------|

requests <- ensure_manifest_columns(read.csv(request_manifest_path, check.names = FALSE, stringsAsFactors = FALSE))
request_species <- canonical_species_name(requests$species_name_canonical)
fallback_species <- canonical_species_name(requests$species_name)
request_species[is.na(request_species)] <- fallback_species[is.na(request_species)]

selected <- rep(TRUE, nrow(requests))
selected <- selected & requests$occurrence_method == "gbif-download"
selected <- selected & !is.na(requests$gbif_download_key) & nzchar(requests$gbif_download_key)
if (length(roles) > 0 && !any(roles %in% c("all", "ALL"))) {
  selected <- selected & requests$species_role %in% roles
}
if (length(species_filter) > 0) {
  selected <- selected & request_species %in% species_filter
}
if (!include_imported) {
  selected <- selected & !(requests$import_status %in% c("cleaned", "already_cleaned"))
}

selected_idx <- which(selected)
selected_idx <- selected_idx[order(requests$run_priority[selected_idx], request_species[selected_idx])]
if (is.finite(max_species)) {
  selected_idx <- head(selected_idx, max_species)
}

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
run_dir <- ensure_dir(file.path(
  repo_root(),
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "gbif_download_fetch_runs",
  timestamp
))
log_dir <- ensure_dir(file.path(run_dir, "logs"))
run_summary_path <- file.path(run_dir, "gbif_download_fetch_summary.csv")
script_occurrences <- file.path(
  repo_root(),
  "scripts",
  "sdms",
  "present",
  "occurrences",
  "01_prepare_gbif_occurrences.R"
)
rows <- vector("list", length(selected_idx))

cat("Selected GBIF download requests:", length(selected_idx), "\n")
if (dry_run) {
  cat("Dry run: ready downloads will not be imported or cleaned.\n")
}

# -----------------------------------------------------------------------------|
# Check status and fetch ready requests ----
# -----------------------------------------------------------------------------|

for (pos in seq_along(selected_idx)) {
  row_idx <- selected_idx[[pos]]
  request <- requests[row_idx, , drop = FALSE]
  species <- request_species[[row_idx]]
  if (is.na(species) || !nzchar(species)) {
    species <- request$species_name[[1]]
  }
  species_safe <- safe_species_name(species)
  download_key <- request$gbif_download_key[[1]]
  checked_status <- request$gbif_status[[1]]
  notes <- NA_character_

  if (refresh_status) {
    status_row <- tryCatch(
      gbif_download_status_row(download_key),
      error = function(err) {
        notes <<- conditionMessage(err)
        data.frame(
          gbif_download_key = download_key,
          gbif_status = NA_character_,
          gbif_doi = NA_character_,
          gbif_download_link = NA_character_,
          gbif_total_records = NA_integer_,
          gbif_size = NA_integer_,
          gbif_created = NA_character_,
          gbif_modified = NA_character_,
          status_checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
          stringsAsFactors = FALSE
        )
      }
    )
    checked_status <- status_row$gbif_status[[1]]
    requests <- update_status_columns(requests, row_idx, status_row)
    write_request_manifest(requests, request_manifest_path)
  }

  cleaned_path <- cleaned_occurrence_path(species)
  summary_path <- occurrence_summary_path(species)
  occurrence_log <- NA_character_
  occurrence_exit_status <- NA_integer_
  import_status <- "skipped_not_ready"

  if (file.exists(cleaned_path) && !redownload_occurrences) {
    import_status <- "already_cleaned"
  } else if (is.na(checked_status) || !checked_status %in% fetch_statuses) {
    import_status <- if (is.na(checked_status)) "status_check_failed" else paste0("skipped_", tolower(checked_status))
  } else if (dry_run) {
    import_status <- "ready_dry_run"
  } else {
    occurrence_log <- file.path(log_dir, paste0(species_safe, "__fetch_clean.log"))
    occurrence_args <- c(
      "--species", species,
      "--method", "gbif-download",
      "--gbif-download-key", download_key,
      "--start-year", as.character(suppressWarnings(as.integer(request$start_year[[1]]))),
      "--end-year", as.character(suppressWarnings(as.integer(request$end_year[[1]]))),
      "--manifest", target_manifest_path,
      "--min-points", as.character(min_points)
    )
    if (!update_target_manifest) {
      occurrence_args <- c(occurrence_args, "--no-manifest-update")
    }
    if (redownload_occurrences) {
      occurrence_args <- c(occurrence_args, "--redownload")
    }

    occurrence_exit_status <- run_rscript(script_occurrences, occurrence_args, occurrence_log)
    import_status <- if (occurrence_exit_status == 0 && file.exists(cleaned_path)) {
      "cleaned"
    } else {
      classify_import_failure(occurrence_log)
    }
  }

  requests$import_status[row_idx] <- import_status
  requests$cleaned_path[row_idx] <- if (file.exists(cleaned_path)) cleaned_path else NA_character_
  requests$occurrence_summary_path[row_idx] <- if (file.exists(summary_path)) summary_path else NA_character_
  if (!is.na(notes)) {
    requests$notes[row_idx] <- notes
  }
  write_request_manifest(requests, request_manifest_path)

  rows[[pos]] <- data.frame(
    species_name = species,
    species_role = request$species_role[[1]],
    run_priority = request$run_priority[[1]],
    gbif_download_key = download_key,
    gbif_status = checked_status,
    import_status = import_status,
    occurrence_exit_status = occurrence_exit_status,
    occurrence_log = occurrence_log,
    cleaned_path = if (file.exists(cleaned_path)) cleaned_path else NA_character_,
    occurrence_summary_path = if (file.exists(summary_path)) summary_path else NA_character_,
    checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    notes = notes,
    stringsAsFactors = FALSE
  )

  write.csv(do.call(rbind, rows[seq_len(pos)]), run_summary_path, row.names = FALSE, na = "")
  cat("[", pos, "/", length(selected_idx), "] ", species, ": ", checked_status, ", ", import_status, "\n", sep = "")
}

summary <- if (length(rows) == 0) {
  data.frame()
} else {
  do.call(rbind, rows)
}
write.csv(summary, run_summary_path, row.names = FALSE, na = "")
write_request_manifest(requests, request_manifest_path)

cat("Wrote fetch run summary:", run_summary_path, "\n")
cat("Updated GBIF request manifest:", request_manifest_path, "\n")
