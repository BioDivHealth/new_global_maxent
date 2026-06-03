#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 04_submit_gbif_download_requests.R ----
# -----------------------------------------------------------------------------|
# Purpose: Submit asynchronous GBIF occurrence-download requests for Chikungunya
#          SDM targets and save the returned download keys immediately.
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
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y")),
    seed_existing_downloads = TRUE,
    refresh_existing_status = TRUE,
    resubmit_existing = FALSE,
    max_new_submissions = 3,
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
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  seed_existing_downloads = TRUE,
  refresh_existing_status = TRUE,
  resubmit_existing = FALSE,
  max_new_submissions = 3,
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

empty_request_manifest <- function() {
  data.frame(
    species_name = character(),
    species_name_canonical = character(),
    species_role = character(),
    sdm_needed_for_disease = character(),
    run_priority = integer(),
    occurrence_method = character(),
    start_year = integer(),
    end_year = integer(),
    taxon_key = integer(),
    gbif_matched_name = character(),
    gbif_download_key = character(),
    request_status = character(),
    gbif_status = character(),
    gbif_doi = character(),
    gbif_download_link = character(),
    gbif_total_records = integer(),
    gbif_size = integer(),
    gbif_created = character(),
    gbif_modified = character(),
    submitted_at = character(),
    status_checked_at = character(),
    import_status = character(),
    cleaned_path = character(),
    occurrence_summary_path = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

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

update_status_columns <- function(requests, row_idx, status_row) {
  for (col in intersect(names(status_row), names(requests))) {
    requests[[col]][row_idx] <- status_row[[col]][[1]]
  }
  requests
}

existing_request_index <- function(requests, species, start_year, end_year) {
  if (nrow(requests) == 0) {
    return(integer())
  }

  species_key <- canonical_species_name(species)
  request_species_key <- canonical_species_name(requests$species_name_canonical)
  if (!"species_name_canonical" %in% names(requests) || all(is.na(request_species_key))) {
    request_species_key <- canonical_species_name(requests$species_name)
  }

  which(
    request_species_key == species_key &
      suppressWarnings(as.integer(requests$start_year)) == start_year &
      suppressWarnings(as.integer(requests$end_year)) == end_year &
      !is.na(requests$gbif_download_key) &
      nzchar(requests$gbif_download_key)
  )
}

seed_existing_gbif_downloads <- function(requests, target_manifest) {
  manifest_paths <- list.files(
    file.path(repo_root(), "sdms", "runs", "chikungunya", "calibration", "occurrences"),
    pattern = "^raw_download_manifest[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  manifest_paths <- manifest_paths[grepl("/gbif-download/raw/raw_download_manifest[.]csv$", manifest_paths)]
  if (length(manifest_paths) == 0) {
    return(requests)
  }

  target_manifest$species_name_canonical <- canonical_species_name(target_manifest$species_name)
  requests <- ensure_manifest_columns(requests)

  for (manifest_path in manifest_paths) {
    raw_manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("species_name", "method", "start_year", "end_year", "gbif_download_key") %in% names(raw_manifest))) {
      next
    }
    if (nrow(raw_manifest) == 0 || raw_manifest$method[[1]] != "gbif-download") {
      next
    }
    if (is.na(raw_manifest$gbif_download_key[[1]]) || !nzchar(raw_manifest$gbif_download_key[[1]])) {
      next
    }

    species <- canonical_species_name(raw_manifest$species_name[[1]])
    start_year_local <- suppressWarnings(as.integer(raw_manifest$start_year[[1]]))
    end_year_local <- suppressWarnings(as.integer(raw_manifest$end_year[[1]]))
    if (length(existing_request_index(requests, species, start_year_local, end_year_local)) > 0) {
      next
    }

    target_idx <- target_manifest$species_name_canonical == species
    target_row <- if (any(target_idx)) target_manifest[which(target_idx)[[1]], , drop = FALSE] else NULL
    species_safe <- safe_species_name(species)
    cleaned_path <- file.path(
      repo_root(),
      "sdms",
      "runs",
      "chikungunya",
      "calibration",
      "occurrences",
      species_safe,
      "gbif-download",
      "cleaned",
      paste0(species_safe, "_cleaned.csv")
    )
    summary_path <- file.path(
      repo_root(),
      "sdms",
      "runs",
      "chikungunya",
      "calibration",
      "occurrences",
      species_safe,
      "gbif-download",
      "occurrence_preparation_summary.csv"
    )

    seeded_row <- data.frame(
      species_name = species,
      species_name_canonical = species,
      species_role = if (!is.null(target_row)) target_row$species_role[[1]] else NA_character_,
      sdm_needed_for_disease = if (!is.null(target_row)) target_row$sdm_needed_for_disease[[1]] else NA_character_,
      run_priority = if (!is.null(target_row)) target_row$run_priority[[1]] else NA_integer_,
      occurrence_method = "gbif-download",
      start_year = start_year_local,
      end_year = end_year_local,
      taxon_key = NA_integer_,
      gbif_matched_name = NA_character_,
      gbif_download_key = raw_manifest$gbif_download_key[[1]],
      request_status = "seeded_from_existing_raw",
      gbif_status = NA_character_,
      gbif_doi = NA_character_,
      gbif_download_link = NA_character_,
      gbif_total_records = suppressWarnings(as.integer(coalesce_scalar(raw_manifest$raw_rows))),
      gbif_size = NA_integer_,
      gbif_created = NA_character_,
      gbif_modified = NA_character_,
      submitted_at = coalesce_scalar(raw_manifest$downloaded_at),
      status_checked_at = NA_character_,
      import_status = if (file.exists(cleaned_path)) "already_cleaned" else "raw_available",
      cleaned_path = if (file.exists(cleaned_path)) cleaned_path else NA_character_,
      occurrence_summary_path = if (file.exists(summary_path)) summary_path else NA_character_,
      notes = paste0("Seeded from ", manifest_path),
      stringsAsFactors = FALSE
    )
    requests <- rbind(ensure_manifest_columns(requests), ensure_manifest_columns(seeded_row))
  }

  ensure_manifest_columns(requests)
}

refresh_gbif_request_statuses <- function(requests) {
  requests <- ensure_manifest_columns(requests)
  if (nrow(requests) == 0) {
    return(requests)
  }

  has_key <- !is.na(requests$gbif_download_key) & nzchar(requests$gbif_download_key)
  already_imported <- requests$import_status %in% c("cleaned", "already_cleaned")
  needs_status <- has_key & !already_imported
  if (!any(needs_status)) {
    return(requests)
  }

  for (row_idx in which(needs_status)) {
    download_key <- requests$gbif_download_key[[row_idx]]
    status_row <- tryCatch(
      gbif_download_status_row(download_key),
      error = function(err) {
        existing_note <- coalesce_scalar(requests$notes[row_idx], default = "")
        if (nzchar(existing_note)) {
          existing_note <- paste(existing_note, "|")
        }
        requests$notes[row_idx] <<- paste0(existing_note, "GBIF status refresh failed: ", conditionMessage(err))
        NULL
      }
    )

    if (!is.null(status_row)) {
      requests <- update_status_columns(requests, row_idx, status_row)
    }
  }

  requests
}

count_active_gbif_downloads <- function(requests) {
  requests <- ensure_manifest_columns(requests)
  if (nrow(requests) == 0) {
    return(0L)
  }

  has_key <- !is.na(requests$gbif_download_key) & nzchar(requests$gbif_download_key)
  already_imported <- requests$import_status %in% c("cleaned", "already_cleaned")
  status <- toupper(trimws(as.character(requests$gbif_status)))
  inactive_status <- status %in% c("SUCCEEDED", "FAILED", "KILLED", "CANCELLED", "CANCELED")
  missing_submitted_status <- (!nzchar(status) | is.na(status)) & requests$request_status == "submitted"

  sum(has_key & !already_imported & (!inactive_status | missing_submitted_status), na.rm = TRUE)
}

existing_request_status_label <- function(request) {
  import_status <- coalesce_scalar(request$import_status, default = "")
  gbif_status <- toupper(trimws(coalesce_scalar(request$gbif_status, default = "")))

  if (import_status %in% c("cleaned", "already_cleaned")) {
    return("already_cleaned")
  }
  if (gbif_status == "SUCCEEDED") {
    return("already_succeeded_ready_to_fetch")
  }
  if (gbif_status %in% c("RUNNING", "PREPARING", "SUBMITTED")) {
    return("already_running")
  }
  if (gbif_status %in% c("FAILED", "KILLED", "CANCELLED", "CANCELED")) {
    return(paste0("already_", tolower(gbif_status)))
  }

  "already_requested_status_unknown"
}

# -----------------------------------------------------------------------------|
# Resolve config ----
# -----------------------------------------------------------------------------|

target_manifest_path <- config_arg("target-manifest-path")
request_manifest_path <- config_arg("request-manifest-path")
roles <- split_arg(config_arg("roles"))
species_filter <- split_arg(config_arg("species-filter"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
max_species <- as.numeric(config_arg("max-species"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
resubmit_existing <- as_logical_arg(config_arg("resubmit-existing")) || has_flag(args, "resubmit-existing")
seed_existing_downloads <- as_logical_arg(config_arg("seed-existing-downloads"))
refresh_existing_status <- as_logical_arg(config_arg("refresh-existing-status"))
max_new_submissions <- as.integer(config_arg("max-new-submissions"))
dry_run <- as_logical_arg(config_arg("dry-run")) || has_flag(args, "dry-run")

if (!file.exists(target_manifest_path)) {
  stop("Missing Chikungunya SDM target manifest: ", target_manifest_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Select targets and load existing request ledger ----
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

requests <- if (file.exists(request_manifest_path)) {
  ensure_manifest_columns(read.csv(request_manifest_path, check.names = FALSE, stringsAsFactors = FALSE))
} else {
  empty_request_manifest()
}
if (seed_existing_downloads) {
  seeded_requests <- seed_existing_gbif_downloads(requests, target_manifest)
  if (nrow(seeded_requests) != nrow(requests)) {
    requests <- seeded_requests
    write_request_manifest(requests, request_manifest_path)
    cat("Seeded existing GBIF downloads:", nrow(requests), "request rows now in manifest.\n")
  }
}
if (refresh_existing_status && !dry_run) {
  requests <- refresh_gbif_request_statuses(requests)
  write_request_manifest(requests, request_manifest_path)
}

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
run_dir <- ensure_dir(file.path(
  repo_root(),
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "gbif_download_request_runs",
  timestamp
))
run_summary_path <- file.path(run_dir, "gbif_download_submit_summary.csv")
rows <- vector("list", nrow(targets))
new_submissions <- 0L
submission_limit_reached <- FALSE
active_downloads <- count_active_gbif_downloads(requests)
available_submission_slots <- max(0L, max_new_submissions - active_downloads)

cat("Selected target species:", nrow(targets), "\n")
cat("Active GBIF downloads in request ledger:", active_downloads, "\n")
cat("Available new GBIF submission slots:", available_submission_slots, "\n")
if (dry_run) {
  cat("Dry run: no GBIF requests will be submitted.\n")
}

# -----------------------------------------------------------------------------|
# Submit requests ----
# -----------------------------------------------------------------------------|

for (i in seq_len(nrow(targets))) {
  target <- targets[i, , drop = FALSE]
  species <- target$species_name_canonical[[1]]
  existing_idx <- existing_request_index(requests, species, start_year, end_year)

  row_status <- "submitted"
  notes <- NA_character_
  row_gbif_status <- NA_character_
  row_import_status <- NA_character_
  row_status_checked_at <- NA_character_
  submitted <- NULL

  if (length(existing_idx) > 0 && !resubmit_existing) {
    existing_row <- requests[existing_idx[[length(existing_idx)]], , drop = FALSE]
    row_status <- existing_request_status_label(existing_row)
    row_gbif_status <- coalesce_scalar(existing_row$gbif_status)
    row_import_status <- coalesce_scalar(existing_row$import_status)
    row_status_checked_at <- coalesce_scalar(existing_row$status_checked_at)
    notes <- paste0(
      "Existing GBIF request key ",
      existing_row$gbif_download_key[[1]],
      "; GBIF status = ",
      coalesce_scalar(existing_row$gbif_status, default = "unknown"),
      "; import status = ",
      coalesce_scalar(existing_row$import_status, default = "unknown")
    )
    submitted <- list(
      gbif_download_key = existing_row$gbif_download_key[[1]],
      taxon_key = suppressWarnings(as.integer(existing_row$taxon_key[[1]])),
      matched_name = existing_row$gbif_matched_name[[1]]
    )
  } else if (dry_run) {
    row_status <- "dry_run"
    submitted <- list(gbif_download_key = NA_character_, taxon_key = NA_integer_, matched_name = NA_character_)
  } else if (submission_limit_reached) {
    row_status <- "skipped_submission_limit"
    notes <- "Skipped because GBIF simultaneous-download limit was reached earlier in this run."
    submitted <- list(gbif_download_key = NA_character_, taxon_key = NA_integer_, matched_name = NA_character_)
  } else if (new_submissions >= available_submission_slots) {
    row_status <- "skipped_submission_limit"
    notes <- paste0(
      "Skipped because available submission slots = ",
      available_submission_slots,
      " after counting active GBIF downloads."
    )
    submitted <- list(gbif_download_key = NA_character_, taxon_key = NA_integer_, matched_name = NA_character_)
  } else {
    submitted <- tryCatch(
      submit_gbif_occurrence_download(species, start_year, end_year),
      error = function(err) {
        error_message <- conditionMessage(err)
        if (grepl("too many simultaneous downloads", error_message, ignore.case = TRUE)) {
          row_status <<- "skipped_gbif_simultaneous_download_limit"
          submission_limit_reached <<- TRUE
        } else {
          row_status <<- "failed"
        }
        notes <<- error_message
        list(gbif_download_key = NA_character_, taxon_key = NA_integer_, matched_name = NA_character_)
      }
    )
    if (row_status == "submitted") {
      new_submissions <- new_submissions + 1L
    }
  }

  request_row <- data.frame(
    species_name = species,
    species_name_canonical = species,
    species_role = target$species_role[[1]],
    sdm_needed_for_disease = target$sdm_needed_for_disease[[1]],
    run_priority = target$run_priority[[1]],
    occurrence_method = "gbif-download",
    start_year = start_year,
    end_year = end_year,
    taxon_key = submitted$taxon_key,
    gbif_matched_name = submitted$matched_name,
    gbif_download_key = submitted$gbif_download_key,
    request_status = row_status,
    gbif_status = row_gbif_status,
    gbif_doi = NA_character_,
    gbif_download_link = NA_character_,
    gbif_total_records = NA_integer_,
    gbif_size = NA_integer_,
    gbif_created = NA_character_,
    gbif_modified = NA_character_,
    submitted_at = if (row_status == "submitted") {
      format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    } else {
      NA_character_
    },
    status_checked_at = row_status_checked_at,
    import_status = row_import_status,
    cleaned_path = NA_character_,
    occurrence_summary_path = NA_character_,
    notes = notes,
    stringsAsFactors = FALSE
  )

  rows[[i]] <- request_row
  if (row_status %in% c("submitted", "failed", "skipped_gbif_simultaneous_download_limit") && !dry_run) {
    requests <- rbind(ensure_manifest_columns(requests), ensure_manifest_columns(request_row))
    write_request_manifest(requests, request_manifest_path)
  }

  write.csv(do.call(rbind, rows[seq_len(i)]), run_summary_path, row.names = FALSE, na = "")
  cat("[", i, "/", nrow(targets), "] ", species, ": ", row_status, "\n", sep = "")
}

summary <- if (length(rows) == 0) {
  data.frame()
} else {
  do.call(rbind, rows)
}
write.csv(summary, run_summary_path, row.names = FALSE, na = "")

if (!dry_run) {
  write_request_manifest(requests, request_manifest_path)
}

cat("Wrote submit run summary:", run_summary_path, "\n")
if (!dry_run) {
  cat("Wrote GBIF request manifest:", request_manifest_path, "\n")
}
