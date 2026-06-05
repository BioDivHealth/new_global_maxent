#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 03_prepare_chikungunya_occurrences_batch.R ----
# -----------------------------------------------------------------------------|
# Purpose: Download and clean occurrence records across an SDM target manifest
#          in one pass. Defaults remain compatible with the Chikungunya target
#          manifest.
#
# For high-volume `gbif-download` vector work, prefer the two-phase submit/fetch
# scripts so GBIF jobs can finish asynchronously.
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
    occurrence_method = "direct-gbif",
    prepare_occurrences = FALSE,
    redownload_occurrences = FALSE,
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y"))
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_batch_config <- list(
  target_manifest_path = file.path(here::here(), "sdms", "runs", "chikungunya", "sdm_target_manifest.csv"),
  occurrence_root = file.path(here::here(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  batch_run_root = file.path(here::here(), "sdms", "runs", "chikungunya", "calibration", "occurrence_batch_runs"),
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  occurrence_method = "direct-gbif",
  prepare_occurrences = FALSE,
  redownload_occurrences = FALSE,
  update_target_manifest = TRUE,
  occurrence_download_attempts = 3,
  occurrence_retry_sleep_seconds = 60,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  min_points = 20
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Config parsing helpers ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
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

run_rscript_with_retries <- function(script, script_args, log_path, attempts, sleep_seconds) {
  attempts <- max(1L, as.integer(attempts))
  sleep_seconds <- max(0, as.numeric(sleep_seconds))
  final_status <- NA_integer_

  for (attempt in seq_len(attempts)) {
    attempt_log_path <- if (attempts == 1L) {
      log_path
    } else {
      sub("[.]log$", paste0("__attempt", attempt, ".log"), log_path)
    }

    cat("  occurrence attempt ", attempt, "/", attempts, ": ", attempt_log_path, "\n", sep = "")
    final_status <- run_rscript(script, script_args, attempt_log_path)

    if (final_status == 0L) {
      if (attempt_log_path != log_path) {
        file.copy(attempt_log_path, log_path, overwrite = TRUE)
      }
      return(final_status)
    }

    if (attempt < attempts && sleep_seconds > 0) {
      Sys.sleep(sleep_seconds)
    }
  }

  final_status
}

# -----------------------------------------------------------------------------|
# Resolve batch settings ----
# -----------------------------------------------------------------------------|

target_manifest_path <- config_arg("target-manifest-path")
occurrence_root <- config_arg("occurrence-root")
batch_run_root <- config_arg("batch-run-root")
roles <- split_arg(config_arg("roles"))
species_filter <- split_arg(config_arg("species-filter"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
max_species <- as.numeric(config_arg("max-species"))
occurrence_method <- config_arg("occurrence-method")
prepare_occurrences <- as_logical_arg(config_arg("prepare-occurrences")) || has_flag(args, "prepare-occurrences")
redownload_occurrences <- as_logical_arg(config_arg("redownload-occurrences")) || has_flag(args, "redownload-occurrences")
update_target_manifest <- as_logical_arg(config_arg("update-target-manifest"))
occurrence_download_attempts <- as.integer(config_arg("occurrence-download-attempts"))
occurrence_retry_sleep_seconds <- as.numeric(config_arg("occurrence-retry-sleep-seconds"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
min_points <- as.integer(config_arg("min-points"))

if (!file.exists(target_manifest_path)) {
  stop("Missing SDM target manifest: ", target_manifest_path, call. = FALSE)
}

if (!occurrence_method %in% c("direct-gbif", "spatial-spp", "gbif-download")) {
  stop("`occurrence_method` must be one of: direct-gbif, spatial-spp, gbif-download", call. = FALSE)
}

if (prepare_occurrences && occurrence_method == "gbif-download") {
  warning(
    "`gbif-download` in this one-pass batch submits, waits, imports, and cleans serially. ",
    "For vector batches, prefer 04_submit_gbif_download_requests.R followed by ",
    "05_fetch_gbif_download_requests.R.",
    call. = FALSE
  )
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
  batch_run_root,
  timestamp
))
log_dir <- ensure_dir(file.path(batch_dir, "logs"))
summary_path <- file.path(batch_dir, "occurrence_batch_summary.csv")
script_occurrences <- file.path(
  repo_root(),
  "scripts",
  "sdms",
  "present",
  "occurrences",
  "01_prepare_one_gbif_species.R"
)

rows <- vector("list", nrow(targets))

# -----------------------------------------------------------------------------|
# Run occurrence preparation ----
# -----------------------------------------------------------------------------|

for (i in seq_len(nrow(targets))) {
  target <- targets[i, , drop = FALSE]
  species <- target$species_name_canonical[[1]]
  species_safe <- safe_species_name(species)
  occurrence_path <- file.path(
    occurrence_root,
    species_safe,
    occurrence_method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )

  occurrence_status <- if (file.exists(occurrence_path)) "already_prepared" else "missing"
  occurrence_log <- NA_character_
  occurrence_exit_status <- NA_integer_

  if (prepare_occurrences && (!file.exists(occurrence_path) || redownload_occurrences)) {
    occurrence_log <- file.path(log_dir, paste0(species_safe, "__occurrences.log"))
    occurrence_args <- c(
      "--species", species,
      "--method", occurrence_method,
      "--start-year", as.character(start_year),
      "--end-year", as.character(end_year),
      "--manifest", target_manifest_path,
      "--occurrence-root", occurrence_root,
      "--min-points", as.character(min_points)
    )
    if (!update_target_manifest) {
      occurrence_args <- c(occurrence_args, "--no-manifest-update")
    }
    if (redownload_occurrences) {
      occurrence_args <- c(occurrence_args, "--redownload")
    }

    occurrence_exit_status <- run_rscript_with_retries(
      script_occurrences,
      occurrence_args,
      occurrence_log,
      occurrence_download_attempts,
      occurrence_retry_sleep_seconds
    )
    occurrence_status <- if (occurrence_exit_status == 0 && file.exists(occurrence_path)) {
      "prepared"
    } else {
      "failed_after_retries"
    }
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
    occurrence_exit_status = occurrence_exit_status,
    occurrence_log = occurrence_log,
    start_year = start_year,
    end_year = end_year,
    prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  write.csv(do.call(rbind, rows[seq_len(i)]), summary_path, row.names = FALSE, na = "")
  cat("[", i, "/", nrow(targets), "] ", species, ": ", occurrence_status, "\n", sep = "")
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

cat("Wrote occurrence batch summary:", summary_path, "\n")
cat("Batch directory:", batch_dir, "\n")
