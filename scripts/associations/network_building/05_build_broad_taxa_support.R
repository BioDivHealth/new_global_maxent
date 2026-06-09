#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 05_build_broad_taxa_support.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run the broad-taxa candidate-strain support stage as one stable
#          entrypoint, with external NCBI metadata refresh kept opt-in.
# Inputs : broad-taxa manual seed table and current WHO analysis-unit tables
# Outputs: broad-taxa candidate strains and optional existing/refreshable NCBI
#          candidate metadata support files
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# 1. Load required libraries and path helpers ----
# -----------------------------------------------------------------------------|

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package `here` is required.", call. = FALSE)
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package `readr` is required.", call. = FALSE)
}

source(here::here("scripts", "associations", "working_inputs.R"))

# -----------------------------------------------------------------------------|
# 2. Parse command-line arguments ----
# -----------------------------------------------------------------------------|

args <- commandArgs(trailingOnly = TRUE)
valid_flags <- c("--refresh-ncbi-metadata", "--help", "-h")
unknown_args <- setdiff(args, valid_flags)
if (length(unknown_args) > 0) {
  stop("Unknown arguments: ", paste(unknown_args, collapse = ", "), call. = FALSE)
}

if (any(args %in% c("--help", "-h"))) {
  cat(
    paste(
      "Broad-taxa support wrapper",
      "",
      "Usage:",
      "  Rscript scripts/associations/network_building/05_build_broad_taxa_support.R",
      "  Rscript scripts/associations/network_building/05_build_broad_taxa_support.R --refresh-ncbi-metadata",
      "",
      "Options:",
      "  --refresh-ncbi-metadata  Also run the NCBI Datasets metadata refresh.",
      "                           Default mode rebuilds candidate strains only",
      "                           and summarizes existing NCBI outputs when",
      "                           present.",
      "  --help, -h               Show this help message.",
      sep = "\n"
    ),
    "\n"
  )
  quit(status = 0)
}

refresh_ncbi_metadata <- "--refresh-ncbi-metadata" %in% args

# -----------------------------------------------------------------------------|
# 3. Define stage runner and broad-taxa outputs ----
# -----------------------------------------------------------------------------|

network_building_script <- function(filename) {
  normalizePath(
    here::here("scripts", "associations", "network_building", filename),
    mustWork = TRUE
  )
}

run_stage <- function(stage_file) {
  stage_path <- network_building_script(stage_file)
  rscript <- normalizePath(file.path(R.home("bin"), "Rscript"), mustWork = TRUE)

  cat("Running broad-taxa support stage:", stage_file, "\n")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("Broad-taxa support stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed broad-taxa support stage:", stage_file, "\n")
}

summarize_csv_output <- function(name, path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Expected output is missing: ", path, call. = FALSE)
    }
    return(NULL)
  }

  data <- readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
  data.frame(
    output = name,
    rows = nrow(data),
    columns = ncol(data),
    path = path,
    check.names = FALSE
  )
}

summarize_text_output <- function(name, path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  data.frame(
    output = name,
    rows = length(readLines(path, warn = FALSE)),
    columns = NA_integer_,
    path = path,
    check.names = FALSE
  )
}

candidate_output <- c(
  who_broad_taxa_candidate_strains = who_diseases_broad_taxa_staged_path(
    "who_broad_taxa_candidate_strains.csv"
  )
)

ncbi_csv_outputs <- c(
  who_broad_taxa_candidate_strains_ncbi_metadata = who_diseases_broad_taxa_staged_path(
    "who_broad_taxa_candidate_strains_ncbi_metadata.csv"
  ),
  who_broad_taxa_candidate_strains_ncbi_enriched = who_diseases_broad_taxa_staged_path(
    "who_broad_taxa_candidate_strains_ncbi_enriched.csv"
  )
)

ncbi_text_outputs <- c(
  who_broad_taxa_candidate_strains_ncbi_raw = who_diseases_broad_taxa_staged_path(
    "who_broad_taxa_candidate_strains_ncbi_raw.jsonl"
  )
)

stages <- c("1_3_WHO_Broad_Taxa_Candidate_Strains.R")
if (refresh_ncbi_metadata) {
  stages <- c(stages, "1_4_NCBI_Broad_Taxa_Candidate_Metadata.R")
}

# -----------------------------------------------------------------------------|
# 4. Run broad-taxa support stages ----
# -----------------------------------------------------------------------------|

if (refresh_ncbi_metadata) {
  cat("Refreshing NCBI candidate metadata after rebuilding candidate strains.\n")
} else {
  cat("Rebuilding candidate strains only; NCBI metadata refresh is opt-in.\n")
}

invisible(lapply(stages, run_stage))

# -----------------------------------------------------------------------------|
# 5. Print output summary ----
# -----------------------------------------------------------------------------|

candidate_summary <- do.call(
  rbind,
  Map(summarize_csv_output, names(candidate_output), unname(candidate_output))
)

ncbi_csv_summary <- do.call(
  rbind,
  Map(
    summarize_csv_output,
    names(ncbi_csv_outputs),
    unname(ncbi_csv_outputs),
    MoreArgs = list(required = FALSE)
  )
)

ncbi_text_summary <- do.call(
  rbind,
  Map(summarize_text_output, names(ncbi_text_outputs), unname(ncbi_text_outputs))
)

output_summary <- do.call(
  rbind,
  Filter(Negate(is.null), list(candidate_summary, ncbi_csv_summary, ncbi_text_summary))
)

if (!refresh_ncbi_metadata) {
  missing_ncbi_outputs <- c(
    ncbi_csv_outputs[!file.exists(ncbi_csv_outputs)],
    ncbi_text_outputs[!file.exists(ncbi_text_outputs)]
  )

  if (length(missing_ncbi_outputs) > 0) {
    cat("Existing NCBI metadata outputs not found:\n")
    cat(paste(unname(missing_ncbi_outputs), collapse = "\n"), "\n")
    cat("Run with --refresh-ncbi-metadata to recreate them when needed.\n")
  } else {
    cat("Reused existing NCBI metadata outputs; default mode did not refresh them.\n")
  }
}

cat("Broad-taxa support wrapper complete. Output summary:\n")
print(output_summary, row.names = FALSE)
