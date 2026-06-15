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
source(here::here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "network_wrapper_helpers.R"
))
require_network_wrapper_packages()
source_network_working_inputs()

# -----------------------------------------------------------------------------|
# 2. Parse command-line arguments ----
# -----------------------------------------------------------------------------|

args <- parse_network_wrapper_args(
  valid_flags = c("--refresh-ncbi-metadata", "--help", "-h"),
  help_text = paste(
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
  )
)

refresh_ncbi_metadata <- "--refresh-ncbi-metadata" %in% args

# -----------------------------------------------------------------------------|
# 3. Define stages and broad-taxa outputs ----
# -----------------------------------------------------------------------------|

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

stages <- c("stages/broad_taxa_support/1_3_WHO_Broad_Taxa_Candidate_Strains.R")
if (refresh_ncbi_metadata) {
  stages <- c(
    stages,
    "stages/broad_taxa_support/1_4_NCBI_Broad_Taxa_Candidate_Metadata.R"
  )
}

# -----------------------------------------------------------------------------|
# 4. Run broad-taxa support stages ----
# -----------------------------------------------------------------------------|

if (refresh_ncbi_metadata) {
  cat("Refreshing NCBI candidate metadata after rebuilding candidate strains.\n")
} else {
  cat("Rebuilding candidate strains only; NCBI metadata refresh is opt-in.\n")
}

invisible(lapply(
  stages,
  run_stage,
  running_label = "broad-taxa support",
  failure_label = "Broad-taxa support"
))

# -----------------------------------------------------------------------------|
# 5. Print output summary ----
# -----------------------------------------------------------------------------|

candidate_summary <- summarize_csv_outputs(candidate_output)
ncbi_csv_summary <- summarize_csv_outputs(ncbi_csv_outputs, required = FALSE)
ncbi_text_summary <- summarize_text_outputs(ncbi_text_outputs)

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
