#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02_build_master_plus_registry.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run the current master-plus source matching and registry stages as
#          one stable entrypoint without changing the underlying stage scripts.
# Inputs : disease-master analysis units, manual name-resolution decisions,
#          local VIRION/CLOVER source tables, and optional transmission rules
# Outputs: master-plus pathogen matches, taxonomy review rows, registry table,
#          and host-query units
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
# 2. Define stage runner and contract outputs ----
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

  cat("Running network-building master-plus registry stage:", stage_file, "\n")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("Master-plus registry stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed network-building master-plus registry stage:", stage_file, "\n")
}

summarize_output <- function(name, path) {
  if (!file.exists(path)) {
    stop("Expected output is missing: ", path, call. = FALSE)
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

master_plus_registry_stages <- c(
  "stages/master_plus_registry/1_2c_Master_Virion_Clover_Matches.R",
  "stages/master_plus_registry/1_2d_Master_WHO_Analysis_Unit_Bridge.R"
)

contract_outputs <- c(
  master_pathogen_virion_clover_candidates = who_diseases_staged_master_expansion_path(
    "master_pathogen_virion_clover_candidates.csv"
  ),
  master_pathogen_virion_clover_matches = who_diseases_staged_master_expansion_path(
    "master_pathogen_virion_clover_matches.csv"
  ),
  master_pathogen_external_taxonomy_review = who_diseases_staged_pathogen_matching_path(
    "master_pathogen_external_taxonomy_review.csv"
  ),
  master_plus_who_analysis_units = who_master_plus_analysis_units_path(),
  master_pathogen_host_query_units = who_diseases_host_query_path(
    "master_pathogen_host_query_units.csv"
  )
)

# -----------------------------------------------------------------------------|
# 3. Run master-plus registry stages ----
# -----------------------------------------------------------------------------|

invisible(lapply(master_plus_registry_stages, run_stage))

# -----------------------------------------------------------------------------|
# 4. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- do.call(
  rbind,
  Map(summarize_output, names(contract_outputs), unname(contract_outputs))
)

cat("Master-plus registry wrapper complete. Contract output summary:\n")
print(output_summary, row.names = FALSE)
