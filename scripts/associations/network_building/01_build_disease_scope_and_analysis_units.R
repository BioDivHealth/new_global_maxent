#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_build_disease_scope_and_analysis_units.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run the current WHO disease-scope and analysis-unit stages as one
#          stable entrypoint without changing the underlying stage scripts.
# Inputs : WHO regional tables, disease/pathogen lookups, and disease master list
# Outputs: WHO pathogen backbone, zoonotic subset, analysis-unit tables, and
#          disease-master expansion review outputs
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

  cat("Running network-building disease-scope stage:", stage_file, "\n")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("Disease-scope stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed network-building disease-scope stage:", stage_file, "\n")
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

disease_scope_stages <- c(
  "1_WHO_Diseases.R",
  "1_1_WHO_Diseases_Zoonotic_Filter.R",
  "1_2_WHO_Pathogen_Analysis_Units.R",
  "1_2b_Disease_Master_Analysis_Units.R"
)

contract_outputs <- c(
  who_pathogens_diseases = who_raw_pathogens_path(),
  who_pathogens_diseases_zoonotic = who_pathogens_diseases_zoonotic_path(),
  who_pathogen_analysis_units = who_pathogen_analysis_units_path(),
  who_pathogen_analysis_units_keep = who_pathogen_analysis_units_keep_path(),
  master_disease_analysis_units = who_master_disease_analysis_units_path(),
  master_disease_name_resolution_review = who_diseases_staged_master_expansion_path(
    "master_disease_name_resolution_review.csv"
  )
)

# -----------------------------------------------------------------------------|
# 3. Run disease-scope stages ----
# -----------------------------------------------------------------------------|

invisible(lapply(disease_scope_stages, run_stage))

# -----------------------------------------------------------------------------|
# 4. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- do.call(
  rbind,
  Map(summarize_output, names(contract_outputs), unname(contract_outputs))
)

cat("Disease-scope wrapper complete. Contract output summary:\n")
print(output_summary, row.names = FALSE)
