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
# 2. Define stages and contract outputs ----
# -----------------------------------------------------------------------------|

disease_scope_stages <- c(
  "stages/disease_scope/1_WHO_Diseases.R",
  "stages/disease_scope/1_1_WHO_Diseases_Zoonotic_Filter.R",
  "stages/disease_scope/1_2_WHO_Pathogen_Analysis_Units.R",
  "stages/disease_scope/1_2b_Disease_Master_Analysis_Units.R"
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

invisible(lapply(
  disease_scope_stages,
  run_stage,
  running_label = "network-building disease-scope",
  failure_label = "Disease-scope"
))

# -----------------------------------------------------------------------------|
# 4. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- summarize_csv_outputs(contract_outputs)

cat("Disease-scope wrapper complete. Contract output summary:\n")
print(output_summary, row.names = FALSE)
