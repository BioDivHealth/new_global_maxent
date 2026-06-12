#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 03_build_master_plus_host_network.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run the current master-plus host-evidence, QA-cleaning, and combined
#          host-network stages as one stable entrypoint without changing the
#          underlying stage scripts.
# Inputs : master-plus host-query units, local VIRION/CLOVER source tables,
#          host taxonomy standardization outputs, and WHO host network
# Outputs: master host-species evidence, QA summary, cleaned host evidence, and
#          master-plus WHO host network
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

master_plus_host_network_stages <- c(
  "stages/master_plus_host_network/1_2e_Master_Host_Species.R",
  "stages/master_plus_host_network/1_2f_Master_Host_Species_QA_Clean.R",
  "stages/master_plus_host_network/4_2_Combine_WHO_Master_Host_Network.R"
)

contract_outputs <- c(
  master_pathogen_host_species = who_master_pathogen_host_species_path(),
  master_pathogen_host_species_summary = who_master_pathogen_host_species_summary_path(),
  master_pathogen_host_species_clean = who_master_pathogen_host_species_clean_path(),
  master_plus_who_host_network = who_network_host_pathogen_path(
    "master_plus_who_host_network.csv"
  )
)

# -----------------------------------------------------------------------------|
# 3. Run master-plus host-network stages ----
# -----------------------------------------------------------------------------|

invisible(lapply(
  master_plus_host_network_stages,
  run_stage,
  running_label = "network-building master-plus host-network",
  failure_label = "Master-plus host-network"
))

# -----------------------------------------------------------------------------|
# 4. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- summarize_csv_outputs(contract_outputs)

cat("Master-plus host-network wrapper complete. Contract output summary:\n")
print(output_summary, row.names = FALSE)
