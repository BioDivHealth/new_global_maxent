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

  cat("Running network-building master-plus host-network stage:", stage_file, "\n")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("Master-plus host-network stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed network-building master-plus host-network stage:", stage_file, "\n")
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

master_plus_host_network_stages <- c(
  "1_2e_Master_Host_Species.R",
  "1_2f_Master_Host_Species_QA_Clean.R",
  "4_2_Combine_WHO_Master_Host_Network.R"
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

invisible(lapply(master_plus_host_network_stages, run_stage))

# -----------------------------------------------------------------------------|
# 4. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- do.call(
  rbind,
  Map(summarize_output, names(contract_outputs), unname(contract_outputs))
)

cat("Master-plus host-network wrapper complete. Contract output summary:\n")
print(output_summary, row.names = FALSE)
