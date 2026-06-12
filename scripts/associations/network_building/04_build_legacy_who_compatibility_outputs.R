#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 04_build_legacy_who_compatibility_outputs.R ----
# -----------------------------------------------------------------------------|
# Purpose: Run the legacy WHO-only CLOVER/VIRION compatibility stages as one
#          stable entrypoint while reusing existing host-taxonomy outputs by
#          default.
# Inputs : WHO pathogen backbone, local CLOVER/VIRION source tables, existing
#          standardized host taxonomy CSVs, and current WHO analysis-unit tables
# Outputs: legacy WHO-only source components, combined raw network, canonical
#          lookup, canonical network, and canonical zoonotic network
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
valid_flags <- c("--refresh-host-taxonomy", "--help", "-h")
unknown_args <- setdiff(args, valid_flags)
if (length(unknown_args) > 0) {
  stop("Unknown arguments: ", paste(unknown_args, collapse = ", "), call. = FALSE)
}

if (any(args %in% c("--help", "-h"))) {
  cat(
    paste(
      "Legacy WHO compatibility network wrapper",
      "",
      "Usage:",
      "  Rscript scripts/associations/network_building/04_build_legacy_who_compatibility_outputs.R",
      "  Rscript scripts/associations/network_building/04_build_legacy_who_compatibility_outputs.R --refresh-host-taxonomy",
      "",
      "Options:",
      "  --refresh-host-taxonomy  Also run the CLOVER and VIRION host taxonomy",
      "                           standardization stages. Default mode reuses",
      "                           existing standardized host-taxonomy CSVs.",
      "  --help, -h               Show this help message.",
      sep = "\n"
    ),
    "\n"
  )
  quit(status = 0)
}

refresh_host_taxonomy <- "--refresh-host-taxonomy" %in% args

# -----------------------------------------------------------------------------|
# 3. Define stage runner and compatibility outputs ----
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

  cat("Running legacy WHO compatibility stage:", stage_file, "\n")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("Legacy WHO compatibility stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed legacy WHO compatibility stage:", stage_file, "\n")
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

required_taxonomy_outputs <- c(
  clover_host_species_standardized = file.path(
    who_clover_dir,
    "clover_host_species_standardized.csv"
  ),
  who_host_species_standardized = file.path(
    who_virion_dir,
    "who_host_species_standardized.csv"
  )
)

compatibility_outputs <- c(
  who_bacteria_clover_taxid = file.path(who_clover_dir, "who_bacteria_clover_taxid.csv"),
  who_bacteria_clover_hosts = file.path(who_clover_dir, "who_bacteria_clover_hosts.csv"),
  who_bacteria_clover_unique_hosts = file.path(
    who_clover_dir,
    "who_bacteria_clover_unique_hosts.csv"
  ),
  clover_host_species_standardized = required_taxonomy_outputs[[
    "clover_host_species_standardized"
  ]],
  clover_who_network = who_network_source_component_path("clover_who_network.csv"),
  who_pathogens_virion_taxid = file.path(who_virion_dir, "who_pathogens_virion_taxid.csv"),
  who_pathogens_virion_hosts_long = file.path(
    who_virion_dir,
    "who_pathogens_virion_hosts_long.csv"
  ),
  who_pathogens_virion_hosts_summary = file.path(
    who_virion_dir,
    "who_pathogens_virion_hosts_summary.csv"
  ),
  who_host_species_standardized = required_taxonomy_outputs[[
    "who_host_species_standardized"
  ]],
  virion_who_network = who_network_source_component_path("virion_who_network.csv"),
  combined_who_network = who_raw_network_path(),
  combined_who_pathogen_canonical_lookup = who_network_canonicalization_path(
    "combined_who_pathogen_canonical_lookup.csv"
  ),
  combined_who_network_canonical = who_canonical_network_path(),
  combined_who_network_canonical_zoonotic = who_canonical_zoonotic_network_path()
)

default_stages <- c(
  "stages/legacy_who_compatibility/2_1_CLOVER.R",
  "stages/legacy_who_compatibility/2_3_CLOVER_Network.R",
  "stages/legacy_who_compatibility/3_1_Match_WHO_Virion.R",
  "stages/legacy_who_compatibility/3_2_WHO_Virion_Hosts.R",
  "stages/legacy_who_compatibility/3_4_VIRION_Networks.R",
  "stages/legacy_who_compatibility/4_CombineNetworks.R",
  "stages/legacy_who_compatibility/4_1_Canonicalize_Combined_WHO_Network.R"
)

taxonomy_refresh_stages <- c(
  "stages/legacy_who_compatibility/2_1_CLOVER.R",
  "stages/legacy_who_compatibility/2_2_CLOVER_Host_Clean.R",
  "stages/legacy_who_compatibility/2_3_CLOVER_Network.R",
  "stages/legacy_who_compatibility/3_1_Match_WHO_Virion.R",
  "stages/legacy_who_compatibility/3_2_WHO_Virion_Hosts.R",
  "stages/legacy_who_compatibility/3_3_Host_Species_Clean.R",
  "stages/legacy_who_compatibility/3_4_VIRION_Networks.R",
  "stages/legacy_who_compatibility/4_CombineNetworks.R",
  "stages/legacy_who_compatibility/4_1_Canonicalize_Combined_WHO_Network.R"
)

# -----------------------------------------------------------------------------|
# 4. Run legacy compatibility stages ----
# -----------------------------------------------------------------------------|

if (!refresh_host_taxonomy) {
  missing_taxonomy <- required_taxonomy_outputs[!file.exists(required_taxonomy_outputs)]
  if (length(missing_taxonomy) > 0) {
    stop(
      "Missing required host-taxonomy outputs for default mode: ",
      paste(missing_taxonomy, collapse = "; "),
      ". Recreate them separately or rerun with --refresh-host-taxonomy.",
      call. = FALSE
    )
  }

  cat("Reusing existing host-taxonomy outputs:\n")
  cat(paste(unname(required_taxonomy_outputs), collapse = "\n"), "\n")
  stages <- default_stages
} else {
  cat("Refreshing host taxonomy before rebuilding compatibility networks.\n")
  stages <- taxonomy_refresh_stages
}

invisible(lapply(stages, run_stage))

# -----------------------------------------------------------------------------|
# 5. Print output summary ----
# -----------------------------------------------------------------------------|

output_summary <- do.call(
  rbind,
  Map(summarize_output, names(compatibility_outputs), unname(compatibility_outputs))
)

cat("Legacy WHO compatibility wrapper complete. Output summary:\n")
print(output_summary, row.names = FALSE)
