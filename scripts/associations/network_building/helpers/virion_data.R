# -----------------------------------------------------------------------------|
# virion_data.R ----
# -----------------------------------------------------------------------------|
# Purpose: Backward-compatible shim for scripts that source this file and expect
#          a loaded `virion_data` object.
#
# For reusable loader functions without automatically loading local VIRION data,
# source `scripts/associations/network_building/helpers/virion_loaders.R`.
# -----------------------------------------------------------------------------|
#
# VIRION: A database of host-virus interactions
# https://github.com/viralemergence/virion
#
# VIRION is a comprehensive database that combines data from:
# - CLOVER (static source)
# - PREDICT (static source)
# - GenBank (dynamic source)
#
# The database contains over 1,162,000 host-virus interactions
# across 9,521 viruses and 3,692 hosts
#
# Data is available via:
# - virionData R package: https://github.com/viralemergence/virionData
# - Zenodo: https://zenodo.org/records/10418723

# Load required libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(vroom)  # For reading compressed CSV files
library(viridis)
library(fs)     # For file system operations
library(kableExtra)  # For nice tables
library(jsonlite)    # For JSON parsing
library(rlang)  # For dynamic column references

source(file.path("scripts", "associations", "working_inputs.R"))

# =============================================================================
# Data loading helpers
# =============================================================================
# library(remotes)
# remotes::install_github("viralemergence/virionData", force = TRUE)
source(file.path(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "virion_loaders.R"
))

# =============================================================================

virion_data <- load_virion_data()
