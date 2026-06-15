# VIRION Data Analysis Script
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
# DATA LOADING FUNCTIONS
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


# I guess next steps will be to use virion_data$taxonomy_virus [Virus & VirusFamily]
# Columns and match ours WHO pathogens to those to then identify hosts
