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

#' Load VIRION data from local files
#' @param data_path Path to VIRION data directory
#' @param files Vector of file names to load (default: all main files)
#' @return List containing loaded data frames
load_virion_data <- function(data_path = virion_source_version_dir,
                            files = c("virion.csv.gz", "edgelist.csv", 
                                     "taxonomy_host.csv", "taxonomy_virus.csv",
                                     "provenance.csv.gz", "detection.csv.gz", 
                                     "temporal.csv.gz")) {
  read_virion_csv <- function(path) {
    tryCatch(
      readr::read_csv(path, show_col_types = FALSE),
      error = function(e) {
        if (grepl("input string .* is invalid", conditionMessage(e), ignore.case = TRUE)) {
          warning(
            "Encountered invalid text encoding while reading ", basename(path),
            "; retrying with Latin1 decoding."
          )
          return(
            readr::read_csv(
              path,
              show_col_types = FALSE,
              locale = readr::locale(encoding = "Latin1")
            )
          )
        }
        stop(e)
      }
    )
  }
  
  # Check if data directory exists
  if (!dir.exists(data_path)) {
    stop("VIRION data directory not found. Please download data from Zenodo or use virionData package.")
  }
  
  # Load each file
  data_list <- list()
  
  for (file in files) {
    file_path <- file.path(data_path, file)
    
    if (file.exists(file_path)) {
      cat("Loading:", file, "\n")
      
      table_name <- gsub("\\.csv(\\.gz)?$", "", file)
      data_list[[table_name]] <- read_virion_csv(file_path)
    } else {
      warning("File not found:", file_path)
    }
  }
  
  return(data_list)
}

#' Load VIRION data using the virionData package
#' @param version Version of VIRION data to load (default: "latest")
#' @param tables Vector of table names to load (default: all main tables)
#' @return List containing loaded data frames

load_virion_package <- function(version = "latest", 
                                tables = c("virion", "edgelist", "taxonomy_host", 
                                          "taxonomy_virus", "provenance", "detection", "temporal")) {
  # Check if virionData package is available
  if (!require(virionData, quietly = TRUE)) {
    cat("virionData package not installed. Installing...\n")
    if (!require(remotes, quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("viralemergence/virionData")
    library(virionData)
  }
  
  # Initialize data list
  virion_data <- list()
  
  read_virion_csv <- function(path) {
    tryCatch(
      readr::read_csv(path, show_col_types = FALSE),
      error = function(e) {
        if (grepl("input string .* is invalid", conditionMessage(e), ignore.case = TRUE)) {
          warning(
            "Encountered invalid text encoding while reading ", basename(path),
            "; retrying with Latin1 decoding."
          )
          return(
            readr::read_csv(
              path,
              show_col_types = FALSE,
              locale = readr::locale(encoding = "Latin1")
            )
          )
        }
        stop(e)
      }
    )
  }
  
  cat("Loading VIRION data tables...\n")
  
  # Get the latest version data using get_versioned_data
  cat("  - Downloading VIRION data from Zenodo...\n")
  data_path <- virionData::get_versioned_data(version = version, dir_path = virion_source_dir)
  
  # Load the CSV files directly
  if ("virion" %in% tables) {
    cat("  - Loading main virion interactions...\n")
    virion_data$virion <- read_virion_csv(file.path(data_path, "virion.csv.gz"))
  }
  
  if ("edgelist" %in% tables) {
    cat("  - Loading edgelist...\n")
    virion_data$edgelist <- readr::read_csv(file.path(data_path, "edgelist.csv"), show_col_types = FALSE)
  }
  
  if ("taxonomy_host" %in% tables) {
    cat("  - Loading host taxonomy...\n")
    virion_data$taxonomy_host <- readr::read_csv(file.path(data_path, "taxonomy_host.csv"), show_col_types = FALSE)
  }
  
  if ("taxonomy_virus" %in% tables) {
    cat("  - Loading virus taxonomy...\n")
    virion_data$taxonomy_virus <- readr::read_csv(file.path(data_path, "taxonomy_virus.csv"), show_col_types = FALSE)
  }
  
  if ("provenance" %in% tables) {
    cat("  - Loading provenance data...\n")
    virion_data$provenance <- read_virion_csv(file.path(data_path, "provenance.csv.gz"))
  }
  
  if ("detection" %in% tables) {
    cat("  - Loading detection methods...\n")
    virion_data$detection <- vroom::vroom(file.path(data_path, "detection.csv.gz"), show_col_types = FALSE)
  }
  
  if ("temporal" %in% tables) {
    cat("  - Loading temporal data...\n")
    virion_data$temporal <- vroom::vroom(file.path(data_path, "temporal.csv.gz"), show_col_types = FALSE)
  }
  
  cat("Data loading complete!\n")
  return(virion_data)
}

#' Get available VIRION data versions and metadata
#' @return List with version information and metadata
get_virion_versions <- function() {
  if (!require(virionData, quietly = TRUE)) {
    stop("virionData package not installed")
  }
  
  # Get available versions using virionData functions
  versions <- virionData::list_deposit_versions()
  summary_info <- virionData::deposit_summary()
  
  return(list(
    versions = versions,
    summary = summary_info
  ))
}


# =============================================================================

virion_data <- load_virion_data()


# I guess next steps will be to use virion_data$taxonomy_virus [Virus & VirusFamily]
# Columns and match ours WHO pathogens to those to then identify hosts
