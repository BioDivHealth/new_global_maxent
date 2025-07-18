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

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load VIRION data from local files
#' @param data_path Path to VIRION data directory
#' @param files Vector of file names to load (default: all main files)
#' @return List containing loaded data frames
load_virion_data <- function(data_path = "data/virion/", 
                            files = c("virion.csv.gz", "edgelist.csv", 
                                     "taxonomy_host.csv", "taxonomy_virus.csv",
                                     "provenance.csv.gz", "detection.csv.gz", 
                                     "temporal.csv.gz")) {
  
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
      
      # Use vroom for compressed files, read_csv for regular files
      if (grepl("\\.gz$", file)) {
        data_list[[gsub("\\.csv\\.gz$", "", file)]] <- vroom(file_path, show_col_types = FALSE)
      } else {
        data_list[[gsub("\\.csv$", "", file)]] <- read_csv(file_path, show_col_types = FALSE)
      }
    } else {
      warning("File not found:", file_path)
    }
  }
  
  return(data_list)
}

#' Load VIRION data using the virionData package
#' @param version Version of VIRION data to load
#' @return List containing loaded data frames
load_virion_package <- function(version = "latest") {
  # Check if virionData package is available
  if (!require(virionData, quietly = TRUE)) {
    cat("virionData package not installed. Installing...\n")
    if (!require(remotes, quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("viralemergence/virionData")
    library(virionData)
  }
  
  # Load data using the package
  virion_data <- list(
    virion = get_virion_data(version = version),
    edgelist = get_edgelist_data(version = version),
    taxonomy_host = get_taxonomy_host_data(version = version),
    taxonomy_virus = get_taxonomy_virus_data(version = version),
    provenance = get_provenance_data(version = version),
    detection = get_detection_data(version = version),
    temporal = get_temporal_data(version = version)
  )
  
  return(virion_data)
}

#' Get available VIRION data versions
#' @return Vector of available versions
get_virion_versions <- function() {
  if (!require(virionData, quietly = TRUE)) {
    stop("virionData package not installed")
  }
  
  # This would depend on the actual virionData package implementation
  # For now, return common versions
  return(c("latest", "0.2.1", "0.2.0"))
}

# =============================================================================
# DATA EXPLORATION FUNCTIONS
# =============================================================================

#' Get summary statistics for VIRION data
#' @param virion_data List of VIRION data frames
#' @return Summary statistics
get_virion_summary <- function(virion_data) {
  
  summary_stats <- list()
  
  # Main virion data summary
  if ("virion" %in% names(virion_data)) {
    virion <- virion_data$virion
    summary_stats$virion <- list(
      total_records = nrow(virion),
      unique_hosts = n_distinct(virion$Host),
      unique_viruses = n_distinct(virion$Virus),
      unique_interactions = n_distinct(paste(virion$Host, virion$Virus, sep = "_")),
      date_range = range(virion$DetectionDate, na.rm = TRUE)
    )
  }
  
  # Taxonomy summaries
  if ("taxonomy_host" %in% names(virion_data)) {
    summary_stats$host_taxonomy <- virion_data$taxonomy_host %>%
      group_by(Class) %>%
      summarise(count = n(), .groups = 'drop') %>%
      arrange(desc(count))
  }
  
  if ("taxonomy_virus" %in% names(virion_data)) {
    summary_stats$virus_taxonomy <- virion_data$taxonomy_virus %>%
      group_by(Family) %>%
      summarise(count = n(), .groups = 'drop') %>%
      arrange(desc(count))
  }
  
  return(summary_stats)
}

#' Find virus families by host class
#' @param virion_data List of VIRION data frames
#' @param host_class Optional host class to filter by
#' @return Data frame with virus family counts by host class
get_virus_families_by_host <- function(virion_data, host_class = NULL) {
  
  # Join virion data with taxonomy
  if (all(c("virion", "taxonomy_host", "taxonomy_virus") %in% names(virion_data))) {
    
    virus_host_data <- virion_data$virion %>%
      left_join(virion_data$taxonomy_host, by = "Host") %>%
      left_join(virion_data$taxonomy_virus, by = "Virus") %>%
      filter(!is.na(Class) & !is.na(Family))
    
    if (!is.null(host_class)) {
      virus_host_data <- virus_host_data %>%
        filter(Class == host_class)
    }
    
    result <- virus_host_data %>%
      group_by(Class, Family) %>%
      summarise(
        interaction_count = n(),
        unique_hosts = n_distinct(Host),
        unique_viruses = n_distinct(Virus),
        .groups = 'drop'
      ) %>%
      arrange(Class, desc(interaction_count))
    
    return(result)
  } else {
    stop("Required data frames not found in virion_data")
  }
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Plot virus family distribution by host class
#' @param virus_host_data Output from get_virus_families_by_host()
#' @param top_n Number of top families to show per class
#' @return ggplot object
plot_virus_families_by_host <- function(virus_host_data, top_n = 10) {
  
  # Get top families per class
  top_families <- virus_host_data %>%
    group_by(Class) %>%
    slice_max(order_by = interaction_count, n = top_n) %>%
    ungroup()
  
  ggplot(top_families, aes(x = reorder(Family, interaction_count), 
                           y = interaction_count, fill = Class)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~Class, scales = "free_y") +
    labs(title = "Top Virus Families by Host Class",
         x = "Virus Family",
         y = "Number of Interactions") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8))
}

#' Plot temporal trends in virus discovery
#' @param virion_data List of VIRION data frames
#' @param by_year Logical, aggregate by year (default: TRUE)
#' @return ggplot object
plot_temporal_trends <- function(virion_data, by_year = TRUE) {
  
  if ("virion" %in% names(virion_data)) {
    
    temporal_data <- virion_data$virion %>%
      filter(!is.na(DetectionDate)) %>%
      mutate(year = lubridate::year(DetectionDate))
    
    if (by_year) {
      temporal_data <- temporal_data %>%
        group_by(year) %>%
        summarise(
          new_interactions = n(),
          new_hosts = n_distinct(Host),
          new_viruses = n_distinct(Virus),
          .groups = 'drop'
        ) %>%
        gather(key = "metric", value = "count", -year)
      
      ggplot(temporal_data, aes(x = year, y = count, color = metric)) +
        geom_line(size = 1) +
        geom_point() +
        labs(title = "Temporal Trends in Virus Discovery",
             x = "Year",
             y = "Count",
             color = "Metric") +
        theme_minimal() +
        scale_color_viridis_d()
    } else {
      # Monthly trends
      temporal_data <- temporal_data %>%
        mutate(month_year = format(DetectionDate, "%Y-%m")) %>%
        group_by(month_year) %>%
        summarise(count = n(), .groups = 'drop') %>%
        mutate(date = as.Date(paste0(month_year, "-01")))
      
      ggplot(temporal_data, aes(x = date, y = count)) +
        geom_line(color = "steelblue") +
        labs(title = "Monthly Virus Discovery Trends",
             x = "Date",
             y = "Number of New Interactions") +
        theme_minimal()
    }
  } else {
    stop("Virion data not found")
  }
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate host-virus network metrics
#' @param virion_data List of VIRION data frames
#' @return Network metrics
calculate_network_metrics <- function(virion_data) {
  
  if ("virion" %in% names(virion_data)) {
    
    # Create interaction matrix
    interactions <- virion_data$virion %>%
      group_by(Host, Virus) %>%
      summarise(interaction_count = n(), .groups = 'drop')
    
    # Calculate basic metrics
    metrics <- list(
      total_interactions = nrow(interactions),
      unique_hosts = n_distinct(interactions$Host),
      unique_viruses = n_distinct(interactions$Virus),
      average_viruses_per_host = n_distinct(interactions$Virus) / n_distinct(interactions$Host),
      average_hosts_per_virus = n_distinct(interactions$Host) / n_distinct(interactions$Virus)
    )
    
    # Host degree distribution
    host_degrees <- interactions %>%
      group_by(Host) %>%
      summarise(degree = n(), .groups = 'drop')
    
    metrics$host_degree_stats <- summary(host_degrees$degree)
    
    # Virus degree distribution
    virus_degrees <- interactions %>%
      group_by(Virus) %>%
      summarise(degree = n(), .groups = 'drop')
    
    metrics$virus_degree_stats <- summary(virus_degrees$degree)
    
    return(metrics)
  } else {
    stop("Virion data not found")
  }
}

#' Find most connected hosts and viruses
#' @param virion_data List of VIRION data frames
#' @param top_n Number of top hosts/viruses to return
#' @return List with top hosts and viruses
find_most_connected <- function(virion_data, top_n = 10) {
  
  if ("virion" %in% names(virion_data)) {
    
    interactions <- virion_data$virion %>%
      group_by(Host, Virus) %>%
      summarise(interaction_count = n(), .groups = 'drop')
    
    # Top hosts
    top_hosts <- interactions %>%
      group_by(Host) %>%
      summarise(
        virus_count = n_distinct(Virus),
        total_interactions = sum(interaction_count),
        .groups = 'drop'
      ) %>%
      arrange(desc(virus_count)) %>%
      head(top_n)
    
    # Top viruses
    top_viruses <- interactions %>%
      group_by(Virus) %>%
      summarise(
        host_count = n_distinct(Host),
        total_interactions = sum(interaction_count),
        .groups = 'drop'
      ) %>%
      arrange(desc(host_count)) %>%
      head(top_n)
    
    return(list(
      top_hosts = top_hosts,
      top_viruses = top_viruses
    ))
  } else {
    stop("Virion data not found")
  }
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Example workflow
example_virion_analysis <- function() {
  
  cat("Loading VIRION data...\n")
  
  # Try to load data (uncomment the method you prefer)
  # Method 1: From local files
  # virion_data <- load_virion_data("data/virion/")
  
  # Method 2: Using virionData package
  # virion_data <- load_virion_package()
  
  # Get summary statistics
  # summary_stats <- get_virion_summary(virion_data)
  # print(summary_stats)
  
  # Find virus families by host class
  # virus_host_data <- get_virus_families_by_host(virion_data)
  # print(head(virus_host_data, 20))
  
  # Create visualizations
  # p1 <- plot_virus_families_by_host(virus_host_data)
  # p2 <- plot_temporal_trends(virion_data)
  
  # Calculate network metrics
  # network_metrics <- calculate_network_metrics(virion_data)
  # print(network_metrics)
  
  # Find most connected species
  # most_connected <- find_most_connected(virion_data)
  # print("Top hosts:")
  # print(most_connected$top_hosts)
  # print("Top viruses:")
  # print(most_connected$top_viruses)
  
  cat("VIRION analysis functions loaded successfully!\n")
  cat("To use, uncomment the example code above and provide your data path.\n")
}

#' Quick start function to demonstrate VIRION capabilities
#' @param use_package Logical, whether to use virionData package (default: TRUE)
#' @return List of analysis results
quick_virion_demo <- function(use_package = TRUE) {
  
  cat("=== VIRION Database Quick Demo ===\n")
  cat("VIRION contains 1,162,000+ host-virus interactions\n")
  cat("across 9,521 viruses and 3,692 hosts\n\n")
  
  if (use_package) {
    cat("Loading data via virionData package...\n")
    virion_data <- load_virion_package()
  } else {
    cat("Loading data from local files...\n")
    virion_data <- load_virion_data()
  }
  
  # Basic summary
  summary_stats <- get_virion_summary(virion_data)
  cat("\n=== Summary Statistics ===\n")
  print(summary_stats$virion)
  
  # Network metrics
  network_metrics <- calculate_network_metrics(virion_data)
  cat("\n=== Network Metrics ===\n")
  cat("Total interactions:", network_metrics$total_interactions, "\n")
  cat("Unique hosts:", network_metrics$unique_hosts, "\n")
  cat("Unique viruses:", network_metrics$unique_viruses, "\n")
  cat("Avg viruses per host:", round(network_metrics$average_viruses_per_host, 2), "\n")
  cat("Avg hosts per virus:", round(network_metrics$average_hosts_per_virus, 2), "\n")
  
  # Most connected species
  most_connected <- find_most_connected(virion_data, top_n = 5)
  cat("\n=== Top 5 Most Connected Hosts ===\n")
  print(most_connected$top_hosts)
  
  cat("\n=== Top 5 Most Connected Viruses ===\n")
  print(most_connected$top_viruses)
  
  return(list(
    summary = summary_stats,
    network = network_metrics,
    most_connected = most_connected
  ))
}

# Run example (commented out by default)
# example_virion_analysis()
# quick_virion_demo()
