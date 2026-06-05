#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 02_extract_local_vector_occurrences.R ----
# -----------------------------------------------------------------------------|
# Purpose: Copy local VectorMap and MapVEu occurrence records for vector
#          species into the SDM occurrence workspace.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# RStudio config: edit this block before sourcing the script ----
# -----------------------------------------------------------------------------|

if (!exists("local_occurrence_config", inherits = FALSE)) {
  local_occurrence_config <- list(
    target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
    output_occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
    output_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "local_vector_occurrence_sources_manifest.csv")
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_local_occurrence_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
  output_occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  output_manifest_path = file.path(
    repo_root(),
    "sdms",
    "runs",
    "vector_sdm_push",
    "local_vector_occurrence_sources_manifest.csv"
  ),
  vectormap_mosquito_path = file.path(
    repo_root(),
    "pathogen_association_data",
    "source_data",
    "vectormap",
    "raw",
    "MosquitoMap2_2627680870621077260.csv"
  ),
  mapveu_sample_path = file.path(
    repo_root(),
    "pathogen_association_data",
    "source_data",
    "mapveu",
    "raw",
    "VBP_MEGA_Sample_subsettedData.txt"
  ),
  mapveu_collection_site_path = file.path(
    repo_root(),
    "pathogen_association_data",
    "source_data",
    "mapveu",
    "raw",
    "VBP_MEGA_Collection site_subsettedData.txt"
  ),
  mapveu_collection_path = file.path(
    repo_root(),
    "pathogen_association_data",
    "source_data",
    "mapveu",
    "raw",
    "VBP_MEGA_Collection_subsettedData.txt"
  )
)

local_occurrence_config <- utils::modifyList(default_local_occurrence_config, local_occurrence_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, local_occurrence_config[[config_key]])
}

local_occurrence_config$target_manifest_path <- config_arg("target-manifest-path")
local_occurrence_config$output_occurrence_root <- config_arg("output-occurrence-root")
local_occurrence_config$output_manifest_path <- config_arg("output-manifest-path")

# -----------------------------------------------------------------------------|
# Normalisation and validation helpers ----
# -----------------------------------------------------------------------------|

canonical_key <- function(x) {
  x <- trimws(tolower(as.character(x)))
  gsub("[[:space:]]+", " ", x)
}

canonical_display <- function(x) {
  x <- canonical_key(x)
  vapply(strsplit(x, " ", fixed = TRUE), function(parts) {
    if (length(parts) == 0 || !nzchar(parts[[1]])) {
      return(NA_character_)
    }
    parts[[1]] <- paste0(toupper(substr(parts[[1]], 1, 1)), substr(parts[[1]], 2, nchar(parts[[1]])))
    paste(parts, collapse = " ")
  }, character(1))
}

require_columns <- function(data, cols, label) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

coord_summary <- function(data, lat_col, lon_col) {
  lat <- suppressWarnings(as.numeric(data[[lat_col]]))
  lon <- suppressWarnings(as.numeric(data[[lon_col]]))
  coord_ok <- !is.na(lat) & !is.na(lon)
  unique_coords <- unique(data.frame(latitude = lat[coord_ok], longitude = lon[coord_ok]))

  list(
    rows = nrow(data),
    coordinate_rows = sum(coord_ok),
    unique_coordinate_rows = nrow(unique_coords)
  )
}

write_species_source <- function(data, species_name, source_method, source_species_col, lat_col, lon_col) {
  # Keep the source table intact, but prepend standard SDM coordinate columns so
  # these files sit next to GBIF raw occurrences with comparable metadata.
  if (nrow(data) == 0) {
    return(NA_character_)
  }

  species_safe <- safe_species_name(species_name)
  output_dir <- ensure_dir(file.path(
    local_occurrence_config$output_occurrence_root,
    species_safe,
    source_method,
    "raw"
  ))
  output_path <- file.path(output_dir, paste0(species_safe, "_", source_method, "_occurrences.csv"))

  source_species <- as.character(data[[source_species_col]])
  lat <- suppressWarnings(as.numeric(data[[lat_col]]))
  lon <- suppressWarnings(as.numeric(data[[lon_col]]))
  standardized <- data.frame(
    source_method = source_method,
    manifest_species_name = species_name,
    source_species_name = source_species,
    decimalLongitude = lon,
    decimalLatitude = lat,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(cbind(standardized, data), output_path, na = "")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}

# -----------------------------------------------------------------------------|
# Targets ----
# -----------------------------------------------------------------------------|

target_manifest <- read.csv(
  local_occurrence_config$target_manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
require_columns(target_manifest, c("species_name", "species_role"), "Target manifest")

target_vectors <- target_manifest[target_manifest$species_role == "vector", , drop = FALSE]
target_vectors$species_key <- canonical_key(target_vectors$species_name)
target_vectors$species_name_display <- canonical_display(target_vectors$species_name)
target_keys <- target_vectors$species_key

# -----------------------------------------------------------------------------|
# Read VectorMap MosquitoMap records ----
# -----------------------------------------------------------------------------|

cat("Reading VectorMap MosquitoMap records...\n")
vectormap <- read.csv(
  local_occurrence_config$vectormap_mosquito_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
)
require_columns(
  vectormap,
  c("ScientificName", "DecimalLatitude", "DecimalLongitude"),
  "VectorMap MosquitoMap"
)
vectormap$species_key <- canonical_key(vectormap$ScientificName)

# -----------------------------------------------------------------------------|
# Read and join MapVEu sample records ----
# -----------------------------------------------------------------------------|

cat("Reading MapVEu sample records...\n")
mapveu_sample <- data.table::fread(
  local_occurrence_config$mapveu_sample_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = FALSE,
  fill = TRUE,
  quote = ""
)
require_columns(
  mapveu_sample,
  c("Sample_ID", "Collection_ID", "Collection_site_ID", "Study_ID", "species [OBI_0001909]"),
  "MapVEu sample"
)
mapveu_sample$species_key <- canonical_key(mapveu_sample[["species [OBI_0001909]"]])
mapveu_sample <- mapveu_sample[mapveu_sample$species_key %in% target_keys, , drop = FALSE]

cat("Reading MapVEu collection-site records...\n")
mapveu_site <- data.table::fread(
  local_occurrence_config$mapveu_collection_site_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = FALSE,
  fill = TRUE,
  quote = ""
)
require_columns(
  mapveu_site,
  c("Collection_site_ID", "Latitude [OBI_0001620]", "Longitude [OBI_0001621]"),
  "MapVEu collection site"
)

cat("Reading MapVEu collection metadata...\n")
mapveu_collection <- data.table::fread(
  local_occurrence_config$mapveu_collection_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = FALSE,
  fill = TRUE,
  quote = ""
)
require_columns(mapveu_collection, c("Collection_ID", "Collection_site_ID", "Study_ID"), "MapVEu collection")
mapveu_collection <- mapveu_collection[
  mapveu_collection$Collection_ID %in% unique(mapveu_sample$Collection_ID),
  ,
  drop = FALSE
]

mapveu <- merge(
  mapveu_sample,
  mapveu_site,
  by = "Collection_site_ID",
  all.x = TRUE,
  suffixes = c("_sample", "_site")
)
mapveu <- merge(
  mapveu,
  mapveu_collection,
  by = "Collection_ID",
  all.x = TRUE,
  suffixes = c("", "_collection")
)

# -----------------------------------------------------------------------------|
# Write per-species source records and manifest rows ----
# -----------------------------------------------------------------------------|

manifest_rows <- list()

for (i in seq_len(nrow(target_vectors))) {
  species_name <- target_vectors$species_name_display[[i]]
  species_key <- target_vectors$species_key[[i]]

  vm_species <- vectormap[vectormap$species_key == species_key, , drop = FALSE]
  vm_path <- write_species_source(
    vm_species,
    species_name,
    "vectormap",
    "ScientificName",
    "DecimalLatitude",
    "DecimalLongitude"
  )
  vm_counts <- coord_summary(vm_species, "DecimalLatitude", "DecimalLongitude")

  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    species_name = species_name,
    source_method = "vectormap",
    source_dataset = "VectorMap MosquitoMap2",
    output_path = vm_path,
    rows = vm_counts$rows,
    coordinate_rows = vm_counts$coordinate_rows,
    unique_coordinate_rows = vm_counts$unique_coordinate_rows,
    source_species_column = "ScientificName",
    latitude_column = "DecimalLatitude",
    longitude_column = "DecimalLongitude",
    source_path = normalizePath(local_occurrence_config$vectormap_mosquito_path, winslash = "/", mustWork = TRUE),
    extracted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  mv_species <- mapveu[mapveu$species_key == species_key, , drop = FALSE]
  mv_path <- write_species_source(
    mv_species,
    species_name,
    "mapveu",
    "species [OBI_0001909]",
    "Latitude [OBI_0001620]",
    "Longitude [OBI_0001621]"
  )
  mv_counts <- coord_summary(mv_species, "Latitude [OBI_0001620]", "Longitude [OBI_0001621]")

  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    species_name = species_name,
    source_method = "mapveu",
    source_dataset = "MapVEu sample + collection site + collection metadata",
    output_path = mv_path,
    rows = mv_counts$rows,
    coordinate_rows = mv_counts$coordinate_rows,
    unique_coordinate_rows = mv_counts$unique_coordinate_rows,
    source_species_column = "species [OBI_0001909]",
    latitude_column = "Latitude [OBI_0001620]",
    longitude_column = "Longitude [OBI_0001621]",
    source_path = paste(
      normalizePath(local_occurrence_config$mapveu_sample_path, winslash = "/", mustWork = TRUE),
      normalizePath(local_occurrence_config$mapveu_collection_site_path, winslash = "/", mustWork = TRUE),
      normalizePath(local_occurrence_config$mapveu_collection_path, winslash = "/", mustWork = TRUE),
      sep = "; "
    ),
    extracted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  cat(
    species_name,
    "- VectorMap rows:", vm_counts$rows,
    "MapVEu rows:", mv_counts$rows,
    "\n"
  )
}

# -----------------------------------------------------------------------------|
# Write combined local-source manifest ----
# -----------------------------------------------------------------------------|

manifest <- do.call(rbind, manifest_rows)
invisible(ensure_dir(dirname(local_occurrence_config$output_manifest_path)))
write.csv(manifest, local_occurrence_config$output_manifest_path, row.names = FALSE, na = "")

cat("Wrote local vector occurrence manifest:", local_occurrence_config$output_manifest_path, "\n")
