#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_build_vector_sdm_target_manifests.R ----
# -----------------------------------------------------------------------------|
# Purpose: Build disease-vector and unique vector-species manifests for the
#          present-day vector SDM push.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# RStudio config: edit this block before sourcing the script ----
# -----------------------------------------------------------------------------|

if (!exists("manifest_config", inherits = FALSE)) {
  manifest_config <- list(
    recommended_next_action = "find_or_build_vector_sdm",
    roles = "vector",
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y")),
    overwrite = TRUE
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_manifest_config <- list(
  pilot_dir = file.path(repo_root(), "pathogen_association_data", "readiness", "disease_modelling_pilot_package"),
  output_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push"),
  recommended_next_action = "find_or_build_vector_sdm",
  roles = "vector",
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  occurrence_method = "combined",
  gbif_occurrence_method = "gbif-download",
  sdm_needed_for_disease = "yes",
  default_run_status = "pending",
  min_points = 20,
  overwrite = TRUE
)

manifest_config <- utils::modifyList(default_manifest_config, manifest_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, manifest_config[[config_key]])
}

# -----------------------------------------------------------------------------|
# Helpers ----
# -----------------------------------------------------------------------------|

require_columns <- function(data, cols, label) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

optional_columns <- function(data, cols) {
  present <- intersect(cols, names(data))
  data[, present, drop = FALSE]
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

collapse_integer <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(sort(unique(x)), collapse = "; ")
}

rank_species <- function(data) {
  evidence_rank <- match(
    tolower(as.character(data$best_evidence_level)),
    c("confirmed", "probable", "candidate", "weak")
  )
  evidence_rank[is.na(evidence_rank)] <- 99L

  competence_rank <- match(
    tolower(as.character(data$vector_competence_status)),
    c("competent", "mixed", "unclear", "not_competent")
  )
  competence_rank[is.na(competence_rank)] <- 99L

  bites_humans_rank <- ifelse(tolower(as.character(data$bites_humans)) %in% c("true", "yes", "1"), 0L, 1L)
  disease_count_rank <- -suppressWarnings(as.integer(data$disease_count))

  order(evidence_rank, competence_rank, bites_humans_rank, disease_count_rank, data$species_name)
}

# -----------------------------------------------------------------------------|
# Resolve paths and config ----
# -----------------------------------------------------------------------------|

pilot_dir <- config_arg("pilot-dir")
output_root <- ensure_dir(config_arg("output-root"))
recommended_next_action <- split_arg(config_arg("recommended-next-action"))
roles <- split_arg(config_arg("roles"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
occurrence_method <- config_arg("occurrence-method")
gbif_occurrence_method <- config_arg("gbif-occurrence-method")
sdm_needed_for_disease <- config_arg("sdm-needed-for-disease")
default_run_status <- config_arg("default-run-status")
min_points <- as.integer(config_arg("min-points"))
overwrite <- as_logical_arg(config_arg("overwrite")) || has_flag(args, "overwrite")

pilot_evidence_summary_path <- file.path(pilot_dir, "pilot_evidence_summary.csv")
pilot_sdm_species_path <- file.path(pilot_dir, "pilot_sdm_species.csv")
pilot_vectors_path <- file.path(pilot_dir, "pilot_vectors.csv")

disease_targets_path <- file.path(output_root, "disease_vector_sdm_targets.csv")
species_targets_path <- file.path(output_root, "vector_species_sdm_targets.csv")
summary_path <- file.path(output_root, "manifest_build_summary.csv")

if (!overwrite && any(file.exists(c(disease_targets_path, species_targets_path, summary_path)))) {
  stop("Output manifests already exist. Set overwrite = TRUE to replace them.", call. = FALSE)
}

for (path in c(pilot_evidence_summary_path, pilot_sdm_species_path, pilot_vectors_path)) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path, call. = FALSE)
  }
}

# -----------------------------------------------------------------------------|
# Load readiness inputs ----
# -----------------------------------------------------------------------------|

evidence_summary <- read.csv(pilot_evidence_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
sdm_species <- read.csv(pilot_sdm_species_path, check.names = FALSE, stringsAsFactors = FALSE)
pilot_vectors <- read.csv(pilot_vectors_path, check.names = FALSE, stringsAsFactors = FALSE)

require_columns(
  evidence_summary,
  c("analysis_unit_id", "readiness_disease_name", "recommended_next_action"),
  "pilot_evidence_summary.csv"
)
require_columns(
  sdm_species,
  c("analysis_unit_id", "readiness_disease_name", "disease_name", "species_name", "species_role", "sdm_needed_for_disease", "sdm_available", "sdm_species"),
  "pilot_sdm_species.csv"
)
require_columns(
  pilot_vectors,
  c("analysis_unit_id", "readiness_disease_name", "disease_name", "species_name"),
  "pilot_vectors.csv"
)

target_diseases <- evidence_summary[
  evidence_summary$recommended_next_action %in% recommended_next_action,
  ,
  drop = FALSE
]

if (nrow(target_diseases) == 0) {
  stop("No diseases matched recommended_next_action: ", paste(recommended_next_action, collapse = ", "), call. = FALSE)
}

target_ids <- unique(target_diseases$analysis_unit_id)

disease_context_cols <- c(
  "analysis_unit_id",
  "readiness_disease_name",
  "recommended_next_action",
  "readiness_blocker",
  "range_limiting_layer",
  "roster_vector_rows",
  "vector_sdm_species_available",
  "direct_vector_species_or_taxa",
  "direct_vector_confirmed_rows",
  "direct_vector_probable_rows",
  "direct_vector_candidate_rows",
  "direct_vector_competent_rows",
  "direct_vector_transmission_yes_rows",
  "direct_vector_natural_infection_yes_rows",
  "genbank_distinct_countries_or_territories",
  "who_don_distinct_countries"
)
disease_context <- optional_columns(target_diseases, disease_context_cols)

sdm_targets <- sdm_species[
  sdm_species$analysis_unit_id %in% target_ids &
    sdm_species$species_role %in% roles,
  ,
  drop = FALSE
]

if (nrow(sdm_targets) == 0) {
  stop("No SDM species rows matched selected diseases and roles.", call. = FALSE)
}

vector_cols <- c(
  "analysis_unit_id",
  "readiness_disease_name",
  "disease_name",
  "species_name",
  "tax_id",
  "vector_group",
  "vector_taxon_rank",
  "vector_join_key",
  "has_disease_vector_evidence",
  "has_host_vector_evidence",
  "has_competence_evidence",
  "best_evidence_level",
  "best_evidence_basis",
  "vector_record_sources",
  "bites_humans",
  "bites_humans_basis",
  "vector_competence_status",
  "transmission_demonstrated",
  "natural_infection_reported",
  "vector_role_hint",
  "uncertainty_reason",
  "taxonomy_caution"
)
vector_context <- optional_columns(pilot_vectors, vector_cols)

disease_targets <- merge(
  sdm_targets,
  disease_context,
  by = c("analysis_unit_id", "readiness_disease_name"),
  all.x = TRUE,
  sort = FALSE
)
disease_targets <- merge(
  disease_targets,
  vector_context,
  by = c("analysis_unit_id", "readiness_disease_name", "disease_name", "species_name"),
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_vector")
)

disease_targets$species_name <- canonical_species_name(disease_targets$species_name)
disease_targets$species_name_canonical <- disease_targets$species_name
disease_targets$species_key <- tolower(disease_targets$species_name_canonical)
disease_targets$run_status <- default_run_status
disease_targets$run_priority <- seq_len(nrow(disease_targets))
disease_targets$occurrence_method <- occurrence_method
disease_targets$gbif_occurrence_method <- gbif_occurrence_method
disease_targets$start_year <- start_year
disease_targets$end_year <- end_year
disease_targets$min_points <- min_points
disease_targets$target_status <- ifelse(
  tolower(as.character(disease_targets$sdm_available)) %in% c("true", "yes", "1"),
  "sdm_already_available",
  "needs_vector_sdm"
)

disease_targets <- disease_targets[
  order(disease_targets$readiness_disease_name, disease_targets$species_name_canonical),
  ,
  drop = FALSE
]
disease_targets$run_priority <- seq_len(nrow(disease_targets))

# -----------------------------------------------------------------------------|
# Collapse to one operational row per vector species ----
# -----------------------------------------------------------------------------|

species_groups <- split(disease_targets, disease_targets$species_key)
species_rows <- lapply(species_groups, function(group) {
  data.frame(
    species_name = first_non_missing(group$species_name_canonical),
    species_name_canonical = first_non_missing(group$species_name_canonical),
    species_role = "vector",
    sdm_needed_for_disease = sdm_needed_for_disease,
    sdm_available = any(tolower(as.character(group$sdm_available)) %in% c("true", "yes", "1"), na.rm = TRUE),
    sdm_species = collapse_unique(group$sdm_species),
    run_status = default_run_status,
    run_priority = NA_integer_,
    occurrence_method = occurrence_method,
    gbif_occurrence_method = gbif_occurrence_method,
    start_year = start_year,
    end_year = end_year,
    min_points = min_points,
    disease_count = length(unique(group$analysis_unit_id)),
    disease_names = collapse_unique(group$readiness_disease_name),
    analysis_unit_ids = collapse_unique(group$analysis_unit_id),
    tax_ids = collapse_integer(group$tax_id),
    vector_group = first_non_missing(group$vector_group),
    vector_taxon_rank = first_non_missing(group$vector_taxon_rank),
    best_evidence_level = first_non_missing(group$best_evidence_level),
    best_evidence_basis = collapse_unique(group$best_evidence_basis),
    vector_record_sources = collapse_unique(group$vector_record_sources),
    bites_humans = first_non_missing(group$bites_humans),
    bites_humans_basis = collapse_unique(group$bites_humans_basis),
    vector_competence_status = first_non_missing(group$vector_competence_status),
    transmission_demonstrated = first_non_missing(group$transmission_demonstrated),
    natural_infection_reported = first_non_missing(group$natural_infection_reported),
    uncertainty_reason = collapse_unique(group$uncertainty_reason),
    taxonomy_caution = collapse_unique(group$taxonomy_caution),
    occurrence_status = "not_prepared",
    occurrence_rows_raw = NA_integer_,
    occurrence_rows_clean = NA_integer_,
    passes_min_points = NA,
    gbif_occurrence_path = NA_character_,
    local_occurrence_path = NA_character_,
    combined_occurrence_path = NA_character_,
    model_status = "not_run",
    model_output_path = NA_character_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
})
species_targets <- do.call(rbind, species_rows)
species_targets <- species_targets[rank_species(species_targets), , drop = FALSE]
species_targets$run_priority <- seq_len(nrow(species_targets))

summary <- data.frame(
  analysis_unit_id = target_diseases$analysis_unit_id,
  readiness_disease_name = target_diseases$readiness_disease_name,
  recommended_next_action = target_diseases$recommended_next_action,
  roster_vector_rows = suppressWarnings(as.integer(target_diseases$roster_vector_rows)),
  disease_vector_target_rows = vapply(
    target_diseases$analysis_unit_id,
    function(id) sum(disease_targets$analysis_unit_id == id, na.rm = TRUE),
    integer(1)
  ),
  unique_vector_species = vapply(
    target_diseases$analysis_unit_id,
    function(id) length(unique(disease_targets$species_key[disease_targets$analysis_unit_id == id])),
    integer(1)
  ),
  vector_sdm_species_available = suppressWarnings(as.integer(target_diseases$vector_sdm_species_available)),
  range_limiting_layer = target_diseases$range_limiting_layer,
  manifest_built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------|
# Write outputs ----
# -----------------------------------------------------------------------------|

write.csv(disease_targets, disease_targets_path, row.names = FALSE, na = "")
write.csv(species_targets, species_targets_path, row.names = FALSE, na = "")
write.csv(summary, summary_path, row.names = FALSE, na = "")

cat("Wrote disease-vector manifest:", disease_targets_path, "\n")
cat("Wrote vector species manifest:", species_targets_path, "\n")
cat("Wrote manifest build summary:", summary_path, "\n")
cat("Disease-vector rows:", nrow(disease_targets), "\n")
cat("Unique vector species:", nrow(species_targets), "\n")
