#!/usr/bin/env Rscript
################################################################################
# 02_build_modelling_evidence_tiers_handoff.R
################################################################################
# Purpose: Build a read-only prototype evidence-tier surface for modelling.
#
# Output : pathogen_association_data/readiness/evidence_tiers/
#
# Notes  : This script derives a compact modelling-readiness handoff surface
#          from the current readiness package. It does not edit canonical
#          evidence, source-check ledgers, role assignments, or SDM manifests.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tibble, tidyr)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "readiness",
  "helpers",
  "disease_modelling_readiness_helpers.R"
))

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_package_dir <- file.path(readiness_dir, "disease_modelling_pilot_package")
output_dir <- file.path(readiness_dir, "evidence_tiers")
chikungunya_delivery_dir <- Sys.getenv(
  "CHIKUNGUNYA_DELIVERY_DIR",
  unset = "/Volumes/LaCie/new_global_maxent/sdms/delivery/chikungunya_vector_sdm_delivery_20260609"
)
delivery_package_dir <- file.path(chikungunya_delivery_dir, "readiness", "disease_modelling_pilot_package")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  repo_hosts = file.path(repo_package_dir, "pilot_hosts.csv"),
  repo_vectors = file.path(repo_package_dir, "pilot_vectors.csv"),
  repo_sdm_species = file.path(repo_package_dir, "pilot_sdm_species.csv"),
  role_modelling_features = role_modelling_features_path(),
  vector_modelling_features = vector_modelling_features_path(),
  delivery_hosts = file.path(delivery_package_dir, "pilot_hosts.csv"),
  delivery_vectors = file.path(delivery_package_dir, "pilot_vectors.csv"),
  delivery_sdm_species = file.path(delivery_package_dir, "pilot_sdm_species.csv"),
  delivery_model_qc = file.path(chikungunya_delivery_dir, "model_qc_summary.csv")
)

# ------------------------------------------------------------------------------|
#      Small Helpers -----------------------------------------------------------|
# ------------------------------------------------------------------------------|
read_optional_csv <- function(path) {
  read_csv_layer(path, required = FALSE)
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(x))
}

delivery_available <- function() {
  all(file.exists(c(paths$delivery_hosts, paths$delivery_vectors, paths$delivery_sdm_species)))
}

build_sdm_lookup <- function(sdm_species) {
  if (is.null(sdm_species)) {
    return(tibble(
      analysis_unit_id = character(),
      species_name = character(),
      species_role = character(),
      sdm_needed_for_disease = character(),
      sdm_available = logical(),
      sdm_species = character()
    ))
  }

  sdm_species %>%
    mutate(
      species_role = clean_text(species_role),
      sdm_available = missing_as_false(sdm_available)
    ) %>%
    select(any_of(c(
      "analysis_unit_id",
      "species_name",
      "species_role",
      "sdm_needed_for_disease",
      "sdm_available",
      "sdm_species"
    ))) %>%
    distinct(analysis_unit_id, species_name, species_role, .keep_all = TRUE)
}

add_sdm_fields <- function(data, sdm_lookup) {
  joined <- data %>%
    left_join(
      sdm_lookup,
      by = c("analysis_unit_id", "species_name", "species_role"),
      suffix = c("", "_sdm")
    )

  sdm_cols <- c(
    "sdm_needed_for_disease",
    "sdm_needed_for_disease_sdm",
    "sdm_available",
    "sdm_available_sdm",
    "sdm_species",
    "sdm_species_sdm"
  )
  for (col in sdm_cols) {
    if (!col %in% names(joined)) {
      joined[[col]] <- NA
    }
  }

  joined %>%
    mutate(
      sdm_needed_for_disease = first_non_empty(sdm_needed_for_disease, sdm_needed_for_disease_sdm),
      sdm_available = missing_as_false(first_non_empty(sdm_available, sdm_available_sdm)),
      sdm_species = first_non_empty(sdm_species, sdm_species_sdm)
    ) %>%
    select(-any_of(c("sdm_needed_for_disease_sdm", "sdm_available_sdm", "sdm_species_sdm")))
}

add_model_quality <- function(data, model_qc = NULL) {
  if (!is.null(model_qc)) {
    qc <- model_qc %>%
      select(any_of(c(
        "species_name",
        "species_role",
        "model_quality",
        "retained_models",
        "min_boyce",
        "max_boyce",
        "min_test_auc",
        "max_test_auc",
        "min_max_tss",
        "max_max_tss",
        "prediction_tif_path"
      ))) %>%
      distinct(species_name, species_role, .keep_all = TRUE)

    data <- data %>%
      left_join(qc, by = c("species_name", "species_role"))
  }

  model_cols <- c(
    "model_quality",
    "retained_models",
    "min_boyce",
    "max_boyce",
    "min_test_auc",
    "max_test_auc",
    "min_max_tss",
    "max_max_tss",
    "prediction_tif_path"
  )
  for (col in model_cols) {
    if (!col %in% names(data)) {
      data[[col]] <- NA
    }
  }

  data %>%
    mutate(
      retained_models = as_numeric_clean(retained_models),
      min_boyce = as_numeric_clean(min_boyce),
      model_quality_raw = clean_text(model_quality),
      model_quality_tier = case_when(
        !sdm_available ~ "no_sdm",
        str_detect(model_quality_raw, regex("^diagnostic", ignore_case = TRUE)) ~ "diagnostic",
        str_detect(model_quality_raw, regex("production", ignore_case = TRUE)) &
          !is.na(min_boyce) & min_boyce >= 0.5 &
          !is.na(retained_models) & retained_models >= 5 ~ "usable_sdm",
        str_detect(model_quality_raw, regex("existing_host_model", ignore_case = TRUE)) ~ "existing_host_model",
        sdm_available ~ "available_unscored",
        TRUE ~ "no_sdm"
      ),
      model_quality_weight = case_when(
        model_quality_tier == "usable_sdm" ~ 1,
        model_quality_tier == "existing_host_model" ~ 0.75,
        model_quality_tier == "available_unscored" ~ 0.5,
        model_quality_tier == "diagnostic" ~ 0.25,
        TRUE ~ 0
      )
    )
}

feature_join_key <- function(disease_name, species_name, tax_id) {
  paste(
    clean_key(disease_name),
    clean_key(species_name),
    clean_key(tax_id),
    sep = "|"
  )
}

build_role_feature_lookup <- function(role_features) {
  feature_cols <- c(
    "role_modelling_feature_id",
    "disease_name",
    "species_name",
    "tax_id",
    "modelling_role_proxy",
    "modelling_role_proxy_basis",
    "modelling_role_proxy_confidence",
    "modelling_role_proxy_rule_id",
    "modelling_role_proxy_needs_review",
    "host_role_bucket",
    "host_role_evidence_basis",
    "host_role_weight",
    "role_proxy_applied",
    "group_proxy_applied",
    "profile_group_proxy",
    "taxonomy_ok",
    "host_detection_tier",
    "host_direct_detection_supported",
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "biological_evidence_tier",
    "tier_rule_id",
    "host_evidence_missingness_reason"
  )
  missing_cols <- setdiff(feature_cols, names(role_features))
  if (length(missing_cols) > 0) {
    stop(
      "Role modelling feature table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  role_features %>%
    filter(species_role == "host") %>%
    select(all_of(feature_cols)) %>%
    mutate(
      .feature_join_key = feature_join_key(
        disease_name,
        species_name,
        tax_id
      )
    ) %>%
    distinct(.feature_join_key, .keep_all = TRUE) %>%
    select(-disease_name, -species_name, -tax_id)
}

build_vector_feature_lookup <- function(vector_features) {
  feature_cols <- c(
    "vector_modelling_feature_id",
    "disease_name",
    "species_name",
    "tax_id",
    "taxonomy_ok",
    "bites_humans_known",
    "bites_humans_true",
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "biological_evidence_tier",
    "tier_rule_id",
    "vector_evidence_missingness_reason"
  )
  missing_cols <- setdiff(feature_cols, names(vector_features))
  if (length(missing_cols) > 0) {
    stop(
      "Vector modelling feature table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  vector_features %>%
    filter(species_role == "vector") %>%
    select(all_of(feature_cols)) %>%
    mutate(
      .feature_join_key = feature_join_key(
        disease_name,
        species_name,
        tax_id
      )
    ) %>%
    distinct(.feature_join_key, .keep_all = TRUE) %>%
    select(-disease_name, -species_name, -tax_id)
}

add_role_feature_fields <- function(data, role_feature_lookup) {
  joined <- data %>%
    mutate(
      .feature_join_key = feature_join_key(
        disease_name,
        species_name,
        tax_id
      )
    ) %>%
    left_join(role_feature_lookup, by = ".feature_join_key")

  unmatched <- joined %>%
    filter(is.na(role_modelling_feature_id)) %>%
    transmute(source_dataset, analysis_unit_id, disease_name, species_name, tax_id)

  if (nrow(unmatched) > 0) {
    preview <- unmatched %>%
      mutate(label = paste(source_dataset, analysis_unit_id, disease_name, species_name, tax_id, sep = " | ")) %>%
      pull(label) %>%
      head(10)

    stop(
      "Role modelling feature table does not cover all host handoff rows: ",
      paste(preview, collapse = "; "),
      call. = FALSE
    )
  }

  joined %>%
    select(-.feature_join_key, -role_modelling_feature_id)
}

add_vector_feature_fields <- function(data, vector_feature_lookup) {
  joined <- data %>%
    mutate(
      .feature_join_key = feature_join_key(
        disease_name,
        species_name,
        tax_id
      )
    ) %>%
    left_join(vector_feature_lookup, by = ".feature_join_key")

  unmatched <- joined %>%
    filter(is.na(vector_modelling_feature_id)) %>%
    transmute(source_dataset, analysis_unit_id, disease_name, species_name, tax_id)

  if (nrow(unmatched) > 0) {
    preview <- unmatched %>%
      mutate(label = paste(source_dataset, analysis_unit_id, disease_name, species_name, tax_id, sep = " | ")) %>%
      pull(label) %>%
      head(10)

    stop(
      "Vector modelling feature table does not cover all vector handoff rows: ",
      paste(preview, collapse = "; "),
      call. = FALSE
    )
  }

  joined %>%
    select(-.feature_join_key, -vector_modelling_feature_id)
}

derive_host_tiers <- function(hosts, sdm_lookup, source_dataset, role_feature_lookup,
                              model_qc = NULL) {
  hosts %>%
    mutate(
      source_dataset = source_dataset,
      species_role = "host"
    ) %>%
    add_sdm_fields(sdm_lookup) %>%
    add_model_quality(model_qc) %>%
    add_role_feature_fields(role_feature_lookup) %>%
    mutate(
      no_spatial_layer = !sdm_available,
      missingness_reason = pmap_chr(
        list(
          host_evidence_missingness_reason,
          if_else(sdm_available, NA_character_, "no_sdm")
        ),
        combine_reasons
      )
    ) %>%
    select(-host_evidence_missingness_reason)
}

derive_vector_tiers <- function(vectors, sdm_lookup, source_dataset, vector_feature_lookup,
                                model_qc = NULL) {
  vectors %>%
    mutate(
      source_dataset = source_dataset,
      species_role = "vector"
    ) %>%
    add_sdm_fields(sdm_lookup) %>%
    add_model_quality(model_qc) %>%
    add_vector_feature_fields(vector_feature_lookup) %>%
    mutate(
      no_spatial_layer = !sdm_available,
      missingness_reason = pmap_chr(
        list(
          vector_evidence_missingness_reason,
          if_else(sdm_available, NA_character_, "no_sdm")
        ),
        combine_reasons
      )
    ) %>%
    select(-vector_evidence_missingness_reason)
}

build_host_role_bucket_counts <- function(tiered_rows) {
  tiered_rows %>%
    filter(species_role == "host") %>%
    group_by(
      source_dataset,
      analysis_unit_id,
      readiness_disease_name,
      host_role_bucket,
      host_role_evidence_basis,
      modelling_role_proxy_confidence,
      modelling_role_proxy_needs_review
    ) %>%
    summarise(rows = n(), .groups = "drop") %>%
    arrange(source_dataset, readiness_disease_name, desc(rows), host_role_bucket)
}

write_table <- function(data, filename, dir = output_dir) {
  path <- file.path(dir, filename)
  write_csv(data, path, na = "")
  path
}

csv_row_count <- function(path) {
  max(length(readLines(path, warn = FALSE)) - 1, 0)
}

csv_column_count <- function(path) {
  length(names(readr::read_csv(path, n_max = 0, show_col_types = FALSE)))
}

# ------------------------------------------------------------------------------|
#      Inputs ------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_hosts <- read_optional_csv(paths$repo_hosts)
repo_vectors <- read_optional_csv(paths$repo_vectors)
repo_sdm <- read_optional_csv(paths$repo_sdm_species)
role_features <- read_csv_layer(paths$role_modelling_features, required = TRUE)
vector_features <- read_csv_layer(paths$vector_modelling_features, required = TRUE)

if (is.null(repo_hosts) || is.null(repo_vectors) || is.null(repo_sdm)) {
  stop("Repository pilot package inputs are missing under: ", repo_package_dir, call. = FALSE)
}

role_feature_lookup <- build_role_feature_lookup(role_features)
vector_feature_lookup <- build_vector_feature_lookup(vector_features)

delivery_hosts <- if (delivery_available()) read_optional_csv(paths$delivery_hosts) else NULL
delivery_vectors <- if (delivery_available()) read_optional_csv(paths$delivery_vectors) else NULL
delivery_sdm <- if (delivery_available()) read_optional_csv(paths$delivery_sdm_species) else NULL
delivery_model_qc <- read_optional_csv(paths$delivery_model_qc)

# ------------------------------------------------------------------------------|
#      Tier Derivation ---------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_sdm_lookup <- build_sdm_lookup(repo_sdm)

repo_tiered <- bind_rows(
  derive_host_tiers(repo_hosts, repo_sdm_lookup, "repo_pilot", role_feature_lookup),
  derive_vector_tiers(repo_vectors, repo_sdm_lookup, "repo_pilot", vector_feature_lookup)
)

delivery_tiered <- tibble()
if (!is.null(delivery_hosts) && !is.null(delivery_vectors) && !is.null(delivery_sdm)) {
  delivery_sdm_lookup <- build_sdm_lookup(delivery_sdm)
  delivery_tiered <- bind_rows(
    derive_host_tiers(
      delivery_hosts,
      delivery_sdm_lookup,
      "chikungunya_delivery",
      role_feature_lookup,
      delivery_model_qc
    ),
    derive_vector_tiers(
      delivery_vectors,
      delivery_sdm_lookup,
      "chikungunya_delivery",
      vector_feature_lookup,
      delivery_model_qc
    )
  )
}

tiered_rows <- bind_rows(repo_tiered, delivery_tiered) %>%
  select(any_of(c(
    "source_dataset",
    "analysis_unit_id",
    "readiness_disease_name",
    "disease_name",
    "species_role",
    "species_name",
    "tax_id",
    "host_class",
    "host_order",
    "host_family",
    "host_detection_method",
    "host_detection_tier",
    "host_direct_detection_supported",
    "host_role_assignment",
    "host_role_confidence",
    "host_role_assignment_status",
    "host_role_needs_manual_review",
    "modelling_role_proxy",
    "modelling_role_proxy_basis",
    "modelling_role_proxy_confidence",
    "modelling_role_proxy_rule_id",
    "modelling_role_proxy_needs_review",
    "host_role_bucket",
    "host_role_evidence_basis",
    "host_role_weight",
    "role_proxy_applied",
    "group_proxy_applied",
    "vector_group",
    "vector_taxon_rank",
    "vector_join_key",
    "best_evidence_level",
    "best_evidence_basis",
    "has_disease_vector_evidence",
    "has_host_vector_evidence",
    "has_competence_evidence",
    "bites_humans",
    "bites_humans_known",
    "bites_humans_true",
    "vector_competence_status",
    "transmission_demonstrated",
    "natural_infection_reported",
    "vector_role_hint",
    "taxonomy_ok",
    "sdm_needed_for_disease",
    "sdm_available",
    "sdm_species",
    "model_quality_raw",
    "model_quality_tier",
    "model_quality_weight",
    "retained_models",
    "min_boyce",
    "max_boyce",
    "min_test_auc",
    "max_test_auc",
    "min_max_tss",
    "max_max_tss",
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "profile_group_proxy",
    "biological_evidence_tier",
    "no_spatial_layer",
    "missingness_reason",
    "tier_rule_id"
  ))) %>%
  arrange(desc(species_role == "vector"), source_dataset, readiness_disease_name, species_name)

host_role_bucket_counts <- build_host_role_bucket_counts(tiered_rows)

# ------------------------------------------------------------------------------|
#      Minimal Writes ----------------------------------------------------------|
# ------------------------------------------------------------------------------|
outputs <- list(
  tiered_species = write_table(tiered_rows, "tiered_species.csv"),
  host_role_bucket_counts = write_table(host_role_bucket_counts, "host_role_bucket_counts.csv")
)

manifest <- enframe(outputs, name = "table_name", value = "path") %>%
  mutate(
    relative_path = stringr::str_remove(path, paste0("^", stringr::fixed(here::here()), "/?")),
    row_count = map_int(path, csv_row_count),
    column_count = map_int(path, csv_column_count),
    description = case_when(
      table_name == "tiered_species" ~ "One row per readiness species/taxon with derived tier/profile fields.",
      table_name == "host_role_bucket_counts" ~ "Counts of broad modelling-facing host-role buckets and evidence bases.",
      TRUE ~ "Derived prototype output."
    )
  ) %>%
  select(table_name, relative_path, row_count, column_count, description)

invisible(write_table(manifest, "manifest.csv"))

readme <- c(
  "# Modelling Evidence Tiers",
  "",
  "Generated by `scripts/associations/readiness/02_build_modelling_evidence_tiers_handoff.R`.",
  "",
  "These files are prototype modelling/readiness outputs. They do not replace",
  "canonical role evidence, vector evidence, source-check decisions, or SDM",
  "manifests.",
  "",
  "Current scope: evidence-tier handoff tables. SDM/model-quality fields are",
  "retained as contextual metadata for later.",
  "",
  "Host rows consume `role_modelling_features.csv`, which owns host role",
  "proxy/bucket fields, evidence-tier flags, biological evidence tiers, and",
  "host evidence missingness reasons before readiness adds SDM availability",
  "and model-quality fields.",
  "",
  "Vector rows consume `vector_modelling_features.csv`, which owns vector",
  "evidence-tier flags and biological evidence tiers before readiness adds SDM",
  "availability and model-quality fields.",
  "",
  "H5N1/H7N9 avian-influenza host rows and West Nile fever bird rows include",
  "derived modelling-role proxies. These proxy roles are for sensitivity",
  "modelling and review triage; they are not source-backed species-level role",
  "assignments.",
  "",
  "## Storage Rules",
  "",
  "- Derived modelling-readiness handoff outputs live in this folder.",
  "- Generated role-curation candidates live under `pathogen_association_data/staged/role_annotation/source_check_candidates/`.",
  "- Human review decisions should live under `pathogen_association_data/manual/`.",
  "- Canonical accepted role evidence and generated role-owned modelling features remain under `pathogen_association_data/evidence/role_annotation/`.",
  "- Spatial admin extracts should later live under `pathogen_association_data/readiness/admin_modelling/`.",
  "",
  "## Main Files",
  "",
  "- `tiered_species.csv`: derived rows for the repository pilot package and, when available, the Chikungunya delivery bundle.",
  "- `host_role_bucket_counts.csv`: counts for compact host-role buckets by evidence basis.",
  "- `manifest.csv`: row/column counts for the readiness files in this folder.",
  "",
  "## Inline Contract",
  "",
  "- `biological_evidence_tier`: `excluded`, `broad`, `supported`, `strong`, or `strict`.",
  "- `model_quality_tier`: `no_sdm`, `diagnostic`, `available_unscored`, `existing_host_model`, or `usable_sdm`.",
  "- `host_role_bucket`: `reservoir_or_amplifying_host`, `dead_end_or_incidental_host`, `susceptible_or_spillover_host`, `host_presence_only`, or `unknown_or_unreviewed`.",
  "- `host_role_evidence_basis`: `exact_source_backed`, `exact_reviewed_needs_review`, `disease_group_proxy`, `weighted_taxonomic_proxy`, or `candidate_only`.",
  "- `host_role_weight`: draft role weight derived from host-role confidence and bucket; it is separate from SDM/model-quality weight.",
  "",
  "Other exploratory counts can be regenerated from `tiered_species.csv` when",
  "needed; they are intentionally not written as separate files in this minimal",
  "handoff surface."
)
writeLines(readme, file.path(output_dir, "README.md"))

message("Wrote modelling evidence-tier prototype outputs to: ", output_dir)
message("Rows in tiered_species.csv: ", nrow(tiered_rows))
