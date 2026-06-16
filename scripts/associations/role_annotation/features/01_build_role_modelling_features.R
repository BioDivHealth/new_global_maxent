#!/usr/bin/env Rscript
################################################################################
# 01_build_role_modelling_features.R
################################################################################
# Purpose: Build modelling-facing role features from the role-annotation roster
#          and explicit proxy rules.
#
# Outputs: pathogen_association_data/evidence/role_annotation/
#            role_modelling_features.csv
#            vector_modelling_features.csv
#          pathogen_association_data/evidence/role_annotation/qa/
#            role_modelling_feature_summary.csv
#            vector_modelling_feature_summary.csv
#
# Notes  : This script does not edit accepted evidence or assignment tables.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tibble)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "association_data_helpers.R"
))
source(here::here(
  "scripts",
  "associations",
  "role_annotation",
  "features",
  "rules",
  "host_proxy_rules.R"
))

# ------------------------------------------------------------------------------|
paths <- list(
  roster = file.path(role_roster_dir, "species_host_vector_roster.csv"),
  role_modelling_features = role_modelling_features_path(),
  feature_summary = role_modelling_feature_summary_path(),
  vector_modelling_features = vector_modelling_features_path(),
  vector_feature_summary = vector_modelling_feature_summary_path()
)

dir.create(dirname(paths$role_modelling_features), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$feature_summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$vector_modelling_features), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$vector_feature_summary), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
host_detection_supported <- function(method) {
  host_detection_tier(method) %in% c("pcr_or_sequencing", "isolation_or_observation")
}

host_detection_tier <- function(method) {
  method <- clean_text(method)
  case_when(
    str_detect(method, regex("PCR|Sequencing", ignore_case = TRUE)) ~ "pcr_or_sequencing",
    str_detect(method, regex("Isolation|Observation", ignore_case = TRUE)) ~ "isolation_or_observation",
    str_detect(method, regex("Antibod", ignore_case = TRUE)) ~ "serology",
    is.na(method) | str_detect(method, regex("not specified", ignore_case = TRUE)) ~ "not_specified",
    TRUE ~ "other_detection"
  )
}

feature_key <- function(disease_name, species_name, tax_id, vector_join_key = NULL) {
  parts <- list(clean_key(disease_name), clean_key(species_name), clean_key(tax_id))
  if (!is.null(vector_join_key)) {
    parts <- c(parts, list(clean_key(vector_join_key)))
  }
  do.call(paste, c(parts, sep = "|"))
}

prepare_host_features <- function(hosts) {
  hosts %>%
    mutate(
      readiness_disease_name = disease_name,
      species_role = "host",
      taxonomy_ok = !missing_as_false(taxonomy_caution),
      host_detection_method = clean_text(host_detection_method),
      host_detection_tier = host_detection_tier(host_detection_method),
      host_direct_detection_supported = host_detection_supported(host_detection_method),
      host_role_assignment = coalesce(clean_text(host_role_assignment), "host_presence_only"),
      host_role_confidence = coalesce(clean_text(host_role_confidence), "low"),
      host_role_assignment_status = coalesce(clean_text(host_role_assignment_status), "candidate_only"),
      host_role_needs_manual_review = missing_as_false(host_role_needs_manual_review),
      host_role_specific = host_role_assignment != "host_presence_only",
      host_role_source_backed = host_role_assignment_status == "draft_source_backed",
      host_role_medium_high = host_role_confidence %in% c("medium", "high")
    ) %>%
    add_host_modelling_proxy() %>%
    mutate(
      profile_broad = taxonomy_ok,
      profile_supported = taxonomy_ok & (
        host_direct_detection_supported |
          (host_role_source_backed & host_role_specific)
      ),
      profile_strong = taxonomy_ok &
        host_role_source_backed &
        host_role_specific &
        host_role_medium_high,
      profile_strict = profile_strong & !host_role_needs_manual_review,
      biological_evidence_tier = max_tier(
        profile_strict,
        profile_strong,
        profile_supported,
        profile_broad
      ),
      tier_rule_id = "host_v0_1",
      host_evidence_missingness_reason = pmap_chr(
        list(
          if_else(taxonomy_ok, NA_character_, "taxonomy_caution"),
          if_else(host_role_specific, NA_character_, "host_role_presence_only"),
          if_else(host_role_needs_manual_review, "host_role_review_needed", NA_character_)
        ),
        combine_reasons
      )
    ) %>%
    select(-host_role_medium_high)
}

prepare_vector_features <- function(vectors) {
  vectors %>%
    mutate(
      species_role = "vector",
      taxonomy_ok = !missing_as_false(taxonomy_caution),
      has_disease_vector_evidence = missing_as_false(has_disease_vector_evidence),
      has_host_vector_evidence = missing_as_false(has_host_vector_evidence),
      has_competence_evidence = missing_as_false(has_competence_evidence),
      best_evidence_level = clean_text(best_evidence_level),
      vector_competence_status = clean_text(vector_competence_status),
      transmission_demonstrated = clean_text(transmission_demonstrated),
      natural_infection_reported = clean_text(natural_infection_reported),
      bites_humans_known = !is.na(clean_text(bites_humans)),
      bites_humans_true = missing_as_false(bites_humans),
      evidence_level_supported = best_evidence_level %in% c("probable", "confirmed"),
      competence_or_transmission_supported =
        vector_competence_status %in% c("competent", "mixed") |
          transmission_demonstrated %in% c("yes", "mixed"),
      profile_broad = taxonomy_ok & has_disease_vector_evidence,
      profile_supported = taxonomy_ok & (
        evidence_level_supported |
          has_competence_evidence |
          has_host_vector_evidence
      ),
      profile_strong = taxonomy_ok & competence_or_transmission_supported,
      profile_strict = profile_strong &
        evidence_level_supported &
        bites_humans_true,
      biological_evidence_tier = max_tier(
        profile_strict,
        profile_strong,
        profile_supported,
        profile_broad
      ),
      tier_rule_id = "vector_v0_1",
      vector_evidence_missingness_reason = pmap_chr(
        list(
          if_else(taxonomy_ok, NA_character_, "taxonomy_caution"),
          if_else(is.na(clean_text(tax_id)), "tax_id_missing", NA_character_),
          if_else(bites_humans_known, NA_character_, "bites_humans_unknown"),
          if_else(!is.na(clean_text(vector_role_hint)), NA_character_, "vector_role_hint_blank"),
          if_else(!is.na(transmission_demonstrated), NA_character_, "transmission_demonstrated_unknown"),
          if_else(!is.na(natural_infection_reported), NA_character_, "natural_infection_unknown")
        ),
        combine_reasons
      )
    )
}

# ------------------------------------------------------------------------------|
#      Inputs ------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
roster <- read_csv_layer(paths$roster, required = TRUE)

host_features <- roster %>%
  filter(species_role == "host") %>%
  prepare_host_features() %>%
  mutate(
    role_modelling_feature_id = feature_key(
      disease_name,
      species_name,
      tax_id
    ),
    feature_rule_version = "host_proxy_rules_v0_1"
  ) %>%
  select(any_of(c(
    "role_modelling_feature_id",
    "disease_name",
    "species_role",
    "species_name",
    "tax_id",
    "host_class",
    "host_order",
    "host_family",
    "host_detection_method",
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
    "profile_group_proxy",
    "feature_rule_version",
    "taxonomy_ok",
    "host_detection_tier",
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "biological_evidence_tier",
    "tier_rule_id",
    "host_evidence_missingness_reason"
  ))) %>%
  arrange(disease_name, species_name)

duplicated_feature_ids <- unique(
  host_features$role_modelling_feature_id[duplicated(host_features$role_modelling_feature_id)]
)
if (length(duplicated_feature_ids) > 0) {
  stop(
    "Role modelling feature IDs are not unique: ",
    paste(head(duplicated_feature_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

feature_summary <- host_features %>%
  count(
    disease_name,
    host_role_bucket,
    host_role_evidence_basis,
    modelling_role_proxy_confidence,
    modelling_role_proxy_needs_review,
    name = "rows"
  ) %>%
  arrange(disease_name, desc(rows), host_role_bucket)

vector_features <- roster %>%
  filter(species_role == "vector") %>%
  prepare_vector_features() %>%
  mutate(
    vector_modelling_feature_id = feature_key(
      disease_name,
      species_name,
      tax_id
    ),
    feature_rule_version = "vector_tier_rules_v0_1"
  ) %>%
  select(any_of(c(
    "vector_modelling_feature_id",
    "disease_name",
    "species_role",
    "species_name",
    "tax_id",
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
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "biological_evidence_tier",
    "tier_rule_id",
    "vector_evidence_missingness_reason",
    "feature_rule_version"
  ))) %>%
  arrange(disease_name, species_name)

duplicated_vector_feature_ids <- unique(
  vector_features$vector_modelling_feature_id[duplicated(vector_features$vector_modelling_feature_id)]
)
if (length(duplicated_vector_feature_ids) > 0) {
  stop(
    "Vector modelling feature IDs are not unique: ",
    paste(head(duplicated_vector_feature_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

vector_feature_summary <- vector_features %>%
  count(
    disease_name,
    biological_evidence_tier,
    best_evidence_level,
    vector_competence_status,
    bites_humans_true,
    name = "rows"
  ) %>%
  arrange(disease_name, biological_evidence_tier, desc(rows))

write_csv(host_features, paths$role_modelling_features, na = "")
write_csv(feature_summary, paths$feature_summary, na = "")
write_csv(vector_features, paths$vector_modelling_features, na = "")
write_csv(vector_feature_summary, paths$vector_feature_summary, na = "")

message("Wrote role modelling features: ", paths$role_modelling_features)
message("Wrote role modelling feature summary: ", paths$feature_summary)
message("Rows in role_modelling_features.csv: ", nrow(host_features))
message("Wrote vector modelling features: ", paths$vector_modelling_features)
message("Wrote vector modelling feature summary: ", paths$vector_feature_summary)
message("Rows in vector_modelling_features.csv: ", nrow(vector_features))
