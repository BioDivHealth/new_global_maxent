#!/usr/bin/env Rscript
################################################################################
# 04_build_role_gap_source_check_candidates.R
################################################################################
# Purpose: Generate a small source-check-shaped queue for high-priority host
#          rows where current role information remains weak or missing.
#
# Output : pathogen_association_data/staged/role_annotation/source_check_candidates/
#            role_gap_source_check_candidates.csv
#
# Notes  : This is a curation input. It does not modify accepted evidence,
#          assignment tables, or the source-check decision ledger.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr, tibble)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "association_data_helpers.R"
))

candidate_queue_columns <- c(
  "candidate_row_id",
  "batch_order",
  "phase",
  "batch_id",
  "disease_name",
  "entity_type",
  "entity_name",
  "role_assignment",
  "assignment_confidence",
  "review_priority",
  "evidence_source_ids",
  "evidence_basis",
  "review_reason",
  "join_note",
  "source_check_note"
)

priority_diseases <- c(
  "Chikungunya fever",
  "Dengue",
  "West Nile fever",
  "Rift Valley fever",
  "Influenza (H5N1 avian influenza)"
)

make_candidate_id <- function(disease_name, entity_name) {
  stable_candidate_id(
    "role_gap",
    "host_roles_v0_1",
    disease_name,
    entity_name,
    "source_check"
  )
}

features <- read_csv_layer(role_modelling_features_path(), required = TRUE)

missing_cols <- setdiff(
  c(
    "disease_name",
    "species_name",
    "host_class",
    "host_order",
    "host_family",
    "host_detection_method",
    "host_direct_detection_supported",
    "host_role_bucket",
    "host_role_evidence_basis",
    "modelling_role_proxy_needs_review"
  ),
  names(features)
)
if (length(missing_cols) > 0) {
  stop(
    "Role modelling features are missing required columns for gap candidates: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

role_gap_candidates <- features %>%
  filter(
    species_role == "host",
    host_role_bucket %in% c("host_presence_only", "unknown_or_unreviewed") |
      modelling_role_proxy_needs_review
  ) %>%
  mutate(
    review_priority_score =
      4 * as.integer(host_direct_detection_supported) +
      3 * as.integer(disease_name %in% priority_diseases) +
      2 * as.integer(host_class %in% c("mammalia", "aves")) +
      2 * as.integer(host_order %in% c("primates", "rodentia")),
    review_priority = case_when(
      review_priority_score >= 11 ~ "P1_role_gap_source_check",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(review_priority)) %>%
  arrange(desc(review_priority_score), disease_name, species_name) %>%
  mutate(
    candidate_row_id = make_candidate_id(disease_name, species_name),
    batch_order = "role_gap",
    phase = "role_gap_source_check",
    batch_id = "role_gap_host_roles_v0_1",
    entity_type = "host",
    entity_name = species_name,
    role_assignment = "host_role_gap_source_check",
    assignment_confidence = "unknown",
    evidence_source_ids = "",
    evidence_basis = paste(
      "current_bucket=", host_role_bucket,
      "; evidence_basis=", host_role_evidence_basis,
      sep = ""
    ),
    review_reason = paste(
      "role_gap_priority_score=", review_priority_score,
      "; host_class=", coalesce(host_class, ""),
      "; host_order=", coalesce(host_order, ""),
      "; detection=", coalesce(host_detection_method, ""),
      sep = ""
    ),
    join_note = "",
    source_check_note = paste(
      "Source-check whether this host has a source-backed biological role",
      "beyond host presence before adding official role evidence or assignment rows."
    )
  ) %>%
  select(all_of(candidate_queue_columns))

duplicated_candidate_ids <- unique(
  role_gap_candidates$candidate_row_id[duplicated(role_gap_candidates$candidate_row_id)]
)
if (length(duplicated_candidate_ids) > 0) {
  stop(
    "Role-gap candidate IDs are not unique: ",
    paste(head(duplicated_candidate_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

dir.create(role_source_check_candidates_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(role_gap_candidates, role_gap_source_check_candidates_path(), na = "")

message("Wrote role-gap source-check candidates: ", role_gap_source_check_candidates_path())
message("Rows in role_gap_source_check_candidates.csv: ", nrow(role_gap_candidates))
