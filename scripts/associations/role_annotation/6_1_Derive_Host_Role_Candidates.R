#!/usr/bin/env Rscript
################################################################################
# 6_1_Derive_Host_Role_Candidates.R
################################################################################
# Purpose: Derive conservative host-role candidate rows from the canonical WHO
#          disease-pathogen-host backbone for the current role-review scope.
#
# Inputs : master-plus compatibility view for the legacy canonical WHO network
#          WHO diseases helper path for who_pathogens_diseases_zoonotic.csv
#
# Outputs: pathogen_association_data/evidence/role_annotation/
#            host_role_candidates.csv
#          pathogen_association_data/evidence/role_annotation/
#            host_role_candidates_summary.csv
#
# Notes  : This script creates candidate-review rows only. It does not assign
#          final reservoir, amplifier, incidental, or dead-end roles.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr, tidyr)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

safe_lower <- function(x) {
  dplyr::if_else(is.na(x), NA_character_, stringr::str_to_lower(clean_text(x)))
}

is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
}

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
zoonotic_path <- who_pathogens_diseases_zoonotic_path()

output_dir <- role_candidates_dir
output_path <- file.path(output_dir, "host_role_candidates.csv")
summary_path <- file.path(output_dir, "host_role_candidates_summary.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Active role scope -------------------------------------------------------|
# ------------------------------------------------------------------------------|
active_scope <- read_csv(zoonotic_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text)) %>%
  filter(
    is_true(in_gibb_etal) | is_true(in_empres_i),
    Pathogens != "Genus Vesiculovirus"
  ) %>%
  transmute(
    disease_name = Disease_name,
    active_source_pathogen = Pathogens,
    in_gibb_etal = is_true(in_gibb_etal),
    in_empres_i = is_true(in_empres_i),
    priority_prototype_status,
    active_scope_reason = paste0(
      "in_gibb_etal=", in_gibb_etal,
      "; in_empres_i=", in_empres_i,
      "; deferred_broad_vesiculovirus=FALSE"
    )
  )

active_diseases <- active_scope %>%
  group_by(disease_name) %>%
  summarise(
    active_source_pathogens = collapse_unique(active_source_pathogen),
    in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
    in_empres_i = any(in_empres_i, na.rm = TRUE),
    priority_prototype_status = collapse_unique(priority_prototype_status),
    active_scope_reason = collapse_unique(active_scope_reason),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------|
#      Candidate derivation ----------------------------------------------------|
# ------------------------------------------------------------------------------|
livestock_species <- c(
  "Bos taurus",
  "Bubalus bubalis",
  "Camelus bactrianus",
  "Camelus dromedarius",
  "Capra hircus",
  "Equus caballus",
  "Equus ferus",
  "Gallus gallus",
  "Ovis aries",
  "Sus scrofa"
)

network <- read_legacy_compatible_master_plus_network() %>%
  mutate(across(where(is.character), clean_text))

host_role_candidates <- network %>%
  inner_join(active_diseases, by = c("Disease_name" = "disease_name")) %>%
  group_by(
    disease_name = Disease_name,
    active_source_pathogens,
    network_pathogen = Pathogen,
    host = Host,
    host_tax_id = HostTaxID,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    active_scope_reason
  ) %>%
  summarise(
    pathogen_tax_ids = collapse_unique(PathogenTaxID),
    host_class = collapse_unique(HostClass),
    host_order = collapse_unique(HostOrder),
    host_family = collapse_unique(HostFamily),
    network_pathogen_raw_examples = collapse_unique(Pathogen_raw_examples),
    disease_name_raw_examples = collapse_unique(Disease_name_raw_examples),
    detection_method_examples = collapse_unique(DetectionMethod),
    main_source_examples = collapse_unique(MainSource),
    source_row_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    host_name_key = safe_lower(host),
    host_class_key = safe_lower(host_class),
    host_order_key = safe_lower(host_order),
    is_human = host_tax_id == 9606 | host_name_key == "homo sapiens",
    is_livestock_like = host %in% livestock_species,
    is_bird = host_class_key == "aves",
    is_rodent = host_order_key == "rodentia",
    is_bat = host_order_key == "chiroptera",
    is_primate = host_order_key == "primates",
    host_role_candidate = "host_present_in_system",
    candidate_basis = paste(
      "host appears in canonical WHO disease-pathogen-host backbone;",
      "role-specific reservoir/amplifier/incidental/dead-end evidence not assessed"
    ),
    role_confidence = "low",
    needs_manual_review = TRUE,
    review_priority = case_when(
      is_livestock_like ~ "high_livestock_triage",
      is_human ~ "medium_human_triage",
      is_bird | is_rodent | is_bat | is_primate ~ "medium_ecology_triage",
      TRUE ~ "standard_review"
    )
  ) %>%
  select(
    disease_name,
    active_source_pathogens,
    network_pathogen,
    pathogen_tax_ids,
    host,
    host_tax_id,
    host_class,
    host_order,
    host_family,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    active_scope_reason,
    network_pathogen_raw_examples,
    disease_name_raw_examples,
    detection_method_examples,
    main_source_examples,
    source_row_count,
    is_human,
    is_livestock_like,
    is_bird,
    is_rodent,
    is_bat,
    is_primate,
    host_role_candidate,
    candidate_basis,
    role_confidence,
    needs_manual_review,
    review_priority
  ) %>%
  arrange(disease_name, review_priority, host_class, host_order, host)

summary_table <- bind_rows(
  host_role_candidates %>%
    count(disease_name, metric = "candidate_rows", name = "row_count") %>%
    mutate(flag_value = NA_character_),
  host_role_candidates %>%
    count(disease_name, metric = "is_human", flag_value = as.character(is_human), name = "row_count"),
  host_role_candidates %>%
    count(disease_name, metric = "is_livestock_like", flag_value = as.character(is_livestock_like), name = "row_count"),
  host_role_candidates %>%
    count(disease_name, metric = "is_bird", flag_value = as.character(is_bird), name = "row_count"),
  host_role_candidates %>%
    count(disease_name, metric = "is_rodent", flag_value = as.character(is_rodent), name = "row_count"),
  host_role_candidates %>%
    count(disease_name, metric = "is_bat", flag_value = as.character(is_bat), name = "row_count")
) %>%
  select(disease_name, metric, flag_value, row_count) %>%
  arrange(disease_name, metric, desc(flag_value))

write_csv(host_role_candidates, output_path, na = "")
write_csv(summary_table, summary_path, na = "")

message("Wrote host role candidates: ", output_path)
message("Wrote host role summary: ", summary_path)
