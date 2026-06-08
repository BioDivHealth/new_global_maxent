#!/usr/bin/env Rscript
################################################################################
# 6_2_Derive_Species_Host_Vector_Roster.R
################################################################################
# Purpose: Build a disease-species roster for collaborator review, covering both
#          vectored and non-vectored diseases.
#
# Inputs : WHO network helper path for combined_who_network_canonical_zoonotic.csv
#          vector_screening_evidence_path(
#            "disease_vector_links_taxonomy_cleaned_competence_annotated.csv"
#          )
#          pathogen_association_data/evidence/host_vector/
#            vector_host_links_join_ready.csv
#          WHO host-vector helper path for:
#            disease_host_vector_links_expanded_competence_annotated.csv
#          who_zoonotic_pathogens_path()
#
# Outputs: pathogen_association_data/evidence/role_annotation/
#            species_host_vector_roster.csv
#          pathogen_association_data/evidence/role_annotation/
#            species_host_vector_roster_summary.csv
#          pathogen_association_data/evidence/role_annotation/
#            species_host_vector_roster.xlsx
#
# Notes  : This is an evidence roster, not a final biological role assignment.
#          Host rows come from the canonical WHO disease-pathogen-host backbone.
#          Vector rows come from the curated disease-vector evidence table, with
#          host-vector observation flags added from the expanded integration
#          table when available.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr, tidyr, writexl)

source(here::here("scripts", "associations", "working_inputs.R"))

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

collapse_semicolon_values <- function(x) {
  x <- clean_text(x)
  x <- unlist(stringr::str_split(stats::na.omit(x), ";"), use.names = FALSE)
  x <- sort(unique(stats::na.omit(clean_text(x))))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
}

first_non_missing <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

collapse_yes_no <- function(x) {
  values <- unique(stats::na.omit(x))

  if (length(values) == 0) {
    return(FALSE)
  }

  any(values)
}

collapse_yes_no_na <- function(x) {
  values <- unique(stats::na.omit(x))

  if (length(values) == 0) {
    return(NA)
  }

  any(values)
}

summarise_bites_humans_basis <- function(x) {
  values <- clean_text(x)
  values <- unlist(stringr::str_split(stats::na.omit(values), ";"), use.names = FALSE)
  values <- clean_text(values)
  values <- sort(unique(stats::na.omit(values)))

  if (length(values) == 0) {
    return(NA_character_)
  }

  if (all(c("blood_meal", "on_host_occurrence") %in% values)) {
    return("both")
  }

  paste(values, collapse = "; ")
}

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
host_vector_dir <- vector_host_outputs_dir
role_dir <- role_roster_dir

network_path <- who_canonical_zoonotic_network_path()
vector_path <- vector_screening_evidence_path(
  "disease_vector_links_taxonomy_cleaned_competence_annotated.csv"
)
host_vector_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
expanded_path <- who_network_host_vector_path(
  "disease_host_vector_links_expanded_competence_annotated.csv"
)
zoonotic_path <- who_pathogens_diseases_zoonotic_path()

output_path <- file.path(role_dir, "species_host_vector_roster.csv")
summary_path <- file.path(role_dir, "species_host_vector_roster_summary.csv")
xlsx_path <- file.path(role_dir, "species_host_vector_roster.xlsx")
host_role_candidates_path <- file.path(role_dir, "host_role_candidates.csv")
host_role_assignments_path <- file.path(role_dir, "host_role_assignments.csv")

dir.create(role_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Disease Scope Metadata --------------------------------------------------|
# ------------------------------------------------------------------------------|
broad_taxon_exclusions <- c(
  "Genus Vesiculovirus",
  "Subgenus Merbecovirus",
  "Subgenus Sarbecovirus"
)

disease_scope_raw <- read_csv(
  zoonotic_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

broad_disease_exclusions <- disease_scope_raw %>%
  filter(Pathogens %in% broad_taxon_exclusions) %>%
  distinct(Disease_name) %>%
  pull(Disease_name)

disease_scope <- disease_scope_raw %>%
  filter(!Pathogens %in% broad_taxon_exclusions) %>%
  group_by(disease_name = Disease_name) %>%
  summarise(
    source_pathogens = collapse_unique(Pathogens),
    in_gibb_etal = any(is_true(in_gibb_etal), na.rm = TRUE),
    in_empres_i = any(is_true(in_empres_i), na.rm = TRUE),
    priority_prototype_status = collapse_unique(priority_prototype_status),
    in_current_role_review_scope = any(
      is_true(in_gibb_etal) | is_true(in_empres_i),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------|
#      Host Role Summaries -----------------------------------------------------|
# ------------------------------------------------------------------------------|
role_join_key <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_to_lower(x)
  stringr::str_replace_all(x, "\\s+", " ")
}

tax_join_key <- function(x) {
  x <- clean_text(x)
  x[is.na(x)] <- ""
  x
}

host_role_assignments <- read_csv(
  host_role_assignments_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    role_disease_key = role_join_key(disease_name),
    role_species_key = role_join_key(host),
    role_tax_id_key = tax_join_key(host_tax_id),
    assignment_needs_manual_review = is_true(needs_manual_review)
  ) %>%
  filter(
    !is.na(role_disease_key),
    !is.na(role_species_key),
    !is.na(host_role_assignment)
  ) %>%
  group_by(role_disease_key, role_species_key, role_tax_id_key) %>%
  summarise(
    host_role_assignment = collapse_unique(host_role_assignment),
    host_role_confidence_from_assignment = collapse_unique(assignment_confidence),
    host_role_assignment_status = collapse_unique(assignment_status),
    host_role_needs_manual_review_from_assignment = any(assignment_needs_manual_review, na.rm = TRUE),
    .groups = "drop"
  )

host_role_candidates <- read_csv(
  host_role_candidates_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    role_disease_key = role_join_key(disease_name),
    role_species_key = role_join_key(host),
    role_tax_id_key = tax_join_key(host_tax_id)
  ) %>%
  filter(
    !is.na(role_disease_key),
    !is.na(role_species_key)
  ) %>%
  group_by(role_disease_key, role_species_key, role_tax_id_key) %>%
  summarise(
    host_role_candidate_confidence = collapse_unique(role_confidence),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------|
#      Host Rows ---------------------------------------------------------------|
# ------------------------------------------------------------------------------|
host_rows <- read_csv(
  network_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(HostTaxID = clean_text(HostTaxID)) %>%
  filter(!Disease_name %in% broad_disease_exclusions) %>%
  group_by(
    disease_name = Disease_name,
    species_name = Host,
    tax_id = HostTaxID
  ) %>%
  summarise(
    source_pathogens_from_network = collapse_unique(Pathogen),
    pathogen_tax_ids = collapse_unique(PathogenTaxID),
    host_class = collapse_unique(HostClass),
    host_order = collapse_unique(HostOrder),
    host_family = collapse_unique(HostFamily),
    host_detection_method = collapse_semicolon_values(DetectionMethod),
    main_source_examples = collapse_unique(MainSource),
    host_source_row_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    appears_as_host = TRUE,
    appears_as_vector = FALSE,
    evidence_source_type = "who_host_backbone",
    vector_group = NA_character_,
    vector_taxon_rank = NA_character_,
    vector_join_key = NA_character_,
    has_disease_vector_evidence = FALSE,
    has_host_vector_evidence = FALSE,
    has_competence_evidence = FALSE,
    disease_vector_evidence_status = NA_character_,
    best_evidence_level = NA_character_,
    best_evidence_basis = NA_character_,
    vector_record_sources = NA_character_,
    vector_competence_status = NA_character_,
    transmission_demonstrated = NA_character_,
    natural_infection_reported = NA_character_,
    vector_role_hint = NA_character_,
    uncertainty_reason = NA_character_,
    vector_host_record_count = NA_real_,
    vector_source_platform_examples = NA_character_,
    vector_country_examples = NA_character_,
    bites_humans = NA,
    bites_humans_basis = NA_character_,
    taxonomy_caution = FALSE,
    review_note = NA_character_
  )

# ------------------------------------------------------------------------------|
#      Vector Rows -------------------------------------------------------------|
# ------------------------------------------------------------------------------|
expanded_vector_support <- read_csv(
  expanded_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    vector_join_key = clean_text(vector_join_key),
    has_disease_vector_evidence = disease_vector_evidence %in% c(TRUE, "TRUE", "true", "True", 1, "1"),
    has_host_vector_evidence = host_vector_evidence %in% c(TRUE, "TRUE", "true", "True", 1, "1"),
    taxonomy_caution = taxonomy_caution %in% c(TRUE, "TRUE", "true", "True", 1, "1")
  ) %>%
  filter(
    !disease_name %in% broad_disease_exclusions,
    !is.na(disease_name),
    !is.na(vector_join_key)
  ) %>%
  group_by(
    disease_name,
    vector_join_key
  ) %>%
  summarise(
    has_host_vector_evidence = collapse_yes_no(has_host_vector_evidence),
    vector_host_record_count = sum(
      suppressWarnings(as.numeric(vector_host_record_count)),
      na.rm = TRUE
    ),
    vector_source_platform_examples = collapse_unique(source_platform_examples),
    vector_country_examples = collapse_unique(country_examples),
    taxonomy_caution = collapse_yes_no(taxonomy_caution),
    review_note = collapse_unique(review_reason_examples),
    .groups = "drop"
  )

human_vector_support <- read_csv(
  host_vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    is_human_host = host_tax_id == "9606" | stringr::str_to_lower(host) == "homo sapiens"
  ) %>%
  filter(!is.na(vector_join_key)) %>%
  group_by(vector_join_key) %>%
  summarise(
    has_any_host_vector_evidence = TRUE,
    bites_humans = any(is_human_host, na.rm = TRUE),
    bites_humans_basis = summarise_bites_humans_basis(
      interaction_type_examples[is_human_host]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    bites_humans = dplyr::if_else(
      has_any_host_vector_evidence,
      bites_humans,
      NA
    )
  )

vector_rows <- read_csv(
  vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    vector_display_name = dplyr::coalesce(vector_species_taxonomy_cleaned, vector_species),
    vector_join_key = clean_text(vector_join_key),
    has_competence_evidence = !is.na(vector_competence_status)
  ) %>%
  filter(
    !disease_name %in% broad_disease_exclusions,
    !is.na(disease_name),
    !is.na(vector_display_name)
  ) %>%
  left_join(
    expanded_vector_support,
    by = c("disease_name", "vector_join_key")
  ) %>%
  left_join(
    human_vector_support,
    by = "vector_join_key"
  ) %>%
  group_by(
    disease_name,
    species_name = vector_display_name,
    tax_id = NA_character_
  ) %>%
  summarise(
    vector_group = collapse_unique(vector_group),
    vector_taxon_rank = collapse_unique(vector_taxon_rank),
    vector_join_key = collapse_unique(vector_join_key),
    has_disease_vector_evidence = TRUE,
    has_host_vector_evidence = any(has_host_vector_evidence, na.rm = TRUE),
    has_competence_evidence = any(has_competence_evidence, na.rm = TRUE),
    disease_vector_evidence_status = "supported_in_disease_vector_table",
    best_evidence_level = collapse_unique(best_evidence_level),
    best_evidence_basis = collapse_unique(best_evidence_basis),
    vector_record_sources = collapse_unique(record_sources),
    vector_competence_status = collapse_unique(vector_competence_status),
    transmission_demonstrated = collapse_unique(transmission_demonstrated),
    natural_infection_reported = collapse_unique(natural_infection_reported),
    vector_role_hint = collapse_unique(vector_role_hint),
    uncertainty_reason = collapse_unique(uncertainty_reason),
    bites_humans = collapse_yes_no_na(bites_humans),
    bites_humans_basis = collapse_unique(bites_humans_basis),
    vector_host_record_count = sum(
      suppressWarnings(as.numeric(vector_host_record_count)),
      na.rm = TRUE
    ),
    vector_source_platform_examples = collapse_unique(vector_source_platform_examples),
    vector_country_examples = collapse_unique(vector_country_examples),
    taxonomy_caution = any(
      review_needed %in% c(TRUE, "TRUE", "true", "True", 1, "1") |
        taxonomy_caution %in% c(TRUE, "TRUE", "true", "True", 1, "1"),
      na.rm = TRUE
    ),
    review_note = collapse_unique(dplyr::coalesce(review_note.x, review_note.y)),
    .groups = "drop"
  ) %>%
  mutate(
    appears_as_host = FALSE,
    appears_as_vector = TRUE,
    evidence_source_type = dplyr::case_when(
      has_disease_vector_evidence & has_host_vector_evidence ~
        "disease_vector_and_host_vector_evidence",
      has_disease_vector_evidence ~ "disease_vector_evidence_only",
      TRUE ~ "vector_evidence_unspecified"
    ),
    source_pathogens_from_network = NA_character_,
    pathogen_tax_ids = NA_character_,
    host_class = NA_character_,
    host_order = NA_character_,
    host_family = NA_character_,
    host_detection_method = NA_character_,
    main_source_examples = NA_character_,
    host_source_row_count = NA_integer_
  )

# ------------------------------------------------------------------------------|
#      Combined Roster ---------------------------------------------------------|
# ------------------------------------------------------------------------------|
all_columns <- union(names(host_rows), names(vector_rows))

add_missing_columns <- function(df, columns) {
  missing_cols <- setdiff(columns, names(df))

  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      df[[col]] <- NA
    }
  }

  df[, columns]
}

roster <- bind_rows(
  add_missing_columns(host_rows, all_columns),
  add_missing_columns(vector_rows, all_columns)
) %>%
  group_by(disease_name, species_name, tax_id) %>%
  summarise(
    appears_as_host = any(appears_as_host, na.rm = TRUE),
    appears_as_vector = any(appears_as_vector, na.rm = TRUE),
    species_role = dplyr::case_when(
      appears_as_vector ~ "vector",
      appears_as_host ~ "host",
      TRUE ~ "unknown"
    ),
    evidence_source_type = collapse_unique(evidence_source_type),
    source_pathogens_from_network = collapse_unique(source_pathogens_from_network),
    pathogen_tax_ids = collapse_unique(pathogen_tax_ids),
    host_class = collapse_unique(host_class),
    host_order = collapse_unique(host_order),
    host_family = collapse_unique(host_family),
    host_detection_method = collapse_semicolon_values(host_detection_method),
    main_source_examples = collapse_unique(main_source_examples),
    host_source_row_count = sum(suppressWarnings(as.integer(host_source_row_count)), na.rm = TRUE),
    vector_group = collapse_unique(vector_group),
    vector_taxon_rank = collapse_unique(vector_taxon_rank),
    vector_join_key = collapse_unique(vector_join_key),
    has_disease_vector_evidence = any(has_disease_vector_evidence, na.rm = TRUE),
    has_host_vector_evidence = any(has_host_vector_evidence, na.rm = TRUE),
    has_competence_evidence = any(has_competence_evidence, na.rm = TRUE),
    best_evidence_level = collapse_unique(best_evidence_level),
    best_evidence_basis = collapse_unique(best_evidence_basis),
    vector_record_sources = collapse_unique(vector_record_sources),
    vector_competence_status = collapse_unique(vector_competence_status),
    transmission_demonstrated = collapse_unique(transmission_demonstrated),
    natural_infection_reported = collapse_unique(natural_infection_reported),
    vector_role_hint = collapse_unique(vector_role_hint),
    uncertainty_reason = collapse_unique(uncertainty_reason),
    bites_humans = collapse_yes_no_na(bites_humans),
    bites_humans_basis = collapse_unique(bites_humans_basis),
    taxonomy_caution = any(taxonomy_caution, na.rm = TRUE),
    .groups = "drop"
) %>%
  left_join(disease_scope, by = "disease_name") %>%
  mutate(
    role_disease_key = role_join_key(disease_name),
    role_species_key = role_join_key(species_name),
    role_tax_id_key = tax_join_key(tax_id)
  ) %>%
  left_join(
    host_role_assignments,
    by = c("role_disease_key", "role_species_key", "role_tax_id_key")
  ) %>%
  left_join(
    host_role_candidates,
    by = c("role_disease_key", "role_species_key", "role_tax_id_key")
  ) %>%
  mutate(
    disease_has_vector_rows = disease_name %in% unique(vector_rows$disease_name),
    has_host_role_assignment = species_role == "host" & !is.na(host_role_assignment),
    host_role_assignment = dplyr::case_when(
      species_role != "host" ~ NA_character_,
      has_host_role_assignment ~ host_role_assignment,
      TRUE ~ "host_presence_only"
    ),
    host_role_confidence = dplyr::case_when(
      species_role != "host" ~ NA_character_,
      has_host_role_assignment ~ host_role_confidence_from_assignment,
      !is.na(host_role_candidate_confidence) ~ host_role_candidate_confidence,
      TRUE ~ "low"
    ),
    host_role_needs_manual_review = dplyr::case_when(
      species_role != "host" ~ NA,
      has_host_role_assignment ~ host_role_needs_manual_review_from_assignment,
      TRUE ~ TRUE
    ),
    host_role_assignment_status = dplyr::case_when(
      species_role != "host" ~ NA_character_,
      has_host_role_assignment ~ host_role_assignment_status,
      TRUE ~ "candidate_only"
    ),
    review_boundary = dplyr::case_when(
      appears_as_host & !appears_as_vector ~
        "host presence only; final host role not assigned",
      appears_as_vector & has_disease_vector_evidence & has_competence_evidence ~
        "vector evidence plus competence annotation; final vector role not assigned",
      appears_as_vector & has_disease_vector_evidence ~
        "curated disease-vector evidence; competence annotation absent",
      appears_as_vector & has_host_vector_evidence ~
        "host-vector observation only for this disease context; not curated disease-vector support",
      TRUE ~ "review needed"
    )
  ) %>%
  select(
    disease_name,
    species_name,
    tax_id,
    species_role,
    disease_has_vector_rows,
    in_current_role_review_scope,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    source_pathogens,
    host_class,
    host_order,
    host_family,
    host_detection_method,
    host_role_assignment,
    host_role_confidence,
    host_role_needs_manual_review,
    host_role_assignment_status,
    vector_group,
    vector_taxon_rank,
    vector_join_key,
    has_disease_vector_evidence,
    has_host_vector_evidence,
    has_competence_evidence,
    best_evidence_level,
    best_evidence_basis,
    vector_record_sources,
    bites_humans,
    bites_humans_basis,
    vector_competence_status,
    transmission_demonstrated,
    natural_infection_reported,
    vector_role_hint,
    uncertainty_reason,
    taxonomy_caution
  ) %>%
  arrange(disease_name, desc(species_role == "vector"), species_role, species_name)

summary_table <- bind_rows(
  roster %>%
    count(disease_name, metric = "roster_rows", name = "row_count") %>%
    mutate(flag_value = NA_character_),
  roster %>%
    count(disease_name, metric = "species_role", flag_value = species_role, name = "row_count"),
  roster %>%
    count(disease_name, metric = "has_disease_vector_evidence", flag_value = as.character(has_disease_vector_evidence), name = "row_count"),
  roster %>%
    count(disease_name, metric = "has_competence_evidence", flag_value = as.character(has_competence_evidence), name = "row_count")
) %>%
  select(disease_name, metric, flag_value, row_count) %>%
  arrange(disease_name, metric, desc(flag_value))

column_dictionary <- tibble::tribble(
  ~column, ~description,
  "disease_name", "WHO disease label used in the current canonical disease-pathogen-host workflow.",
  "species_name", "Host species name or vector taxon name represented by this row.",
  "tax_id", "NCBI TaxID where available. Host rows usually have TaxIDs; vector rows are currently name-keyed and usually blank.",
  "species_role", "Whether this row represents a host or vector in the current evidence roster.",
  "disease_has_vector_rows", "TRUE if the disease has one or more curated vector rows in this roster; FALSE for diseases represented only by host rows.",
  "in_current_role_review_scope", "TRUE if the disease/pathogen source row is currently in the role-review scope after excluding broad genus/subgenus placeholders.",
  "in_gibb_etal", "TRUE if the source pathogen/disease row is marked as included in the Gibb et al. priority set.",
  "in_empres_i", "TRUE if the source pathogen/disease row is marked as included in the EMPRES-i set.",
  "priority_prototype_status", "WHO source status carried from the upstream disease table: priority, prototype, both, or none.",
  "source_pathogens", "Upstream pathogen or analysis-unit names contributing to this disease after broad genus/subgenus exclusions.",
  "host_class", "Host taxonomic class for host rows, when available.",
  "host_order", "Host taxonomic order for host rows, when available.",
  "host_family", "Host taxonomic family for host rows, when available.",
  "host_detection_method", "For host rows, harmonised host-pathogen detection method from the canonical WHO backbone, including CLOVER and VIRION evidence such as PCR/Sequencing and Isolation/Observation.",
  "host_role_assignment", "Compact reviewed or fallback host role used for modelling handoff. Source-backed rows come from host_role_assignments.csv; unassigned host rows are marked host_presence_only.",
  "host_role_confidence", "Confidence for host_role_assignment. Source-backed rows use assignment_confidence; fallback host_presence_only rows use host_role_candidates.csv confidence where available.",
  "host_role_needs_manual_review", "TRUE when the host role assignment remains review-flagged or is only a candidate/presence fallback.",
  "host_role_assignment_status", "Status for host_role_assignment, such as draft_source_backed or candidate_only.",
  "vector_group", "Broad vector group for vector rows, such as mosquito, tick, flea, or midge.",
  "vector_taxon_rank", "Taxonomic grain of the vector name after cleanup, such as species, genus, or infraspecific.",
  "vector_join_key", "Normalized vector name key used to join disease-vector, host-vector, and competence evidence.",
  "has_disease_vector_evidence", "TRUE if the vector appears in the curated disease-vector table for this disease.",
  "has_host_vector_evidence", "TRUE if this vector also has VectorMap/MapVEu host-vector evidence in the integrated host-vector layer.",
  "has_competence_evidence", "TRUE if this disease-vector pair has a joined vector competence annotation.",
  "best_evidence_level", "Best curated disease-vector evidence level, such as confirmed, probable, candidate, or poor/unsupported.",
  "best_evidence_basis", "Short basis for the best disease-vector evidence level, such as field, lab, field+lab, review, or EFSA-derived evidence.",
  "vector_record_sources", "Whether vector evidence came from the literature review, EFSA, or both.",
  "bites_humans", "For vector rows, TRUE if VectorMap/MapVEu host-vector evidence links this vector to Homo sapiens; FALSE if host-vector evidence exists but not to humans; blank if no host-vector evidence is available for that vector.",
  "bites_humans_basis", "Evidence basis for bites_humans when TRUE: blood_meal, on_host_occurrence, or both.",
  "vector_competence_status", "Joined competence status for the disease-vector pair, such as competent, mixed, not_competent, or unclear.",
  "transmission_demonstrated", "Whether transmission was demonstrated in the competence evidence where extractable.",
  "natural_infection_reported", "Whether natural infection was reported in the competence evidence where extractable.",
  "vector_role_hint", "Source-language hint for vector role, such as primary_vector, bridge_vector, or sylvatic_vector, where captured.",
  "uncertainty_reason", "Compact caveat flags from competence extraction, such as field_detection_only, temperature_dependent, or no_transmission_demonstrated.",
  "taxonomy_caution", "TRUE when vector taxonomy cleanup or matching raised a review caution."
)

write_csv(roster, output_path, na = "")
write_csv(summary_table, summary_path, na = "")
write_xlsx(
  list(
    roster = roster,
    column_descriptions = column_dictionary
  ),
  xlsx_path
)

message("Wrote species host/vector roster: ", output_path)
message("Wrote species host/vector roster summary: ", summary_path)
message("Wrote species host/vector roster workbook: ", xlsx_path)
