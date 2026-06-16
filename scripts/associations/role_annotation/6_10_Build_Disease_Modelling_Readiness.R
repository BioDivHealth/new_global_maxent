#!/usr/bin/env Rscript
################################################################################
# 6_10_Build_Disease_Modelling_Readiness.R
################################################################################
# Purpose: Join disease master-list modelling rules to live evidence summaries.
#
# Output : pathogen_association_data/readiness/
#            disease_modelling_readiness.csv
#            disease_modelling_pilot.csv
#            disease_modelling_pilot_package/
#            disease_modelling_pilot_package.rds
#            disease_modelling_pilot_package.xlsx, when writexl is available
#
# Notes  : This is a planning/readiness surface. It does not assign biological
#          roles, infer vector evidence, or expand any upstream evidence layer.
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

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
source(here::here("scripts", "associations", "role_annotation", "disease_modelling_readiness_helpers.R"))

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
role_dir <- role_annotation_dir
qa_dir <- role_qa_dir
genbank_legacy_dir <- genbank_simple_legacy_dir
sdm_dir <- here::here("sdms")

dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(readiness_dir, recursive = TRUE, showWarnings = FALSE)

genbank_readiness_summary_path <- file.path(
  genbank_simple_evidence_dir,
  "genbank_readiness_disease_country_summary_standardized.csv"
)
genbank_readiness_summary_path <- prefer_existing_path(
  genbank_readiness_summary_path,
  file.path(
    genbank_legacy_dir,
    "genbank_readiness_disease_country_summary_standardized.csv"
  )
)
genbank_standard_summary_path <- file.path(
  genbank_simple_evidence_dir,
  "genbank_disease_country_summary_standardized.csv"
)
genbank_standard_summary_path <- prefer_existing_path(
  genbank_standard_summary_path,
  file.path(genbank_legacy_dir, "genbank_disease_country_summary_standardized.csv")
)
genbank_summary_source <- if (file.exists(genbank_readiness_summary_path)) {
  "readiness_combined"
} else {
  "standard"
}

paths <- list(
  master = who_master_plus_analysis_units_path(),
  rules = who_diseases_transmission_rules_path(
    "master_plus_who_transmission_rules_manual_reviewed_v2.csv"
  ),
  master_disease_units = who_master_disease_analysis_units_path(),
  disease_evidence_readiness = file.path(qa_dir, "disease_evidence_readiness.csv"),
  vector_evidence_readiness = file.path(qa_dir, "vector_evidence_readiness_by_disease.csv"),
  species_host_vector_roster = file.path(role_dir, "species_host_vector_roster.csv"),
  who_don_modelling_ready = prefer_existing_path(
    file.path(who_don_v2_final_dir, "who_don_modelling_ready.csv"),
    file.path(who_don_v2_legacy_dir, "final", "who_don_modelling_ready.csv")
  ),
  genbank_disease_country_summary = if (genbank_summary_source == "readiness_combined") {
    genbank_readiness_summary_path
  } else {
    genbank_standard_summary_path
  },
  accessible_sdm_species = file.path(sdm_dir, "outputs", "catalog", "accessible_sdm_species.csv"),
  sdm_projections = file.path(sdm_dir, "outputs", "projections", "projection_manifest.csv"),
  sdm_comparisons = file.path(sdm_dir, "outputs", "comparisons", "comparison_manifest.csv"),
  virion_taxid = file.path(who_virion_dir, "who_pathogens_virion_taxid.csv"),
  clover_taxid = file.path(who_clover_dir, "who_bacteria_clover_taxid.csv")
)

# ------------------------------------------------------------------------------|
#      Inputs ------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
master <- read_csv_layer(paths$master, required = TRUE)
rules <- read_csv_layer(paths$rules, required = TRUE)
master_disease_units <- read_csv_layer(paths$master_disease_units)
disease_qa <- read_csv_layer(paths$disease_evidence_readiness)
vector_qa <- read_csv_layer(paths$vector_evidence_readiness)
roster <- read_csv_layer(paths$species_host_vector_roster)
who_don <- read_csv_layer(paths$who_don_modelling_ready)
genbank <- read_csv_layer(paths$genbank_disease_country_summary)
accessible_sdm_species <- read_csv_layer(paths$accessible_sdm_species)
sdm_projections <- read_csv_layer(paths$sdm_projections)
sdm_comparisons <- read_csv_layer(paths$sdm_comparisons)
virion_taxid <- read_csv_layer(paths$virion_taxid)
clover_taxid <- read_csv_layer(paths$clover_taxid)

if (!"analysis_unit_id" %in% names(master)) {
  stop("`analysis_unit_id` is missing from the master analysis-unit table.", call. = FALSE)
}

pathogen_taxid_lookup <- build_combined_pathogen_taxid_lookup(virion_taxid, clover_taxid)

# ------------------------------------------------------------------------------|
#      Master Readiness Rows ---------------------------------------------------|
# ------------------------------------------------------------------------------|
rule_cols <- c(
  "vectored_status",
  "generalist_status",
  "transmission_complexity",
  "guild",
  "host_sdm_needed",
  "vector_sdm_needed",
  "host_range_rule",
  "vector_range_rule",
  "range_limiting_layer",
  "transmission_rule_notes",
  "transmission_rule_review_status",
  "modelling_scope_status",
  "modelling_scope_reason"
)

rules_for_join <- rules %>%
  select(any_of(c("analysis_unit_id", rule_cols))) %>%
  distinct(analysis_unit_id, .keep_all = TRUE) %>%
  rename_with(~ paste0(.x, "_rules"), all_of(intersect(rule_cols, names(.))))

base <- master %>%
  mutate(.row_order = row_number()) %>%
  left_join(rules_for_join, by = "analysis_unit_id")

for (col in rule_cols) {
  base[[col]] <- coalesce_column(base, paste0(col, "_rules"), col)
}

base <- base %>%
  select(-any_of(paste0(rule_cols, "_rules"))) %>%
  mutate(
    readiness_disease_name = first_non_empty(
      source_disease_name,
      disease_master_name,
      analysis_unit_label,
      analysis_unit
    )
  )

if (!is.null(master_disease_units) && "master_row" %in% names(master_disease_units)) {
  master_flags <- master_disease_units %>%
    select(any_of(c("master_row", "in_master_who"))) %>%
    mutate(master_row = as.character(master_row)) %>%
    rename(in_master_who_from_master_units = in_master_who) %>%
    distinct(master_row, .keep_all = TRUE)

  base <- base %>%
    mutate(master_row = as.character(master_row)) %>%
    left_join(master_flags, by = "master_row") %>%
    mutate(
      in_master_who = case_when(
        is_true(in_master_who_from_master_units) ~ TRUE,
        !is.na(in_master_who_from_master_units) ~ FALSE,
        include_as_analysis_unit == "yes_existing_who" ~ TRUE,
        row_type == "source_row" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    select(-any_of("in_master_who_from_master_units"))
} else {
  base <- base %>%
    mutate(
      in_master_who = include_as_analysis_unit == "yes_existing_who" | row_type == "source_row"
    )
}

candidate_fields <- c("source_disease_name", "disease_master_name", "analysis_unit_label", "analysis_unit")

# ------------------------------------------------------------------------------|
#      Evidence Summaries ------------------------------------------------------|
# ------------------------------------------------------------------------------|
disease_qa_summary <- if (is.null(disease_qa)) {
  empty_disease_summary()
} else {
  disease_qa %>%
    select(any_of(c(
      "disease_name",
      "host_candidate_rows",
      "roster_host_rows",
      "roster_vector_rows",
      "roster_vectors_with_host_vector_evidence",
      "roster_vectors_with_human_biting",
      "roster_taxonomy_caution_rows",
      "host_role_evidence_rows",
      "host_role_evidence_substantive_rows",
      "vector_role_evidence_rows",
      "vector_role_evidence_substantive_rows",
      "host_role_assignment_rows",
      "host_role_substantive_assignment_rows",
      "vector_role_assignment_rows",
      "vector_role_substantive_assignment_rows",
      "has_substantive_host_role_evidence",
      "has_substantive_vector_role_evidence",
      "role_assignment_status"
    )))
}

vector_summary <- if (is.null(vector_qa)) {
  empty_disease_summary()
} else {
  vector_qa %>%
    transmute(
      disease_name,
      direct_disease_vector_rows = vector_rows,
      direct_vector_species_or_taxa = vector_species_or_taxa,
      direct_vector_confirmed_rows = confirmed_rows,
      direct_vector_probable_rows = probable_rows,
      direct_vector_candidate_rows = candidate_rows,
      direct_vector_competence_joined_rows = competence_joined_rows,
      direct_vector_competent_rows = competent_rows,
      direct_vector_mixed_competence_rows = mixed_competence_rows,
      direct_vector_not_competent_rows = not_competent_rows,
      direct_vector_unclear_competence_rows = unclear_competence_rows,
      direct_vector_transmission_yes_rows = transmission_yes_rows,
      direct_vector_natural_infection_yes_rows = natural_infection_yes_rows,
      direct_vector_taxonomy_review_needed_rows = taxonomy_review_needed_rows
    )
}

genbank_summary <- if (is.null(genbank)) {
  empty_disease_summary()
} else {
  genbank %>%
    mutate(disease_name = Disease_name) %>%
    tidyr::separate_rows(disease_name, sep = ";\\s*") %>%
    mutate(disease_name = clean_text(disease_name)) %>%
    filter(!is.na(disease_name)) %>%
    group_by(disease_name) %>%
    summarise(
      genbank_country_rows = dplyr::n(),
      genbank_distinct_countries_or_territories = dplyr::n_distinct(country_standardized),
      genbank_records_with_country = sum(records_with_country, na.rm = TRUE),
      .groups = "drop"
    )
}

who_don_summary <- if (is.null(who_don)) {
  empty_disease_summary()
} else {
  who_don %>%
    group_by(disease_name = disease_label_standard) %>%
    summarise(
      who_don_focal_country_rows = dplyr::n(),
      who_don_distinct_records = dplyr::n_distinct(record_key),
      who_don_distinct_countries = dplyr::n_distinct(country_standard),
      who_don_high_confidence_rows = sum(scope_confidence == "high", na.rm = TRUE),
      who_don_manual_review_rows = sum(is_true(needs_review), na.rm = TRUE),
      .groups = "drop"
    )
}

sdm_summary <- if (is.null(roster)) {
  empty_disease_summary()
} else {
  accessible_species <- if (is.null(accessible_sdm_species) || !"species" %in% names(accessible_sdm_species)) {
    character()
  } else {
    unique(clean_key(accessible_sdm_species$species))
  }
  projection_species <- available_species(sdm_projections)
  comparison_species <- available_species(sdm_comparisons)

  roster %>%
    mutate(
      species_key = clean_key(species_name),
      species_role = clean_key(species_role)
    ) %>%
    filter(!is.na(disease_name), !is.na(species_key)) %>%
    group_by(disease_name) %>%
    summarise(
      host_sdm_species_available = dplyr::n_distinct(
        species_name[species_role == "host" & species_key %in% accessible_species],
        na.rm = TRUE
      ),
      vector_sdm_species_available = dplyr::n_distinct(
        species_name[species_role == "vector" & species_key %in% accessible_species],
        na.rm = TRUE
      ),
      host_sdm_projection_species_available = dplyr::n_distinct(
        species_name[species_role == "host" & species_key %in% projection_species],
        na.rm = TRUE
      ),
      vector_sdm_projection_species_available = dplyr::n_distinct(
        species_name[species_role == "vector" & species_key %in% projection_species],
        na.rm = TRUE
      ),
      host_sdm_comparison_species_available = dplyr::n_distinct(
        species_name[species_role == "host" & species_key %in% comparison_species],
        na.rm = TRUE
      ),
      vector_sdm_comparison_species_available = dplyr::n_distinct(
        species_name[species_role == "vector" & species_key %in% comparison_species],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------|
#      Join Evidence -----------------------------------------------------------|
# ------------------------------------------------------------------------------|
readiness <- base %>%
  join_evidence_by_names(disease_qa_summary, "evidence_qa", candidate_fields) %>%
  join_evidence_by_names(vector_summary, "direct_vector", candidate_fields) %>%
  join_evidence_by_names(genbank_summary, "genbank", candidate_fields) %>%
  join_evidence_by_names(who_don_summary, "who_don", candidate_fields) %>%
  join_evidence_by_names(sdm_summary, "sdm", candidate_fields)

# ------------------------------------------------------------------------------|
#      Derived Status Fields ---------------------------------------------------|
# ------------------------------------------------------------------------------|
count_cols <- c(
  "host_candidate_rows",
  "roster_host_rows",
  "roster_vector_rows",
  "roster_vectors_with_host_vector_evidence",
  "roster_vectors_with_human_biting",
  "roster_taxonomy_caution_rows",
  "host_role_evidence_rows",
  "host_role_evidence_substantive_rows",
  "vector_role_evidence_rows",
  "vector_role_evidence_substantive_rows",
  "host_role_assignment_rows",
  "host_role_substantive_assignment_rows",
  "vector_role_assignment_rows",
  "vector_role_substantive_assignment_rows",
  "direct_disease_vector_rows",
  "direct_vector_species_or_taxa",
  "direct_vector_confirmed_rows",
  "direct_vector_probable_rows",
  "direct_vector_candidate_rows",
  "direct_vector_competence_joined_rows",
  "direct_vector_competent_rows",
  "direct_vector_mixed_competence_rows",
  "direct_vector_not_competent_rows",
  "direct_vector_unclear_competence_rows",
  "direct_vector_transmission_yes_rows",
  "direct_vector_natural_infection_yes_rows",
  "direct_vector_taxonomy_review_needed_rows",
  "genbank_country_rows",
  "genbank_distinct_countries_or_territories",
  "genbank_records_with_country",
  "who_don_focal_country_rows",
  "who_don_distinct_records",
  "who_don_distinct_countries",
  "who_don_high_confidence_rows",
  "who_don_manual_review_rows",
  "host_sdm_species_available",
  "vector_sdm_species_available",
  "host_sdm_projection_species_available",
  "vector_sdm_projection_species_available",
  "host_sdm_comparison_species_available",
  "vector_sdm_comparison_species_available"
)

flag_cols <- c(
  "has_substantive_host_role_evidence",
  "has_substantive_vector_role_evidence"
)

readiness <- readiness %>%
  add_missing_count_cols(count_cols) %>%
  add_missing_flag_cols(flag_cols) %>%
  mutate(
    has_substantive_host_role_evidence =
      has_substantive_host_role_evidence | host_role_evidence_substantive_rows > 0,
    has_substantive_vector_role_evidence =
      has_substantive_vector_role_evidence | vector_role_evidence_substantive_rows > 0,
    has_direct_vector_evidence = direct_disease_vector_rows > 0,
    has_genbank_country_evidence = genbank_country_rows > 0,
    has_who_don_focal_evidence = who_don_focal_country_rows > 0,
    has_any_role_evidence_or_assignment =
      has_substantive_host_role_evidence |
        has_substantive_vector_role_evidence |
        host_role_substantive_assignment_rows > 0 |
        vector_role_substantive_assignment_rows > 0,
    role_assignment_status = case_when(
      !is.na(role_assignment_status) ~ role_assignment_status,
      host_role_substantive_assignment_rows > 0 | vector_role_substantive_assignment_rows > 0 ~
        "reviewed_assignment_present",
      has_substantive_host_role_evidence | has_substantive_vector_role_evidence ~
        "evidence_present_no_assignment",
      TRUE ~ "no_role_evidence_or_assignment"
    ),
    evidence_join_status = case_when(
      if_any(
        all_of(c(
          "evidence_qa_join_status",
          "direct_vector_join_status",
          "genbank_join_status",
          "who_don_join_status",
          "sdm_join_status"
        )),
        ~ .x == "matched_primary_name"
      ) ~ "matched_primary_name",
      if_any(
        all_of(c(
          "evidence_qa_join_status",
          "direct_vector_join_status",
          "genbank_join_status",
          "who_don_join_status",
          "sdm_join_status"
        )),
        ~ .x == "matched_alternate_name"
      ) ~ "matched_alternate_name",
      if_any(
        all_of(c(
          "evidence_qa_join_status",
          "direct_vector_join_status",
          "genbank_join_status",
          "who_don_join_status",
          "sdm_join_status"
        )),
        ~ .x == "ambiguous_not_joined"
      ) ~ "ambiguous_not_joined",
      TRUE ~ "unmatched_no_evidence"
    ),
    direct_vector_evidence_status = case_when(
      has_direct_vector_evidence ~ "direct_vector_evidence_present",
      vectored_status == "vectored" ~ "conceptual_vectored_no_direct_vector_evidence",
      vectored_status == "non_vectored" ~ "not_expected_non_vectored",
      TRUE ~ "mixed_or_uncertain_review"
    ),
    country_evidence_status = case_when(
      has_who_don_focal_evidence & has_genbank_country_evidence ~ "who_don_and_genbank",
      has_who_don_focal_evidence ~ "who_don_only",
      has_genbank_country_evidence ~ "genbank_only",
      TRUE ~ "no_country_evidence"
    ),
    host_sdm_required = host_sdm_needed == "yes",
    vector_sdm_required = vector_sdm_needed == "yes",
    host_sdm_optional = host_sdm_needed %in% c("optional", "uncertain"),
    vector_sdm_optional = vector_sdm_needed %in% c("optional", "uncertain"),
    sdm_availability_status = case_when(
      modelling_scope_status != "include" | range_limiting_layer == "not_model_ready" ~
        "sdm_not_required_or_not_model_ready",
      host_sdm_required & host_sdm_species_available == 0 ~ "missing_required_host_sdm",
      vector_sdm_required & vector_sdm_species_available == 0 ~ "missing_required_vector_sdm",
      (host_sdm_optional & host_sdm_species_available == 0) |
        (vector_sdm_optional & vector_sdm_species_available == 0) ~ "optional_sdms_missing",
      !host_sdm_required & !vector_sdm_required & !host_sdm_optional & !vector_sdm_optional ~
        "sdm_not_required_or_not_model_ready",
      TRUE ~ "required_sdms_available"
    ),
    readiness_blocker = case_when(
      modelling_scope_status != "include" ~ "scope_not_include",
      transmission_rule_review_status != "reviewed" ~ "transmission_rule_not_reviewed",
      vectored_status %in% c("vectored", "mixed_or_indirect") & !has_direct_vector_evidence ~
        "missing_direct_vector_evidence",
      host_sdm_required & host_sdm_species_available == 0 ~ "missing_required_host_sdm",
      vector_sdm_required & has_direct_vector_evidence & vector_sdm_species_available == 0 ~
        "missing_required_vector_sdm",
      !has_any_role_evidence_or_assignment ~ "missing_role_evidence",
      !has_who_don_focal_evidence & !has_genbank_country_evidence ~ "missing_country_evidence",
      TRUE ~ "none"
    ),
    recommended_next_action = case_when(
      modelling_scope_status != "include" ~ "review_or_defer_scope",
      transmission_rule_review_status != "reviewed" ~ "review_transmission_rule",
      vectored_status %in% c("vectored", "mixed_or_indirect") & !has_direct_vector_evidence ~
        "curate_direct_vector_evidence",
      host_sdm_required & host_sdm_species_available == 0 ~ "find_or_build_host_sdm",
      vector_sdm_required & has_direct_vector_evidence & vector_sdm_species_available == 0 ~
        "find_or_build_vector_sdm",
      !has_any_role_evidence_or_assignment ~ "curate_role_evidence",
      !has_who_don_focal_evidence & !has_genbank_country_evidence ~ "review_country_evidence",
      TRUE ~ "ready_for_model_spec_review"
    )
  ) %>%
  select(-host_sdm_required, -vector_sdm_required, -host_sdm_optional, -vector_sdm_optional)

# ------------------------------------------------------------------------------|
#      Output Ordering And Validation -----------------------------------------|
# ------------------------------------------------------------------------------|
output_specs <- readiness_output_specs()
core_cols <- output_specs$core
readiness_cols <- output_specs$readiness
evidence_cols <- output_specs$evidence
sdm_cols <- output_specs$sdm
decision_cols <- output_specs$decision
match_cols <- output_specs$match
slim_cols <- output_specs$slim

readiness_full <- readiness %>%
  arrange(.row_order) %>%
  select(
    any_of(c(core_cols, readiness_cols, evidence_cols, sdm_cols, decision_cols, match_cols)),
    everything()
  ) %>%
  mutate(
    pathogen_species_name = case_when(
      is_species_like_name(source_msl39_viral_name) ~ source_msl39_viral_name,
      is_species_like_name(source_pathogen) ~ source_pathogen,
      is_species_like_name(matched_pathogen_names) ~ matched_pathogen_names,
      TRUE ~ NA_character_
    ),
    pathogen_species_name = format_taxon_name(pathogen_species_name),
    pathogen_species_name_source = case_when(
      is_species_like_name(source_msl39_viral_name) ~ "source_msl39_viral_name",
      is_species_like_name(source_pathogen) ~ "source_pathogen",
      is_species_like_name(matched_pathogen_names) ~ "matched_pathogen_names",
      TRUE ~ "unresolved_or_aggregate"
    ),
    .pathogen_taxid_lookup = purrr::pmap_chr(
      list(pathogen_species_name, source_msl39_viral_name, source_pathogen, matched_pathogen_names),
      ~ lookup_pathogen_taxid(pathogen_taxid_lookup, ...)
    ),
    .pathogen_taxid_lookup_source = purrr::pmap_chr(
      list(pathogen_species_name, source_msl39_viral_name, source_pathogen, matched_pathogen_names),
      ~ lookup_pathogen_taxid(pathogen_taxid_lookup, ..., return_source = TRUE)
    ),
    pathogen_taxid = first_non_empty(matched_taxids, .pathogen_taxid_lookup),
    pathogen_taxid_source = case_when(
      !is.na(matched_taxids) & matched_taxids != "" & !is.na(preferred_match_source) ~
        paste("matched_taxids", preferred_match_source, sep = ":"),
      !is.na(matched_taxids) & matched_taxids != "" ~ "matched_taxids",
      !is.na(.pathogen_taxid_lookup) & .pathogen_taxid_lookup != "" ~ .pathogen_taxid_lookup_source,
      TRUE ~ "unresolved"
    )
  ) %>%
  select(-any_of(c(".pathogen_taxid_lookup", ".pathogen_taxid_lookup_source"))) %>%
  select(-any_of(".row_order"))

required_output_cols <- c(core_cols, readiness_cols, evidence_cols, sdm_cols, decision_cols)
missing_required_cols <- setdiff(required_output_cols, names(readiness_full))
if (length(missing_required_cols) > 0) {
  stop("Readiness output is missing required columns: ", paste(missing_required_cols, collapse = ", "), call. = FALSE)
}

missing_slim_cols <- setdiff(slim_cols, names(readiness_full))
if (length(missing_slim_cols) > 0) {
  stop("Slim readiness output is missing columns: ", paste(missing_slim_cols, collapse = ", "), call. = FALSE)
}

held_analysis_unit_ids <- readiness_full %>%
  filter(include_as_analysis_unit == "hold" | modelling_scope_status == "exclude_hold") %>%
  pull(analysis_unit_id)

readiness_slim <- readiness_full %>%
  filter(!analysis_unit_id %in% held_analysis_unit_ids) %>%
  select(all_of(slim_cols))

pilot_cols <- c(
  "pilot_subset",
  "pilot_priority",
  "pilot_next_step",
  slim_cols,
  "modelling_scope_reason"
)

readiness_pilot <- readiness_full %>%
  filter(
    !analysis_unit_id %in% held_analysis_unit_ids,
    is_true(in_master_who)
  ) %>%
  mutate(
    pilot_subset = case_when(
      modelling_scope_status == "include" ~ "who_core",
      modelling_scope_status == "review_before_modelling" ~ "who_scope_review",
      modelling_scope_status == "defer_broad_or_aggregate_unit" ~ "who_defer_aggregate",
      TRUE ~ "who_other_review"
    ),
    pilot_priority = case_when(
      pilot_subset == "who_core" ~ 1L,
      pilot_subset == "who_scope_review" ~ 2L,
      pilot_subset == "who_defer_aggregate" ~ 3L,
      TRUE ~ 4L
    ),
    pilot_next_step = case_when(
      pilot_subset == "who_scope_review" ~ "decide_scope_before_modelling",
      pilot_subset == "who_defer_aggregate" ~ "defer_or_split_aggregate_unit",
      TRUE ~ recommended_next_action
    )
  ) %>%
  arrange(pilot_priority, analysis_unit_id) %>%
  select(all_of(pilot_cols))

validate_readiness_outputs(
  readiness_full = readiness_full,
  readiness_slim = readiness_slim,
  readiness_pilot = readiness_pilot,
  master = master,
  slim_cols = slim_cols,
  pilot_cols = pilot_cols,
  held_analysis_unit_ids = held_analysis_unit_ids
)

# ------------------------------------------------------------------------------|
#      Pilot Package Tables ----------------------------------------------------|
# ------------------------------------------------------------------------------|
package_context_cols <- c(
  "analysis_unit_id",
  "readiness_disease_name",
  "pathogen_species_name",
  "pathogen_taxid",
  "pilot_subset",
  "pilot_priority",
  "pilot_next_step",
  "modelling_scope_status",
  "recommended_next_action",
  "readiness_blocker",
  "host_sdm_needed",
  "vector_sdm_needed",
  "direct_vector_evidence_status",
  "country_evidence_status",
  "sdm_availability_status"
)

package_child_context_cols <- c(
  "analysis_unit_id",
  "readiness_disease_name"
)

package_match_cols <- names(readiness_full)[
  stringr::str_ends(names(readiness_full), "_matched_disease_name")
]

pilot_context <- readiness_pilot %>%
  select(any_of(package_context_cols)) %>%
  left_join(
    readiness_full %>% select(analysis_unit_id, all_of(package_match_cols)),
    by = "analysis_unit_id"
  )

pilot_roster <- join_pilot_layer(
  roster,
  disease_col = "disease_name",
  layers = c("evidence_qa", "sdm"),
  pilot_context = pilot_context,
  context_cols = package_context_cols
)

pilot_hosts <- pilot_roster %>%
  filter(clean_key(species_role) == "host") %>%
  arrange(pilot_priority, analysis_unit_id, species_name) %>%
  select(any_of(c(
    package_child_context_cols,
    "disease_name",
    "species_name",
    "tax_id",
    "host_class",
    "host_order",
    "host_family",
    "host_detection_method",
    "host_role_assignment",
    "host_role_confidence",
    "host_role_needs_manual_review",
    "host_role_assignment_status",
    "source_pathogens",
    "in_current_role_review_scope",
    "in_gibb_etal",
    "in_empres_i",
    "taxonomy_caution"
  )))

pilot_vectors <- pilot_roster %>%
  filter(clean_key(species_role) == "vector") %>%
  arrange(pilot_priority, analysis_unit_id, species_name) %>%
  select(any_of(c(
    package_child_context_cols,
    "disease_name",
    "species_name",
    "tax_id",
    "vector_group",
    "vector_taxon_rank",
    "vector_join_key",
    "vector_role_assignment",
    "vector_role_confidence",
    "vector_role_needs_manual_review",
    "vector_role_assignment_status",
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
  )))

genbank_country_rows <- join_pilot_layer(
  genbank,
  disease_col = "Disease_name",
  layers = "genbank",
  pilot_context = pilot_context,
  context_cols = package_context_cols
)

pilot_genbank_countries <- genbank_country_rows %>%
  transmute(
    across(any_of(package_child_context_cols)),
    country_source = "genbank",
    source_disease_name = Disease_name,
    country = country_standardized,
    country_status,
    evidence_rows = NA_integer_,
    records_with_country = suppressWarnings(as.numeric(records_with_country)),
    distinct_records = NA_integer_,
    high_confidence_rows = NA_integer_,
    manual_review_rows = NA_integer_,
    source_methods = NA_character_,
    claim_types = NA_character_,
    pathogens,
    target_ids,
    latest_publication_datetime_utc = NA_character_
  )

who_don_country_summary <- if (is.null(who_don)) {
  tibble(disease_label_standard = character(), country_standard = character())
} else {
  who_don %>%
    group_by(disease_label_standard, country_standard) %>%
    summarise(
      evidence_rows = dplyr::n(),
      distinct_records = dplyr::n_distinct(record_key),
      high_confidence_rows = sum(scope_confidence == "high", na.rm = TRUE),
      manual_review_rows = sum(is_true(needs_review), na.rm = TRUE),
      source_methods = collapse_unique(source_method),
      claim_types = collapse_unique(claim_type),
      latest_publication_datetime_utc = collapse_latest(publication_datetime_utc),
      .groups = "drop"
    )
}

who_don_country_rows <- join_pilot_layer(
  who_don_country_summary,
  disease_col = "disease_label_standard",
  layers = "who_don",
  pilot_context = pilot_context,
  context_cols = package_context_cols
)

pilot_who_don_countries <- who_don_country_rows %>%
  transmute(
    across(any_of(package_child_context_cols)),
    country_source = "who_don",
    source_disease_name = disease_label_standard,
    country = country_standard,
    country_status = "reported_by_who_don",
    evidence_rows,
    records_with_country = NA_real_,
    distinct_records,
    high_confidence_rows,
    manual_review_rows,
    source_methods,
    claim_types,
    pathogens = NA_character_,
    target_ids = NA_character_,
    latest_publication_datetime_utc
  )

pilot_countries <- bind_rows(pilot_genbank_countries, pilot_who_don_countries) %>%
  left_join(
    readiness_pilot %>% transmute(analysis_unit_id, .pilot_sort_order = pilot_priority),
    by = "analysis_unit_id"
  ) %>%
  arrange(.pilot_sort_order, analysis_unit_id, country_source, country) %>%
  select(-.pilot_sort_order)

sdm_species_summary <- if (is.null(accessible_sdm_species) || !"species" %in% names(accessible_sdm_species)) {
  tibble(species_key = character())
} else {
  accessible_sdm_species %>%
    mutate(species_key = clean_key(species)) %>%
    group_by(species_key) %>%
    summarise(
      sdm_species = collapse_unique(species),
      sdm_model_file_count = sum(suppressWarnings(as.numeric(model_file_count)), na.rm = TRUE),
      sdm_model_paths = collapse_unique(model_paths),
      sdm_has_top_level_model = any(is_true(has_top_level_model), na.rm = TRUE),
      sdm_has_nested_model = any(is_true(has_nested_model), na.rm = TRUE),
      sdm_nested_contexts = collapse_unique(nested_contexts),
      .groups = "drop"
    )
}

sdm_projection_summary <- if (is.null(sdm_projections)) {
  tibble(species_key = character())
} else {
  sdm_projections %>%
    mutate(
      species_key = clean_key(species),
      scenario = paste(gcm, ssp, time_slice, sep = ":")
    ) %>%
    group_by(species_key) %>%
    summarise(
      sdm_projection_runs = dplyr::n(),
      sdm_projection_success_runs = sum(status == "success", na.rm = TRUE),
      sdm_projection_scenarios = collapse_unique(scenario),
      .groups = "drop"
    )
}

sdm_comparison_summary <- if (is.null(sdm_comparisons)) {
  tibble(species_key = character())
} else {
  sdm_comparisons %>%
    mutate(
      species_key = clean_key(species),
      comparison_scenario = paste(gcm, ssp, time_slice, mode, sep = ":")
    ) %>%
    group_by(species_key) %>%
    summarise(
      sdm_comparison_runs = dplyr::n(),
      sdm_comparison_success_runs = sum(status == "success", na.rm = TRUE),
      sdm_comparison_scenarios = collapse_unique(comparison_scenario),
      .groups = "drop"
    )
}

pilot_sdm_species <- pilot_roster %>%
  mutate(
    .pilot_sort_order = pilot_priority,
    species_key = clean_key(species_name),
    species_role_clean = clean_key(species_role),
    sdm_needed_for_disease = case_when(
      species_role_clean == "host" ~ host_sdm_needed,
      species_role_clean == "vector" ~ vector_sdm_needed,
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(
    across(any_of(c(
      ".pilot_sort_order",
      package_child_context_cols,
      "disease_name",
      "species_name",
      "tax_id",
      "species_role",
      "sdm_needed_for_disease",
      "host_role_assignment",
      "host_role_confidence",
      "host_role_needs_manual_review",
      "host_role_assignment_status",
      "vector_role_assignment",
      "vector_role_confidence",
      "vector_role_needs_manual_review",
      "vector_role_assignment_status",
      "species_key"
    )))
  ) %>%
  left_join(sdm_species_summary, by = "species_key") %>%
  left_join(sdm_projection_summary, by = "species_key") %>%
  left_join(sdm_comparison_summary, by = "species_key") %>%
  mutate(
    sdm_available = !is.na(sdm_species),
    sdm_projection_available = coalesce(sdm_projection_success_runs, 0L) > 0,
    sdm_comparison_available = coalesce(sdm_comparison_success_runs, 0L) > 0
  ) %>%
  arrange(.pilot_sort_order, analysis_unit_id, species_role, species_name) %>%
  select(
    any_of(package_child_context_cols),
    disease_name,
    species_name,
    tax_id,
    species_role,
    sdm_needed_for_disease,
    host_role_assignment,
    host_role_confidence,
    host_role_needs_manual_review,
    host_role_assignment_status,
    vector_role_assignment,
    vector_role_confidence,
    vector_role_needs_manual_review,
    vector_role_assignment_status,
    sdm_available,
    sdm_species
  )

pilot_evidence_summary <- readiness_full %>%
  inner_join(
    readiness_pilot %>% select(analysis_unit_id, pilot_subset, pilot_priority, pilot_next_step),
    by = "analysis_unit_id"
  ) %>%
  arrange(pilot_priority, analysis_unit_id) %>%
  select(any_of(c(
    package_context_cols,
    "vectored_status",
    "generalist_status",
    "transmission_complexity",
    "guild",
    "range_limiting_layer",
    "role_assignment_status",
    "evidence_join_status",
    count_cols,
    flag_cols,
    "has_direct_vector_evidence",
    "has_genbank_country_evidence",
    "has_who_don_focal_evidence",
    package_match_cols
  )))

pilot_package_data_tables <- list(
  disease_modelling_pilot = readiness_pilot,
  pilot_hosts = pilot_hosts,
  pilot_vectors = pilot_vectors,
  pilot_countries = pilot_countries,
  pilot_sdm_species = pilot_sdm_species,
  pilot_evidence_summary = pilot_evidence_summary
)

package_generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
pilot_package_manifest <- build_pilot_package_manifest(
  data_tables = pilot_package_data_tables,
  generated_at_utc = package_generated_at_utc,
  sources = pilot_package_source_descriptions()
)

pilot_package_tables <- c(
  list(manifest = pilot_package_manifest),
  pilot_package_data_tables
)

validate_pilot_package_tables(
  data_tables = pilot_package_data_tables,
  pilot_ids = readiness_pilot$analysis_unit_id
)

# ------------------------------------------------------------------------------|
#      Write Output And Console Summary ---------------------------------------|
# ------------------------------------------------------------------------------|
output_path <- file.path(readiness_dir, "disease_modelling_readiness.csv")
pilot_output_path <- file.path(readiness_dir, "disease_modelling_pilot.csv")
full_output_path <- file.path(readiness_dir, "disease_modelling_readiness_full.csv")
pilot_package_dir <- file.path(readiness_dir, "disease_modelling_pilot_package")
pilot_package_rds_path <- file.path(readiness_dir, "disease_modelling_pilot_package.rds")
pilot_package_xlsx_path <- file.path(readiness_dir, "disease_modelling_pilot_package.xlsx")
readr::write_csv(readiness_slim, output_path, na = "")
readr::write_csv(readiness_pilot, pilot_output_path, na = "")
readr::write_csv(readiness_full, full_output_path, na = "")
dir.create(pilot_package_dir, recursive = TRUE, showWarnings = FALSE)
purrr::iwalk(pilot_package_tables, function(table, table_name) {
  readr::write_csv(table, file.path(pilot_package_dir, paste0(table_name, ".csv")), na = "")
})
writeLines(
  pilot_package_readme_lines(),
  con = file.path(pilot_package_dir, "README.md")
)
saveRDS(pilot_package_tables, pilot_package_rds_path)

pilot_package_xlsx_written <- FALSE
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(pilot_package_tables, pilot_package_xlsx_path)
  pilot_package_xlsx_written <- TRUE
} else {
  warning("Package `writexl` is not available; skipping pilot package XLSX output.", call. = FALSE)
}

message("Wrote disease modelling readiness: ", output_path)
message("Wrote disease modelling pilot: ", pilot_output_path)
message("Wrote full disease modelling readiness audit table: ", full_output_path)
message("Wrote disease modelling pilot package folder: ", pilot_package_dir)
message("Wrote disease modelling pilot package RDS: ", pilot_package_rds_path)
if (pilot_package_xlsx_written) {
  message("Wrote disease modelling pilot package XLSX: ", pilot_package_xlsx_path)
}
message("Rows written: ", nrow(readiness_slim))
message("Slim columns written: ", ncol(readiness_slim))
message("Pilot rows written: ", nrow(readiness_pilot))
message("Pilot columns written: ", ncol(readiness_pilot))
message("Full columns written: ", ncol(readiness_full))
message("Pilot package table rows:")
print(pilot_package_manifest %>% select(table_name, rows, columns))
message("Recommended next actions:")
print(readiness_slim %>% count(recommended_next_action, name = "rows") %>% arrange(desc(rows)))
message("Pilot subsets:")
print(readiness_pilot %>% count(pilot_subset, pilot_priority, name = "rows") %>% arrange(pilot_priority))
message("Direct vector evidence statuses:")
print(readiness_slim %>% count(direct_vector_evidence_status, name = "rows") %>% arrange(desc(rows)))

message(
  "Direct vector evidence disease rows: ",
  sum(readiness_slim$has_direct_vector_evidence, na.rm = TRUE)
)
message(
  "Non-WHO master-list rows with direct vector evidence: ",
  sum(!is_true(readiness_slim$in_master_who) & readiness_slim$has_direct_vector_evidence, na.rm = TRUE)
)
message(
  "Summed WHO DON focal rows in readiness output: ",
  sum(readiness_full$who_don_focal_country_rows, na.rm = TRUE)
)
message(
  "Summed GenBank disease-country rows in readiness output: ",
  sum(readiness_full$genbank_country_rows, na.rm = TRUE)
)
message(
  "GenBank country summary source: ",
  genbank_summary_source,
  " (",
  paths$genbank_disease_country_summary,
  ")"
)

report_unmatched(
  "Evidence QA",
  disease_qa_summary$disease_name,
  readiness_full$evidence_qa_matched_disease_name
)
report_unmatched(
  "Direct vector evidence",
  vector_summary$disease_name,
  readiness_full$direct_vector_matched_disease_name
)
report_unmatched(
  "GenBank",
  genbank_summary$disease_name,
  readiness_full$genbank_matched_disease_name
)
report_unmatched(
  "WHO DON",
  who_don_summary$disease_name,
  readiness_full$who_don_matched_disease_name
)
report_unmatched(
  "SDM roster",
  sdm_summary$disease_name,
  readiness_full$sdm_matched_disease_name
)
