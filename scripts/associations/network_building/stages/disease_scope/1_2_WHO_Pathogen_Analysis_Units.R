# ------------------------------------------------------------------------------
# 1_2_WHO_Pathogen_Analysis_Units.R
# ------------------------------------------------------------------------------
# Purpose: Derive a curation-ready analysis-unit table from the zoonotic WHO
#          shortlist by:
#          1. keeping already specific rows as analysis units
#          2. flagging broad placeholder rows that are too coarse to model
#          3. expanding selected broad taxa into narrower candidate units using
#             a transparent decision framework
#
# Input  : who_zoonotic_pathogens_path()
# Output : who_pathogen_analysis_units_path()
#          who_pathogen_analysis_units_keep_path()
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "disease_scope_helpers.R"
))

who_provenance_cols <- disease_scope_provenance_cols()

classify_source_scope <- function(pathogen) {
  pathogen <- disease_scope_clean_text(pathogen)

  dplyr::case_when(
    pathogen %in% c("Subgenus Sarbecovirus", "Subgenus Merbecovirus") ~ "broad_subgenus",
    pathogen == "Genus Vesiculovirus" ~ "broad_genus",
    pathogen == "Salmonella enterica non typhoidal serovars" ~ "broad_serovar_group",
    pathogen == "Orthoreovirus mammalis" ~ "broad_species_group",
    stringr::str_detect(pathogen, "^Alphainfluenzavirus influenzae \\(") ~ "subtype_level",
    TRUE ~ "specific_unit"
  )
}

classify_analysis_rank <- function(pathogen) {
  scope <- classify_source_scope(pathogen)

  dplyr::case_when(
    scope == "subtype_level" ~ "subtype",
    scope == "specific_unit" ~ "species_or_equivalent",
    scope == "broad_serovar_group" ~ "serovar_group",
    scope == "broad_species_group" ~ "species_group",
    scope == "broad_subgenus" ~ "subgenus",
    scope == "broad_genus" ~ "genus",
    TRUE ~ "review"
  )
}

infer_transmission_context <- function(pathogen, disease_name) {
  pathogen <- disease_scope_clean_text(pathogen)
  disease_name <- disease_scope_clean_text(disease_name)

  vector_associated <- c(
    "Orthonairovirus haemorrhagiae",
    "Phlebovirus riftense",
    "Orthoflavivirus denguei",
    "Orthoflavivirus flavi",
    "Orthoflavivirus zikaense",
    "Orthoflavivirus nilense",
    "Orthoflavivirus encephalitidis",
    "Alphavirus chikungunya",
    "Alphavirus venezuelan",
    "Orthobunyavirus oropoucheense",
    "Bandavirus dabieense",
    "Genus Vesiculovirus"
  )

  dplyr::case_when(
    pathogen %in% vector_associated ~ "vector_associated_or_mixed",
    disease_name %in% c(
      "Plague",
      "Mpox (Monkeypox)"
    ) ~ "mixed_or_context_dependent",
    TRUE ~ "host_linked_non_vector"
  )
}

input_path <- who_pathogens_diseases_zoonotic_path()
output_path <- who_pathogen_analysis_units_path()
output_keep_path <- who_pathogen_analysis_units_keep_path()

split_source_pathogens <- c(
  "Subgenus Sarbecovirus",
  "Subgenus Merbecovirus",
  "Genus Vesiculovirus"
)

review_broad_pathogens <- c(
  "Salmonella enterica non typhoidal serovars",
  "Orthoreovirus mammalis"
)

who_zoonotic <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), disease_scope_clean_text))

base_units <- who_zoonotic %>%
  transmute(
    row_type = "source_row",
    family = Family,
    pheic_risk = `PHEIC risk`,
    source_pathogen = Pathogens,
    source_previous_name = previous_name,
    source_msl39_viral_name = msl39_viral_name,
    source_disease_name = Disease_name,
    is_priority_pathogen,
    is_prototype_pathogen,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    region_africa,
    region_americas,
    region_europe,
    region_mediterranean,
    region_se_asia,
    region_western_pacific,
    source_unit_scope = classify_source_scope(Pathogens),
    analysis_unit = Pathogens,
    analysis_unit_label = dplyr::coalesce(previous_name, msl39_viral_name, Pathogens),
    analysis_unit_rank = classify_analysis_rank(Pathogens),
    analysis_decision = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "drop_broad_placeholder",
      Pathogens == "Alphainfluenzavirus influenzae (H5N1)" ~ "keep",
      stringr::str_detect(Pathogens, "^Alphainfluenzavirus influenzae \\(") ~ "drop_not_prioritized",
      Pathogens %in% review_broad_pathogens ~ "review_split_needed",
      TRUE ~ "keep"
    ),
    decision_rule_trigger = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "broad_taxon_replaced_by_candidate_units",
      Pathogens == "Alphainfluenzavirus influenzae (H5N1)" ~ "influenza_keep_h5n1_only",
      stringr::str_detect(Pathogens, "^Alphainfluenzavirus influenzae \\(") ~ "influenza_drop_currently_human_to_human_or_not_prioritized",
      Pathogens == "Salmonella enterica non typhoidal serovars" ~ "broad_serovar_group",
      Pathogens == "Orthoreovirus mammalis" ~ "mixed_species_group",
      TRUE ~ "specific_source_row"
    ),
    transmission_context = infer_transmission_context(Pathogens, Disease_name),
    human_infection_status = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "broad_placeholder",
      stringr::str_detect(Pathogens, "^Alphainfluenzavirus influenzae \\(") ~ "yes",
      TRUE ~ "yes"
    ),
    host_link_status = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "review",
      TRUE ~ "yes_or_expected_from_source"
    ),
    vector_data_status = dplyr::case_when(
      infer_transmission_context(Pathogens, Disease_name) == "vector_associated_or_mixed" ~ "yes_or_expected",
      Pathogens %in% split_source_pathogens ~ "review",
      TRUE ~ "not_central_or_unknown"
    ),
    amplifier_data_status = dplyr::case_when(
      Pathogens %in% c("Genus Vesiculovirus", "Subgenus Sarbecovirus", "Subgenus Merbecovirus") ~ "review",
      TRUE ~ "unknown"
    ),
    example_members = NA_character_,
    rationale = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "Broad placeholder row retained for provenance only; use narrower candidate rows below.",
      Pathogens == "Alphainfluenzavirus influenzae (H5N1)" ~ "Retained as the only active influenza analysis unit from the current WHO shortlist, following the decision to keep H5N1 but drop the other current influenza subtype rows.",
      stringr::str_detect(Pathogens, "^Alphainfluenzavirus influenzae \\(") ~ "Dropped from the active keep set following the decision to remove the other current influenza subtype rows and only keep H5N1 plus an added H7N9 unit.",
      Pathogens == "Salmonella enterica non typhoidal serovars" ~ "Broad serovar group kept for now but flagged as too coarse for final ecological analysis.",
      Pathogens == "Orthoreovirus mammalis" ~ "Broad mixed-unit row kept for now but flagged as too coarse for final ecological analysis.",
      TRUE ~ "Specific zoonotic source row retained as the current analysis unit."
    ),
    notes = dplyr::case_when(
      Pathogens %in% split_source_pathogens ~ "Do not use this broad placeholder as the final modelling or host/vector/amplifier unit.",
      Pathogens == "Alphainfluenzavirus influenzae (H5N1)" ~ "Kept as the only influenza subtype from the current WHO shortlist in the active analysis-unit set.",
      stringr::str_detect(Pathogens, "^Alphainfluenzavirus influenzae \\(") ~ "Excluded from the active analysis-unit keep file under the current influenza decision rule.",
      TRUE ~ NA_character_
    )
  )

split_candidates <- tibble::tribble(
  ~source_pathogen, ~analysis_unit, ~analysis_unit_label, ~analysis_unit_rank, ~analysis_decision, ~decision_rule_trigger, ~transmission_context, ~human_infection_status, ~host_link_status, ~vector_data_status, ~amplifier_data_status, ~example_members, ~rationale, ~notes,
  "Subgenus Sarbecovirus", "Severe acute respiratory syndrome-related coronavirus (SARS-CoV-1)", "SARS-CoV-1", "species", "keep", "clear_human_infection_and_host_link", "host_linked_non_vector", "yes", "yes", "no", "unknown", NA_character_, "Clear human-pathogenic sarbecovirus with well-established host-link relevance.", "Retained as a distinct analysis unit rather than folded into the broader sarbecovirus placeholder.",
  "Subgenus Sarbecovirus", "Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)", "SARS-CoV-2", "species", "keep", "clear_human_infection_and_host_link", "host_linked_non_vector", "yes", "yes", "no", "unknown", NA_character_, "Clear human-pathogenic sarbecovirus with substantial host-link evidence.", "Retained as a distinct analysis unit rather than folded into the broader sarbecovirus placeholder.",
  "Subgenus Sarbecovirus", "SARS-like bat sarbecoviruses", "SARS-like bat sarbecoviruses", "species_group", "keep", "host_link_and_zoonotic_potential_group", "host_linked_non_vector", "zoonotic_potential", "yes", "no", "unknown", "RaTG13; RmYN02; WIV1; SHC014; Rs3367; RsSHC014; Rs4028; BtKY72; PDF-2370", "Retained as an active wildlife sarbecovirus analysis group because host-link data and zoonotic concern exist even where direct human infection is unclear or absent.", "Treat the listed viruses as supporting examples inside this group rather than splitting further for now.",
  "Subgenus Merbecovirus", "Middle East respiratory syndrome-related coronavirus (MERS-CoV)", "MERS-CoV", "species", "keep", "clear_human_infection_and_host_link", "host_linked_non_vector", "yes", "yes", "no", "yes", NA_character_, "Clear human-pathogenic merbecovirus with established host and amplifier relevance.", "Retained as a distinct analysis unit rather than folded into the broader merbecovirus placeholder.",
  "Subgenus Merbecovirus", "MERS-like bat merbecoviruses", "MERS-like bat merbecoviruses", "species_group", "keep", "host_link_and_zoonotic_potential_group", "host_linked_non_vector", "zoonotic_potential", "yes", "no", "unknown", "HKU4; HKU5; PDF-2180; BtVs-BetaCoV/SC2013; NeoCoV", "Retained as an active wildlife merbecovirus analysis group because host-link data and zoonotic concern exist even where direct human infection is unclear or absent.", "Treat the listed viruses as supporting examples inside this group rather than splitting further for now.",
  "Genus Vesiculovirus", "Vesicular stomatitis Indiana virus", "VSIV", "species", "keep", "clear_human_infection_and_vector_data", "vector_associated_or_mixed", "yes", "yes", "yes", "yes", NA_character_, "Clear zoonotic vesiculovirus with host, vector, and amplifier relevance.", "Retained as a distinct analysis unit rather than folded into the broader vesiculovirus placeholder.",
  "Genus Vesiculovirus", "Vesicular stomatitis New Jersey virus", "VSNJV", "species", "keep", "clear_human_infection_and_vector_data", "vector_associated_or_mixed", "yes", "yes", "yes", "yes", NA_character_, "Clear zoonotic vesiculovirus with host, vector, and amplifier relevance.", "Retained as a distinct analysis unit rather than folded into the broader vesiculovirus placeholder.",
  "Genus Vesiculovirus", "Vesicular stomatitis Alagoas virus", "VSAV", "species", "drop_not_prioritized", "supporting_example_only_for_current_vesiculovirus_scope", "vector_associated_or_mixed", "limited_or_unclear", "yes", "yes", "yes", NA_character_, "Not retained as an active analytic unit in the current vesiculovirus scope; track as a supporting example only.", "Do not split further into a standalone active unit for now.",
  "Genus Vesiculovirus", "Cocal virus", "Cocal virus", "species", "drop_not_prioritized", "supporting_example_only_for_current_vesiculovirus_scope", "vector_associated_or_mixed", "limited_or_unclear", "yes", "partial_or_unknown", "unknown", NA_character_, "Not retained as an active analytic unit in the current vesiculovirus scope; track as a supporting example only.", "Do not split further into a standalone active unit for now.",
  "Genus Vesiculovirus", "Maraba virus", "Maraba virus", "species", "drop_not_prioritized", "supporting_example_only_for_current_vesiculovirus_scope", "vector_associated_or_mixed", "documented_but_limited", "partial_or_unknown", "partial_or_unknown", "unknown", NA_character_, "Not retained as an active analytic unit in the current vesiculovirus scope; track as a supporting example only.", "Do not split further into a standalone active unit for now."
)

expanded_units <- who_zoonotic %>%
  transmute(
    family = Family,
    pheic_risk = `PHEIC risk`,
    source_pathogen = Pathogens,
    source_previous_name = previous_name,
    source_msl39_viral_name = msl39_viral_name,
    source_disease_name = Disease_name,
    is_priority_pathogen,
    is_prototype_pathogen,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    region_africa,
    region_americas,
    region_europe,
    region_mediterranean,
    region_se_asia,
    region_western_pacific,
    source_unit_scope = classify_source_scope(Pathogens)
  ) %>%
  inner_join(split_candidates, by = "source_pathogen") %>%
  mutate(row_type = "split_candidate") %>%
  select(
    row_type,
    family,
    pheic_risk,
    source_pathogen,
    source_previous_name,
    source_msl39_viral_name,
    source_disease_name,
    all_of(who_provenance_cols),
    source_unit_scope,
    analysis_unit,
    analysis_unit_label,
    analysis_unit_rank,
    analysis_decision,
    decision_rule_trigger,
    transmission_context,
    human_infection_status,
    host_link_status,
    vector_data_status,
    amplifier_data_status,
    example_members,
    rationale,
    notes
  )

manual_additions <- tibble::tribble(
  ~row_type, ~family, ~pheic_risk, ~source_pathogen, ~source_previous_name, ~source_msl39_viral_name, ~source_disease_name, ~is_priority_pathogen, ~is_prototype_pathogen, ~in_gibb_etal, ~in_empres_i, ~priority_prototype_status, ~region_africa, ~region_americas, ~region_europe, ~region_mediterranean, ~region_se_asia, ~region_western_pacific, ~source_unit_scope, ~analysis_unit, ~analysis_unit_label, ~analysis_unit_rank, ~analysis_decision, ~decision_rule_trigger, ~transmission_context, ~human_infection_status, ~host_link_status, ~vector_data_status, ~amplifier_data_status, ~example_members, ~rationale, ~notes,
  "manual_addition", "Orthomyxoviridae", "High", "Alphainfluenzavirus influenzae (H7N9)", "Influenza A", "Alphainfluenzavirus influenzae", "Influenza (H7N9 avian influenza)", TRUE, FALSE, FALSE, FALSE, "priority", "none", "none", "none", "none", "none", "none", "subtype_level", "Alphainfluenzavirus influenzae (H7N9)", "Influenza A (H7N9)", "subtype", "keep", "manual_addition_from_decision_framework", "host_linked_non_vector", "yes", "yes_or_expected_from_source", "not_central_or_unknown", "unknown", NA_character_, "Added as an active influenza analysis unit following the decision to keep H5N1 and add H7N9 while dropping the other current influenza subtype rows.", "Manual addition requested in the decision framework even though H7N9 is not present as a source row in the current WHO zoonotic shortlist."
)

analysis_units <- bind_rows(base_units, expanded_units, manual_additions) %>%
  mutate(
    row_sort = dplyr::case_when(
      row_type == "source_row" ~ 1L,
      row_type == "split_candidate" ~ 2L,
      row_type == "manual_addition" ~ 3L,
      TRUE ~ 4L
    ),
    decision_sort = dplyr::case_when(
      analysis_decision == "keep" ~ 1L,
      analysis_decision == "review" ~ 2L,
      analysis_decision == "review_split_needed" ~ 3L,
      analysis_decision == "drop_broad_placeholder" ~ 4L,
      analysis_decision == "drop_not_prioritized" ~ 5L,
      TRUE ~ 6L
    )
  ) %>%
  arrange(source_disease_name, source_pathogen, row_sort, decision_sort, analysis_unit) %>%
  select(-row_sort, -decision_sort)

analysis_units_keep <- analysis_units %>%
  filter(analysis_decision == "keep") %>%
  transmute(
    family,
    pheic_risk,
    source_pathogen,
    source_previous_name,
    source_msl39_viral_name,
    source_disease_name,
    is_priority_pathogen,
    is_prototype_pathogen,
    in_gibb_etal,
    in_empres_i,
    priority_prototype_status,
    region_africa,
    region_americas,
    region_europe,
    region_mediterranean,
    region_se_asia,
    region_western_pacific,
    analysis_unit,
    analysis_unit_label,
    analysis_unit_rank,
    transmission_context,
    human_infection_status,
    host_link_status,
    vector_data_status,
    amplifier_data_status,
    example_members,
    do_not_split_further_yet = TRUE,
    keep_reason = decision_rule_trigger,
    rationale,
    notes
  ) %>%
  distinct(analysis_unit, .keep_all = TRUE) %>%
  arrange(family, source_disease_name, analysis_unit)

write_csv(analysis_units, output_path, na = "")
write_csv(analysis_units_keep, output_keep_path, na = "")

cat("Source zoonotic rows:", nrow(who_zoonotic), "\n")
cat("Derived analysis-unit rows written:", nrow(analysis_units), "\n")
cat("Rows kept as analysis units:", sum(analysis_units$analysis_decision == "keep"), "\n")
cat("Rows flagged for review:", sum(analysis_units$analysis_decision %in% c("review", "review_split_needed")), "\n")
cat("Broad placeholder rows retained for provenance only:", sum(analysis_units$analysis_decision == "drop_broad_placeholder"), "\n")
cat("Active keep rows written:", nrow(analysis_units_keep), "\n")
cat("Wrote analysis-unit table to", output_path, "\n")
cat("Wrote active keep table to", output_keep_path, "\n")
