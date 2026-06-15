# ------------------------------------------------------------------------------
# 1_2b_Disease_Master_Analysis_Units.R
# ------------------------------------------------------------------------------
# Purpose: Combine the disease master list with the current WHO analysis-unit
#          table without overwriting either source.
#
#          This script creates an additive scaffold for "all analysis units":
#          - rows already represented in the WHO analysis-unit table inherit the
#            resolved pathogen/analysis-unit fields from that table
#          - rows that are new to the WHO shortlist are retained as disease-level
#            review rows, preserving master-list guild/tier/source metadata for later
#            pathogen-name resolution
#
# Input  : dr/disease_master_list_v2.xlsx
#          who_pathogen_analysis_units_path()
# Output : who_master_disease_analysis_units_path()
#          who_diseases_staged_master_expansion_path(
#            "master_disease_name_resolution_review.csv"
#          )
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, readxl, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "disease_scope_helpers.R"
))

input_master_path <- here("dr", "disease_master_list_v2.xlsx")
analysis_units_path <- who_pathogen_analysis_units_path()
output_path <- who_master_disease_analysis_units_path()
review_path <- who_diseases_staged_master_expansion_path(
  "master_disease_name_resolution_review.csv"
)

master <- readxl::read_excel(input_master_path, sheet = "Disease Master List") %>%
  disease_scope_standardize_master_cols()

other_disease_rows <- master %>%
  filter(
    !is_section_header,
    str_detect(disease_master_name, regex("\\bOther\\b", ignore_case = TRUE))
  )

master_disease_rows <- master %>%
  filter(!is_section_header) %>%
  filter(!str_detect(disease_master_name, regex("\\bOther\\b", ignore_case = TRUE))) %>%
  select(
    master_row,
    disease_master_name,
    disease_master_key,
    pathogen_family_master,
    in_master_who,
    in_master_gibb,
    in_master_empres_i,
    in_master_atlas,
    master_gbif_checked,
    master_guild,
    master_livestock_amplified,
    master_tier,
    master_key_host_vector,
    master_notes
  )

analysis_units <- readr::read_csv(
  analysis_units_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), disease_scope_clean_text))

existing_index <- disease_scope_build_existing_unit_index(analysis_units)

master_with_lookup <- master_disease_rows %>%
  left_join(
    disease_scope_master_existing_aliases() %>%
      select(disease_master_key, existing_lookup_name, existing_lookup_key),
    by = "disease_master_key"
  ) %>%
  mutate(primary_lookup_key = coalesce(existing_lookup_key, disease_master_key))

matched_rows <- master_with_lookup %>%
  left_join(
    existing_index,
    by = c("primary_lookup_key" = "match_key"),
    relationship = "many-to-many"
  ) %>%
  group_by(master_row) %>%
  arrange(
    desc(analysis_decision == "keep"),
    match_field,
    analysis_unit,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

combined_units <- matched_rows %>%
  mutate(
    master_list_source = "disease_master_list_v2",
    name_resolution_status = case_when(
      !is.na(unit_row) & !is.na(existing_lookup_key) ~ "matched_existing_unit_by_manual_alias",
      !is.na(unit_row) ~ "matched_existing_unit_by_exact_key",
      TRUE ~ "needs_pathogen_name_resolution"
    ),
    combined_row_type = case_when(
      !is.na(unit_row) ~ "existing_who_analysis_unit",
      TRUE ~ "master_list_review_unit"
    ),
    family = coalesce(family, pathogen_family_master),
    source_disease_name = coalesce(source_disease_name, disease_master_name),
    source_pathogen = source_pathogen,
    source_unit_scope = coalesce(source_unit_scope, "disease_level_review"),
    analysis_unit = coalesce(analysis_unit, disease_master_name),
    analysis_unit_label = coalesce(analysis_unit_label, disease_master_name),
    analysis_unit_rank = coalesce(analysis_unit_rank, "disease_or_syndrome_review"),
    analysis_decision = coalesce(analysis_decision, "review_name_resolution"),
    decision_rule_trigger = coalesce(decision_rule_trigger, "master_list_new_or_unmatched"),
    transmission_context = coalesce(transmission_context, master_guild),
    amplifier_data_status = case_when(
      !is.na(amplifier_data_status) ~ amplifier_data_status,
      master_livestock_amplified ~ "livestock_amplifier_flag_from_master_list",
      TRUE ~ "unknown"
    ),
    rationale = coalesce(
      rationale,
      "Disease appears in the disease master list but has not yet been resolved to a curated pathogen analysis unit."
    ),
    notes = coalesce(notes, master_notes)
  ) %>%
  select(
    master_list_source,
    master_row,
    disease_master_name,
    pathogen_family_master,
    in_master_who,
    in_master_gibb,
    in_master_empres_i,
    in_master_atlas,
    master_gbif_checked,
    master_guild,
    master_livestock_amplified,
    master_tier,
    master_key_host_vector,
    master_notes,
    name_resolution_status,
    combined_row_type,
    existing_lookup_name,
    match_field,
    row_type,
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
  ) %>%
  arrange(
    dplyr::case_when(
      master_tier == "1" ~ 1L,
      master_tier == "1/2" ~ 2L,
      master_tier == "1?" ~ 3L,
      master_tier == "2" ~ 4L,
      master_tier == "2/3" ~ 5L,
      master_tier == "3" ~ 6L,
      TRUE ~ 7L
    ),
    master_row
  )

review_rows <- combined_units %>%
  filter(name_resolution_status == "needs_pathogen_name_resolution") %>%
  select(
    master_row,
    disease_master_name,
    pathogen_family_master,
    in_master_who,
    in_master_gibb,
    in_master_empres_i,
    in_master_atlas,
    master_guild,
    master_livestock_amplified,
    master_tier,
    master_key_host_vector,
    master_notes,
    suggested_next_step = decision_rule_trigger
  )

readr::write_csv(combined_units, output_path, na = "")
readr::write_csv(review_rows, review_path, na = "")

cat("Disease master rows:", nrow(master_disease_rows), "\n")
cat("Dropped ambiguous 'Other' disease rows:", nrow(other_disease_rows), "\n")
if (nrow(other_disease_rows) > 0) {
  print(other_disease_rows %>% select(master_row, disease_master_name), n = Inf)
}
cat("Rows matched to existing WHO analysis units:", sum(combined_units$combined_row_type == "existing_who_analysis_unit"), "\n")
cat("Rows needing pathogen-name resolution:", nrow(review_rows), "\n")
cat("Wrote combined master analysis-unit scaffold to:\n")
cat(output_path, "\n")
cat("Wrote name-resolution review rows to:\n")
cat(review_path, "\n")
