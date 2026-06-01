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

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

clean_key <- function(x) {
  x %>%
    clean_text() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("&", " and ") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

flag_from_mark <- function(x) {
  x <- clean_text(x)
  !is.na(x)
}

standardize_master_cols <- function(master_raw) {
  master_raw %>%
    rename(
      disease_master_name = Disease,
      pathogen_family_master = `Pathogen family`,
      master_who_flag = WHO,
      master_gibb_flag = Gibb,
      master_empres_i_flag = `EMPRES-i`,
      master_atlas_flag = Atlas,
      master_gbif_check = `GBIF\ncheck`,
      master_guild = Guild,
      master_livestock_amplified_flag = `Livestock\namplified`,
      master_tier = Tier,
      master_key_host_vector = `Key host/vector`,
      master_notes = Notes
    ) %>%
    mutate(
      master_row = row_number(),
      across(where(is.character), clean_text),
      is_section_header = is.na(pathogen_family_master) &
        is.na(master_who_flag) &
        is.na(master_gibb_flag) &
        is.na(master_empres_i_flag) &
        is.na(master_atlas_flag) &
        is.na(master_gbif_check) &
        is.na(master_guild) &
        is.na(master_livestock_amplified_flag) &
        is.na(master_tier) &
        is.na(master_key_host_vector) &
        is.na(master_notes),
      in_master_who = flag_from_mark(master_who_flag),
      in_master_gibb = flag_from_mark(master_gibb_flag),
      in_master_empres_i = flag_from_mark(master_empres_i_flag),
      in_master_atlas = flag_from_mark(master_atlas_flag),
      master_gbif_checked = flag_from_mark(master_gbif_check),
      master_livestock_amplified = flag_from_mark(master_livestock_amplified_flag),
      disease_master_key = clean_key(disease_master_name)
    )
}

build_existing_unit_index <- function(analysis_units) {
  index_fields <- analysis_units %>%
    mutate(unit_row = row_number()) %>%
    select(
      unit_row,
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
    )

  bind_rows(
    index_fields %>%
      transmute(unit_row, match_field = "source_disease_name", match_key = clean_key(source_disease_name)),
    index_fields %>%
      transmute(unit_row, match_field = "source_pathogen", match_key = clean_key(source_pathogen)),
    index_fields %>%
      transmute(unit_row, match_field = "source_previous_name", match_key = clean_key(source_previous_name)),
    index_fields %>%
      transmute(unit_row, match_field = "source_msl39_viral_name", match_key = clean_key(source_msl39_viral_name)),
    index_fields %>%
      transmute(unit_row, match_field = "analysis_unit", match_key = clean_key(analysis_unit)),
    index_fields %>%
      transmute(unit_row, match_field = "analysis_unit_label", match_key = clean_key(analysis_unit_label))
  ) %>%
    filter(!is.na(match_key)) %>%
    distinct(match_key, unit_row, .keep_all = TRUE) %>%
    left_join(index_fields, by = "unit_row")
}

# Explicit aliases for common acronyms and short disease labels in the disease master list.
# These map only into already curated WHO analysis-unit labels; unmatched rows are
# deliberately left for review instead of inventing pathogen names here.
master_existing_aliases <- tibble::tribble(
  ~disease_master_name, ~existing_lookup_name,
  "CCHF", "Crimean-Congo hemorrhagic fever",
  "Rift Valley fever", "Rift Valley fever",
  "Ebola", "Ebola virus disease",
  "Marburg", "Marburg virus disease",
  "Lassa", "Lassa fever",
  "Nipah", "Nipah virus disease",
  "Hendra", "Hendra virus disease",
  "MERS", "MERS-CoV",
  "Mpox", "Mpox (Monkeypox)",
  "Avian influenza (H5N1)", "Alphainfluenzavirus influenzae (H5N1)",
  "Plague", "Plague",
  "Oropouche", "Oropouche fever",
  "HCPS (hantaviruses)", "Hantavirus pulmonary syndrome",
  "Argentine HF (Junin)", "Argentine hemorrhagic fever",
  "Sarbecoviruses (SARS-like)", "Subgenus Sarbecovirus",
  "Tick-borne encephalitis", "Tick-borne encephalitis",
  "SFTS (Bandavirus)", "Severe fever with thrombocytopenia syndrome (SFTS)",
  "VEE", "Venezuelan equine encephalitis",
  "Borna disease", "Borna disease (encephalitis)",
  "Lujo HF", "Lujo hemorrhagic fever",
  "Hepatitis E (zoonotic)", "Hepatitis E",
  "Dengue", "Dengue",
  "Zika", "Zika virus disease",
  "Chikungunya", "Chikungunya fever",
  "Yellow fever", "Yellow fever",
  "West Nile", "West Nile fever"
) %>%
  mutate(
    disease_master_key = clean_key(disease_master_name),
    existing_lookup_key = clean_key(existing_lookup_name)
  )

input_master_path <- here("dr", "disease_master_list_v2.xlsx")
analysis_units_path <- who_pathogen_analysis_units_path()
output_path <- who_master_disease_analysis_units_path()
review_path <- who_diseases_staged_master_expansion_path(
  "master_disease_name_resolution_review.csv"
)

master <- readxl::read_excel(input_master_path, sheet = "Disease Master List") %>%
  standardize_master_cols()

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
  mutate(across(where(is.character), clean_text))

existing_index <- build_existing_unit_index(analysis_units)

master_with_lookup <- master_disease_rows %>%
  left_join(
    master_existing_aliases %>% select(disease_master_key, existing_lookup_name, existing_lookup_key),
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
