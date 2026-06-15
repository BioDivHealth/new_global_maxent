# -----------------------------------------------------------------------------|
# disease_scope_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for WHO disease-scope stage scripts.
# -----------------------------------------------------------------------------|

disease_scope_provenance_cols <- function() {
  c(
    "is_priority_pathogen",
    "is_prototype_pathogen",
    "in_gibb_etal",
    "in_empres_i",
    "priority_prototype_status",
    "region_africa",
    "region_americas",
    "region_europe",
    "region_mediterranean",
    "region_se_asia",
    "region_western_pacific"
  )
}

disease_scope_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

disease_scope_clean_key <- function(x) {
  key <- disease_scope_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

disease_scope_flag_from_mark <- function(x) {
  x <- disease_scope_clean_text(x)
  !is.na(x)
}

disease_scope_first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

disease_scope_standardize_master_cols <- function(master_raw) {
  master_raw %>%
    dplyr::rename(
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
    dplyr::mutate(
      master_row = dplyr::row_number(),
      dplyr::across(where(is.character), disease_scope_clean_text),
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
      in_master_who = disease_scope_flag_from_mark(master_who_flag),
      in_master_gibb = disease_scope_flag_from_mark(master_gibb_flag),
      in_master_empres_i = disease_scope_flag_from_mark(master_empres_i_flag),
      in_master_atlas = disease_scope_flag_from_mark(master_atlas_flag),
      master_gbif_checked = disease_scope_flag_from_mark(master_gbif_check),
      master_livestock_amplified = disease_scope_flag_from_mark(master_livestock_amplified_flag),
      disease_master_key = disease_scope_clean_key(disease_master_name)
    )
}

disease_scope_build_existing_unit_index <- function(analysis_units) {
  index_fields <- analysis_units %>%
    dplyr::mutate(unit_row = dplyr::row_number()) %>%
    dplyr::select(
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

  dplyr::bind_rows(
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "source_disease_name", match_key = disease_scope_clean_key(source_disease_name)),
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "source_pathogen", match_key = disease_scope_clean_key(source_pathogen)),
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "source_previous_name", match_key = disease_scope_clean_key(source_previous_name)),
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "source_msl39_viral_name", match_key = disease_scope_clean_key(source_msl39_viral_name)),
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "analysis_unit", match_key = disease_scope_clean_key(analysis_unit)),
    index_fields %>%
      dplyr::transmute(unit_row, match_field = "analysis_unit_label", match_key = disease_scope_clean_key(analysis_unit_label))
  ) %>%
    dplyr::filter(!is.na(match_key)) %>%
    dplyr::distinct(match_key, unit_row, .keep_all = TRUE) %>%
    dplyr::left_join(index_fields, by = "unit_row")
}

disease_scope_master_existing_aliases <- function() {
  tibble::tribble(
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
    dplyr::mutate(
      disease_master_key = disease_scope_clean_key(disease_master_name),
      existing_lookup_key = disease_scope_clean_key(existing_lookup_name)
    )
}
