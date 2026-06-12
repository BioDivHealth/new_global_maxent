# ------------------------------------------------------------------------------
# 1_2d_Master_WHO_Analysis_Unit_Bridge.R
# ------------------------------------------------------------------------------
# Purpose: Build compact additive tables that let disease-master analysis units
#          be used alongside the existing WHO analysis-unit tables.
#
# Inputs : who_master_disease_analysis_units_path()
#          who_diseases_name_resolution_path(
#            "master_disease_name_resolution_manual.csv"
#          )
#          who_diseases_staged_master_expansion_path(
#            "master_pathogen_virion_clover_matches.csv"
#          )
#          who_diseases_transmission_rules_path(
#            "master_plus_who_transmission_rules_manual.csv"
#          ) (optional)
#          who_pathogen_analysis_units_path()
#
# Outputs: who_master_plus_analysis_units_path()
#          who_diseases_host_query_path("master_pathogen_host_query_units.csv")
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

source(here("scripts", "associations", "working_inputs.R"))

master_units_path <- who_master_disease_analysis_units_path()
manual_path <- who_diseases_name_resolution_path(
  "master_disease_name_resolution_manual.csv"
)
matches_path <- who_diseases_staged_master_expansion_path(
  "master_pathogen_virion_clover_matches.csv"
)
who_units_path <- who_pathogen_analysis_units_path()
transmission_rules_path <- who_diseases_transmission_rules_path(
  "master_plus_who_transmission_rules_manual.csv"
)

combined_output_path <- who_master_plus_analysis_units_path()
host_query_output_path <- who_diseases_host_query_path(
  "master_pathogen_host_query_units.csv"
)

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

clean_key <- function(x) {
  x %>%
    clean_text() %>%
    str_to_lower() %>%
    str_replace_all("&", " and ") %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
}

coalesce_chr <- function(...) {
  coalesce(!!!map(list(...), as.character))
}

pick_preferred <- function(preferred_source, virion_value, clover_value) {
  case_when(
    preferred_source == "virion" ~ as.character(virion_value),
    preferred_source == "clover" ~ as.character(clover_value),
    TRUE ~ NA_character_
  )
}

transmission_rule_columns <- c(
  "analysis_unit_id",
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

required_paths <- c(master_units_path, manual_path, matches_path, who_units_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = "; "))
}

master_units <- read_csv(master_units_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

manual_units <- read_csv(manual_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    across(where(is.character), clean_text),
    analysis_unit_id = paste0("master_", master_row)
  ) %>%
  select(
    master_row,
    manual_analysis_unit_id = analysis_unit_id,
    manual_resolved_disease_name = resolved_disease_name,
    manual_resolved_pathogen_name = resolved_pathogen_name,
    manual_resolved_pathogen_rank = resolved_pathogen_rank,
    manual_resolved_taxid = resolved_taxid,
    pathogen_aliases,
    include_as_analysis_unit,
    split_group,
    resolution_source,
    resolution_notes
  )

master_matches <- read_csv(matches_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    across(where(is.character), clean_text),
    across(
      c(
        clover_matched_taxids,
        virion_matched_taxids,
        clover_best_match_type,
        virion_best_match_type,
        clover_match_status,
        virion_match_status,
        preferred_match_source,
        overall_match_status,
        preferred_source_match_status,
        match_review_notes
      ),
      as.character
    )
  ) %>%
  select(
    analysis_unit_id,
    master_row,
    match_resolved_disease_name = resolved_disease_name,
    match_resolved_pathogen_name = resolved_pathogen_name,
    match_resolved_pathogen_rank = resolved_pathogen_rank,
    clover_matched_pathogen_names,
    virion_matched_pathogen_names,
    clover_matched_taxids,
    virion_matched_taxids,
    clover_matched_families,
    virion_matched_families,
    clover_matched_source_types,
    virion_matched_source_types,
    clover_best_match_type,
    virion_best_match_type,
    clover_match_status,
    virion_match_status,
    preferred_match_source,
    overall_match_status,
    preferred_source_match_status,
    match_review_flag,
    shared_species_proxy_flag,
    match_review_notes
  )

who_units <- read_csv(who_units_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    across(where(is.character), clean_text),
    who_unit_row = row_number(),
    who_analysis_unit_key = clean_key(analysis_unit),
    who_analysis_unit_label_key = clean_key(analysis_unit_label),
    who_source_disease_key = clean_key(source_disease_name),
    who_source_pathogen_key = clean_key(source_pathogen)
  )

transmission_rules <- tibble(
  analysis_unit_id = character(),
  vectored_status = character(),
  generalist_status = character(),
  transmission_complexity = character(),
  guild = character(),
  host_sdm_needed = character(),
  vector_sdm_needed = character(),
  host_range_rule = character(),
  vector_range_rule = character(),
  range_limiting_layer = character(),
  transmission_rule_notes = character(),
  transmission_rule_review_status = character(),
  modelling_scope_status = character(),
  modelling_scope_reason = character()
)

if (file.exists(transmission_rules_path)) {
  transmission_rules_raw <- read_csv(transmission_rules_path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text))

  missing_transmission_cols <- setdiff(transmission_rule_columns, names(transmission_rules_raw))
  if (length(missing_transmission_cols) > 0) {
    stop(
      "master_plus_who_transmission_rules_manual.csv missing required columns: ",
      paste(missing_transmission_cols, collapse = ", ")
    )
  }

  transmission_rules <- transmission_rules_raw %>%
    select(all_of(transmission_rule_columns)) %>%
    distinct(analysis_unit_id, .keep_all = TRUE)
}

bridge <- master_units %>%
  left_join(manual_units, by = "master_row") %>%
  left_join(master_matches, by = "master_row") %>%
  mutate(
    analysis_unit_id = coalesce(analysis_unit_id, manual_analysis_unit_id, paste0("master_", master_row)),
    bridge_source = "disease_master_list",
    bridge_row_type = case_when(
      combined_row_type == "existing_who_analysis_unit" & !is.na(manual_resolved_pathogen_name) ~ "master_and_who_manually_resolved",
      combined_row_type == "existing_who_analysis_unit" ~ "master_and_who_existing_resolution",
      !is.na(manual_resolved_pathogen_name) ~ "master_manual_resolution",
      TRUE ~ "master_unresolved"
    ),
    bridge_to_existing_who = combined_row_type == "existing_who_analysis_unit",
    active_master_analysis_unit = include_as_analysis_unit == "yes",
    resolved_disease_name_final = coalesce_chr(
      match_resolved_disease_name,
      manual_resolved_disease_name,
      source_disease_name,
      disease_master_name
    ),
    resolved_pathogen_name_final = coalesce_chr(
      match_resolved_pathogen_name,
      manual_resolved_pathogen_name,
      analysis_unit,
      source_pathogen
    ),
    resolved_pathogen_rank_final = coalesce_chr(
      match_resolved_pathogen_rank,
      manual_resolved_pathogen_rank,
      analysis_unit_rank
    ),
    include_status_final = case_when(
      !is.na(include_as_analysis_unit) ~ include_as_analysis_unit,
      analysis_decision == "keep" ~ "yes_existing_who",
      analysis_decision == "review_name_resolution" ~ "review",
      TRUE ~ "not_reviewed"
    ),
    host_query_source = preferred_match_source,
    host_query_pathogen_names = pick_preferred(
      preferred_match_source,
      virion_matched_pathogen_names,
      clover_matched_pathogen_names
    ),
    host_query_taxids = pick_preferred(
      preferred_match_source,
      virion_matched_taxids,
      clover_matched_taxids
    ),
    host_query_include_default = active_master_analysis_unit &
      preferred_source_match_status == "preferred_source_matched" &
      !coalesce(match_review_flag, FALSE) &
      !coalesce(shared_species_proxy_flag, FALSE),
    host_query_bucket = case_when(
      include_status_final == "hold" ~ "hold",
      include_status_final == "review" ~ "manual_review",
      !active_master_analysis_unit ~ "inactive_or_not_reviewed",
      coalesce(shared_species_proxy_flag, FALSE) ~ "shared_species_proxy_review",
      coalesce(match_review_flag, FALSE) ~ "match_review",
      preferred_source_match_status == "preferred_source_matched" ~ "default_clean",
      overall_match_status == "matched_or_candidate" ~ "fallback_source_review",
      TRUE ~ "unmatched"
    )
  ) %>%
  select(
    bridge_source,
    bridge_row_type,
    bridge_to_existing_who,
    analysis_unit_id,
    master_row,
    disease_master_name,
    resolved_disease_name_final,
    resolved_pathogen_name_final,
    resolved_pathogen_rank_final,
    include_status_final,
    active_master_analysis_unit,
    split_group,
    pathogen_aliases,
    pathogen_family_master,
    master_tier,
    master_guild,
    master_livestock_amplified,
    master_key_host_vector,
    in_master_who,
    in_master_gibb,
    in_master_empres_i,
    in_master_atlas,
    master_gbif_checked,
    master_notes,
    name_resolution_status,
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
    notes,
    resolution_source,
    resolution_notes,
    clover_matched_pathogen_names,
    virion_matched_pathogen_names,
    clover_matched_taxids,
    virion_matched_taxids,
    clover_matched_families,
    virion_matched_families,
    clover_matched_source_types,
    virion_matched_source_types,
    clover_best_match_type,
    virion_best_match_type,
    clover_match_status,
    virion_match_status,
    preferred_match_source,
    overall_match_status,
    preferred_source_match_status,
    match_review_flag,
    shared_species_proxy_flag,
    match_review_notes,
    host_query_include_default,
    host_query_bucket,
    host_query_source,
    host_query_pathogen_names,
    host_query_taxids
  ) %>%
  arrange(master_row)

master_keys_for_who_overlap <- bridge %>%
  transmute(
    analysis_unit_key = clean_key(analysis_unit),
    analysis_unit_label_key = clean_key(analysis_unit_label),
    source_disease_key = clean_key(source_disease_name),
    source_pathogen_key = clean_key(source_pathogen)
  )

who_only_units <- who_units %>%
  filter(
    !who_analysis_unit_key %in% master_keys_for_who_overlap$analysis_unit_key,
    !who_analysis_unit_label_key %in% master_keys_for_who_overlap$analysis_unit_label_key,
    !who_source_disease_key %in% master_keys_for_who_overlap$source_disease_key,
    !who_source_pathogen_key %in% master_keys_for_who_overlap$source_pathogen_key
  ) %>%
  transmute(
    bridge_source = "who_pathogen_analysis_units",
    bridge_row_type = "who_only_existing_analysis_unit",
    bridge_to_existing_who = TRUE,
    analysis_unit_id = paste0("who_", who_unit_row),
    master_row = NA_integer_,
    disease_master_name = NA_character_,
    resolved_disease_name_final = coalesce_chr(source_disease_name, analysis_unit_label),
    resolved_pathogen_name_final = coalesce_chr(analysis_unit, source_pathogen),
    resolved_pathogen_rank_final = analysis_unit_rank,
    include_status_final = if_else(analysis_decision == "keep", "yes_existing_who", "review"),
    active_master_analysis_unit = FALSE,
    split_group = NA_character_,
    pathogen_aliases = NA_character_,
    pathogen_family_master = NA_character_,
    master_tier = NA_character_,
    master_guild = NA_character_,
    master_livestock_amplified = NA,
    master_key_host_vector = NA_character_,
    in_master_who = FALSE,
    in_master_gibb = FALSE,
    in_master_empres_i = FALSE,
    in_master_atlas = FALSE,
    master_gbif_checked = FALSE,
    master_notes = NA_character_,
    name_resolution_status = "who_only_existing_unit",
    existing_lookup_name = NA_character_,
    match_field = NA_character_,
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
    notes,
    resolution_source = "existing_who_analysis_unit",
    resolution_notes = NA_character_,
    clover_matched_pathogen_names = NA_character_,
    virion_matched_pathogen_names = NA_character_,
    clover_matched_taxids = NA_character_,
    virion_matched_taxids = NA_character_,
    clover_matched_families = NA_character_,
    virion_matched_families = NA_character_,
    clover_matched_source_types = NA_character_,
    virion_matched_source_types = NA_character_,
    clover_best_match_type = NA_character_,
    virion_best_match_type = NA_character_,
    clover_match_status = NA_character_,
    virion_match_status = NA_character_,
    preferred_match_source = NA_character_,
    overall_match_status = NA_character_,
    preferred_source_match_status = NA_character_,
    match_review_flag = NA,
    shared_species_proxy_flag = NA,
    match_review_notes = NA_character_,
    host_query_include_default = FALSE,
    host_query_bucket = "who_only_needs_master_match",
    host_query_source = NA_character_,
    host_query_pathogen_names = NA_character_,
    host_query_taxids = NA_character_
  )

combined_units <- bind_rows(bridge, who_only_units) %>%
  left_join(transmission_rules, by = "analysis_unit_id") %>%
  arrange(
    bridge_source != "disease_master_list",
    master_row,
    resolved_pathogen_name_final
  )

combined_units_compact <- combined_units %>%
  transmute(
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
    notes,
    vectored_status,
    generalist_status,
    transmission_complexity,
    guild,
    host_sdm_needed,
    vector_sdm_needed,
    host_range_rule,
    vector_range_rule,
    range_limiting_layer,
    transmission_rule_notes,
    transmission_rule_review_status,
    modelling_scope_status,
    modelling_scope_reason,
    analysis_unit_id,
    master_row,
    disease_master_name,
    master_tier,
    master_guild,
    master_livestock_amplified,
    master_key_host_vector,
    include_as_analysis_unit = include_status_final,
    preferred_match_source,
    matched_pathogen_names = host_query_pathogen_names,
    matched_taxids = host_query_taxids,
    match_review_flag = coalesce(match_review_flag, FALSE),
    shared_species_proxy_flag = coalesce(shared_species_proxy_flag, FALSE),
    match_review_notes,
    host_query_bucket
  )

host_query_units <- combined_units %>%
  filter(
    bridge_source == "disease_master_list",
    include_status_final == "yes"
  ) %>%
  transmute(
    analysis_unit_id,
    master_row,
    disease_master_name,
    resolved_disease_name = resolved_disease_name_final,
    resolved_pathogen_name = resolved_pathogen_name_final,
    resolved_pathogen_rank = resolved_pathogen_rank_final,
    preferred_match_source,
    host_query_include_default,
    host_query_bucket,
    host_query_source,
    host_query_pathogen_names,
    host_query_taxids,
    match_review_flag = coalesce(match_review_flag, FALSE),
    shared_species_proxy_flag = coalesce(shared_species_proxy_flag, FALSE),
    match_review_notes,
    split_group,
    pathogen_family_master,
    master_tier,
    master_guild,
    master_livestock_amplified,
    master_key_host_vector
  ) %>%
  arrange(host_query_bucket, master_row)

stopifnot(nrow(bridge) == nrow(master_units))
stopifnot(!anyDuplicated(bridge$analysis_unit_id))
stopifnot(nrow(host_query_units) == sum(bridge$include_status_final == "yes", na.rm = TRUE))
stopifnot(nrow(combined_units_compact) == nrow(combined_units))

write_csv(combined_units_compact, combined_output_path, na = "")
write_csv(host_query_units, host_query_output_path, na = "")

cat("Disease-master rows:", nrow(bridge), "\n")
cat("WHO-only rows appended:", nrow(who_only_units), "\n")
cat("Compact combined master + WHO rows:", nrow(combined_units_compact), "\n")
cat("Host-query rows:", nrow(host_query_units), "\n")
cat("Default clean host-query rows:", sum(host_query_units$host_query_include_default), "\n")
cat("Host-query buckets:\n")
print(count(host_query_units, host_query_bucket), n = Inf)
cat("Transmission rule rows joined:", nrow(transmission_rules), "\n")
cat("Wrote combined output:", combined_output_path, "\n")
cat("Wrote host-query output:", host_query_output_path, "\n")
