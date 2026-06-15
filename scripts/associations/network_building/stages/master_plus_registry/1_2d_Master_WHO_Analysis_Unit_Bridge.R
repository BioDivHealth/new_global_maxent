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
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_registry_helpers.R"
))

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
  mutate(across(where(is.character), registry_clean_text))

manual_units <- read_csv(manual_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    across(where(is.character), registry_clean_text),
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
    across(where(is.character), registry_clean_text),
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
    across(where(is.character), registry_clean_text),
    who_unit_row = row_number(),
    who_analysis_unit_key = registry_clean_key(analysis_unit),
    who_analysis_unit_label_key = registry_clean_key(analysis_unit_label),
    who_source_disease_key = registry_clean_key(source_disease_name),
    who_source_pathogen_key = registry_clean_key(source_pathogen)
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
    mutate(across(where(is.character), registry_clean_text))

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

bridge <- registry_build_master_bridge(master_units, manual_units, master_matches)
who_only_units <- registry_build_who_only_units(who_units, bridge)

combined_units <- bind_rows(bridge, who_only_units) %>%
  left_join(transmission_rules, by = "analysis_unit_id") %>%
  arrange(
    bridge_source != "disease_master_list",
    master_row,
    resolved_pathogen_name_final
  )

combined_units_compact <- registry_compact_analysis_units(combined_units)
host_query_units <- registry_build_host_query_units(combined_units)

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
