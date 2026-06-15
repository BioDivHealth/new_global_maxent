# ------------------------------------------------------------------------------
# 1_2c_Master_Virion_Clover_Matches.R
# ------------------------------------------------------------------------------
# Purpose: Match active, concrete disease-master analysis units to local VIRION
#          and CLOVER pathogen taxonomies.
#
# Inputs : who_diseases_name_resolution_path(
#            "master_disease_name_resolution_manual.csv"
#          )
#          local VIRION and CLOVER source tables
#
# Outputs: who_diseases_staged_master_expansion_path(
#            "master_pathogen_virion_clover_candidates.csv"
#          )
#          who_diseases_staged_master_expansion_path(
#            "master_pathogen_virion_clover_matches.csv"
#          )
# ------------------------------------------------------------------------------

# ------------------------------| Load libraries |------------------------------
library(tidyverse)
library(here)
library(stringdist)
library(magrittr)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_registry_helpers.R"
))

# ------------------------------| Helper paths |-------------------------------
manual_path <- who_diseases_name_resolution_path(
  "master_disease_name_resolution_manual.csv"
)
alias_path <- who_diseases_pathogen_matching_manual_path(
  "master_pathogen_aliases.csv"
)
candidate_output_path <- who_diseases_staged_master_expansion_path(
  "master_pathogen_virion_clover_candidates.csv"
)
match_output_path <- who_diseases_staged_master_expansion_path(
  "master_pathogen_virion_clover_matches.csv"
)
external_review_output_path <- who_diseases_staged_pathogen_matching_path(
  "master_pathogen_external_taxonomy_review.csv"
)

clover_dir <- file.path(
  clover_source_dir,
  "clover", "clover_1.0_allpathogens"
)

clover_paths <- file.path(
  clover_dir,
  c(
    "CLOVER_1.0_Bacteria_AssociationsFlatFile.csv",
    "CLOVER_1.0_Viruses_AssociationsFlatFile.csv",
    "CLOVER_1.0_HelminthProtozoaFungi_AssociationsFlatFile.csv"
  )
)

# ------------------------------| Helpers |------------------------------------
rank_in_scope <- c("species", "species_complex", "subspecies")
include_states_in_scope <- "yes"

external_taxonomy_review <- tribble(
  ~resolved_pathogen_name, ~external_source, ~external_taxid, ~external_taxon_name, ~external_rank, ~external_parent_taxon, ~external_source_url, ~external_review_notes,
  "Alkhumra hemorrhagic fever virus", "NCBI Taxonomy", "172148", "Alkhumra hemorrhagic fever virus", "no rank", "Orthoflavivirus kyasanurense / Kyasanur Forest disease virus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=172148&mode=Info", "Found in NCBI but not in the local VIRION/CLOVER tables used here. Also standardizes spelling from Alkhurma to Alkhumra.",
  "Rocio virus", "NCBI Taxonomy", "64315", "Rocio virus", "no rank", "Orthoflavivirus ilheusense / Ilheus virus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=64315&mode=Info", "Found in NCBI but not in the local VIRION/CLOVER tables used here; NCBI places it under Orthoflavivirus ilheusense."
) %>%
  mutate(query_key = registry_normalize_name(resolved_pathogen_name))

manual_aliases <- read_csv(alias_path, show_col_types = FALSE, na = c("", "NA")) %>%
  rename(
    source = alias_source,
    source_name = alias_name
  ) %>%
  filter(!is.na(resolved_pathogen_name), !is.na(source), !is.na(source_name)) %>%
  mutate(
    source = str_to_lower(source),
    query_key = registry_normalize_name(resolved_pathogen_name),
    source_key = registry_normalize_name(source_name),
    alias_review_flag = alias_type %in% c(
      "shared_species_proxy",
      "species_proxy",
      "narrow_local_match",
      "parent_species_match"
    )
  )

unexpected_alias_sources <- setdiff(unique(manual_aliases$source), c("virion", "clover"))
if (length(unexpected_alias_sources) > 0) {
  stop(
    "Unexpected alias_source values in ",
    alias_path,
    ": ",
    paste(unexpected_alias_sources, collapse = ", ")
  )
}

# ------------------------------| Load query rows |----------------------------
manual_units <- read_csv(manual_path, show_col_types = FALSE, na = c("", "NA")) %>%
  filter(
    resolved_pathogen_rank %in% rank_in_scope,
    include_as_analysis_unit %in% include_states_in_scope
  ) %>%
  mutate(
    analysis_unit_id = paste0("master_", master_row),
    query_key = registry_normalize_name(resolved_pathogen_name),
    preferred_match_source = if_else(
      str_detect(str_to_lower(pathogen_family_master), "viridae$|virus|lyssa|hanta|arena|flavi|toga|paramyxo|peribunya|reo|pox"),
      "virion",
      "clover"
    )
  )

stopifnot(nrow(manual_units) > 0)

# ------------------------------| Load VIRION taxonomy |-----------------------
if (!exists("virion_data")) {
  source(file.path("scripts", "associations", "network_building", "helpers", "virion_loaders.R"))
  virion_data <- load_virion_data()
}

virion_taxonomy <- virion_data$taxonomy_virus %>%
  transmute(
    source = "virion",
    source_pathogen_name = Virus,
    source_taxid = VirusTaxID,
    source_family = VirusFamily,
    source_type = "virus"
  ) %>%
  distinct()

# ------------------------------| Load CLOVER taxonomy |-----------------------
missing_clover_files <- clover_paths[!file.exists(clover_paths)]
if (length(missing_clover_files) > 0) {
  stop(
    "Missing CLOVER input files: ",
    paste(missing_clover_files, collapse = "; ")
  )
}

clover_taxonomy <- map_dfr(
  clover_paths,
  ~ read_csv(.x, show_col_types = FALSE, na = c("", "NA"))
) %>%
  filter(!is.na(Pathogen), Pathogen != "") %>%
  transmute(
    source = "clover",
    source_pathogen_name = Pathogen,
    source_taxid = PathogenTaxID,
    source_family = PathogenFamily,
    source_type = PathogenType
  ) %>%
  distinct()

# ------------------------------| Match |--------------------------------------
source_taxonomies <- list(
  virion = virion_taxonomy,
  clover = clover_taxonomy
)

all_candidates <- imap_dfr(
  source_taxonomies,
  ~ registry_make_source_matches(
    manual_units,
    .x,
    source_name = .y,
    manual_aliases = manual_aliases,
    max_dist = 0.08
  )
) %>%
  mutate(
    match_status = case_when(
      match_type %in% c("exact", "manual_alias") ~ "accepted_candidate",
      match_type == "fuzzy_candidate" & match_distance <= 0.03 ~ "strong_review_candidate",
      match_type == "fuzzy_candidate" ~ "review_candidate",
      TRUE ~ "review_candidate"
    )
  ) %>%
  arrange(master_row, source, match_status, match_distance, source_pathogen_name)

best_matches <- all_candidates %>%
  group_by(analysis_unit_id, source) %>%
  arrange(
    match(match_type, c("exact", "manual_alias", "fuzzy_candidate")),
    match_distance,
    source_pathogen_name,
    .by_group = TRUE
  ) %>%
  summarise(
    master_row = first(master_row),
    disease_master_name = first(disease_master_name),
    resolved_disease_name = first(resolved_disease_name),
    resolved_pathogen_name = first(resolved_pathogen_name),
    resolved_pathogen_rank = first(resolved_pathogen_rank),
    include_as_analysis_unit = first(include_as_analysis_unit),
    split_group = first(split_group),
    matched_pathogen_names = registry_collapse_unique(source_pathogen_name),
    matched_taxids = registry_collapse_unique(source_taxid),
    matched_families = registry_collapse_unique(source_family),
    matched_source_types = registry_collapse_unique(source_type),
    alias_types = registry_collapse_unique(alias_type),
    alias_notes = registry_collapse_unique(alias_notes),
    alias_review_flag = any(alias_review_flag, na.rm = TRUE),
    best_match_type = first(match_type),
    best_match_distance = first(match_distance),
    match_status = first(match_status),
    candidate_count = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(
      matched_pathogen_names, matched_taxids, matched_families,
      matched_source_types, best_match_type, best_match_distance,
      match_status, candidate_count, alias_types, alias_notes,
      alias_review_flag
    ),
    names_glue = "{source}_{.value}"
  ) %>%
  right_join(
    manual_units %>%
      select(
        analysis_unit_id, master_row, disease_master_name, resolved_disease_name,
        resolved_pathogen_name, resolved_pathogen_rank, include_as_analysis_unit,
        split_group, preferred_match_source
      ),
    by = c(
      "analysis_unit_id", "master_row", "disease_master_name",
      "resolved_disease_name", "resolved_pathogen_name",
      "resolved_pathogen_rank", "include_as_analysis_unit", "split_group"
    )
  ) %>%
  mutate(
    overall_match_status = case_when(
      !is.na(virion_matched_taxids) | !is.na(clover_matched_taxids) ~ "matched_or_candidate",
      TRUE ~ "unmatched"
    ),
    preferred_source_match_status = case_when(
      preferred_match_source == "virion" & !is.na(virion_matched_taxids) ~ "preferred_source_matched",
      preferred_match_source == "clover" & !is.na(clover_matched_taxids) ~ "preferred_source_matched",
      overall_match_status == "matched_or_candidate" ~ "fallback_source_matched",
      TRUE ~ "unmatched"
    ),
    match_review_flag = coalesce(virion_alias_review_flag, FALSE) | coalesce(clover_alias_review_flag, FALSE),
    shared_species_proxy_flag = str_detect(
      coalesce(paste(virion_alias_types, clover_alias_types, sep = "; "), ""),
      "shared_species_proxy"
    )
  ) %>%
  rowwise() %>%
  mutate(
    match_review_notes = registry_collapse_unique(c(virion_alias_notes, clover_alias_notes))
  ) %>%
  ungroup() %>%
  arrange(master_row)

# ------------------------------| Save outputs |-------------------------------
write_csv(all_candidates, candidate_output_path, na = "")
write_csv(best_matches, match_output_path, na = "")

external_review <- best_matches %>%
  filter(overall_match_status == "unmatched") %>%
  mutate(query_key = registry_normalize_name(resolved_pathogen_name)) %>%
  left_join(
    external_taxonomy_review %>% select(-resolved_pathogen_name),
    by = "query_key"
  ) %>%
  select(
    master_row, disease_master_name, resolved_disease_name, resolved_pathogen_name,
    resolved_pathogen_rank, include_as_analysis_unit, split_group,
    external_source, external_taxid, external_taxon_name, external_rank,
    external_parent_taxon, external_source_url, external_review_notes
  )

write_csv(external_review, external_review_output_path, na = "")

# ------------------------------| Console summary |----------------------------
cat("Master units in scope:", nrow(manual_units), "\n")
cat("Candidate rows written:", nrow(all_candidates), "\n")
cat("Units with VIRION candidates:", sum(!is.na(best_matches$virion_matched_taxids)), "\n")
cat("Units with CLOVER candidates:", sum(!is.na(best_matches$clover_matched_taxids)), "\n")
cat("Units still unmatched:", sum(best_matches$overall_match_status == "unmatched"), "\n")
cat("Candidate output:", candidate_output_path, "\n")
cat("Best-match output:", match_output_path, "\n")
cat("External review output:", external_review_output_path, "\n")
