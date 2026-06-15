# ------------------------------------------------------------------------------
# 1_3_WHO_Broad_Taxa_Candidate_Strains.R
# ------------------------------------------------------------------------------
# Purpose: Build a curation inventory of ICTV-supported candidate strains and
#          exemplar viruses for the broad taxa currently under review.
#
# Input  : who_diseases_broad_taxa_manual_path(
#            "who_broad_taxa_candidate_strains_seed.csv"
#          ) plus the current zoonotic WHO analysis-unit shortlist.
# Output : who_diseases_broad_taxa_staged_path(
#            "who_broad_taxa_candidate_strains.csv"
#          )
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "broad_taxa_support_helpers.R"
))

candidate_seed_path <- who_diseases_broad_taxa_manual_path(
  "who_broad_taxa_candidate_strains_seed.csv"
)
output_path <- who_diseases_broad_taxa_staged_path(
  "who_broad_taxa_candidate_strains.csv"
)
analysis_units_keep_path <- who_pathogen_analysis_units_keep_path()

if (!file.exists(candidate_seed_path)) {
  stop("Candidate strain seed table not found: ", candidate_seed_path)
}

analysis_units_keep <- read_csv(
  analysis_units_keep_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), broad_taxa_clean_text))

candidate_strains <- read_csv(
  candidate_seed_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), broad_taxa_clean_text))

candidate_strains <- candidate_strains %>%
  mutate(
    source_rationale = decision_reason,
    source_notes = notes
  ) %>%
  select(-notes) %>%
  left_join(
    analysis_units_keep %>%
      transmute(
        analysis_unit,
        keep_reason,
        keep_rationale = rationale,
        keep_notes = notes,
        do_not_split_further_yet
      ),
    by = c("proposed_active_unit" = "analysis_unit")
  ) %>%
  mutate(
    keep_reason = dplyr::coalesce(keep_reason, decision_reason),
    keep_rationale = dplyr::coalesce(keep_rationale, source_rationale),
    keep_notes = dplyr::coalesce(keep_notes, source_notes),
    support_status = dplyr::case_when(
      decision == "keep_active" ~ "active_unit",
      decision == "keep_supporting_example" ~ "supporting_example",
      TRUE ~ "review_only"
    ),
    ncbi_followup_hint = dplyr::case_when(
      decision == "keep_active" ~ paste0(virus_name, " | ", isolate, " | ", accession),
      TRUE ~ paste0(virus_name, " | ", isolate)
    )
  ) %>%
  select(
    broad_group,
    ictv_genus,
    ictv_subgenus,
    ictv_species,
    virus_name,
    isolate,
    accession,
    available_sequence,
    abbrev,
    candidate_role,
    proposed_active_unit,
    decision,
    decision_reason,
    support_status,
    ncbi_followup_hint,
    ictv_taxon_anchor,
    ictv_taxon_source,
    keep_reason,
    keep_rationale,
    keep_notes,
    source_rationale,
    source_notes
  ) %>%
  arrange(broad_group, desc(decision == "keep_active"), ictv_species, virus_name, isolate)

write_csv(candidate_strains, output_path, na = "")

cat("Candidate strain rows written:", nrow(candidate_strains), "\n")
cat("Active-unit rows:", sum(candidate_strains$decision == "keep_active"), "\n")
cat("Supporting-example rows:", sum(candidate_strains$decision == "keep_supporting_example"), "\n")
cat("Review-only rows:", sum(candidate_strains$decision == "review_only"), "\n")
cat("Wrote candidate strain table to", output_path, "\n")
