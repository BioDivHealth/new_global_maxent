# ------------------------------------------------------------------------------
# 5_11_Host_Vector_Join_QA.R
# ------------------------------------------------------------------------------
# Purpose: Audit host-vector join readiness, integrated row coverage, blocked
#          rows, and taxonomy cautions across the disease-level and
#          pathogen-level host-vector-pathogen outputs.
#
# Inputs : WHO network helper path for master_plus_who_host_network.csv,
#          filtered to legacy canonical zoonotic associations
#          vector_screening_evidence_path(
#            "disease_vector_links_taxonomy_cleaned.csv"
#          )
#          vector_screening_evidence_path("pathogen_vector_links_filled.csv")
#          pathogen_association_data/evidence/host_vector/
#          vector_host_links_join_ready.csv
#          vector_host_links_join_blocked.csv
#          WHO host-vector helper paths for:
#          disease_host_vector_links.csv
#          pathogen_host_vector_links.csv
# Outputs: WHO network QA helper paths:
#          host_vector_join_qa_summary.csv
#          host_vector_join_missing_host_tax_id.csv
#          host_vector_join_unmatched_disease_vectors.csv
#          host_vector_join_unmatched_pathogen_vectors.csv
#          host_vector_join_taxonomy_caution_rows.csv
#          host_vector_join_disease_coverage.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))
source(here("scripts", "associations", "host_vector_integration", "host_vector_join_helpers.R"))

host_vector_dir <- vector_host_outputs_dir

who_path <- who_network_host_pathogen_path("master_plus_who_host_network.csv")
disease_vector_path <- vector_screening_evidence_path(
  "disease_vector_links_taxonomy_cleaned.csv"
)
pathogen_vector_path <- vector_screening_evidence_path("pathogen_vector_links_filled.csv")
host_vector_join_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
host_vector_blocked_path <- file.path(host_vector_dir, "vector_host_links_join_blocked.csv")
disease_output_path <- who_network_host_vector_path("disease_host_vector_links.csv")
pathogen_output_path <- who_network_host_vector_path("pathogen_host_vector_links.csv")

summary_path <- who_network_qa_path("host_vector_join_qa_summary.csv")
missing_taxid_path <- who_network_qa_path("host_vector_join_missing_host_tax_id.csv")
unmatched_disease_path <- who_network_qa_path("host_vector_join_unmatched_disease_vectors.csv")
unmatched_pathogen_path <- who_network_qa_path("host_vector_join_unmatched_pathogen_vectors.csv")
taxonomy_caution_path <- who_network_qa_path("host_vector_join_taxonomy_caution_rows.csv")
disease_coverage_path <- who_network_qa_path("host_vector_join_disease_coverage.csv")

who_network <- read_clean_csv(who_path) %>%
  filter_legacy_compatible_host_network()
disease_vectors <- read_clean_csv(disease_vector_path)
pathogen_vectors <- read_clean_csv(pathogen_vector_path)
host_vector_join <- read_clean_csv(host_vector_join_path)
host_vector_blocked <- read_clean_csv(host_vector_blocked_path)
disease_output <- read_clean_csv(
  disease_output_path,
  col_types = cols(
    review_reason_examples = col_character()
  )
)
pathogen_output <- read_clean_csv(
  pathogen_output_path,
  col_types = cols(
    review_reason_examples = col_character()
  )
)

host_vector_keys <- prepare_host_vector_keys(host_vector_join)

missing_host_tax_id <- host_vector_blocked %>%
  filter(stringr::str_detect(coalesce(block_reason, ""), "missing_host_tax_id"))

unmatched_disease_vectors <- disease_vectors %>%
  mutate(vector_join_key = normalize_vector_key(vector_species_taxonomy_cleaned)) %>%
  filter(!is.na(vector_join_key)) %>%
  anti_join(host_vector_keys, by = "vector_join_key") %>%
  arrange(disease_name, vector_species_taxonomy_cleaned)

unmatched_pathogen_vectors <- pathogen_vectors %>%
  mutate(vector_join_key = normalize_vector_key(candidate_vector_species)) %>%
  filter(!is.na(vector_join_key), assignment_basis != "no_disease_vector_match") %>%
  anti_join(host_vector_keys, by = "vector_join_key") %>%
  arrange(disease_name, pathogen, candidate_vector_species)

taxonomy_caution_rows <- bind_rows(
  disease_output %>% mutate(output_level = "disease"),
  pathogen_output %>% mutate(output_level = "pathogen")
) %>%
  filter(taxonomy_caution %in% TRUE)

disease_vector_coverage <- disease_vectors %>%
  mutate(vector_join_key = normalize_vector_key(vector_species_taxonomy_cleaned)) %>%
  group_by(disease_name) %>%
  summarise(
    total_disease_vector_rows = n(),
    total_distinct_disease_vectors = n_distinct(vector_join_key),
    disease_vectors_with_host_overlap = n_distinct(vector_join_key[vector_join_key %in% host_vector_keys$vector_join_key]),
    .groups = "drop"
  ) %>%
  left_join(
    disease_output %>%
      group_by(disease_name) %>%
      summarise(
        final_disease_host_vector_rows = n(),
        final_distinct_hosts = n_distinct(host_tax_id),
        final_distinct_vectors = n_distinct(vector_join_key),
        taxonomy_caution_rows = sum(taxonomy_caution %in% TRUE),
        .groups = "drop"
      ),
    by = "disease_name"
  ) %>%
  mutate(across(starts_with("final_"), ~ coalesce(.x, 0L))) %>%
  mutate(taxonomy_caution_rows = coalesce(taxonomy_caution_rows, 0L)) %>%
  arrange(desc(final_disease_host_vector_rows), disease_name)

qa_summary <- tibble::tribble(
  ~metric, ~value,
  "who_network_rows", as.character(nrow(who_network)),
  "who_distinct_diseases", as.character(n_distinct(who_network$Disease_name)),
  "who_distinct_pathogens", as.character(n_distinct(who_network$PathogenTaxID)),
  "who_distinct_hosts", as.character(n_distinct(who_network$HostTaxID)),
  "host_vector_join_ready_rows", as.character(nrow(host_vector_join)),
  "host_vector_join_blocked_rows", as.character(nrow(host_vector_blocked)),
  "host_vector_blocked_missing_host_tax_id_rows", as.character(nrow(missing_host_tax_id)),
  "host_vector_join_distinct_hosts", as.character(n_distinct(host_vector_join$host_tax_id)),
  "host_vector_join_distinct_vectors", as.character(n_distinct(host_vector_join$vector_join_key)),
  "disease_vector_rows", as.character(nrow(disease_vectors)),
  "disease_vector_unmatched_rows", as.character(nrow(unmatched_disease_vectors)),
  "pathogen_vector_rows", as.character(nrow(pathogen_vectors)),
  "pathogen_vector_unmatched_rows", as.character(nrow(unmatched_pathogen_vectors)),
  "disease_host_vector_rows", as.character(nrow(disease_output)),
  "pathogen_host_vector_rows", as.character(nrow(pathogen_output)),
  "diseases_with_disease_level_matches", as.character(n_distinct(disease_output$disease_name)),
  "diseases_with_pathogen_level_matches", as.character(n_distinct(pathogen_output$disease_name)),
  "taxonomy_caution_rows_total", as.character(nrow(taxonomy_caution_rows))
)

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
write_csv(qa_summary, summary_path, na = "")
write_csv(missing_host_tax_id, missing_taxid_path, na = "")
write_csv(unmatched_disease_vectors, unmatched_disease_path, na = "")
write_csv(unmatched_pathogen_vectors, unmatched_pathogen_path, na = "")
write_csv(taxonomy_caution_rows, taxonomy_caution_path, na = "")
write_csv(disease_vector_coverage, disease_coverage_path, na = "")

cat("QA summary rows written:", nrow(qa_summary), "\n")
cat("Missing host-taxid blocked rows:", nrow(missing_host_tax_id), "\n")
cat("Unmatched disease-vector rows:", nrow(unmatched_disease_vectors), "\n")
cat("Unmatched pathogen-vector rows:", nrow(unmatched_pathogen_vectors), "\n")
cat("Taxonomy caution rows:", nrow(taxonomy_caution_rows), "\n")
cat("Disease coverage rows:", nrow(disease_vector_coverage), "\n")
cat("Wrote QA summary to", summary_path, "\n")
