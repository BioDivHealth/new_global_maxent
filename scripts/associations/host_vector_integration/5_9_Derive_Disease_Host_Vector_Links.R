# ------------------------------------------------------------------------------
# 5_9_Derive_Disease_Host_Vector_Links.R
# ------------------------------------------------------------------------------
# Purpose: Derive one row per disease-host-vector by joining the WHO
#          disease-host network to the canonical disease-vector table and the
#          observational host-vector join table.
#
# Inputs : WHO network helper path for master_plus_who_host_network.csv,
#          filtered to legacy canonical zoonotic associations
#          vector_screening_evidence_path(
#            "disease_vector_links_taxonomy_cleaned.csv"
#          )
#          pathogen_association_data/evidence/host_vector/
#          vector_host_links_join_ready.csv
# Output : WHO host-vector helper path for disease_host_vector_links.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))
source(here("scripts", "associations", "host_vector_integration", "host_vector_join_helpers.R"))

host_vector_dir <- vector_host_outputs_dir

who_path <- who_network_host_pathogen_path("master_plus_who_host_network.csv")
disease_vector_path <- vector_screening_evidence_path(
  "disease_vector_links_taxonomy_cleaned.csv"
)
host_vector_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
output_path <- who_network_host_vector_path("disease_host_vector_links.csv")

who_network <- read_clean_csv(who_path) %>%
  filter_legacy_compatible_host_network()
disease_vectors <- read_clean_csv(disease_vector_path)
host_vectors <- read_clean_csv(host_vector_path)

disease_host_network <- prepare_disease_host_network(who_network)
disease_vector_joinable <- prepare_disease_vector_joinable(disease_vectors)
host_vector_joinable <- prepare_host_vector_joinable(host_vectors)

disease_host_vector_links <- disease_host_network %>%
  inner_join(
    disease_vector_joinable,
    by = "disease_name_join",
    relationship = "many-to-many"
  ) %>%
  inner_join(
    host_vector_joinable,
    by = c("host_tax_id", "vector_join_key"),
    relationship = "many-to-many"
  ) %>%
  transmute(
    disease_name = Disease_name,
    disease_name_clean,
    host = host,
    host_tax_id,
    host_class,
    host_order,
    host_family,
    pathogen_count_in_disease_host_network,
    pathogen_examples,
    detection_method_examples,
    main_source_examples,
    vector_species,
    vector_group,
    best_evidence_level,
    best_evidence_basis,
    record_sources,
    supporting_row_count,
    disease_vector_taxon_rank,
    disease_vector_review_needed,
    host_vector_species = hv_vector_species,
    vector_taxon_rank = hv_vector_taxon_rank,
    vector_species_needs_review = hv_vector_species_needs_review,
    vector_host_record_count = hv_record_count,
    vector_name_taxonomy_examples = hv_vector_name_taxonomy_examples,
    source_platform_examples = hv_source_platform_examples,
    source_dataset_examples = hv_source_dataset_examples,
    interaction_type_examples = hv_interaction_type_examples,
    country_examples = hv_country_examples,
    review_reason_examples = hv_review_reason_examples,
    vector_join_key,
    link_type = "disease_host_vector_derived",
    vector_join_match_type = "exact_normalized_vector_name",
    taxonomy_caution = dplyr::if_else(
      coalesce(hv_vector_species_needs_review, FALSE) |
        coalesce(disease_vector_review_needed, FALSE) |
        coalesce(hv_vector_taxon_rank, "species") != "species" |
        coalesce(disease_vector_taxon_rank, "species") != "species",
      TRUE,
      FALSE
    )
  ) %>%
  arrange(disease_name, host, vector_species)

duplicate_key_count <- disease_host_vector_links %>%
  count(disease_name, host_tax_id, vector_join_key, name = "n") %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate disease_name + host_tax_id + vector_join_key rows found in disease-host-vector output")
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(disease_host_vector_links, output_path, na = "")

cat("WHO disease-host rows used:", nrow(disease_host_network), "\n")
cat("Disease-vector rows used:", nrow(disease_vector_joinable), "\n")
cat("Host-vector join rows used:", nrow(host_vector_joinable), "\n")
cat("Disease-host-vector rows written:", nrow(disease_host_vector_links), "\n")
cat("Distinct diseases in output:", n_distinct(disease_host_vector_links$disease_name), "\n")
cat("Taxonomy caution rows:", sum(disease_host_vector_links$taxonomy_caution %in% TRUE), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Wrote disease-host-vector links to", output_path, "\n")
