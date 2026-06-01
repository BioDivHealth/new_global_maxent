# ------------------------------------------------------------------------------
# 5_9b_Derive_Disease_Host_Vector_Links_Expanded.R
# ------------------------------------------------------------------------------
# Purpose: Build an expanded disease-host-vector table for the screened disease
#          subset by keeping all observed host-vector combinations for WHO hosts
#          and then marking whether curated disease-vector evidence is present.
#
# Inputs : WHO network helper path for combined_who_network_canonical_zoonotic.csv
#          vector_screening_evidence_path(
#            "disease_vector_links_taxonomy_cleaned.csv"
#          )
#          pathogen_association_data/evidence/host_vector/
#          vector_host_links_join_ready.csv
# Output : WHO host-vector helper paths for:
#          disease_host_vector_links_expanded.csv
#          disease_host_vector_links_expanded_summary.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

host_vector_dir <- vector_host_outputs_dir

who_path <- who_working_network_path()
disease_vector_path <- vector_screening_evidence_path(
  "disease_vector_links_taxonomy_cleaned.csv"
)
host_vector_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
expanded_path <- who_network_host_vector_path("disease_host_vector_links_expanded.csv")
summary_path <- who_network_host_vector_path("disease_host_vector_links_expanded_summary.csv")

who_network <- read_csv(
  who_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

disease_vectors <- read_csv(
  disease_vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

host_vectors <- read_csv(
  host_vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    host_tax_id = clean_text(host_tax_id),
    vector_join_key = normalize_vector_key(vector_join_key)
  )

screened_diseases <- disease_vectors %>%
  filter(!is.na(disease_name)) %>%
  distinct(disease_name) %>%
  pull(disease_name)

disease_host_network <- who_network %>%
  filter(
    !is.na(Disease_name),
    Disease_name %in% screened_diseases,
    !is.na(HostTaxID),
    !is.na(Host)
  ) %>%
  mutate(
    disease_name_join = normalize_name_for_match(Disease_name),
    host_tax_id = clean_text(HostTaxID)
  ) %>%
  group_by(disease_name_join, Disease_name, host_tax_id) %>%
  summarise(
    host = first_non_missing(Host),
    host_class = first_non_missing(HostClass),
    host_order = first_non_missing(HostOrder),
    host_family = first_non_missing(HostFamily),
    pathogen_count_in_disease_host_network = n_distinct(PathogenTaxID),
    pathogen_examples = collapse_unique(Pathogen),
    detection_method_examples = collapse_unique(DetectionMethod),
    main_source_examples = collapse_unique(MainSource),
    .groups = "drop"
  )

disease_vector_joinable <- disease_vectors %>%
  filter(!is.na(disease_name), !is.na(vector_species_taxonomy_cleaned)) %>%
  mutate(
    disease_name_join = normalize_name_for_match(disease_name),
    vector_join_key = normalize_vector_key(vector_species_taxonomy_cleaned)
  ) %>%
  group_by(disease_name_join, vector_join_key) %>%
  summarise(
    disease_name_clean = first_non_missing(disease_name_clean),
    vector_species = first_non_missing(vector_species_taxonomy_cleaned),
    vector_group = collapse_unique(vector_group),
    best_evidence_level = first_non_missing(best_evidence_level),
    best_evidence_basis = first_non_missing(best_evidence_basis),
    record_sources = collapse_unique(record_sources),
    supporting_row_count = sum(suppressWarnings(as.integer(supporting_row_count)), na.rm = TRUE),
    disease_vector_taxon_rank = first_non_missing(vector_taxon_rank),
    disease_vector_review_needed = any(review_needed %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

host_vector_joinable <- host_vectors %>%
  filter(!is.na(host_tax_id), !is.na(vector_join_key)) %>%
  rename(
    hv_host = host,
    hv_host_class = host_class,
    hv_host_order = host_order,
    hv_host_family = host_family,
    hv_vector_species = vector_species,
    hv_vector_taxon_rank = vector_taxon_rank,
    hv_vector_species_needs_review = vector_species_needs_review,
    hv_vector_name_taxonomy_examples = vector_name_taxonomy_examples,
    hv_source_platform_examples = source_platform_examples,
    hv_source_dataset_examples = source_dataset_examples,
    hv_interaction_type_examples = interaction_type_examples,
    hv_country_examples = country_examples,
    hv_review_reason_examples = review_reason_examples,
    hv_record_count = record_count
  )

disease_host_vector_links_expanded <- disease_host_network %>%
  inner_join(
    host_vector_joinable,
    by = "host_tax_id",
    relationship = "many-to-many"
  ) %>%
  left_join(
    disease_vector_joinable,
    by = c("disease_name_join", "vector_join_key"),
    relationship = "many-to-many"
  ) %>%
  transmute(
    disease_name = Disease_name,
    disease_name_clean = coalesce(disease_name_clean, normalize_name_for_match(Disease_name)),
    host = host,
    host_tax_id,
    host_class,
    host_order,
    host_family,
    pathogen_count_in_disease_host_network,
    pathogen_examples,
    detection_method_examples,
    main_source_examples,
    vector_species = coalesce(vector_species, hv_vector_species),
    vector_group,
    best_evidence_level,
    best_evidence_basis,
    record_sources,
    supporting_row_count = if_else(is.na(supporting_row_count), NA_real_, as.numeric(supporting_row_count)),
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
    host_vector_evidence = TRUE,
    disease_vector_evidence = !is.na(best_evidence_level),
    disease_vector_evidence_status = if_else(
      !is.na(best_evidence_level),
      "supported_in_disease_vector_table",
      "not_supported_in_disease_vector_table"
    ),
    link_status = if_else(
      !is.na(best_evidence_level),
      "confirmed_by_both",
      "host_vector_only_candidate"
    ),
    screening_scope = "screened_disease_subset",
    link_type = "disease_host_vector_expanded",
    vector_join_match_type = "exact_normalized_vector_name",
    taxonomy_caution = dplyr::if_else(
      coalesce(hv_vector_species_needs_review, FALSE) |
        coalesce(disease_vector_review_needed, FALSE) |
        coalesce(hv_vector_taxon_rank, "species") != "species" |
        (!is.na(disease_vector_taxon_rank) & coalesce(disease_vector_taxon_rank, "species") != "species"),
      TRUE,
      FALSE
    )
  ) %>%
  arrange(disease_name, host, vector_species)

duplicate_key_count <- disease_host_vector_links_expanded %>%
  count(disease_name, host_tax_id, vector_join_key, name = "n") %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate disease_name + host_tax_id + vector_join_key rows found in expanded disease-host-vector output")
}

expanded_summary <- disease_host_vector_links_expanded %>%
  group_by(disease_name) %>%
  summarise(
    expanded_rows = n(),
    confirmed_by_both_rows = sum(link_status == "confirmed_by_both"),
    host_vector_only_candidate_rows = sum(link_status == "host_vector_only_candidate"),
    distinct_hosts = n_distinct(host_tax_id),
    distinct_vectors = n_distinct(vector_join_key),
    taxonomy_caution_rows = sum(taxonomy_caution %in% TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(expanded_rows), disease_name)

dir.create(dirname(expanded_path), recursive = TRUE, showWarnings = FALSE)
write_csv(disease_host_vector_links_expanded, expanded_path, na = "")
write_csv(expanded_summary, summary_path, na = "")

cat("WHO disease-host rows used:", nrow(disease_host_network), "\n")
cat("Disease-vector rows used:", nrow(disease_vector_joinable), "\n")
cat("Host-vector join rows used:", nrow(host_vector_joinable), "\n")
cat("Expanded disease-host-vector rows written:", nrow(disease_host_vector_links_expanded), "\n")
cat("Distinct diseases in output:", n_distinct(disease_host_vector_links_expanded$disease_name), "\n")
cat("Confirmed-by-both rows:", sum(disease_host_vector_links_expanded$link_status == "confirmed_by_both"), "\n")
cat("Host-vector-only candidate rows:", sum(disease_host_vector_links_expanded$link_status == "host_vector_only_candidate"), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Wrote expanded disease-host-vector table to", expanded_path, "\n")
cat("Wrote expanded summary table to", summary_path, "\n")
