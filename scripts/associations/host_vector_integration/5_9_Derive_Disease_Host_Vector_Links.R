# ------------------------------------------------------------------------------
# 5_9_Derive_Disease_Host_Vector_Links.R
# ------------------------------------------------------------------------------
# Purpose: Derive one row per disease-host-vector by joining the WHO
#          disease-host network to the canonical disease-vector table and the
#          observational host-vector join table.
#
# Inputs : pathogen_association_data/WHO/networks/
#          combined_who_network_canonical_zoonotic.csv
#          pathogen_association_data/WHO/vector_screening/
#          disease_vector_links_taxonomy_cleaned.csv
#          pathogen_association_data/vector_host/outputs/
#          vector_host_links_join_ready.csv
# Output : pathogen_association_data/WHO/networks/
#          disease_host_vector_links.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

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

normalize_name_for_match <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "haemorrh", "hemorrh")
  x <- stringr::str_replace_all(x, "&", " and ")
  x <- stringr::str_replace_all(x, "[/]", " ")
  x <- stringr::str_replace_all(x, "[-–—]", " ")
  x <- stringr::str_replace_all(x, "[()\\[\\],.;:*'`\"]", " ")
  x <- stringr::str_replace_all(x, "\\bviruses\\b", "virus")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

normalize_vector_key <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_to_lower(x)
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

first_non_missing <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

networks_dir <- here("pathogen_association_data", "WHO", "networks")
vector_dir <- here("pathogen_association_data", "WHO", "vector_screening")
host_vector_dir <- here("pathogen_association_data", "vector_host", "outputs")
vector_output_dir <- file.path(vector_dir, "outputs")

who_path <- who_working_network_path()
disease_vector_path <- file.path(vector_output_dir, "disease_vector_links_taxonomy_cleaned.csv")
host_vector_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
output_path <- file.path(networks_dir, "disease_host_vector_links.csv")

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

disease_host_network <- who_network %>%
  filter(!is.na(Disease_name), !is.na(HostTaxID), !is.na(Host)) %>%
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

write_csv(disease_host_vector_links, output_path, na = "")

cat("WHO disease-host rows used:", nrow(disease_host_network), "\n")
cat("Disease-vector rows used:", nrow(disease_vector_joinable), "\n")
cat("Host-vector join rows used:", nrow(host_vector_joinable), "\n")
cat("Disease-host-vector rows written:", nrow(disease_host_vector_links), "\n")
cat("Distinct diseases in output:", n_distinct(disease_host_vector_links$disease_name), "\n")
cat("Taxonomy caution rows:", sum(disease_host_vector_links$taxonomy_caution %in% TRUE), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Wrote disease-host-vector links to", output_path, "\n")
