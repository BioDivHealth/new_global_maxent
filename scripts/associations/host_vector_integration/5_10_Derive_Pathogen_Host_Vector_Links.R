# ------------------------------------------------------------------------------
# 5_10_Derive_Pathogen_Host_Vector_Links.R
# ------------------------------------------------------------------------------
# Purpose: Derive one row per disease-pathogen-host-vector by joining the WHO
#          pathogen-host network to the pathogen-vector backfill output and the
#          observational host-vector join table.
#
# Inputs : WHO network helper path for combined_who_network_canonical_zoonotic.csv
#          vector_screening_evidence_path("pathogen_vector_links_filled.csv")
#          pathogen_association_data/evidence/host_vector/
#          vector_host_links_join_ready.csv
# Output : WHO host-vector helper path for pathogen_host_vector_links.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

host_vector_dir <- vector_host_outputs_dir

who_path <- who_working_network_path()
pathogen_vector_path <- vector_screening_evidence_path("pathogen_vector_links_filled.csv")
host_vector_path <- file.path(host_vector_dir, "vector_host_links_join_ready.csv")
output_path <- who_network_host_vector_path("pathogen_host_vector_links.csv")

who_network <- read_csv(
  who_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

pathogen_vectors <- read_csv(
  pathogen_vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(pathogen_tax_id = clean_text(pathogen_tax_id))

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

who_pathogen_host <- who_network %>%
  filter(!is.na(Disease_name), !is.na(Pathogen), !is.na(PathogenTaxID), !is.na(Host), !is.na(HostTaxID)) %>%
  transmute(
    disease_name = Disease_name,
    disease_name_join = normalize_name_for_match(Disease_name),
    pathogen = Pathogen,
    pathogen_tax_id = clean_text(PathogenTaxID),
    pathogen_type = PathogenType,
    pathogen_family = PathogenFamily,
    pathogen_genus = PathogenGenus,
    host = Host,
    host_tax_id = clean_text(HostTaxID),
    host_class = HostClass,
    host_order = HostOrder,
    host_family = HostFamily,
    detection_method = DetectionMethod,
    main_source = MainSource
  ) %>%
  group_by(disease_name_join, disease_name, pathogen, pathogen_tax_id, host_tax_id) %>%
  summarise(
    pathogen_type = first_non_missing(pathogen_type),
    pathogen_family = first_non_missing(pathogen_family),
    pathogen_genus = first_non_missing(pathogen_genus),
    host = first_non_missing(host),
    host_class = first_non_missing(host_class),
    host_order = first_non_missing(host_order),
    host_family = first_non_missing(host_family),
    detection_method = collapse_unique(detection_method),
    main_source = collapse_unique(main_source),
    .groups = "drop"
  )

pathogen_vector_joinable <- pathogen_vectors %>%
  filter(
    !is.na(disease_name),
    !is.na(pathogen),
    !is.na(pathogen_tax_id),
    !is.na(candidate_vector_species),
    assignment_basis != "no_disease_vector_match"
  ) %>%
  mutate(
    disease_name_join = normalize_name_for_match(disease_name),
    vector_join_key = normalize_vector_key(candidate_vector_species)
  ) %>%
  dplyr::select(
    pv_disease_name = disease_name,
    pv_disease_name_clean = disease_name_clean,
    disease_name_join,
    pathogen,
    pathogen_tax_id,
    pv_pathogen_type = pathogen_type,
    pv_pathogen_family = pathogen_family,
    pv_pathogen_genus = pathogen_genus,
    pheic_risk,
    candidate_vector_species,
    candidate_vector_group,
    evidence_strength,
    vector_evidence_basis,
    vector_record_sources,
    vector_supporting_row_count,
    assignment_basis,
    vector_join_key
  ) %>%
  distinct()

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

pathogen_host_vector_links <- who_pathogen_host %>%
  inner_join(
    pathogen_vector_joinable,
    by = c("disease_name_join", "pathogen", "pathogen_tax_id"),
    relationship = "many-to-many"
  ) %>%
  inner_join(
    host_vector_joinable,
    by = c("host_tax_id", "vector_join_key"),
    relationship = "many-to-many"
  ) %>%
  transmute(
    disease_name = disease_name,
    disease_name_clean = pv_disease_name_clean,
    pathogen,
    pathogen_tax_id,
    pathogen_type = coalesce(pv_pathogen_type, pathogen_type),
    pathogen_family = coalesce(pv_pathogen_family, pathogen_family),
    pathogen_genus = coalesce(pv_pathogen_genus, pathogen_genus),
    pheic_risk,
    host,
    host_tax_id,
    host_class,
    host_order,
    host_family,
    detection_method,
    main_source,
    vector_species = candidate_vector_species,
    vector_group = candidate_vector_group,
    evidence_strength,
    vector_evidence_basis,
    vector_record_sources,
    vector_supporting_row_count,
    assignment_basis,
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
    link_type = "pathogen_host_vector_derived",
    vector_join_match_type = "exact_normalized_vector_name",
    taxonomy_caution = dplyr::if_else(
      coalesce(hv_vector_species_needs_review, FALSE) |
        coalesce(hv_vector_taxon_rank, "species") != "species",
      TRUE,
      FALSE
    )
  ) %>%
  arrange(disease_name, pathogen, host, vector_species)

duplicate_key_count <- pathogen_host_vector_links %>%
  count(disease_name, pathogen_tax_id, host_tax_id, vector_join_key, name = "n") %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate disease_name + pathogen_tax_id + host_tax_id + vector_join_key rows found in pathogen-host-vector output")
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(pathogen_host_vector_links, output_path, na = "")

cat("WHO pathogen-host rows used:", nrow(who_pathogen_host), "\n")
cat("Pathogen-vector rows used:", nrow(pathogen_vector_joinable), "\n")
cat("Host-vector join rows used:", nrow(host_vector_joinable), "\n")
cat("Pathogen-host-vector rows written:", nrow(pathogen_host_vector_links), "\n")
cat("Distinct diseases in output:", n_distinct(pathogen_host_vector_links$disease_name), "\n")
cat("Distinct pathogens in output:", n_distinct(pathogen_host_vector_links$pathogen_tax_id), "\n")
cat("Taxonomy caution rows:", sum(pathogen_host_vector_links$taxonomy_caution %in% TRUE), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Wrote pathogen-host-vector links to", output_path, "\n")
