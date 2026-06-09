# ------------------------------------------------------------------------------
# host_vector_join_helpers.R
# ------------------------------------------------------------------------------
# Purpose: Shared read and join-table preparation helpers for the host-vector
#          integration scripts.
# ------------------------------------------------------------------------------

read_clean_csv <- function(path, col_types = NULL) {
  args <- list(
    file = path,
    show_col_types = FALSE,
    na = c("", "NA")
  )

  if (!is.null(col_types)) {
    args$col_types <- col_types
  }

  do.call(readr::read_csv, args) %>%
    mutate(across(where(is.character), clean_text))
}

filter_legacy_compatible_host_network <- function(who_network) {
  legacy_flag <- "in_legacy_canonical_zoonotic_pathogen_host"

  if (!legacy_flag %in% names(who_network)) {
    stop(
      "Master-plus host network is missing required legacy compatibility flag: ",
      legacy_flag,
      call. = FALSE
    )
  }

  who_network %>%
    filter(.data[[legacy_flag]] %in% TRUE)
}

prepare_disease_host_network <- function(who_network, screened_diseases = NULL) {
  out <- who_network %>%
    filter(!is.na(Disease_name), !is.na(HostTaxID), !is.na(Host))

  if (!is.null(screened_diseases)) {
    out <- out %>%
      filter(Disease_name %in% screened_diseases)
  }

  out %>%
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
}

prepare_who_pathogen_host_network <- function(who_network) {
  who_network %>%
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
}

prepare_disease_vector_joinable <- function(disease_vectors) {
  disease_vectors %>%
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
}

prepare_pathogen_vector_joinable <- function(pathogen_vectors) {
  pathogen_vectors %>%
    mutate(pathogen_tax_id = clean_text(pathogen_tax_id)) %>%
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
}

prepare_host_vector_joinable <- function(host_vectors) {
  host_vectors %>%
    mutate(
      host_tax_id = clean_text(host_tax_id),
      vector_join_key = normalize_vector_key(vector_join_key)
    ) %>%
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
}

prepare_host_vector_keys <- function(host_vector_join) {
  host_vector_join %>%
    mutate(vector_join_key = normalize_vector_key(vector_join_key)) %>%
    filter(!is.na(vector_join_key)) %>%
    distinct(vector_join_key)
}
