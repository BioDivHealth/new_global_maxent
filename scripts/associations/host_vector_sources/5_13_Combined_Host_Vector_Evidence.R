# ------------------------------------------------------------------------------
# 5_13_Combined_Host_Vector_Evidence.R
# ------------------------------------------------------------------------------
# Purpose: Combine the analysis-ready VectorMap and MapVEu host-vector evidence
#          tables into one record-level evidence file plus one deduplicated
#          host-vector summary, while preserving source provenance.
#
# Inputs : pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_analysis_ready.csv
#          pathogen_association_data/staged/mapveu/outputs/
#          mapveu_vector_host_links_analysis_ready.csv
# Outputs: pathogen_association_data/evidence/host_vector/
#          vector_host_links_analysis_ready.csv
#          pathogen_association_data/evidence/host_vector/
#          vector_host_links_analysis_summary.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
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

add_missing_columns <- function(df, all_columns) {
  missing_cols <- setdiff(all_columns, names(df))

  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      df[[col]] <- NA
    }
  }

  df[, all_columns]
}

vectormap_path <- file.path(
  vectormap_outputs_dir,
  "vectormap_vector_host_links_analysis_ready.csv"
)
mapveu_path <- file.path(
  mapveu_outputs_dir,
  "mapveu_vector_host_links_analysis_ready.csv"
)

output_dir <- vector_host_outputs_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

combined_ready_path <- file.path(output_dir, "vector_host_links_analysis_ready.csv")
combined_summary_path <- file.path(output_dir, "vector_host_links_analysis_summary.csv")

who_hosts <- read_legacy_compatible_master_plus_network() %>%
  mutate(across(where(is.character), clean_text)) %>%
  transmute(
    matched_who_host = clean_text(Host),
    who_matched_who_host_tax_id = clean_text(HostTaxID),
    who_matched_who_host_class = clean_text(HostClass),
    who_matched_who_host_order = clean_text(HostOrder),
    who_matched_who_host_family = clean_text(HostFamily)
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  group_by(matched_who_host) %>%
  summarise(
    who_matched_who_host_tax_id = first_non_missing(who_matched_who_host_tax_id),
    who_matched_who_host_class = first_non_missing(who_matched_who_host_class),
    who_matched_who_host_order = first_non_missing(who_matched_who_host_order),
    who_matched_who_host_family = first_non_missing(who_matched_who_host_family),
    .groups = "drop"
  )

vectormap <- read_csv(
  vectormap_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    vectormap_row_id = as.integer(vectormap_row_id),
    record_id = as.character(record_id),
    matched_who_host_tax_id = as.character(matched_who_host_tax_id),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    earliest_date_collected = as.character(earliest_date_collected),
    latest_date_collected = as.character(latest_date_collected),
    host_prevalence_percent = NA_real_,
    source_platform = "VectorMap",
    source_record_id = clean_text(record_id),
    source_record_type = "vector_occurrence_or_blood_meal_record",
    source_study_name = NA_character_,
    pubmed_id = NA_character_,
    doi = NA_character_,
    host_match_method = NA_character_,
    vector_species_source = NA_character_,
    review_needed = NA,
    review_reason = NA_character_,
    blood_meal_assay_id = NA_character_,
    sample_id = NA_character_,
    collection_id = NA_character_,
    collection_site_id = NA_character_,
    study_id = NA_character_,
    host_organism_raw = NA_character_,
    vector_species_raw = NA_character_,
    vector_species_assay_raw = NA_character_,
    vector_species_sample_raw = NA_character_,
    host_presence = NA_character_,
    host_prevalence_percent = NA_real_,
    collection_host_organism = NA_character_,
    collection_device = NA_character_,
    collection_start_date = earliest_date_collected,
    collection_end_date = latest_date_collected
  )

mapveu <- read_csv(
  mapveu_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    mapveu_row_id = as.integer(mapveu_row_id),
    blood_meal_assay_id = as.character(blood_meal_assay_id),
    sample_id = as.character(sample_id),
    collection_id = as.character(collection_id),
    collection_site_id = as.character(collection_site_id),
    study_id = as.character(study_id),
    matched_who_host_tax_id = as.character(matched_who_host_tax_id),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    collection_start_date = as.character(collection_start_date),
    collection_end_date = as.character(collection_end_date),
    host_prevalence_percent = suppressWarnings(as.numeric(host_prevalence_percent)),
    source_platform = "MapVEu",
    source_record_id = clean_text(blood_meal_assay_id),
    source_record_type = "blood_meal_assay_record",
    record_id = NA_character_,
    basis_of_record = NA_character_,
    source_citation = NA_character_,
    vector_species_collapse_method = NA_character_,
    vector_family = NA_character_,
    vector_genus = NA_character_,
    vector_species = NA_character_,
    vector_scientific_name = NA_character_,
    state_province = NA_character_,
    county = NA_character_,
    locality = NA_character_,
    earliest_date_collected = collection_start_date,
    latest_date_collected = collection_end_date,
    collecting_method = NA_character_,
    associated_pathogen = NA_character_,
    associated_parasite = NA_character_
  )

all_columns <- union(names(vectormap), names(mapveu))

combined_ready <- bind_rows(
  add_missing_columns(vectormap, all_columns),
  add_missing_columns(mapveu, all_columns)
) %>%
  left_join(who_hosts, by = "matched_who_host") %>%
  mutate(
    matched_who_host_tax_id = dplyr::coalesce(
      matched_who_host_tax_id,
      who_matched_who_host_tax_id
    ),
    matched_who_host_class = dplyr::coalesce(
      matched_who_host_class,
      who_matched_who_host_class
    ),
    matched_who_host_order = dplyr::coalesce(
      matched_who_host_order,
      who_matched_who_host_order
    ),
    matched_who_host_family = dplyr::coalesce(
      matched_who_host_family,
      who_matched_who_host_family
    )
  ) %>%
  transmute(
    source_platform,
    source_dataset,
    source_record_id,
    source_record_type,
    vectormap_row_id,
    mapveu_row_id,
    record_id,
    blood_meal_assay_id,
    sample_id,
    collection_id,
    collection_site_id,
    study_id,
    source_study_name,
    pubmed_id,
    doi,
    interaction_type,
    basis_of_record,
    source_citation,
    matched_who_host,
    matched_who_host_tax_id,
    matched_who_host_class,
    matched_who_host_order,
    matched_who_host_family,
    host_match_method,
    vector_species_analysis,
    vector_species_needs_review,
    vector_name_taxonomy_cleaned,
    vector_taxon_rank,
    vector_species_source,
    vector_species_collapse_method,
    review_needed,
    review_reason,
    vector_family,
    vector_genus,
    vector_species,
    vector_scientific_name,
    vector_species_raw,
    vector_species_assay_raw,
    vector_species_sample_raw,
    host_organism_raw,
    host_presence,
    host_prevalence_percent,
    collection_host_organism,
    collection_device,
    country,
    state_province,
    county,
    locality,
    latitude,
    longitude,
    earliest_date_collected,
    latest_date_collected,
    collection_start_date,
    collection_end_date,
    collecting_method,
    associated_pathogen,
    associated_parasite
  )

if (any(is.na(combined_ready$matched_who_host))) {
  stop("Combined analysis-ready table contains missing matched_who_host values")
}

if (any(is.na(combined_ready$vector_species_analysis))) {
  stop("Combined analysis-ready table contains missing vector_species_analysis values")
}

combined_summary <- combined_ready %>%
  group_by(
    matched_who_host,
    vector_species_analysis,
    vector_species_needs_review
  ) %>%
  summarise(
    matched_who_host_tax_id = first_non_missing(matched_who_host_tax_id),
    matched_who_host_class = first_non_missing(matched_who_host_class),
    matched_who_host_order = first_non_missing(matched_who_host_order),
    matched_who_host_family = first_non_missing(matched_who_host_family),
    vector_name_taxonomy_examples = collapse_unique(vector_name_taxonomy_cleaned),
    source_platform_examples = collapse_unique(source_platform),
    source_dataset_examples = collapse_unique(source_dataset),
    interaction_type_examples = collapse_unique(interaction_type),
    country_examples = collapse_unique(country),
    review_reason_examples = collapse_unique(review_reason),
    record_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(matched_who_host, vector_species_analysis)

duplicate_summary_keys <- combined_summary %>%
  count(matched_who_host, vector_species_analysis, vector_species_needs_review) %>%
  filter(n > 1)

if (nrow(duplicate_summary_keys) > 0) {
  stop("Combined summary has duplicate host-vector-review keys")
}

write_csv(combined_ready, combined_ready_path, na = "")
write_csv(combined_summary, combined_summary_path, na = "")

cat("VectorMap analysis-ready rows:", nrow(vectormap), "\n")
cat("MapVEu analysis-ready rows:", nrow(mapveu), "\n")
cat("Combined analysis-ready rows:", nrow(combined_ready), "\n")
cat("Combined unique hosts:", n_distinct(combined_ready$matched_who_host), "\n")
cat("Combined unique vectors:", n_distinct(combined_ready$vector_species_analysis), "\n")
cat("Combined summary pairs:", nrow(combined_summary), "\n")
cat("Wrote combined analysis-ready table to", combined_ready_path, "\n")
cat("Wrote combined summary table to", combined_summary_path, "\n")
