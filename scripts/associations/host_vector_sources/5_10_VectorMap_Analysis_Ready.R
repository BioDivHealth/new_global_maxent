# ------------------------------------------------------------------------------
# 5_10_VectorMap_Analysis_Ready.R
# ------------------------------------------------------------------------------
# Purpose: Build analysis-ready VectorMap host-vector outputs from the cleaned
#          WHO-host-filtered table. One output stays at the original evidence
#          record grain, and the second collapses to unique host-vector pairs.
#
# Input  : pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_who_vector_cleaned.csv
# Outputs: pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_analysis_ready.csv
#          pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_analysis_summary.csv
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

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

outputs_dir <- vectormap_outputs_dir

input_path <- file.path(outputs_dir, "vectormap_vector_host_links_who_vector_cleaned.csv")
analysis_ready_path <- file.path(outputs_dir, "vectormap_vector_host_links_analysis_ready.csv")
analysis_summary_path <- file.path(outputs_dir, "vectormap_vector_host_links_analysis_summary.csv")

vectormap_cleaned <- read_csv(
  input_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_cols <- c(
  "vectormap_row_id",
  "record_id",
  "source_dataset",
  "interaction_type",
  "basis_of_record",
  "source_citation",
  "matched_who_host",
  "matched_who_host_tax_id",
  "matched_who_host_class",
  "matched_who_host_order",
  "matched_who_host_family",
  "vector_species_collapsed_for_analysis",
  "vector_name_taxonomy_cleaned",
  "vector_species_collapse_method",
  "vector_taxon_rank",
  "vector_family",
  "vector_genus",
  "vector_species",
  "vector_scientific_name",
  "country",
  "state_province",
  "county",
  "locality",
  "latitude",
  "longitude",
  "earliest_date_collected",
  "latest_date_collected",
  "collecting_method",
  "associated_pathogen",
  "associated_parasite"
)

missing_cols <- setdiff(required_cols, names(vectormap_cleaned))
if (length(missing_cols) > 0) {
  stop(
    "vectormap_vector_host_links_who_vector_cleaned.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

excluded_genus_only_n <- vectormap_cleaned %>%
  filter(vector_species_collapse_method == "exclude_genus_only") %>%
  nrow()

excluded_missing_name_n <- vectormap_cleaned %>%
  filter(vector_species_collapse_method == "missing_vector_name") %>%
  nrow()

analysis_ready <- vectormap_cleaned %>%
  filter(
    !vector_species_collapse_method %in% c("exclude_genus_only", "missing_vector_name"),
    !is.na(vector_species_collapsed_for_analysis)
  ) %>%
  transmute(
    vectormap_row_id,
    record_id,
    source_dataset,
    interaction_type,
    basis_of_record,
    source_citation,
    matched_who_host,
    matched_who_host_tax_id,
    matched_who_host_class,
    matched_who_host_order,
    matched_who_host_family,
    vector_species_analysis = vector_species_collapsed_for_analysis,
    vector_species_needs_review = vector_species_collapse_method == "collapse_to_binomial_base_present",
    vector_name_taxonomy_cleaned,
    vector_species_collapse_method,
    vector_taxon_rank,
    vector_family,
    vector_genus,
    vector_species,
    vector_scientific_name,
    country,
    state_province,
    county,
    locality,
    latitude,
    longitude,
    earliest_date_collected,
    latest_date_collected,
    collecting_method,
    associated_pathogen,
    associated_parasite
  )

analysis_summary <- analysis_ready %>%
  group_by(
    matched_who_host,
    matched_who_host_tax_id,
    matched_who_host_class,
    matched_who_host_order,
    matched_who_host_family,
    vector_species_analysis,
    vector_species_needs_review
  ) %>%
  summarise(
    vector_name_taxonomy_examples = collapse_unique(vector_name_taxonomy_cleaned),
    interaction_type_examples = collapse_unique(interaction_type),
    source_dataset_examples = collapse_unique(source_dataset),
    country_examples = collapse_unique(country),
    associated_pathogen_examples = collapse_unique(associated_pathogen),
    associated_parasite_examples = collapse_unique(associated_parasite),
    record_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(matched_who_host, vector_species_analysis)

expected_analysis_ready_rows <- nrow(vectormap_cleaned) - excluded_genus_only_n - excluded_missing_name_n
if (nrow(analysis_ready) != expected_analysis_ready_rows) {
  stop("Unexpected row count in analysis-ready VectorMap output")
}

if (any(is.na(analysis_ready$vector_species_analysis))) {
  stop("analysis_ready output contains missing vector_species_analysis values")
}

write_csv(analysis_ready, analysis_ready_path, na = "")
write_csv(analysis_summary, analysis_summary_path, na = "")

cat("Input cleaned VectorMap rows:", nrow(vectormap_cleaned), "\n")
cat("Excluded genus-only rows:", excluded_genus_only_n, "\n")
cat("Excluded missing-name rows:", excluded_missing_name_n, "\n")
cat(
  "Retained review-flagged rows:",
  sum(analysis_ready$vector_species_needs_review, na.rm = TRUE),
  "\n"
)
cat("Analysis-ready record rows:", nrow(analysis_ready), "\n")
cat("Analysis summary pairs:", nrow(analysis_summary), "\n")
cat("Wrote record-level analysis table to", analysis_ready_path, "\n")
cat("Wrote host-vector summary table to", analysis_summary_path, "\n")
