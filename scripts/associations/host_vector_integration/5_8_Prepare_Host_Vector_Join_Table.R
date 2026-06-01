# ------------------------------------------------------------------------------
# 5_8_Prepare_Host_Vector_Join_Table.R
# ------------------------------------------------------------------------------
# Purpose: Build one canonical host-vector join table from the combined
#          MapVEu/VectorMap analysis-ready evidence, while preserving blocked
#          non-joinable rows for QA.
#
# Input  : pathogen_association_data/evidence/host_vector/
#          vector_host_links_analysis_ready.csv
# Outputs: pathogen_association_data/evidence/host_vector/
#          vector_host_links_join_ready.csv
#          vector_host_links_join_blocked.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

collapse_rank <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  if (length(x) == 1) {
    return(x[[1]])
  }

  "mixed"
}

vector_host_dir <- vector_host_outputs_dir
input_path <- file.path(vector_host_dir, "vector_host_links_analysis_ready.csv")
join_ready_path <- file.path(vector_host_dir, "vector_host_links_join_ready.csv")
blocked_path <- file.path(vector_host_dir, "vector_host_links_join_blocked.csv")

host_vector_ready <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    host_tax_id = clean_text(matched_who_host_tax_id),
    vector_join_key = normalize_vector_key(vector_species_analysis)
  )

blocked_rows <- host_vector_ready %>%
  mutate(
    block_reason = case_when(
      is.na(host_tax_id) & is.na(vector_join_key) ~ "missing_host_tax_id_and_vector_join_key",
      is.na(host_tax_id) ~ "missing_host_tax_id",
      is.na(vector_join_key) ~ "missing_vector_join_key",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(block_reason))

join_ready <- host_vector_ready %>%
  filter(!is.na(host_tax_id), !is.na(vector_join_key)) %>%
  group_by(host_tax_id, vector_join_key) %>%
  summarise(
    host = first_non_missing(matched_who_host),
    host_class = first_non_missing(matched_who_host_class),
    host_order = first_non_missing(matched_who_host_order),
    host_family = first_non_missing(matched_who_host_family),
    vector_species = first_non_missing(vector_species_analysis),
    vector_taxon_rank = collapse_rank(vector_taxon_rank),
    vector_species_needs_review = any(vector_species_needs_review %in% TRUE, na.rm = TRUE),
    vector_name_taxonomy_examples = collapse_unique(vector_name_taxonomy_cleaned),
    source_platform_examples = collapse_unique(source_platform),
    source_dataset_examples = collapse_unique(source_dataset),
    interaction_type_examples = collapse_unique(interaction_type),
    country_examples = collapse_unique(country),
    review_reason_examples = collapse_unique(review_reason),
    record_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(host, vector_species)

duplicate_key_count <- join_ready %>%
  count(host_tax_id, vector_join_key, name = "n") %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate host_tax_id + vector_join_key rows found in join-ready output")
}

write_csv(join_ready, join_ready_path, na = "")
write_csv(blocked_rows, blocked_path, na = "")

cat("Input host-vector record rows:", nrow(host_vector_ready), "\n")
cat("Join-ready host-vector rows written:", nrow(join_ready), "\n")
cat("Blocked host-vector rows written:", nrow(blocked_rows), "\n")
cat("Distinct hosts in join-ready output:", n_distinct(join_ready$host_tax_id), "\n")
cat("Distinct vectors in join-ready output:", n_distinct(join_ready$vector_join_key), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Wrote join-ready host-vector table to", join_ready_path, "\n")
cat("Wrote blocked host-vector rows to", blocked_path, "\n")
