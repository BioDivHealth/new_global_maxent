# ------------------------------------------------------------------------------
# 5_5_Consolidate_Vector_Evidence.R
# ------------------------------------------------------------------------------
# Purpose: Collapse the standardized disease-vector evidence table into a
#          canonical WHO disease-vector table with one row per disease-vector
#          pair while preserving provenance summaries.
#
# Input  : vector_table_with_efsa_standardized.csv
# Outputs: disease_vector_links.csv
#          disease_vector_link_gaps.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tidyr)

source(here("scripts", "associations", "working_inputs.R"))

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
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

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

pick_most_frequent_label <- function(x) {
  x <- clean_text(x)
  x <- stats::na.omit(x)

  if (length(x) == 0) {
    return(NA_character_)
  }

  counts <- sort(table(x), decreasing = TRUE)
  top_count <- unname(counts[[1]])
  candidates <- sort(names(counts)[counts == top_count])
  candidates[[1]]
}

rank_evidence_level <- function(x) {
  x <- stringr::str_to_lower(clean_text(x))
  dplyr::case_when(
    x == "confirmed" ~ 1L,
    x == "probable" ~ 2L,
    x == "candidate" ~ 3L,
    x == "poor/unsupported" ~ 4L,
    TRUE ~ 5L
  )
}

vector_output_dir <- vector_screening_staged_outputs_dir
dir.create(vector_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_screening_qa_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- vector_screening_efsa_staged_path("vector_table_with_efsa_standardized.csv")
scaffold_path <- vector_screening_staged_path("pathogen_vector_links.csv")
output_path <- file.path(vector_output_dir, "disease_vector_links.csv")
gap_output_path <- file.path(vector_screening_qa_dir, "disease_vector_link_gaps.csv")

vector_table <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

scaffold <- read_csv(
  scaffold_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_cols <- c(
  "matched_network_disease_name",
  "disease_name_clean",
  "v_species",
  "v_species_clean",
  "v_group",
  "v_group_clean",
  "evidence_level",
  "evidence_basis",
  "location",
  "source",
  "notes",
  "record_source"
)

missing_cols <- setdiff(required_cols, names(vector_table))
if (length(missing_cols) > 0) {
  stop(
    "vector_table_with_efsa_standardized.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

vector_table_grouped <- vector_table %>%
  filter(
    !is.na(matched_network_disease_name),
    !is.na(v_species_clean),
    !is.na(v_group_clean)
  ) %>%
  mutate(
    disease_name = matched_network_disease_name,
    disease_name_clean = normalize_name_for_match(matched_network_disease_name),
    evidence_rank = rank_evidence_level(evidence_level),
    record_source_rank = dplyr::case_when(
      record_source == "lit_review" ~ 1L,
      record_source == "EFSA" ~ 2L,
      TRUE ~ 3L
    ),
    has_location = if_else(!is.na(location), 0L, 1L),
    source_alpha = dplyr::coalesce(source, "zzzz")
  )

best_evidence_rows <- vector_table_grouped %>%
  arrange(
    disease_name_clean,
    v_species_clean,
    v_group_clean,
    evidence_rank,
    record_source_rank,
    has_location,
    source_alpha
  ) %>%
  group_by(disease_name_clean, v_species_clean, v_group_clean) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    disease_name_clean,
    vector_species_clean = v_species_clean,
    vector_group_clean = v_group_clean,
    best_evidence_level = evidence_level,
    best_evidence_basis = evidence_basis
  )

disease_vector_links <- vector_table_grouped %>%
  group_by(disease_name_clean, v_species_clean, v_group_clean) %>%
  summarise(
    disease_name = dplyr::first(disease_name),
    vector_species = pick_most_frequent_label(v_species),
    vector_group = pick_most_frequent_label(v_group),
    record_sources = collapse_unique(record_source),
    source_count = dplyr::n_distinct(record_source, na.rm = TRUE),
    supporting_row_count = dplyr::n(),
    locations_summary = collapse_unique(location),
    sources_summary = collapse_unique(source),
    notes_summary = collapse_unique(notes),
    has_lit_review = any(record_source == "lit_review", na.rm = TRUE),
    has_efsa = any(record_source == "EFSA", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    best_evidence_rows,
    by = c(
      "disease_name_clean",
      "v_species_clean" = "vector_species_clean",
      "v_group_clean" = "vector_group_clean"
    )
  ) %>%
  transmute(
    disease_name,
    disease_name_clean,
    vector_species,
    vector_species_clean = v_species_clean,
    vector_group,
    vector_group_clean = v_group_clean,
    best_evidence_level,
    best_evidence_basis,
    record_sources,
    source_count,
    supporting_row_count,
    locations_summary,
    sources_summary,
    notes_summary,
    has_lit_review,
    has_efsa
  ) %>%
  arrange(disease_name, vector_group, vector_species)

duplicate_key_count <- disease_vector_links %>%
  count(disease_name_clean, vector_species_clean, vector_group_clean, name = "n") %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate keys found in disease_vector_links.csv output")
}

scaffold_counts <- scaffold %>%
  mutate(disease_name_clean_join = normalize_name_for_match(disease_name)) %>%
  group_by(disease_name, disease_name_clean_join) %>%
  summarise(
    in_scaffold = TRUE,
    n_pathogen_rows = dplyr::n(),
    .groups = "drop"
  )

canonical_counts <- disease_vector_links %>%
  group_by(disease_name, disease_name_clean) %>%
  summarise(
    in_canonical_vectors = TRUE,
    n_vector_rows = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(disease_name_clean_join = disease_name_clean)

gap_table <- full_join(
  scaffold_counts,
  canonical_counts,
  by = "disease_name_clean_join",
  suffix = c("_scaffold", "_canonical")
) %>%
  mutate(
    disease_name = dplyr::coalesce(disease_name_scaffold, disease_name_canonical),
    in_scaffold = dplyr::coalesce(in_scaffold, FALSE),
    in_canonical_vectors = dplyr::coalesce(in_canonical_vectors, FALSE),
    n_pathogen_rows = dplyr::coalesce(n_pathogen_rows, 0L),
    n_vector_rows = dplyr::coalesce(n_vector_rows, 0L),
    gap_status = dplyr::case_when(
      in_scaffold & in_canonical_vectors ~ "matched",
      in_scaffold & !in_canonical_vectors ~ "missing_canonical_vectors",
      !in_scaffold & in_canonical_vectors ~ "canonical_only",
      TRUE ~ "unknown"
    )
  ) %>%
  select(
    disease_name,
    in_scaffold,
    in_canonical_vectors,
    n_pathogen_rows,
    n_vector_rows,
    gap_status
  ) %>%
  arrange(desc(gap_status == "missing_canonical_vectors"), disease_name)

write_csv(disease_vector_links, output_path, na = "")
write_csv(gap_table, gap_output_path, na = "")

unresolved_disease_count <- sum(gap_table$gap_status == "missing_canonical_vectors")

cat("Input evidence rows:", nrow(vector_table_grouped), "\n")
cat("Canonical disease-vector rows written:", nrow(disease_vector_links), "\n")
cat("Distinct canonical diseases:", n_distinct(disease_vector_links$disease_name), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Gap rows written:", nrow(gap_table), "\n")
cat("Unresolved scaffold diseases:", unresolved_disease_count, "\n")
cat("Wrote canonical table to", output_path, "\n")
cat("Wrote gap table to", gap_output_path, "\n")
