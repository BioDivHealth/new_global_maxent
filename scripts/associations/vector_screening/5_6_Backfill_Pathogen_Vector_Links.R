# ------------------------------------------------------------------------------
# 5_6_Backfill_Pathogen_Vector_Links.R
# ------------------------------------------------------------------------------
# Purpose: Expand the pathogen-vector scaffold into one row per
#          disease-pathogen-vector by inheriting canonical disease-level vector
#          evidence to every pathogen row for that disease.
#
# Inputs : pathogen_vector_links.csv
#          disease_vector_links.csv
# Output : pathogen_vector_links_filled.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

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

extract_vector_genus <- function(x) {
  x <- clean_text(x)
  has_binomial <- !is.na(x) & stringr::str_detect(x, "^[A-Z][a-z]+\\s+[a-z]")
  out <- rep(NA_character_, length(x))
  out[has_binomial] <- stringr::word(x[has_binomial], 1)
  out
}

scaffold_path <- vector_screening_staged_path("pathogen_vector_links.csv")
canonical_path <- vector_screening_staged_path("disease_vector_links.csv")
output_path <- file.path(vector_screening_evidence_dir, "pathogen_vector_links_filled.csv")
dir.create(vector_screening_evidence_dir, recursive = TRUE, showWarnings = FALSE)

scaffold <- read_csv(
  scaffold_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

canonical_vectors <- read_csv(
  canonical_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_scaffold_cols <- c(
  "disease_name",
  "pathogen",
  "candidate_vector_species",
  "candidate_vector_genus",
  "candidate_vector_family",
  "candidate_vector_group",
  "vector_status",
  "evidence_type",
  "evidence_strength"
)

missing_scaffold_cols <- setdiff(required_scaffold_cols, names(scaffold))
if (length(missing_scaffold_cols) > 0) {
  stop(
    "pathogen_vector_links.csv is missing required columns: ",
    paste(missing_scaffold_cols, collapse = ", ")
  )
}

required_canonical_cols <- c(
  "disease_name",
  "disease_name_clean",
  "vector_species",
  "vector_group",
  "best_evidence_level",
  "best_evidence_basis",
  "record_sources",
  "supporting_row_count"
)

missing_canonical_cols <- setdiff(required_canonical_cols, names(canonical_vectors))
if (length(missing_canonical_cols) > 0) {
  stop(
    "disease_vector_links.csv is missing required columns: ",
    paste(missing_canonical_cols, collapse = ", ")
  )
}

scaffold_joinable <- scaffold %>%
  mutate(disease_name_join = normalize_name_for_match(disease_name))

canonical_joinable <- canonical_vectors %>%
  mutate(disease_name_join = normalize_name_for_match(disease_name))

matched_rows <- scaffold_joinable %>%
  inner_join(
    canonical_joinable %>%
      dplyr::select(
        disease_name_join,
        vector_species,
        vector_group,
        best_evidence_level,
        best_evidence_basis,
        record_sources,
        supporting_row_count
      ),
    by = "disease_name_join",
    relationship = "many-to-many"
  ) %>%
  mutate(
    candidate_vector_species = vector_species,
    candidate_vector_genus = extract_vector_genus(vector_species),
    candidate_vector_family = NA_character_,
    candidate_vector_group = vector_group,
    vector_status = "inherited_from_disease_level",
    evidence_type = "disease_vector_evidence",
    evidence_strength = best_evidence_level,
    assignment_basis = "disease_level_inherited",
    vector_evidence_basis = best_evidence_basis,
    vector_record_sources = record_sources,
    vector_supporting_row_count = supporting_row_count
  ) %>%
  dplyr::select(-disease_name_join, -vector_species, -vector_group, -best_evidence_level,
         -best_evidence_basis, -record_sources, -supporting_row_count)

unmatched_rows <- scaffold_joinable %>%
  anti_join(
    canonical_joinable %>% distinct(disease_name_join),
    by = "disease_name_join"
  ) %>%
  mutate(
    candidate_vector_species = NA_character_,
    candidate_vector_genus = NA_character_,
    candidate_vector_family = NA_character_,
    candidate_vector_group = NA_character_,
    vector_status = NA_character_,
    evidence_type = NA_character_,
    evidence_strength = NA_character_,
    assignment_basis = "no_disease_vector_match",
    vector_evidence_basis = NA_character_,
    vector_record_sources = NA_character_,
    vector_supporting_row_count = NA_integer_
  ) %>%
  dplyr::select(-disease_name_join)

pathogen_vector_links_filled <- bind_rows(matched_rows, unmatched_rows) %>%
  arrange(disease_name, pathogen, candidate_vector_group, candidate_vector_species)

duplicate_key_count <- pathogen_vector_links_filled %>%
  mutate(
    vector_species_key = dplyr::coalesce(candidate_vector_species, "__NA__"),
    assignment_key = dplyr::coalesce(assignment_basis, "__NA__")
  ) %>%
  count(
    disease_name,
    pathogen,
    pathogen_tax_id,
    vector_species_key,
    assignment_key,
    name = "n"
  ) %>%
  filter(n > 1) %>%
  nrow()

if (duplicate_key_count > 0) {
  stop("Duplicate keys found in pathogen_vector_links_filled.csv output")
}

filled_vector_rows <- sum(
  !is.na(pathogen_vector_links_filled$candidate_vector_species) &
    pathogen_vector_links_filled$candidate_vector_species != ""
)

unresolved_disease_count <- pathogen_vector_links_filled %>%
  filter(assignment_basis == "no_disease_vector_match") %>%
  distinct(disease_name) %>%
  nrow()

write_csv(pathogen_vector_links_filled, output_path, na = "")

cat("Input scaffold rows:", nrow(scaffold), "\n")
cat("Canonical disease-vector rows used:", nrow(canonical_vectors), "\n")
cat("Filled pathogen-vector rows written:", nrow(pathogen_vector_links_filled), "\n")
cat("Distinct diseases in output:", n_distinct(pathogen_vector_links_filled$disease_name), "\n")
cat("Duplicate key check passed:", duplicate_key_count == 0, "\n")
cat("Rows with assigned vectors:", filled_vector_rows, "\n")
cat("Unresolved scaffold diseases:", unresolved_disease_count, "\n")
cat("Wrote filled scaffold output to", output_path, "\n")
