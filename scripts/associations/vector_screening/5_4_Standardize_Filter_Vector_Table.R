# ------------------------------------------------------------------------------
# 5_4_Standardize_Filter_Vector_Table.R
# ------------------------------------------------------------------------------
# Purpose: Standardize disease/vector naming in the combined EFSA + lit-review
#          dataset, then keep only rows with diseases present in the combined
#          WHO network.
#
# Inputs : vector_table_with_efsa.csv
#          (from 5_3_Combine_LitReview_EFSA_Vector_Table.R)
#          master_plus_who_host_network.csv compatibility slice
#
# Output : vector_table_with_efsa_standardized.csv
#          vector_table_with_efsa_unmatched_diseases.csv
# ------------------------------------------------------------------------------

# ------------------------------| Load libraries |------------------------------
library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

# ------------------------------| Helper functions |----------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

# Aligns with 5_2 conventions while staying conservative.
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

remove_generic_disease_suffix <- function(x) {
  x <- normalize_name_for_match(x)
  x <- stringr::str_replace(
    x,
    stringr::regex(
      "(?:\\s+(virus|disease|fever|syndrome|infection))+$",
      ignore_case = TRUE
    ),
    ""
  )
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

extract_parenthetical_alias <- function(x) {
  x <- clean_text(x)
  alias <- stringr::str_match(x, "\\(([^)]+)\\)")[, 2]
  normalize_name_for_match(alias)
}

normalize_vector_group <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "-", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

# ------------------------------| Define paths |--------------------------------
vector_output_dir <- vector_screening_efsa_outputs_dir
vector_input_path <- vector_screening_efsa_staged_path("vector_table_with_efsa.csv")
dir.create(vector_output_dir, recursive = TRUE, showWarnings = FALSE)

output_standardized_path <- file.path(
  vector_output_dir,
  "vector_table_with_efsa_standardized.csv"
)
output_unmatched_path <- file.path(
  vector_output_dir,
  "vector_table_with_efsa_unmatched_diseases.csv"
)

# ------------------------------| Load datasets |-------------------------------
vector_table <- read_csv(
  vector_input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

combined_network <- read_legacy_compatible_master_plus_network() %>%
  mutate(across(where(is.character), clean_text))

# ------------------------------| Validate required columns |--------------------
required_vector_cols <- c("disease", "v_species", "v_group", "record_source")
missing_vector_cols <- setdiff(required_vector_cols, names(vector_table))
if (length(missing_vector_cols) > 0) {
  stop(
    "vector_table_with_efsa.csv is missing required columns: ",
    paste(missing_vector_cols, collapse = ", ")
  )
}

if (!("Disease_name" %in% names(combined_network))) {
  stop("The master-plus compatibility network is missing required column: Disease_name")
}

# ------------------------------| Prepare disease mapping |----------------------
# Keep manual mappings small and explicit for high-value pathogen-style labels
# that should resolve to disease labels in the WHO network.
manual_disease_map <- tibble::tribble(
  ~source_name_clean, ~network_disease_name,
  "sftsv", "Severe fever with thrombocytopenia syndrome (SFTS)"
)

manual_disease_map <- manual_disease_map %>%
  mutate(
    source_name_clean = normalize_name_for_match(source_name_clean),
    network_disease_name = clean_text(network_disease_name)
  )

# ------------------------------| Standardize names |---------------------------
vector_table_standardized <- vector_table %>%
  mutate(
    disease_name_clean = normalize_name_for_match(disease),
    disease_name_stem = remove_generic_disease_suffix(disease),
    v_species_clean = normalize_name_for_match(v_species),
    v_group_clean = normalize_vector_group(v_group)
  )

network_diseases <- combined_network %>%
  transmute(
    Disease_name,
    disease_name_clean = normalize_name_for_match(Disease_name),
    disease_name_stem = remove_generic_disease_suffix(Disease_name),
    disease_name_alias = extract_parenthetical_alias(Disease_name)
  ) %>%
  filter(!is.na(disease_name_clean)) %>%
  distinct()

network_disease_aliases <- bind_rows(
  network_diseases %>%
    transmute(
      network_disease_name = Disease_name,
      alias_clean = disease_name_clean,
      alias_type = "disease_name_exact"
    ),
  network_diseases %>%
    transmute(
      network_disease_name = Disease_name,
      alias_clean = disease_name_stem,
      alias_type = "disease_name_stem"
    ),
  network_diseases %>%
    transmute(
      network_disease_name = Disease_name,
      alias_clean = disease_name_alias,
      alias_type = "disease_name_alias"
    )
) %>%
  filter(!is.na(alias_clean)) %>%
  distinct()

network_disease_aliases_unique <- network_disease_aliases %>%
  group_by(alias_clean) %>%
  mutate(alias_candidate_count = dplyr::n_distinct(network_disease_name)) %>%
  ungroup() %>%
  filter(alias_candidate_count == 1) %>%
  distinct(alias_clean, network_disease_name, .keep_all = TRUE)

vector_table_standardized <- vector_table_standardized %>%
  left_join(
    manual_disease_map,
    by = c("disease_name_clean" = "source_name_clean")
  ) %>%
  left_join(
    network_disease_aliases_unique %>%
      transmute(
        disease_name_clean = alias_clean,
        alias_match_network_disease_name = network_disease_name,
        alias_match_type = alias_type
      ),
    by = c("disease_name_clean" = "disease_name_clean")
  ) %>%
  left_join(
    network_disease_aliases_unique %>%
      transmute(
        disease_name_stem = alias_clean,
        stem_match_network_disease_name = network_disease_name,
        stem_match_type = alias_type
      ),
    by = c("disease_name_stem" = "disease_name_stem")
  ) %>%
  mutate(
    matched_network_disease_name = dplyr::coalesce(
      network_disease_name,
      alias_match_network_disease_name,
      stem_match_network_disease_name
    ),
    disease_match_method = dplyr::case_when(
      !is.na(network_disease_name) ~ "manual_map",
      !is.na(alias_match_network_disease_name) ~ paste0("alias_", alias_match_type),
      !is.na(stem_match_network_disease_name) ~ paste0("stem_", stem_match_type),
      TRUE ~ NA_character_
    )
  )

# ------------------------------| Filter to in-network diseases |---------------
vector_table_filtered <- vector_table_standardized %>%
  filter(!is.na(matched_network_disease_name))

# One row per unmatched disease key with examples and source composition.
unmatched_diseases <- vector_table_standardized %>%
  filter(!is.na(disease_name_clean), is.na(matched_network_disease_name)) %>%
  group_by(disease_name_clean) %>%
  summarise(
    disease_examples = paste(unique(stats::na.omit(disease)), collapse = " | "),
    disease_name_stem = paste(unique(stats::na.omit(disease_name_stem)), collapse = " | "),
    n_rows = n(),
    record_sources = paste(sort(unique(stats::na.omit(record_source))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rows), disease_name_clean)

# ------------------------------| Write outputs |-------------------------------
write_csv(vector_table_filtered, output_standardized_path, na = "")
write_csv(unmatched_diseases, output_unmatched_path, na = "")

# ------------------------------| Console summary |-----------------------------
rows_total <- nrow(vector_table_standardized)
rows_matched <- nrow(vector_table_filtered)
rows_unmatched <- rows_total - rows_matched

raw_disease_labels_total <- n_distinct(vector_table_standardized$disease_name_clean, na.rm = TRUE)
raw_disease_labels_matched <- n_distinct(vector_table_filtered$disease_name_clean, na.rm = TRUE)
raw_disease_labels_unmatched <- nrow(unmatched_diseases)
network_diseases_matched <- n_distinct(vector_table_filtered$matched_network_disease_name, na.rm = TRUE)

source_counts_before <- vector_table_standardized %>%
  mutate(record_source = coalesce(record_source, "unknown")) %>%
  count(record_source, name = "n_rows_before")

source_counts_after <- vector_table_filtered %>%
  mutate(record_source = coalesce(record_source, "unknown")) %>%
  count(record_source, name = "n_rows_after")

source_counts_summary <- source_counts_before %>%
  full_join(source_counts_after, by = "record_source") %>%
  mutate(
    n_rows_before = coalesce(n_rows_before, 0L),
    n_rows_after = coalesce(n_rows_after, 0L)
  ) %>%
  arrange(record_source)

cat("Rows in vector table (input):", rows_total, "\n")
cat("Rows retained after disease filter:", rows_matched, "\n")
cat("Rows removed (disease not in network):", rows_unmatched, "\n")
cat("Distinct raw cleaned disease labels in input:", raw_disease_labels_total, "\n")
cat("Distinct raw cleaned disease labels retained:", raw_disease_labels_matched, "\n")
cat("Distinct matched network diseases:", network_diseases_matched, "\n")
cat("Distinct raw cleaned disease labels unmatched:", raw_disease_labels_unmatched, "\n")
cat("Wrote standardized filtered table to", output_standardized_path, "\n")
cat("Wrote unmatched disease summary to", output_unmatched_path, "\n")
cat("\nRows by record_source (before -> after filter):\n")
for (i in seq_len(nrow(source_counts_summary))) {
  cat(
    " - ",
    source_counts_summary$record_source[[i]],
    ": ",
    source_counts_summary$n_rows_before[[i]],
    " -> ",
    source_counts_summary$n_rows_after[[i]],
    "\n",
    sep = ""
  )
}
