# ------------------------------------------------------------------------------
# 5_9_VectorMap_Vector_Name_Cleanup.R
# ------------------------------------------------------------------------------
# Purpose: Add a conservative vector taxonomy-cleanup layer on top of the
#          WHO-host-filtered VectorMap host-vector table while preserving the
#          original VectorMap fields and routing unresolved names to review.
#
# Input  : pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_who_filtered.csv
# Outputs: pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_host_links_who_vector_cleaned.csv
#          pathogen_association_data/staged/vectormap/outputs/
#          vectormap_vector_taxonomy_review_needed.csv
#          pathogen_association_data/manual/vectormap/
#          vectormap_vector_taxonomy_manual_map.csv (seeded if absent)
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

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

normalize_vector_key <- function(x) {
  x <- clean_text(x)

  transliterated <- suppressWarnings(iconv(x, from = "", to = "ASCII//TRANSLIT"))
  transliterated[is.na(transliterated)] <- x[is.na(transliterated)]

  transliterated %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9() .-]+", " ") %>%
    stringr::str_squish() %>%
    dplyr::na_if("")
}

count_vector_words <- function(x) {
  cleaned <- clean_text(x)
  dplyr::if_else(
    is.na(cleaned),
    0L,
    stringr::str_count(cleaned, " ") + 1L
  )
}

extract_vector_binomial <- function(x) {
  cleaned <- clean_text(x)
  parsed <- lapply(cleaned, function(value) {
    if (is.na(value)) {
      return(NA_character_)
    }

    tokens <- stringr::str_split(value, "\\s+", simplify = TRUE)
    tokens <- tokens[tokens != ""]

    if (length(tokens) < 2) {
      return(NA_character_)
    }

    paste(tokens[[1]], tokens[[2]])
  })

  dplyr::na_if(unlist(parsed, use.names = FALSE), "")
}

extract_vector_infraspecific_epithet <- function(x) {
  cleaned <- clean_text(x)
  parsed <- lapply(cleaned, function(value) {
    if (is.na(value)) {
      return(NA_character_)
    }

    tokens <- stringr::str_split(value, "\\s+", simplify = TRUE)
    tokens <- tokens[tokens != ""]

    if (length(tokens) < 3) {
      return(NA_character_)
    }

    tokens[[3]]
  })

  dplyr::na_if(unlist(parsed, use.names = FALSE), "")
}

classify_vector_taxon_rank <- function(x) {
  x <- clean_text(x)
  x_key <- normalize_vector_key(x)

  dplyr::case_when(
    is.na(x_key) ~ NA_character_,
    stringr::str_detect(x_key, "\\bsp\\.?\\b|\\bspp\\.?\\b") ~ "genus_only",
    stringr::str_detect(x_key, "complex") ~ "complex",
    stringr::str_detect(x_key, "\\bgroup\\b") ~ "group",
    stringr::str_detect(x_key, "sensu lato|\\bs\\.?l\\.?\\b") ~ "sensu_lato",
    stringr::str_detect(x_key, " x ") ~ "hybrid",
    stringr::str_detect(x_key, "\\bbiotype\\b|\\bform\\b|\\bsubsp\\.?\\b|\\bssp\\.?\\b|\\bvar\\.?\\b|\\bsubspecies\\b") ~ "infraspecific",
    count_vector_words(x) == 1L ~ "genus_only",
    count_vector_words(x) == 2L ~ "species",
    TRUE ~ "infraspecific"
  )
}

append_method <- function(methods, new_method) {
  if (length(methods) == 0) {
    return(new_method)
  }

  paste(c(methods, new_method), collapse = "; ")
}

apply_rule_cleanup <- function(x) {
  parsed <- lapply(clean_text(x), function(value) {
    if (is.na(value)) {
      return(list(
        vector_name_rule_cleaned = NA_character_,
        vector_cleanup_method = "missing"
      ))
    }

    cleaned <- stringr::str_replace_all(value, "_", " ")
    cleaned <- stringr::str_squish(cleaned)
    methods <- character(0)

    tokens <- stringr::str_split(cleaned, "\\s+", simplify = TRUE)
    tokens <- tokens[tokens != ""]

    if (length(tokens) >= 2 && stringr::str_to_lower(tokens[[1]]) == stringr::str_to_lower(tokens[[2]])) {
      tokens <- tokens[-2]
      cleaned <- paste(tokens, collapse = " ")
      methods <- append_method(methods, "rule_repeated_genus")
    }

    cleaned_key <- stringr::str_to_lower(cleaned)

    if (stringr::str_detect(cleaned, "^[A-Za-z-]+ \\([A-Za-z-]+\\)\\s+")) {
      cleaned <- stringr::str_replace(
        cleaned,
        "^([A-Za-z-]+) \\([A-Za-z-]+\\)\\s+",
        "\\1 "
      )
      methods <- append_method(methods, "rule_drop_parenthetical_subgenus")
      cleaned_key <- stringr::str_to_lower(cleaned)
    }

    if (stringr::str_detect(cleaned_key, "^aedes (ochlerotatus|neomelaniconion)\\b")) {
      cleaned <- stringr::str_replace(
        cleaned,
        regex("^aedes (ochlerotatus|neomelaniconion)\\s+", ignore_case = TRUE),
        "Aedes "
      )
      methods <- append_method(methods, "rule_drop_subgenus_token")
      cleaned_key <- stringr::str_to_lower(cleaned)
    }

    if (stringr::str_detect(cleaned_key, "^culex melanoconion\\b")) {
      cleaned <- stringr::str_replace(
        cleaned,
        regex("^culex melanoconion\\s+", ignore_case = TRUE),
        "Culex "
      )
      methods <- append_method(methods, "rule_drop_subgenus_token")
    }

    tokens <- stringr::str_split(cleaned, "\\s+", simplify = TRUE)
    tokens <- tokens[tokens != ""]

    if (length(tokens) >= 3) {
      author_start <- which(
        seq_along(tokens) >= 3 &
          (
            stringr::str_detect(tokens, "[A-Z]") |
              stringr::str_detect(tokens, "^[12][0-9]{3}$")
          )
      )

      if (length(author_start) > 0) {
        cleaned <- paste(tokens[seq_len(author_start[[1]] - 1)], collapse = " ")
        methods <- append_method(methods, "rule_strip_authorship_suffix")
      }
    }

    cleaned <- stringr::str_squish(cleaned)

    list(
      vector_name_rule_cleaned = cleaned,
      vector_cleanup_method = if (length(methods) == 0) "no_change" else methods
    )
  })

  tibble(
    vector_name_rule_cleaned = vapply(parsed, `[[`, character(1), "vector_name_rule_cleaned"),
    vector_cleanup_method = vapply(parsed, `[[`, character(1), "vector_cleanup_method")
  ) %>%
    mutate(across(everything(), ~ dplyr::na_if(.x, "")))
}

seed_manual_map <- function(path) {
  seeded_map <- tibble(
    source_name = character(),
    canonical_name = character(),
    vector_taxon_rank = character(),
    cleanup_method = character(),
    notes = character()
  )

  write_csv(seeded_map, path, na = "")
}

outputs_dir <- vectormap_outputs_dir
manual_dir <- vectormap_manual_dir

input_path <- file.path(outputs_dir, "vectormap_vector_host_links_who_filtered.csv")
output_path <- file.path(outputs_dir, "vectormap_vector_host_links_who_vector_cleaned.csv")
review_path <- file.path(outputs_dir, "vectormap_vector_taxonomy_review_needed.csv")
manual_map_path <- file.path(manual_dir, "vectormap_vector_taxonomy_manual_map.csv")

dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manual_map_path)) {
  seed_manual_map(manual_map_path)
}

vectormap_links <- read_csv(
  input_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

manual_map <- read_csv(
  manual_map_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_input_cols <- c(
  "source_dataset",
  "vector_family",
  "vector_genus",
  "vector_species",
  "vector_scientific_name"
)
missing_input_cols <- setdiff(required_input_cols, names(vectormap_links))
if (length(missing_input_cols) > 0) {
  stop(
    "vectormap_vector_host_links_who_filtered.csv is missing required columns: ",
    paste(missing_input_cols, collapse = ", ")
  )
}

required_map_cols <- c(
  "source_name",
  "canonical_name",
  "vector_taxon_rank",
  "cleanup_method",
  "notes"
)
missing_map_cols <- setdiff(required_map_cols, names(manual_map))
if (length(missing_map_cols) > 0) {
  stop(
    "vectormap_vector_taxonomy_manual_map.csv is missing required columns: ",
    paste(missing_map_cols, collapse = ", ")
  )
}

manual_map <- manual_map %>%
  mutate(source_name = normalize_vector_key(source_name)) %>%
  distinct(source_name, .keep_all = TRUE)

vectormap_links <- vectormap_links %>%
  mutate(
    vector_name_working = dplyr::case_when(
      !is.na(vector_scientific_name) ~ vector_scientific_name,
      !is.na(vector_genus) & !is.na(vector_species) ~ stringr::str_squish(paste(vector_genus, vector_species)),
      !is.na(vector_genus) ~ vector_genus,
      TRUE ~ NA_character_
    ),
    vector_name_clean = clean_text(vector_name_working)
  )

rule_cleaned <- apply_rule_cleanup(vectormap_links$vector_name_clean)

taxonomy_cleaned <- vectormap_links %>%
  bind_cols(rule_cleaned) %>%
  mutate(
    vector_name_map_key = normalize_vector_key(vector_name_rule_cleaned)
  ) %>%
  left_join(
    manual_map %>%
      transmute(
        source_name,
        manual_canonical_name = canonical_name,
        manual_vector_taxon_rank = vector_taxon_rank,
        manual_cleanup_method = cleanup_method,
        manual_notes = notes
      ),
    by = c("vector_name_map_key" = "source_name")
  ) %>%
  mutate(
    vector_name_taxonomy_cleaned = dplyr::coalesce(
      manual_canonical_name,
      vector_name_rule_cleaned
    ),
    vector_taxon_rank = dplyr::coalesce(
      manual_vector_taxon_rank,
      classify_vector_taxon_rank(vector_name_taxonomy_cleaned)
    ),
    vector_cleanup_method = dplyr::case_when(
      !is.na(manual_canonical_name) ~ manual_cleanup_method,
      TRUE ~ vector_cleanup_method
    ),
    vector_review_note = dplyr::case_when(
      is.na(vector_name_taxonomy_cleaned) ~ "Missing vector name after fallback from scientific name to genus/species.",
      stringr::str_detect(stringr::str_to_lower(vector_name_taxonomy_cleaned), "^stegomyia\\b") ~ "Old genus retained after cleanup; add a manual mapping decision.",
      stringr::str_detect(stringr::str_to_lower(vector_name_taxonomy_cleaned), "\\bsp\\.?\\b|\\bspp\\.?\\b") ~ "Genus-level sp./spp. label retained; do not force to species automatically.",
      vector_taxon_rank == "genus_only" ~ "Genus-only label retained; decide later whether genus-level matching is acceptable.",
      vector_taxon_rank == "complex" ~ "Species-complex label retained; do not collapse to species automatically.",
      vector_taxon_rank == "group" ~ "Group-level label retained; manual decision needed if species-level joins are required.",
      vector_taxon_rank == "sensu_lato" ~ "Sensu lato label retained; manual decision needed if species-level joins are required.",
      vector_taxon_rank == "hybrid" ~ "Hybrid label retained; do not collapse to a single species automatically.",
      vector_taxon_rank == "infraspecific" ~ "Infraspecific label retained; decide later whether species-level collapsing is appropriate.",
      !is.na(manual_notes) ~ manual_notes,
      TRUE ~ NA_character_
    ),
    vector_review_needed = !is.na(vector_review_note)
  ) %>%
  select(
    everything(),
    -vector_name_rule_cleaned,
    -vector_name_map_key,
    -manual_canonical_name,
    -manual_vector_taxon_rank,
    -manual_cleanup_method,
    -manual_notes
  )

binomial_presence_lookup <- taxonomy_cleaned %>%
  mutate(vector_species_binomial = extract_vector_binomial(vector_name_taxonomy_cleaned)) %>%
  filter(!is.na(vector_species_binomial)) %>%
  count(vector_species_binomial, name = "binomial_total_rows")

taxonomy_cleaned <- taxonomy_cleaned %>%
  mutate(
    vector_species_binomial = extract_vector_binomial(vector_name_taxonomy_cleaned),
    vector_infraspecific_epithet = extract_vector_infraspecific_epithet(vector_name_taxonomy_cleaned)
  ) %>%
  left_join(binomial_presence_lookup, by = "vector_species_binomial") %>%
  group_by(vector_name_taxonomy_cleaned) %>%
  mutate(current_label_rows = dplyr::n()) %>%
  ungroup() %>%
  mutate(
    binomial_present_elsewhere = !is.na(vector_species_binomial) &
      !is.na(binomial_total_rows) &
      (
        vector_taxon_rank != "infraspecific" |
          (binomial_total_rows - current_label_rows) > 0
      ),
    vector_species_collapse_method = dplyr::case_when(
      is.na(vector_name_taxonomy_cleaned) ~ "missing_vector_name",
      vector_taxon_rank == "species" ~ "already_species",
      vector_taxon_rank %in% c("genus_only", "complex", "group", "sensu_lato", "hybrid") ~ paste0("exclude_", vector_taxon_rank),
      vector_taxon_rank == "infraspecific" &
        !is.na(vector_infraspecific_epithet) &
        !is.na(vector_species_binomial) &
        stringr::str_to_lower(vector_infraspecific_epithet) ==
          stringr::str_to_lower(stringr::word(vector_species_binomial, 2)) ~
        "auto_drop_repeated_species_epithet",
      vector_taxon_rank == "infraspecific" &
        !is.na(vector_species_binomial) &
        binomial_present_elsewhere ~
        "collapse_to_binomial_base_present",
      vector_taxon_rank == "infraspecific" &
        !is.na(vector_species_binomial) ~
        "review_infraspecific_binomial_absent",
      TRUE ~ "no_species_level_decision"
    ),
    vector_species_collapsed_for_analysis = dplyr::case_when(
      vector_species_collapse_method == "already_species" ~ vector_name_taxonomy_cleaned,
      vector_species_collapse_method %in% c(
        "auto_drop_repeated_species_epithet",
        "collapse_to_binomial_base_present"
      ) ~ vector_species_binomial,
      TRUE ~ NA_character_
    ),
    vector_review_note = dplyr::case_when(
      is.na(vector_name_taxonomy_cleaned) ~ "Missing vector name after fallback from scientific name to genus/species.",
      stringr::str_detect(stringr::str_to_lower(vector_name_taxonomy_cleaned), "^stegomyia\\b") ~ "Old genus retained after cleanup; add a manual mapping decision.",
      stringr::str_detect(stringr::str_to_lower(vector_name_taxonomy_cleaned), "\\bsp\\.?\\b|\\bspp\\.?\\b") ~ "Genus-level sp./spp. label retained; do not force to species automatically.",
      vector_taxon_rank == "genus_only" ~ "Genus-only label retained; decide later whether genus-level matching is acceptable.",
      vector_taxon_rank == "complex" ~ "Species-complex label retained; do not collapse to species automatically.",
      vector_taxon_rank == "group" ~ "Group-level label retained; manual decision needed if species-level joins are required.",
      vector_taxon_rank == "sensu_lato" ~ "Sensu lato label retained; manual decision needed if species-level joins are required.",
      vector_taxon_rank == "hybrid" ~ "Hybrid label retained; do not collapse to a single species automatically.",
      vector_species_collapse_method == "auto_drop_repeated_species_epithet" ~ NA_character_,
      vector_species_collapse_method == "collapse_to_binomial_base_present" ~ "Infraspecific label collapsed to species because the binomial base already exists elsewhere in the cleaned VectorMap table.",
      vector_species_collapse_method == "review_infraspecific_binomial_absent" ~ "Infraspecific label retained; binomial base is not otherwise present in the cleaned VectorMap table, so species-level collapsing still needs review.",
      TRUE ~ NA_character_
    ),
    vector_review_needed = !is.na(vector_review_note)
  ) %>%
  select(-binomial_total_rows, -binomial_present_elsewhere, -current_label_rows)

review_table <- taxonomy_cleaned %>%
  filter(vector_review_needed) %>%
  group_by(
    vector_name_working,
    vector_name_clean,
    vector_name_taxonomy_cleaned,
    vector_species_binomial,
    vector_infraspecific_epithet,
    vector_taxon_rank,
    vector_cleanup_method,
    vector_species_collapse_method,
    vector_species_collapsed_for_analysis,
    vector_review_note
  ) %>%
  summarise(
    source_dataset_examples = collapse_unique(source_dataset),
    vector_family_examples = collapse_unique(vector_family),
    vector_genus_examples = collapse_unique(vector_genus),
    row_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(desc(row_count), vector_name_clean)

if (nrow(taxonomy_cleaned) != nrow(vectormap_links)) {
  stop("Row count changed during VectorMap vector-name cleanup")
}

write_csv(taxonomy_cleaned, output_path, na = "")
write_csv(review_table, review_path, na = "")

cat("Input VectorMap host-vector rows:", nrow(vectormap_links), "\n")
cat("Rows in cleaned output:", nrow(taxonomy_cleaned), "\n")
cat("Manual vector map entries:", nrow(manual_map), "\n")
cat(
  "Rows changed by cleanup:",
  sum(
    !is.na(taxonomy_cleaned$vector_name_clean) &
      !is.na(taxonomy_cleaned$vector_name_taxonomy_cleaned) &
      taxonomy_cleaned$vector_name_clean != taxonomy_cleaned$vector_name_taxonomy_cleaned
  ),
  "\n"
)
cat(
  "Rows with species-level collapsed vector names:",
  sum(!is.na(taxonomy_cleaned$vector_species_collapsed_for_analysis)),
  "\n"
)
cat("Review-needed vector labels:", nrow(review_table), "\n")
cat("Wrote cleaned VectorMap table to", output_path, "\n")
cat("Wrote VectorMap vector review table to", review_path, "\n")
