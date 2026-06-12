# ------------------------------------------------------------------------------
# 5_11_MapVEu_Analysis_Ready.R
# ------------------------------------------------------------------------------
# Purpose: Turn the staged MapVEu blood-meal table into analysis-ready,
#          VectorMap-compatible host-vector outputs.
#
# Input  : pathogen_association_data/staged/mapveu/outputs/
#          mapveu_vector_host_links_raw.csv
#          master_plus_who_host_network.csv compatibility slice
# Outputs: pathogen_association_data/staged/mapveu/outputs/
#          mapveu_vector_host_links_analysis_ready.csv
#          pathogen_association_data/staged/mapveu/outputs/
#          mapveu_vector_host_links_analysis_summary.csv
#          pathogen_association_data/manual/mapveu/
#          mapveu_host_manual_crosswalk.csv (seeded if absent)
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

normalize_key <- function(x) {
  x <- clean_text(x)

  transliterated <- suppressWarnings(iconv(x, from = "", to = "ASCII//TRANSLIT"))
  transliterated[is.na(transliterated)] <- x[is.na(transliterated)]

  transliterated %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9 ]+", " ") %>%
    stringr::str_squish() %>%
    dplyr::na_if("")
}

count_words <- function(x) {
  key <- normalize_key(x)
  dplyr::if_else(is.na(key), 0L, stringr::str_count(key, " ") + 1L)
}

is_non_actionable_host <- function(x) {
  key <- normalize_key(x)
  word_count <- count_words(x)

  is.na(key) |
    key == "" |
    word_count < 2L |
    stringr::str_detect(
      key,
      "\\b(sp|sp\\.|spp|spp\\.|species complex|complex|group|subgroup|unknown|unidentified|organism|organisms|multi species|multi species collection|multi|mixed|other|vertebrates)\\b"
    ) |
    stringr::str_detect(key, "\\band\\b|/")
}

combine_reason_pair <- function(x, y) {
  collapse_unique(c(x, y))
}

seed_host_crosswalk <- function(path) {
  seeded <- tibble::tribble(
    ~raw_host_label, ~source_field, ~matched_who_host, ~include_in_final, ~note,
    "Canis lupus familiaris", "host_organism_raw", "Canis lupus", TRUE, "Conservative domestic synonym mapped to WHO host."
  )

  write_csv(seeded, path, na = "")
}

clean_vector_single <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

classify_vector_taxon_rank <- function(x) {
  key <- normalize_key(x)
  dplyr::case_when(
    is.na(key) ~ NA_character_,
    stringr::str_detect(key, "\\bspecies complex\\b|\\bcomplex\\b") ~ "complex",
    stringr::str_detect(key, "\\bgroup\\b|\\bsubgroup\\b") ~ "group",
    stringr::str_detect(key, "sensu lato") ~ "sensu_lato",
    stringr::str_detect(key, " x ") ~ "hybrid",
    stringr::str_detect(key, "\\bno\\.?\\s*[0-9]+\\b") ~ "infraspecific",
    count_words(x) == 1L ~ "genus_only",
    count_words(x) == 2L ~ "species",
    TRUE ~ "infraspecific"
  )
}

parse_vector_field <- function(x) {
  values <- clean_text(x)

  parsed <- lapply(values, function(value) {
    if (is.na(value)) {
      return(list(
        vector_name_taxonomy_cleaned = NA_character_,
        vector_species_analysis = NA_character_,
        vector_taxon_rank = NA_character_,
        vector_species_needs_review = TRUE,
        vector_review_reason_extra = "missing_vector_name"
      ))
    }

    split_values <- stringr::str_split(value, "\\s*\\|\\s*", simplify = FALSE)[[1]]
    split_values <- clean_vector_single(split_values)
    split_values <- split_values[!is.na(split_values)]
    split_values <- unique(split_values)

    if (length(split_values) == 0) {
      return(list(
        vector_name_taxonomy_cleaned = NA_character_,
        vector_species_analysis = NA_character_,
        vector_taxon_rank = NA_character_,
        vector_species_needs_review = TRUE,
        vector_review_reason_extra = "missing_vector_name"
      ))
    }

    if (length(split_values) == 1) {
      cleaned_value <- split_values[[1]]
      rank <- classify_vector_taxon_rank(cleaned_value)
      needs_review <- !is.na(rank) && rank != "species"

      return(list(
        vector_name_taxonomy_cleaned = cleaned_value,
        vector_species_analysis = cleaned_value,
        vector_taxon_rank = rank,
        vector_species_needs_review = needs_review,
        vector_review_reason_extra = if (needs_review) paste0("vector_", rank) else NA_character_
      ))
    }

    return(list(
      vector_name_taxonomy_cleaned = paste(split_values, collapse = " | "),
      vector_species_analysis = NA_character_,
      vector_taxon_rank = "multiple_labels",
      vector_species_needs_review = TRUE,
      vector_review_reason_extra = "multiple_vector_labels"
    ))
  })

  tibble(
    vector_name_taxonomy_cleaned = vapply(parsed, `[[`, character(1), "vector_name_taxonomy_cleaned"),
    vector_species_analysis = vapply(parsed, `[[`, character(1), "vector_species_analysis"),
    vector_taxon_rank = vapply(parsed, `[[`, character(1), "vector_taxon_rank"),
    vector_species_needs_review = vapply(parsed, `[[`, logical(1), "vector_species_needs_review"),
    vector_review_reason_extra = vapply(parsed, `[[`, character(1), "vector_review_reason_extra")
  ) %>%
    mutate(across(c(vector_name_taxonomy_cleaned, vector_species_analysis, vector_taxon_rank, vector_review_reason_extra), ~ dplyr::na_if(.x, "")))
}

outputs_dir <- mapveu_outputs_dir
manual_dir <- mapveu_manual_dir

input_path <- file.path(outputs_dir, "mapveu_vector_host_links_raw.csv")
host_crosswalk_path <- file.path(manual_dir, "mapveu_host_manual_crosswalk.csv")
analysis_ready_path <- file.path(outputs_dir, "mapveu_vector_host_links_analysis_ready.csv")
analysis_summary_path <- file.path(outputs_dir, "mapveu_vector_host_links_analysis_summary.csv")

dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(host_crosswalk_path)) {
  seed_host_crosswalk(host_crosswalk_path)
}

mapveu_raw <- read_csv(
  input_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(mapveu_row_id = row_number())

who_hosts <- read_legacy_compatible_master_plus_network() %>%
  transmute(
    matched_who_host = clean_text(Host),
    matched_who_host_tax_id = clean_text(HostTaxID),
    matched_who_host_class = clean_text(HostClass),
    matched_who_host_order = clean_text(HostOrder),
    matched_who_host_family = clean_text(HostFamily)
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  group_by(matched_who_host) %>%
  summarise(
    matched_who_host_tax_id = first_non_missing(matched_who_host_tax_id),
    matched_who_host_class = first_non_missing(matched_who_host_class),
    matched_who_host_order = first_non_missing(matched_who_host_order),
    matched_who_host_family = first_non_missing(matched_who_host_family),
    .groups = "drop"
  )

who_lookup <- who_hosts %>%
  mutate(who_host_key = normalize_key(matched_who_host)) %>%
  distinct(who_host_key, .keep_all = TRUE)

host_crosswalk <- read_csv(
  host_crosswalk_path,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
) %>%
  transmute(
    raw_host_label = clean_text(raw_host_label),
    raw_host_key = normalize_key(raw_host_label),
    source_field = clean_text(source_field),
    matched_who_host = clean_text(matched_who_host),
    include_in_final = as.logical(include_in_final),
    note = clean_text(note)
  ) %>%
  filter(
    source_field == "host_organism_raw",
    !is.na(raw_host_key),
    !is.na(matched_who_host),
    include_in_final %in% TRUE
  )

mapveu_hosts <- mapveu_raw %>%
  mutate(
    host_key = normalize_key(host_organism_raw),
    host_non_actionable = is_non_actionable_host(host_organism_raw)
  ) %>%
  left_join(
    who_lookup,
    by = c("host_key" = "who_host_key")
  ) %>%
  left_join(
    host_crosswalk %>%
      select(raw_host_key, crosswalk_who_host = matched_who_host),
    by = c("host_key" = "raw_host_key")
  ) %>%
  mutate(
    matched_who_host = dplyr::coalesce(matched_who_host, crosswalk_who_host),
    host_match_method = dplyr::case_when(
      !is.na(crosswalk_who_host) ~ "manual_crosswalk",
      !is.na(matched_who_host) ~ "exact_host",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-crosswalk_who_host)

host_filtered <- mapveu_hosts %>%
  filter(!is.na(matched_who_host))

excluded_broad_host_rows <- mapveu_hosts %>%
  filter(is.na(matched_who_host), host_non_actionable) %>%
  nrow()

excluded_other_unmatched_rows <- mapveu_hosts %>%
  filter(is.na(matched_who_host), !host_non_actionable) %>%
  nrow()

vector_cleaned <- host_filtered %>%
  bind_cols(parse_vector_field(host_filtered$vector_species_raw)) %>%
  mutate(
    review_needed = dplyr::coalesce(review_needed, FALSE),
    vector_species_needs_review = dplyr::coalesce(vector_species_needs_review, FALSE) | review_needed,
    review_reason = mapply(
      combine_reason_pair,
      review_reason,
      vector_review_reason_extra,
      USE.NAMES = FALSE
    )
  )

analysis_ready <- vector_cleaned %>%
  filter(!is.na(vector_species_analysis)) %>%
  left_join(who_hosts, by = "matched_who_host", suffix = c("", "_who")) %>%
  mutate(
    matched_who_host_tax_id = dplyr::coalesce(
      matched_who_host_tax_id,
      matched_who_host_tax_id_who
    ),
    matched_who_host_class = dplyr::coalesce(
      matched_who_host_class,
      matched_who_host_class_who
    ),
    matched_who_host_order = dplyr::coalesce(
      matched_who_host_order,
      matched_who_host_order_who
    ),
    matched_who_host_family = dplyr::coalesce(
      matched_who_host_family,
      matched_who_host_family_who
    ),
    source_dataset = "MapVEu",
    interaction_type = "blood_meal"
  ) %>%
  transmute(
    mapveu_row_id,
    blood_meal_assay_id,
    sample_id,
    collection_id,
    collection_site_id,
    study_id,
    source_dataset,
    interaction_type,
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
    review_needed,
    review_reason,
    country,
    latitude,
    longitude,
    collection_start_date,
    collection_end_date,
    source_study_name = study_name,
    pubmed_id,
    doi,
    host_organism_raw,
    vector_species_raw,
    vector_species_assay_raw,
    vector_species_sample_raw,
    host_presence,
    host_prevalence_percent,
    collection_host_organism,
    collection_device
  )

if (any(!analysis_ready$matched_who_host %in% who_hosts$matched_who_host)) {
  stop("analysis_ready contains hosts not present in the canonical zoonotic WHO network")
}

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
    source_dataset_examples = collapse_unique(source_dataset),
    country_examples = collapse_unique(country),
    review_reason_examples = collapse_unique(review_reason),
    record_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(matched_who_host, vector_species_analysis)

duplicate_summary_keys <- analysis_summary %>%
  count(matched_who_host, vector_species_analysis, vector_species_needs_review) %>%
  filter(n > 1)

if (nrow(duplicate_summary_keys) > 0) {
  stop("analysis_summary has duplicate host-vector-review keys")
}

write_csv(analysis_ready, analysis_ready_path, na = "")
write_csv(analysis_summary, analysis_summary_path, na = "")

cat("Raw MapVEu rows:", nrow(mapveu_raw), "\n")
cat("Rows retained after WHO-host filtering:", nrow(host_filtered), "\n")
cat("Rows excluded for broad/non-species hosts:", excluded_broad_host_rows, "\n")
cat("Rows excluded for other unmatched hosts:", excluded_other_unmatched_rows, "\n")
cat(
  "Rows retained with vector_species_needs_review = TRUE:",
  sum(analysis_ready$vector_species_needs_review, na.rm = TRUE),
  "\n"
)
cat("Analysis-ready MapVEu rows:", nrow(analysis_ready), "\n")
cat("Analysis summary pairs:", nrow(analysis_summary), "\n")
cat("Wrote analysis-ready MapVEu table to", analysis_ready_path, "\n")
cat("Wrote MapVEu summary table to", analysis_summary_path, "\n")
