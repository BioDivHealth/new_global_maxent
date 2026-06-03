# ------------------------------------------------------------------------------
# 5_5c_Vector_Taxonomy_Package_Review.R
# ------------------------------------------------------------------------------
# Purpose: Audit all cleaned vector names using the existing synonym-retrieval
#          helper plus lightweight package lookups, then write both a full
#          official-name check table and a review-focused suggestions table.
#
# Input  : disease_vector_links_taxonomy_cleaned.csv
# Outputs: vector_taxonomy_official_name_checks.csv
#          vector_taxonomy_package_suggestions.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, rgbif, stringr, taxize, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "New_functions", "get_synonyms.R"))
iucn_redlist_key <- Sys.getenv("IUCN_REDLIST_KEY", unset = Sys.getenv("IUCN_API_KEY", unset = ""))
if (nzchar(iucn_redlist_key)) {
  options(iucn_redlist_key = iucn_redlist_key)
}

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
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

safe_name_backbone <- function(name) {
  tryCatch(
    rgbif::name_backbone(name = name),
    error = function(e) tibble::tibble(
      usageKey = NA_character_,
      scientificName = NA_character_,
      canonicalName = NA_character_,
      rank = NA_character_,
      status = NA_character_,
      matchType = NA_character_,
      note = paste("rgbif error:", conditionMessage(e))
    )
  )
}

safe_itis_lookup <- function(name) {
  tsn_res <- tryCatch(
    taxize::get_tsn_(
      sci_com = name,
      searchtype = "scientific",
      messages = FALSE
    ),
    error = function(e) e
  )

  if (inherits(tsn_res, "error")) {
    return(tibble::tibble(
      tsn = NA_character_,
      scientificName = NA_character_,
      nameUsage = NA_character_,
      acceptedName = NA_character_,
      note = paste("taxize error:", conditionMessage(tsn_res))
    ))
  }

  if (!is.list(tsn_res) || length(tsn_res) == 0 || is.null(tsn_res[[1]]) || nrow(tsn_res[[1]]) == 0) {
    return(tibble::tibble(
      tsn = NA_character_,
      scientificName = NA_character_,
      nameUsage = NA_character_,
      acceptedName = NA_character_,
      note = "taxize returned no ITIS matches"
    ))
  }

  first_hit <- as_tibble(tsn_res[[1]]) %>% slice(1)
  accepted_name <- NA_character_

  if ("tsn" %in% names(first_hit) && !is.na(first_hit$tsn[[1]])) {
    accepted_res <- tryCatch(
      taxize::itis_acceptname(first_hit$tsn[[1]], silent = TRUE),
      error = function(e) NULL
    )

    if (!is.null(accepted_res) && is.data.frame(accepted_res) && nrow(accepted_res) > 0) {
      accepted_name <- accepted_res$acceptedName[[1]]
    }
  }

  tibble::tibble(
    tsn = if ("tsn" %in% names(first_hit)) as.character(first_hit$tsn[[1]]) else NA_character_,
    scientificName = if ("scientificName" %in% names(first_hit)) first_hit$scientificName[[1]] else NA_character_,
    nameUsage = if ("nameUsage" %in% names(first_hit)) first_hit$nameUsage[[1]] else NA_character_,
    acceptedName = accepted_name,
    note = NA_character_
  )
}

empty_itis_lookup <- function() {
  tibble::tibble(
    tsn = NA_character_,
    scientificName = NA_character_,
    nameUsage = NA_character_,
    acceptedName = NA_character_,
    note = "ITIS lookup skipped by VECTOR_AUDIT_SKIP_ITIS"
  )
}

safe_retrieve_syns <- function(name, n_times = 3, Gbif = TRUE, Skip_ITIS = FALSE) {
  lookup <- NULL
  lookup_messages <- character()

  tryCatch(
    {
      lookup_messages <- utils::capture.output({
        lookup <- retrieve_syns_new(
          name,
          n_times = n_times,
          Gbif = Gbif,
          Skip_ITIS = Skip_ITIS
        )
      })
    },
    error = function(e) {
      lookup <<- list(error_message = conditionMessage(e))
    }
  )

  if (!is.null(lookup$error_message)) {
    return(tibble::tibble(
      submitted_name = name,
      official_correct_name = NA_character_,
      official_taxon_level = NA_character_,
      official_synonyms = NA_character_,
      official_iucn_name = NA_character_,
      helper_iucn_name = NA_character_,
      helper_itis_name = NA_character_,
      helper_gbif_name = NA_character_,
      helper_note = paste("retrieve_syns_new error:", lookup$error_message)
    ))
  }

  tax_dat <- lookup$TaxDat

  tibble::tibble(
    submitted_name = lookup$Submitted_name,
    official_correct_name = clean_text(lookup$correct_name),
    official_taxon_level = clean_text(lookup$taxon_level),
    official_synonyms = collapse_unique(lookup$Spp_syn),
    official_iucn_name = collapse_unique(lookup$IUCN_spp),
    helper_iucn_name = if (!is.null(tax_dat) && "IUCN_name" %in% names(tax_dat)) collapse_unique(tax_dat$IUCN_name) else NA_character_,
    helper_itis_name = if (!is.null(tax_dat) && "ITIS_name" %in% names(tax_dat)) collapse_unique(tax_dat$ITIS_name) else NA_character_,
    helper_gbif_name = if (!is.null(tax_dat) && "GBIF_name" %in% names(tax_dat)) collapse_unique(tax_dat$GBIF_name) else NA_character_,
    helper_note = collapse_unique(lookup_messages)
  )
}

input_path <- vector_screening_evidence_path("disease_vector_links_taxonomy_cleaned.csv")
full_output_path <- file.path(
  vector_screening_taxonomy_review_dir,
  "vector_taxonomy_official_name_checks.csv"
)
review_output_path <- file.path(
  vector_screening_taxonomy_review_dir,
  "vector_taxonomy_package_suggestions.csv"
)
dir.create(vector_screening_taxonomy_review_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_cleaned <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_cols <- c(
  "vector_species_original",
  "vector_species_clean",
  "vector_species_taxonomy_cleaned",
  "review_needed",
  "review_note",
  "disease_name"
)
missing_cols <- setdiff(required_cols, names(taxonomy_cleaned))
if (length(missing_cols) > 0) {
  stop(
    "disease_vector_links_taxonomy_cleaned.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

all_names <- taxonomy_cleaned %>%
  filter(!is.na(vector_species_taxonomy_cleaned)) %>%
  group_by(vector_species_taxonomy_cleaned) %>%
  summarise(
    vector_species_original_examples = collapse_unique(vector_species_original),
    vector_species_clean_examples = collapse_unique(vector_species_clean),
    disease_examples = collapse_unique(disease_name),
    review_needed_any = any(review_needed %in% TRUE, na.rm = TRUE),
    review_notes = collapse_unique(review_note),
    .groups = "drop"
  ) %>%
  arrange(vector_species_taxonomy_cleaned)

if (file.exists(full_output_path)) {
  existing_output <- read_csv(
    full_output_path,
    show_col_types = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(across(where(is.character), clean_text))
} else {
  existing_output <- tibble::tibble()
}

completed_names <- if ("vector_species_taxonomy_cleaned" %in% names(existing_output)) {
  unique(stats::na.omit(existing_output$vector_species_taxonomy_cleaned))
} else {
  character()
}

names_to_process <- all_names %>%
  filter(!(vector_species_taxonomy_cleaned %in% completed_names))

limit_env <- suppressWarnings(as.integer(Sys.getenv("VECTOR_NAME_LIMIT", "")))
if (!is.na(limit_env) && limit_env > 0) {
  names_to_process <- names_to_process %>% slice_head(n = limit_env)
}

skip_itis <- tolower(Sys.getenv("VECTOR_AUDIT_SKIP_ITIS", "")) %in% c(
  "1", "true", "yes", "y", "on"
)

results <- existing_output

cat("VECTOR_AUDIT_SKIP_ITIS:", skip_itis, "\n")

for (i in seq_len(nrow(names_to_process))) {
  target_name <- names_to_process$vector_species_taxonomy_cleaned[[i]]
  cat("Checking vector name", i, "of", nrow(names_to_process), ":", target_name, "\n")

  gbif_res <- safe_name_backbone(target_name)
  itis_res <- if (skip_itis) empty_itis_lookup() else safe_itis_lookup(target_name)
  helper_res <- safe_retrieve_syns(
    target_name,
    n_times = 3,
    Gbif = TRUE,
    Skip_ITIS = skip_itis
  )

  row_result <- tibble::tibble(
    vector_species_taxonomy_cleaned = target_name,
    vector_species_original_examples = names_to_process$vector_species_original_examples[[i]],
    vector_species_clean_examples = names_to_process$vector_species_clean_examples[[i]],
    disease_examples = names_to_process$disease_examples[[i]],
    review_needed_any = names_to_process$review_needed_any[[i]],
    review_notes = names_to_process$review_notes[[i]],
    official_correct_name = helper_res$official_correct_name[[1]],
    official_taxon_level = helper_res$official_taxon_level[[1]],
    official_synonyms = helper_res$official_synonyms[[1]],
    official_iucn_name = helper_res$official_iucn_name[[1]],
    official_name_status = dplyr::case_when(
      is.na(helper_res$official_correct_name[[1]]) ~ "no_official_name_returned",
      stringr::str_to_lower(helper_res$official_correct_name[[1]]) == stringr::str_to_lower(target_name) ~ "exact_match",
      TRUE ~ "suggested_change"
    ),
    helper_iucn_name = helper_res$helper_iucn_name[[1]],
    helper_itis_name = helper_res$helper_itis_name[[1]],
    helper_gbif_name = helper_res$helper_gbif_name[[1]],
    helper_note = helper_res$helper_note[[1]],
    gbif_usage_key = if ("usageKey" %in% names(gbif_res)) as.character(gbif_res$usageKey[[1]]) else NA_character_,
    gbif_scientific_name = if ("scientificName" %in% names(gbif_res)) gbif_res$scientificName[[1]] else NA_character_,
    gbif_canonical_name = if ("canonicalName" %in% names(gbif_res)) gbif_res$canonicalName[[1]] else NA_character_,
    gbif_rank = if ("rank" %in% names(gbif_res)) gbif_res$rank[[1]] else NA_character_,
    gbif_status = if ("status" %in% names(gbif_res)) gbif_res$status[[1]] else NA_character_,
    gbif_match_type = if ("matchType" %in% names(gbif_res)) gbif_res$matchType[[1]] else NA_character_,
    gbif_note = if ("note" %in% names(gbif_res)) gbif_res$note[[1]] else NA_character_,
    itis_tsn = itis_res$tsn[[1]],
    itis_scientific_name = itis_res$scientificName[[1]],
    itis_name_usage = itis_res$nameUsage[[1]],
    itis_accepted_name = itis_res$acceptedName[[1]],
    itis_note = itis_res$note[[1]]
  )

  results <- bind_rows(results, row_result)
  write_csv(results, full_output_path, na = "")
}

results <- results %>%
  arrange(vector_species_taxonomy_cleaned)

write_csv(results, full_output_path, na = "")

review_output <- results %>%
  filter(
    review_needed_any %in% TRUE |
      official_name_status != "exact_match" |
      gbif_match_type %in% c("VARIANT", "HIGHERRANK") |
      !is.na(gbif_note) |
      !is.na(itis_note)
  )

write_csv(review_output, review_output_path, na = "")

cat("Distinct cleaned vector names audited:", nrow(results), "\n")
cat(
  "Exact official matches:",
  sum(results$official_name_status == "exact_match", na.rm = TRUE),
  "\n"
)
cat(
  "Suggested changes from helper:",
  sum(results$official_name_status == "suggested_change", na.rm = TRUE),
  "\n"
)
cat(
  "Names with no official name returned:",
  sum(results$official_name_status == "no_official_name_returned", na.rm = TRUE),
  "\n"
)
cat("Review rows written:", nrow(review_output), "\n")
cat("Wrote full audit to", full_output_path, "\n")
cat("Wrote review suggestions to", review_output_path, "\n")
