# ------------------------------------------------------------------------------
# 5_5b_Vector_Name_Cleanup.R
# ------------------------------------------------------------------------------
# Purpose: Add a conservative taxonomy-cleanup layer on top of the canonical
#          disease-vector table produced by 5_5 while preserving the original
#          names and routing unresolved cases to a review file.
#
# Input  : disease_vector_links.csv
# Outputs: disease_vector_links_taxonomy_cleaned.csv
#          vector_taxonomy_review_needed.csv
#          vector_taxonomy_manual_map.csv (seeded if absent)
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

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

classify_vector_taxon_rank <- function(x) {
  x <- clean_text(x)
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    stringr::str_detect(x, "\\bspp\\b") ~ "genus_only",
    stringr::str_detect(x, "complex") ~ "complex",
    stringr::str_detect(x, "group") ~ "group",
    stringr::str_detect(x, "sensu lato|\\bs l\\b") ~ "sensu_lato",
    stringr::str_detect(x, " x ") ~ "hybrid",
    stringr::str_detect(x, "biotype|form") ~ "infraspecific",
    stringr::str_detect(x, "\\bs s\\b") ~ "infraspecific",
    stringr::str_count(x, "\\S+") == 1 ~ "genus_only",
    stringr::str_count(x, "\\S+") == 2 ~ "species",
    TRUE ~ "infraspecific"
  )
}

apply_rule_cleanup <- function(x) {
  x <- clean_text(x)
  method <- rep("no_change", length(x))
  cleaned <- x

  repeated_genus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^([a-z]+) \\1\\b")
  if (any(repeated_genus)) {
    cleaned[repeated_genus] <- stringr::str_replace(
      cleaned[repeated_genus],
      "^([a-z]+) \\1\\s+",
      "\\1 "
    )
    method[repeated_genus] <- "rule_repeated_genus"
  }

  aedes_subgenus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^aedes (ochlerotatus|neomelaniconion)\\b")
  if (any(aedes_subgenus)) {
    cleaned[aedes_subgenus] <- stringr::str_replace(
      cleaned[aedes_subgenus],
      "^aedes (ochlerotatus|neomelaniconion)\\s+",
      "aedes "
    )
    method[aedes_subgenus] <- "rule_drop_subgenus_token"
  }

  culex_subgenus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^culex melanoconion\\b")
  if (any(culex_subgenus)) {
    cleaned[culex_subgenus] <- stringr::str_replace(
      cleaned[culex_subgenus],
      "^culex melanoconion\\s+",
      "culex "
    )
    method[culex_subgenus] <- "rule_drop_subgenus_token"
  }

  tibble::tibble(
    vector_species_rule_cleaned = cleaned,
    taxonomy_cleanup_method = method
  )
}

seed_manual_map <- function(path) {
  seeded_map <- tibble::tribble(
    ~source_name, ~canonical_name, ~vector_taxon_rank, ~cleanup_method, ~notes,
    "stegomyia albopicta", "aedes albopictus", "species", "manual_map", "Old genus usage normalized to the Aedes spelling used elsewhere in the table.",
    "aedes flavicolis", "aedes flavicollis", "species", "manual_map", "Likely spelling variant aligned to the canonical species spelling already present in the table.",
    "aedes luteocephalis", "aedes luteocephalus", "species", "manual_map", "Likely spelling variant aligned to the canonical species spelling already present in the table.",
    "haemagogus spegazzinni", "haemagogus spegazzinii", "species", "manual_map", "Likely spelling variant aligned to the canonical species spelling already present in the table."
  )

  write_csv(seeded_map, path, na = "")
}

input_path <- vector_screening_staged_path("disease_vector_links.csv")
output_path <- file.path(
  vector_screening_evidence_dir,
  "disease_vector_links_taxonomy_cleaned.csv"
)
review_path <- file.path(
  vector_screening_taxonomy_review_dir,
  "vector_taxonomy_review_needed.csv"
)
manual_map_path <- vector_screening_taxonomy_manual_path("vector_taxonomy_manual_map.csv")
dir.create(vector_screening_evidence_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_screening_taxonomy_review_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_screening_taxonomy_manual_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manual_map_path)) {
  seed_manual_map(manual_map_path)
}

disease_vector_links <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

manual_map <- read_csv(
  manual_map_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_input_cols <- c("disease_name", "vector_species", "vector_species_clean")
missing_input_cols <- setdiff(required_input_cols, names(disease_vector_links))
if (length(missing_input_cols) > 0) {
  stop(
    "disease_vector_links.csv is missing required columns: ",
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
    "vector_taxonomy_manual_map.csv is missing required columns: ",
    paste(missing_map_cols, collapse = ", ")
  )
}

manual_map <- manual_map %>%
  mutate(source_name = stringr::str_to_lower(source_name)) %>%
  distinct(source_name, .keep_all = TRUE)

rule_cleaned <- apply_rule_cleanup(disease_vector_links$vector_species_clean)

suspect_names <- c(
  "hyalomma onatoli",
  "verrallina lineato"
)

taxonomy_cleaned <- disease_vector_links %>%
  bind_cols(rule_cleaned) %>%
  mutate(
    vector_species_map_key = stringr::str_to_lower(vector_species_rule_cleaned)
  ) %>%
  left_join(
    manual_map %>%
      transmute(
        source_name = stringr::str_to_lower(source_name),
        manual_canonical_name = canonical_name,
        manual_vector_taxon_rank = vector_taxon_rank,
        manual_cleanup_method = cleanup_method,
        manual_notes = notes
      ),
    by = c("vector_species_map_key" = "source_name")
  ) %>%
  mutate(
    vector_species_original = vector_species,
    vector_species_taxonomy_cleaned = dplyr::coalesce(
      manual_canonical_name,
      vector_species_rule_cleaned
    ),
    vector_taxon_rank = dplyr::coalesce(
      manual_vector_taxon_rank,
      classify_vector_taxon_rank(vector_species_taxonomy_cleaned)
    ),
    taxonomy_cleanup_method = dplyr::case_when(
      !is.na(manual_canonical_name) ~ manual_cleanup_method,
      TRUE ~ taxonomy_cleanup_method
    ),
    review_needed = dplyr::case_when(
      vector_species_clean %in% suspect_names ~ TRUE,
      stringr::str_detect(vector_species_taxonomy_cleaned, "^stegomyia\\b") ~ TRUE,
      TRUE ~ FALSE
    ),
    review_note = dplyr::case_when(
      vector_species_clean %in% suspect_names ~ "Likely spelling or taxonomy issue; manual decision needed.",
      stringr::str_detect(vector_species_taxonomy_cleaned, "^stegomyia\\b") ~ "Old genus retained after rules; add a manual mapping decision.",
      !is.na(manual_notes) ~ manual_notes,
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    everything(),
    -vector_species_map_key,
    -manual_canonical_name,
    -manual_vector_taxon_rank,
    -manual_cleanup_method,
    -manual_notes
  )

review_table <- taxonomy_cleaned %>%
  filter(review_needed) %>%
  transmute(
    disease_name,
    vector_species_original,
    vector_species_clean,
    vector_species_taxonomy_cleaned,
    vector_taxon_rank,
    taxonomy_cleanup_method,
    review_note
  ) %>%
  group_by(
    vector_species_original,
    vector_species_clean,
    vector_species_taxonomy_cleaned,
    vector_taxon_rank,
    taxonomy_cleanup_method,
    review_note
  ) %>%
  summarise(
    disease_examples = collapse_unique(disease_name),
    .groups = "drop"
  ) %>%
  arrange(vector_species_clean)

if (nrow(taxonomy_cleaned) != nrow(disease_vector_links)) {
  stop("Row count changed during taxonomy cleanup")
}

write_csv(taxonomy_cleaned, output_path, na = "")
write_csv(review_table, review_path, na = "")

cat("Input disease-vector rows:", nrow(disease_vector_links), "\n")
cat("Rows in cleaned output:", nrow(taxonomy_cleaned), "\n")
cat("Manual map entries:", nrow(manual_map), "\n")
cat(
  "Rows changed by cleanup:",
  sum(
    taxonomy_cleaned$vector_species_clean != taxonomy_cleaned$vector_species_taxonomy_cleaned,
    na.rm = TRUE
  ),
  "\n"
)
cat(
  "Rows flagged for review:",
  sum(taxonomy_cleaned$review_needed, na.rm = TRUE),
  "\n"
)
cat("Wrote cleaned table to", output_path, "\n")
cat("Wrote review table to", review_path, "\n")
cat("Manual map path:", manual_map_path, "\n")

# Check with taxonomic packages----------------------
sp = taxonomy_cleaned$vector_species[1]

p_load(here, rgbif, taxize, raster, dismo, 
       doParallel, rJava, XML, rgbif, Hmisc, readr, 
       stringr, purrr, dplyr, tidyr, magrittr, tidyverse)

source(here("scripts", "New_functions", "get_synonyms.R"))
iucn_redlist_key <- Sys.getenv("IUCN_REDLIST_KEY", unset = Sys.getenv("IUCN_API_KEY", unset = ""))
if (nzchar(iucn_redlist_key)) {
  options(iucn_redlist_key = iucn_redlist_key)
}

# Helper function from 0_SpList.R -----------------------------------------
collapse_vals <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x)])
  paste(x, collapse = sep)
}
retrieve_syns_new(sp,  
                  n_times=10,
                  Gbif=TRUE)
