# ------------------------------------------------------------------------------
# Build a first-pass MapVEu vector-host table from blood meal assay downloads
# ------------------------------------------------------------------------------

library(pacman)
p_load(data.table, dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))

# Normalize free-text fields and convert blanks to NA.
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

# MapVEu stores many values as JSON-like arrays such as ["Aedes vexans"].
parse_list_field <- function(x) {
  x <- clean_text(x)
  x <- str_replace_all(x, '^\\["', "")
  x <- str_replace_all(x, '"\\]$', "")
  x <- str_replace_all(x, '"\\s*,\\s*"', " | ")
  x <- str_replace_all(x, "^\\[|\\]$", "")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_unique_values <- function(x) {
  x <- clean_text(x)
  x <- unique(stats::na.omit(x))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(sort(x), collapse = " | ")
}

# The study export has a few unquoted line breaks inside titles, so we repair
# lines until they reach the expected column count before parsing.
read_irregular_tsv <- function(path, expected_cols) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  if (length(lines) == 0) {
    stop("File is empty: ", path)
  }

  repaired_lines <- lines[1]
  buffer <- character()

  for (line in lines[-1]) {
    candidate <- paste(c(buffer, line), collapse = " ")
    field_count <- length(strsplit(candidate, "\t", fixed = TRUE)[[1]])

    if (field_count < expected_cols) {
      buffer <- c(buffer, line)
    } else {
      repaired_lines <- c(repaired_lines, candidate)
      buffer <- character()
    }
  }

  if (length(buffer) > 0) {
    repaired_lines <- c(repaired_lines, paste(buffer, collapse = " "))
  }

  split_lines <- strsplit(repaired_lines, "\t", fixed = TRUE)
  field_counts <- vapply(split_lines, length, integer(1))

  if (any(field_counts != expected_cols)) {
    stop("Irregular TSV repair failed for: ", path)
  }

  header <- split_lines[[1]]
  values <- split_lines[-1]
  value_matrix <- do.call(rbind, values)
  colnames(value_matrix) <- header

  as.data.frame(value_matrix, stringsAsFactors = FALSE, check.names = FALSE)
}

output_dir <- mapveu_outputs_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

blood_meal_path <- file.path(mapveu_dir, "VBP_MEGA_Blood meal assay_subsettedData.txt")
sample_path <- file.path(mapveu_dir, "VBP_MEGA_Sample_subsettedData.txt")
species_assay_path <- file.path(mapveu_dir, "VBP_MEGA_Species identification assay_subsettedData.txt")
collection_path <- file.path(mapveu_dir, "VBP_MEGA_Collection_subsettedData.txt")
collection_site_path <- file.path(mapveu_dir, "VBP_MEGA_Collection site_subsettedData.txt")
study_path <- file.path(mapveu_dir, "VBP_MEGA_Study_subsettedData.txt")

output_links_path <- file.path(output_dir, "mapveu_vector_host_links_raw.csv")
output_review_path <- file.path(output_dir, "mapveu_vector_host_links_review.csv")

# Read only the columns needed for first-pass joins and review.
blood_meal <- fread(
  blood_meal_path,
  sep = "\t",
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = c(
    "Blood_meal_assay_ID",
    "Sample_ID",
    "Collection_ID",
    "Collection_site_ID",
    "Study_ID",
    "protocol [OBI_0000272]",
    "blood meal host organism [OBI_0002995]",
    "blood meal host presence [OBI_0002994]",
    "blood meal host prevalence (percent) [OBI_0002993]",
    "number of input specimens to quantitative assay [POPBIO_8000018]"
  )
) %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    blood_meal_protocol = parse_list_field(`protocol [OBI_0000272]`),
    host_organism_raw = parse_list_field(`blood meal host organism [OBI_0002995]`),
    host_presence = clean_text(`blood meal host presence [OBI_0002994]`),
    host_prevalence_percent = `blood meal host prevalence (percent) [OBI_0002993]`,
    n_input_specimens = `number of input specimens to quantitative assay [POPBIO_8000018]`
  ) %>%
  select(
    blood_meal_assay_id = Blood_meal_assay_ID,
    sample_id = Sample_ID,
    collection_id = Collection_ID,
    collection_site_id = Collection_site_ID,
    study_id = Study_ID,
    blood_meal_protocol,
    host_organism_raw,
    host_presence,
    host_prevalence_percent,
    n_input_specimens
  )

blood_sample_ids <- unique(blood_meal$sample_id)
blood_collection_ids <- unique(blood_meal$collection_id)
blood_collection_site_ids <- unique(blood_meal$collection_site_id)
blood_study_ids <- unique(blood_meal$study_id)

sample_tbl <- fread(
  sample_path,
  sep = "\t",
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = c(
    "Sample_ID",
    "Collection_ID",
    "Collection_site_ID",
    "Study_ID",
    "Sample type [EUPATH_0000611]",
    "biological sex [PATO_0000047]",
    "female insect feeding status [EUPATH_0043227]",
    "life cycle stage [UBERON_0000105]",
    "species [OBI_0001909]",
    "species qualifier [IAO_0000078]",
    "specimen count [EUPATH_0043155]"
  )
) %>%
  .[`Sample_ID` %in% blood_sample_ids] %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    sample_type = clean_text(`Sample type [EUPATH_0000611]`),
    biological_sex = clean_text(`biological sex [PATO_0000047]`),
    female_insect_feeding_status = clean_text(`female insect feeding status [EUPATH_0043227]`),
    life_cycle_stage = clean_text(`life cycle stage [UBERON_0000105]`),
    vector_species_sample_raw = parse_list_field(`species [OBI_0001909]`),
    vector_species_sample_qualifier = clean_text(`species qualifier [IAO_0000078]`),
    specimen_count = `specimen count [EUPATH_0043155]`
  ) %>%
  select(
    sample_id = Sample_ID,
    collection_id = Collection_ID,
    collection_site_id = Collection_site_ID,
    study_id = Study_ID,
    sample_type,
    biological_sex,
    female_insect_feeding_status,
    life_cycle_stage,
    vector_species_sample_raw,
    vector_species_sample_qualifier,
    specimen_count
  )

species_assay <- fread(
  species_assay_path,
  sep = "\t",
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = c(
    "Species_identification_assay_ID",
    "Sample_ID",
    "Collection_ID",
    "Collection_site_ID",
    "Study_ID",
    "protocol [OBI_0000272]",
    "organism identification datum [EUPATH_0043194]"
  )
) %>%
  .[`Sample_ID` %in% blood_sample_ids] %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    species_assay_protocol = parse_list_field(`protocol [OBI_0000272]`),
    vector_species_assay_raw = parse_list_field(`organism identification datum [EUPATH_0043194]`)
  ) %>%
  group_by(
    sample_id = Sample_ID,
    collection_id = Collection_ID,
    collection_site_id = Collection_site_ID,
    study_id = Study_ID
  ) %>%
  summarise(
    species_assay_id = collapse_unique_values(Species_identification_assay_ID),
    species_assay_protocol = collapse_unique_values(species_assay_protocol),
    vector_species_assay_raw = collapse_unique_values(vector_species_assay_raw),
    species_assay_rows = n(),
    species_assay_unique_species = n_distinct(stats::na.omit(vector_species_assay_raw)),
    .groups = "drop"
  )

collection_tbl <- fread(
  collection_path,
  sep = "\t",
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = c(
    "Collection_ID",
    "Collection_site_ID",
    "Study_ID",
    "protocol [OBI_0000272]",
    "host organism [EUPATH_0000591]",
    "device [OBI_0000968]",
    "specimen collection start date [EUPATH_0043256]",
    "specimen collection end date [EUPATH_0043257]"
  )
) %>%
  .[`Collection_ID` %in% blood_collection_ids] %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    collection_protocol = parse_list_field(`protocol [OBI_0000272]`),
    collection_host_organism = parse_list_field(`host organism [EUPATH_0000591]`),
    collection_device = parse_list_field(`device [OBI_0000968]`),
    collection_start_date = clean_text(`specimen collection start date [EUPATH_0043256]`),
    collection_end_date = clean_text(`specimen collection end date [EUPATH_0043257]`)
  ) %>%
  select(
    collection_id = Collection_ID,
    collection_site_id = Collection_site_ID,
    study_id = Study_ID,
    collection_protocol,
    collection_host_organism,
    collection_device,
    collection_start_date,
    collection_end_date
  )

collection_site_tbl <- fread(
  collection_site_path,
  sep = "\t",
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = c(
    "Collection_site_ID",
    "Study_ID",
    "continent [GAZ_00000013]",
    "Latitude [OBI_0001620]",
    "Longitude [OBI_0001621]",
    "country [OBI_0001627]",
    "Administrative region, level 1 [ENVO_00000005]",
    "Administrative region, level 2 [ENVO_00000006]",
    "town [POPBIO_8000015]"
  )
) %>%
  .[`Collection_site_ID` %in% blood_collection_site_ids] %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  transmute(
    collection_site_id = Collection_site_ID,
    study_id = Study_ID,
    continent = clean_text(`continent [GAZ_00000013]`),
    latitude = `Latitude [OBI_0001620]`,
    longitude = `Longitude [OBI_0001621]`,
    country = clean_text(`country [OBI_0001627]`),
    admin_region_level_1 = clean_text(`Administrative region, level 1 [ENVO_00000005]`),
    admin_region_level_2 = clean_text(`Administrative region, level 2 [ENVO_00000006]`),
    town = clean_text(`town [POPBIO_8000015]`)
  )

study_tbl <- tryCatch(
  {
    read_irregular_tsv(study_path, expected_cols = 8) %>%
      as_tibble() %>%
      select(
        "Study_ID",
        "Study name [OBI_0001622]",
        "PubMed ID [OBI_0001617]",
        "DOI [OBI_0002110]"
      ) %>%
      filter(Study_ID %in% blood_study_ids) %>%
      mutate(across(where(is.character), clean_text)) %>%
      transmute(
        study_id = Study_ID,
        study_name = clean_text(`Study name [OBI_0001622]`),
        pubmed_id = clean_text(`PubMed ID [OBI_0001617]`),
        doi = clean_text(`DOI [OBI_0002110]`)
      )
  },
  error = function(e) {
    message("Study metadata could not be parsed cleanly; continuing without study_name/pubmed_id/doi.")
    tibble(
      study_id = character(),
      study_name = character(),
      pubmed_id = character(),
      doi = character()
    )
  }
)

# Join blood-meal rows to vector identity and supporting metadata.
vector_host_links_raw <- blood_meal %>%
  left_join(sample_tbl, by = c("sample_id", "collection_id", "collection_site_id", "study_id")) %>%
  left_join(species_assay, by = c("sample_id", "collection_id", "collection_site_id", "study_id")) %>%
  left_join(collection_tbl, by = c("collection_id", "collection_site_id", "study_id")) %>%
  left_join(collection_site_tbl, by = c("collection_site_id", "study_id")) %>%
  left_join(study_tbl, by = "study_id") %>%
  mutate(
    vector_species_raw = coalesce(vector_species_assay_raw, vector_species_sample_raw),
    vector_species_source = case_when(
      !is.na(vector_species_assay_raw) &
        !is.na(vector_species_sample_raw) &
        vector_species_assay_raw == vector_species_sample_raw ~ "species_assay_and_sample_match",
      !is.na(vector_species_assay_raw) &
        !is.na(vector_species_sample_raw) &
        vector_species_assay_raw != vector_species_sample_raw ~ "species_assay_preferred_sample_differs",
      !is.na(vector_species_assay_raw) ~ "species_assay_only",
      !is.na(vector_species_sample_raw) ~ "sample_only",
      TRUE ~ "missing"
    ),
    review_needed = case_when(
      is.na(vector_species_raw) ~ TRUE,
      species_assay_unique_species > 1 ~ TRUE,
      vector_species_source == "species_assay_preferred_sample_differs" ~ TRUE,
      TRUE ~ FALSE
    ),
    review_reason = case_when(
      is.na(vector_species_raw) ~ "missing_vector_species",
      species_assay_unique_species > 1 ~ "multiple_species_assay_values",
      vector_species_source == "species_assay_preferred_sample_differs" ~ "sample_species_differs_from_species_assay",
      TRUE ~ NA_character_
    ),
    source_dataset = "MapVEu"
  ) %>%
  select(
    source_dataset,
    blood_meal_assay_id,
    sample_id,
    collection_id,
    collection_site_id,
    study_id,
    study_name,
    pubmed_id,
    doi,
    vector_species_raw,
    vector_species_assay_raw,
    vector_species_sample_raw,
    vector_species_sample_qualifier,
    vector_species_source,
    species_assay_rows,
    species_assay_unique_species,
    host_organism_raw,
    host_presence,
    host_prevalence_percent,
    n_input_specimens,
    sample_type,
    biological_sex,
    female_insect_feeding_status,
    life_cycle_stage,
    specimen_count,
    collection_host_organism,
    collection_device,
    collection_protocol,
    collection_start_date,
    collection_end_date,
    continent,
    country,
    admin_region_level_1,
    admin_region_level_2,
    town,
    latitude,
    longitude,
    review_needed,
    review_reason
  )

review_tbl <- vector_host_links_raw %>%
  filter(review_needed) %>%
  select(
    blood_meal_assay_id,
    sample_id,
    study_id,
    vector_species_raw,
    vector_species_assay_raw,
    vector_species_sample_raw,
    host_organism_raw,
    review_reason
  )

write_csv(vector_host_links_raw, output_links_path, na = "")
write_csv(review_tbl, output_review_path, na = "")

cat("Blood-meal rows:", nrow(blood_meal), "\n")
cat("Rows written:", nrow(vector_host_links_raw), "\n")
cat("Rows flagged for review:", nrow(review_tbl), "\n")
cat("Unique vector species labels:", n_distinct(vector_host_links_raw$vector_species_raw, na.rm = TRUE), "\n")
cat("Unique host labels:", n_distinct(vector_host_links_raw$host_organism_raw, na.rm = TRUE), "\n")
cat("Wrote raw links to", output_links_path, "\n")
cat("Wrote review rows to", output_review_path, "\n")
