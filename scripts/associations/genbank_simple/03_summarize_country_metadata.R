# ------------------------------------------------------------------------------
# 03_summarize_country_metadata.R
# ------------------------------------------------------------------------------
# Purpose: Bind GenBank-simple checkpoint records and summarize pathogen-country
#          and disease-country coverage.
# Inputs : genbank_simple_manifest.csv
#          pathogen_runs/search_logs/*.csv
#          pathogen_runs/country_records/*.csv
# Outputs: genbank_country_records.csv
#          genbank_pathogen_country_summary.csv
#          genbank_disease_country_summary.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, stringr, tibble)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))

output_dir <- here("pathogen_association_data", "WHO", "genbank_simple")
manifest_path <- file.path(output_dir, "genbank_simple_manifest.csv")
search_log_dir <- file.path(output_dir, "pathogen_runs", "search_logs")
country_record_dir <- file.path(output_dir, "pathogen_runs", "country_records")

manifest <- read_csv(manifest_path, show_col_types = FALSE, na = c("", "NA"))

record_paths <- list.files(country_record_dir, pattern = "\\.csv$", full.names = TRUE)
log_paths <- list.files(search_log_dir, pattern = "\\.csv$", full.names = TRUE)

empty_country_records <- tibble(
  target_id = character(),
  Pathogens = character(),
  Disease_name = character(),
  PathogenTaxID = character(),
  query_used = character(),
  source_db = character(),
  accession_version = character(),
  primary_accession = character(),
  definition = character(),
  organism = character(),
  taxonomy = character(),
  sequence_length = character(),
  country_raw = character(),
  geo_loc_name_raw = character(),
  country = character(),
  lat_lon = character(),
  collection_date = character(),
  host = character(),
  isolate = character(),
  strain = character(),
  isolate_source = character(),
  db_xref = character(),
  source_path = character(),
  source_file = character(),
  accession_key = character()
)

final_country_record_cols <- c(
  names(empty_country_records),
  "in_gibb_etal",
  "in_empres_i",
  "network_pathogen_type",
  "network_zoonotic_status",
  "network_canonicalization_status"
)

country_records <- if (length(record_paths) == 0) {
  empty_country_records
} else {
  purrr::map_dfr(record_paths, function(path) {
    read_csv(
      path,
      col_types = cols(.default = col_character()),
      na = c("", "NA")
    ) %>%
      mutate(source_path = path, source_file = basename(path))
  })
} %>%
  bind_rows(empty_country_records) %>%
  select(all_of(names(empty_country_records)))

search_logs <- if (length(log_paths) == 0) {
  tibble()
} else {
  purrr::map_dfr(log_paths, function(path) {
    read_csv(
      path,
      col_types = cols(.default = col_character()),
      na = c("", "NA")
    ) %>%
      mutate(source_path = path, source_file = basename(path))
  })
}

if (nrow(country_records) > 0) {
  country_records <- country_records %>%
    mutate(
      country = clean_text(country),
      accession_key = dplyr::coalesce(accession_version, primary_accession)
    ) %>%
    distinct(target_id, accession_key, .keep_all = TRUE) %>%
    left_join(
      manifest %>%
        select(
          target_id,
          in_gibb_etal,
          in_empres_i,
          network_pathogen_type,
          network_zoonotic_status,
          network_canonicalization_status
        ),
      by = "target_id"
    )
} else {
  country_records <- country_records %>%
    mutate(
      in_gibb_etal = logical(),
      in_empres_i = logical(),
      network_pathogen_type = character(),
      network_zoonotic_status = character(),
      network_canonicalization_status = character()
    )
}

country_records <- country_records %>%
  select(all_of(final_country_record_cols))

pathogen_country_summary <- if (nrow(country_records) == 0) {
  tibble()
} else {
  country_records %>%
    filter(!is.na(country)) %>%
    group_by(target_id, Pathogens, Disease_name, country) %>%
    summarise(
      records_with_country = n(),
      accessions = collapse_unique(accession_key),
      organisms = collapse_unique(organism),
      hosts = collapse_unique(host),
      in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
      in_empres_i = any(in_empres_i, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Pathogens, Disease_name, country)
}

disease_country_summary <- if (nrow(pathogen_country_summary) == 0) {
  tibble()
} else {
  pathogen_country_summary %>%
    group_by(Disease_name, country) %>%
    summarise(
      records_with_country = sum(records_with_country, na.rm = TRUE),
      pathogens = collapse_unique(Pathogens),
      target_ids = collapse_unique(target_id),
      in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
      in_empres_i = any(in_empres_i, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Disease_name, country)
}

write_csv(country_records, file.path(output_dir, "genbank_country_records.csv"))
write_csv(pathogen_country_summary, file.path(output_dir, "genbank_pathogen_country_summary.csv"))
write_csv(disease_country_summary, file.path(output_dir, "genbank_disease_country_summary.csv"))

if (nrow(search_logs) > 0) {
  write_csv(search_logs, file.path(output_dir, "genbank_search_logs.csv"))
}

message("Wrote country records: ", nrow(country_records))
message("Wrote pathogen-country rows: ", nrow(pathogen_country_summary))
message("Wrote disease-country rows: ", nrow(disease_country_summary))
