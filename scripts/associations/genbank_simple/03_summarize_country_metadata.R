# ------------------------------------------------------------------------------|
#      03_summarize_country_metadata.R ----------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Bind GenBank-simple checkpoint records and summarize pathogen-country
#          and disease-country coverage.
# Inputs : genbank_simple_manifest.csv
#          pathogen_runs/search_logs/*.csv
#          pathogen_runs/country_records/*.csv
#          Optional `GENBANK_SIMPLE_SUMMARY_KIND=readiness_combined` binds the
#          old 19-target run and expanded readiness run against the readiness
#          manifest, writing readiness outputs under the top-level, `qa/`, and
#          `intermediate/` GenBank-simple folders.
# Outputs: genbank_country_records.csv
#          genbank_pathogen_country_summary.csv
#          genbank_disease_country_summary.csv
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, purrr, readr, stringr, tibble)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))

# ------------------------------------------------------------------------------|
#      Resolve run mode and paths ---------------------------------------------
# ------------------------------------------------------------------------------|
output_dir <- here("pathogen_association_data", "WHO", "genbank_simple")
summary_kind <- Sys.getenv("GENBANK_SIMPLE_SUMMARY_KIND", unset = "standard") %>%
  clean_text() %>%
  stringr::str_to_lower()

summary_kind <- case_when(
  summary_kind %in% c("standard", "simple", "current") ~ "standard",
  summary_kind %in% c("readiness", "readiness_combined", "expanded_readiness") ~
    "readiness_combined",
  TRUE ~ NA_character_
)

if (is.na(summary_kind)) {
  stop(
    "GENBANK_SIMPLE_SUMMARY_KIND must be `standard` or `readiness_combined`.",
    call. = FALSE
  )
}

manifest_path <- file.path(
  output_dir,
  if_else(
    summary_kind == "readiness_combined",
    "genbank_simple_readiness_manifest.csv",
    "genbank_simple_manifest.csv"
  )
)

standard_search_log_dir <- file.path(output_dir, "pathogen_runs", "search_logs")
standard_country_record_dir <- file.path(output_dir, "pathogen_runs", "country_records")
readiness_search_log_dir <- file.path(output_dir, "pathogen_runs_readiness", "search_logs")
readiness_country_record_dir <- file.path(output_dir, "pathogen_runs_readiness", "country_records")

search_log_dirs <- if (summary_kind == "readiness_combined") {
  c(standard_search_log_dir, readiness_search_log_dir)
} else {
  standard_search_log_dir
}

country_record_dirs <- if (summary_kind == "readiness_combined") {
  c(standard_country_record_dir, readiness_country_record_dir)
} else {
  standard_country_record_dir
}

output_prefix <- if_else(summary_kind == "readiness_combined", "genbank_readiness_", "genbank_")

manifest <- read_csv(manifest_path, show_col_types = FALSE, na = c("", "NA"))

# ------------------------------------------------------------------------------|
#      Discover checkpoint files ----------------------------------------------
# ------------------------------------------------------------------------------|
record_paths <- unlist(lapply(
  country_record_dirs,
  list.files,
  pattern = "\\.csv$",
  full.names = TRUE
))
log_paths <- unlist(lapply(
  search_log_dirs,
  list.files,
  pattern = "\\.csv$",
  full.names = TRUE
))

# ------------------------------------------------------------------------------|
#      Define combined country-record schema ----------------------------------
# ------------------------------------------------------------------------------|
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
  "network_canonicalization_status",
  "analysis_unit_ids",
  "readiness_row_count",
  "manifest_status",
  "manifest_status_reason",
  "qa_flags",
  "genbank_run_source"
)

# ------------------------------------------------------------------------------|
#      Bind per-target records and logs ---------------------------------------
# ------------------------------------------------------------------------------|
country_records <- if (length(record_paths) == 0) {
  empty_country_records
} else {
  purrr::map_dfr(record_paths, function(path) {
    read_csv(
      path,
      col_types = cols(.default = col_character()),
      na = c("", "NA")
    ) %>%
      mutate(
        source_path = path,
        source_file = basename(path),
        genbank_run_source = if_else(
          str_detect(path, fixed("pathogen_runs_readiness")),
          "readiness_run",
          "standard_19_target_run"
        )
      )
  })
} %>%
  bind_rows(empty_country_records) %>%
  mutate(genbank_run_source = clean_text(genbank_run_source)) %>%
  select(any_of(c(names(empty_country_records), "genbank_run_source")))

search_logs <- if (length(log_paths) == 0) {
  tibble()
} else {
  purrr::map_dfr(log_paths, function(path) {
    read_csv(
      path,
      col_types = cols(.default = col_character()),
      na = c("", "NA")
    ) %>%
      mutate(
        source_path = path,
        source_file = basename(path),
        genbank_run_source = if_else(
          str_detect(path, fixed("pathogen_runs_readiness")),
          "readiness_run",
          "standard_19_target_run"
        )
      )
  })
}

# ------------------------------------------------------------------------------|
#      Build manifest lookup for standard and readiness runs ------------------
# ------------------------------------------------------------------------------|
manifest_join <- if (summary_kind == "readiness_combined") {
  manifest %>%
    transmute(
      target_id,
      readiness_target_id = target_id,
      current_target_id = clean_text(current_target_id),
      readiness_Pathogens = query_pathogen_label,
      readiness_Disease_name = readiness_disease_names,
      readiness_PathogenTaxID = pathogen_taxid,
      source_db,
      query_used,
      analysis_unit_ids,
      readiness_row_count,
      manifest_status,
      manifest_status_reason,
      qa_flags,
      in_gibb_etal = NA,
      in_empres_i = NA,
      network_pathogen_type = NA_character_,
      network_zoonotic_status = NA_character_,
      network_canonicalization_status = NA_character_
    ) %>%
    tidyr::pivot_longer(
      cols = c(target_id, current_target_id),
      names_to = "target_id_source",
      values_to = "record_target_id"
    ) %>%
    filter(!is.na(record_target_id)) %>%
    distinct(record_target_id, .keep_all = TRUE)
} else {
  manifest %>%
    transmute(
      record_target_id = target_id,
      readiness_target_id = target_id,
      readiness_Pathogens = Pathogens,
      readiness_Disease_name = Disease_name,
      readiness_PathogenTaxID = PathogenTaxID,
      in_gibb_etal,
      in_empres_i,
      network_pathogen_type,
      network_zoonotic_status,
      network_canonicalization_status,
      analysis_unit_ids = NA_character_,
      readiness_row_count = NA,
      manifest_status = NA_character_,
      manifest_status_reason = NA_character_,
      qa_flags = NA_character_
  )
}

# ------------------------------------------------------------------------------|
#      Attach manifest metadata to country records ----------------------------
# ------------------------------------------------------------------------------|
if (nrow(country_records) > 0) {
  country_records <- country_records %>%
    mutate(
      country = clean_text(country),
      accession_key = dplyr::coalesce(accession_version, primary_accession)
    ) %>%
    distinct(target_id, accession_key, .keep_all = TRUE) %>%
    left_join(manifest_join, by = c("target_id" = "record_target_id")) %>%
    mutate(
      target_id = dplyr::coalesce(readiness_target_id, target_id),
      Pathogens = dplyr::coalesce(readiness_Pathogens, Pathogens),
      Disease_name = dplyr::coalesce(readiness_Disease_name, Disease_name),
      PathogenTaxID = dplyr::coalesce(readiness_PathogenTaxID, PathogenTaxID)
    )
} else {
  country_records <- country_records %>%
    mutate(
      in_gibb_etal = logical(),
      in_empres_i = logical(),
      network_pathogen_type = character(),
      network_zoonotic_status = character(),
      network_canonicalization_status = character(),
      analysis_unit_ids = character(),
      readiness_row_count = integer(),
      manifest_status = character(),
      manifest_status_reason = character(),
      qa_flags = character(),
      genbank_run_source = character()
    )
}

country_records <- country_records %>%
  select(any_of(final_country_record_cols))

# ------------------------------------------------------------------------------|
#      Summarize pathogen-country and disease-country evidence ----------------
# ------------------------------------------------------------------------------|
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

# ------------------------------------------------------------------------------|
#      Write outputs -----------------------------------------------------------
# ------------------------------------------------------------------------------|
write_csv(
  country_records,
  genbank_simple_file_path(output_dir, paste0(output_prefix, "country_records.csv"), create_parent = TRUE)
)
write_csv(
  pathogen_country_summary,
  genbank_simple_file_path(output_dir, paste0(output_prefix, "pathogen_country_summary.csv"), create_parent = TRUE)
)
write_csv(
  disease_country_summary,
  genbank_simple_file_path(output_dir, paste0(output_prefix, "disease_country_summary.csv"), create_parent = TRUE)
)

if (nrow(search_logs) > 0) {
  write_csv(
    search_logs,
    genbank_simple_file_path(output_dir, paste0(output_prefix, "search_logs.csv"), create_parent = TRUE)
  )
}

message("Summary kind: ", summary_kind)
message("Wrote country records: ", nrow(country_records))
message("Wrote pathogen-country rows: ", nrow(pathogen_country_summary))
message("Wrote disease-country rows: ", nrow(disease_country_summary))
