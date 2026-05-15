# ------------------------------------------------------------------------------|
#      04_quality_checks.R -----------------------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Write a compact QA summary for GenBank-simple manifest and retrieved
#          checkpoint outputs.
# Inputs : genbank_simple_manifest.csv
#          excluded_targets.csv
#          genbank_search_logs.csv
#          genbank_country_records.csv
#          Optional `GENBANK_SIMPLE_SUMMARY_KIND=readiness_combined` reads
#          readiness-prefixed combined outputs from the top-level or
#          reorganized `qa/` and `intermediate/` folders.
# Outputs: genbank_simple_qa_summary.csv
#          genbank_simple_target_qa.csv
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))

# ------------------------------------------------------------------------------|
#      Resolve run mode and output files --------------------------------------
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

manifest_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_simple_readiness_manifest.csv",
  "genbank_simple_manifest.csv"
)
search_log_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_readiness_search_logs.csv",
  "genbank_search_logs.csv"
)
country_records_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_readiness_country_records.csv",
  "genbank_country_records.csv"
)
qa_summary_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_readiness_qa_summary.csv",
  "genbank_simple_qa_summary.csv"
)
target_qa_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_readiness_target_qa.csv",
  "genbank_simple_target_qa.csv"
)

# ------------------------------------------------------------------------------|
#      Load manifest, logs, and country records -------------------------------
# ------------------------------------------------------------------------------|
read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }

  read_csv(path, show_col_types = FALSE, na = c("", "NA"))
}

manifest <- read_optional_csv(file.path(output_dir, manifest_file))
excluded_targets <- read_optional_csv(file.path(output_dir, "excluded_targets.csv"))
search_logs <- read_optional_csv(genbank_simple_existing_file_path(output_dir, search_log_file))
country_records <- read_optional_csv(genbank_simple_existing_file_path(output_dir, country_records_file))

# ------------------------------------------------------------------------------|
#      Normalize empty or missing log inputs ----------------------------------
# ------------------------------------------------------------------------------|
empty_log_cols <- c(
  "target_id",
  "status",
  "records_found",
  "ids_collected",
  "records_parsed",
  "countries_observed",
  "note"
)

if (!all(empty_log_cols %in% names(search_logs))) {
  search_logs <- tibble(
    target_id = character(),
    status = character(),
    records_found = integer(),
    ids_collected = integer(),
    records_parsed = integer(),
    countries_observed = integer(),
    note = character()
  )
}

# ------------------------------------------------------------------------------|
#      Build per-target QA table ----------------------------------------------
# ------------------------------------------------------------------------------|
manifest_for_qa <- if (summary_kind == "readiness_combined") {
  manifest %>%
    transmute(
      target_id,
      current_target_id = clean_text(current_target_id),
      Pathogens = query_pathogen_label,
      Disease_name = readiness_disease_names,
      query_used,
      in_gibb_etal = NA,
      in_empres_i = NA,
      manifest_status,
      qa_flags
    )
} else {
  manifest %>%
    transmute(
      target_id,
      current_target_id = NA_character_,
      Pathogens,
      Disease_name,
      query_used,
      in_gibb_etal,
      in_empres_i,
      manifest_status = NA_character_,
      qa_flags = NA_character_
    )
}

search_logs_for_join <- if (summary_kind == "readiness_combined" && nrow(search_logs) > 0) {
  manifest_for_qa %>%
    filter(!is.na(current_target_id)) %>%
    select(readiness_target_id = target_id, current_target_id) %>%
    right_join(search_logs, by = c("current_target_id" = "target_id")) %>%
    mutate(target_id = dplyr::coalesce(readiness_target_id, current_target_id)) %>%
    select(-readiness_target_id, -current_target_id)
} else {
  search_logs
}

target_qa <- manifest_for_qa %>%
  select(-current_target_id) %>%
  left_join(
    search_logs_for_join %>%
      mutate(
        status = if_else(
          status == "skipped_records_found_exceeds_limit",
          status,
          status
        )
      ) %>%
      select(
        target_id,
        status,
        records_found,
        ids_collected,
        records_parsed,
        countries_observed,
        note
      ),
    by = "target_id"
  ) %>%
  mutate(
    has_retrieval_log = !is.na(status),
    records_found = suppressWarnings(as.integer(records_found)),
    records_parsed = suppressWarnings(as.integer(records_parsed)),
    countries_observed = suppressWarnings(as.integer(countries_observed)),
    qa_flag = case_when(
      !has_retrieval_log ~ "not_run",
      status == "skipped_records_found_exceeds_limit" ~ "skipped_records_found_exceeds_limit",
      status == "no_records" | dplyr::coalesce(records_found, 0L) == 0L ~ "zero_records",
      status != "success" ~ "retrieval_failed",
      dplyr::coalesce(records_parsed, 0L) > 0L & dplyr::coalesce(countries_observed, 0L) == 0L ~ "records_without_country",
      TRUE ~ "ok"
    )
  ) %>%
  arrange(qa_flag, Pathogens, Disease_name)

# ------------------------------------------------------------------------------|
#      Summarize run-level QA metrics -----------------------------------------
# ------------------------------------------------------------------------------|
country_record_rows <- nrow(country_records)
unique_countries <- if ("country" %in% names(country_records)) {
  dplyr::n_distinct(country_records$country, na.rm = TRUE)
} else {
  0L
}
unique_pathogens_with_country <- if (all(c("Pathogens", "country") %in% names(country_records))) {
  dplyr::n_distinct(country_records$Pathogens[!is.na(country_records$country)], na.rm = TRUE)
} else {
  0L
}

qa_summary <- tibble(
  metric = c(
    "manifest_targets",
    "excluded_targets",
    "excluded_coronavirus_scope_deferred",
    "excluded_broad_or_unwanted_influenza",
    "targets_not_run",
    "targets_zero_records",
    "targets_retrieval_failed",
    "targets_skipped_records_found_exceeds_limit",
    "targets_records_without_country",
    "targets_ok",
    "country_record_rows",
    "unique_countries",
    "unique_pathogens_with_country"
  ),
  value = c(
    nrow(manifest),
    nrow(excluded_targets),
    sum(excluded_targets$exclusion_reason == "coronavirus_scope_deferred", na.rm = TRUE),
    sum(excluded_targets$exclusion_reason == "broad_or_unwanted_influenza", na.rm = TRUE),
    sum(target_qa$qa_flag == "not_run", na.rm = TRUE),
    sum(target_qa$qa_flag == "zero_records", na.rm = TRUE),
    sum(target_qa$qa_flag == "retrieval_failed", na.rm = TRUE),
    sum(target_qa$qa_flag == "skipped_records_found_exceeds_limit", na.rm = TRUE),
    sum(target_qa$qa_flag == "records_without_country", na.rm = TRUE),
    sum(target_qa$qa_flag == "ok", na.rm = TRUE),
    country_record_rows,
    unique_countries,
    unique_pathogens_with_country
  )
)

# ------------------------------------------------------------------------------|
#      Write outputs -----------------------------------------------------------
# ------------------------------------------------------------------------------|
write_csv(
  qa_summary,
  genbank_simple_file_path(output_dir, qa_summary_file, create_parent = TRUE)
)
write_csv(
  target_qa,
  genbank_simple_file_path(output_dir, target_qa_file, create_parent = TRUE)
)

message("Summary kind: ", summary_kind)
message("Wrote QA summary rows: ", nrow(qa_summary))
message("Wrote target QA rows: ", nrow(target_qa))
