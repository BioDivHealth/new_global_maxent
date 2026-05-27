library(readr)
library(dplyr)
library(jsonlite)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_paths.R"))

v2_read_csv <- function(path, required_cols = character()) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", path, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  x
}

v2_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
  invisible(path)
}

v2_verbose_qa_enabled <- function() {
  identical(Sys.getenv("WHO_DON_V2_VERBOSE_QA"), "1")
}

v2_write_stage_diagnostic <- function(x, filename) {
  if (!v2_verbose_qa_enabled()) {
    return(invisible(NA_character_))
  }
  v2_write_csv(x, who_don_v2_qa_archive_dir("stage_diagnostics", filename))
}

v2_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  invisible(path)
}

v2_required_final_cols <- c(
  "record_key", "DonId", "record_id", "Title", "publication_datetime_utc",
  "article_url", "country_standard", "disease_label_clean",
  "disease_label_raw", "association_scope", "final_country_role",
  "final_event_country_flag", "final_event_confidence",
  "don_country_report_scope", "strict_focal_event_country_flag"
)

v2_clean_final_path <- function() {
  reference_path <- who_don_v2_reference_dir("who_don_clean_final_seed.csv")
  if (file.exists(reference_path)) {
    return(reference_path)
  }
  who_don_clean_output_dir("final", "who_don_country_disease_event_focal_scope_evidence_final.csv")
}

v2_clean_modelling_path <- function() {
  reference_path <- who_don_v2_reference_dir("who_don_clean_modelling_seed.csv")
  if (file.exists(reference_path)) {
    return(reference_path)
  }
  who_don_clean_output_dir("final", "who_don_country_disease_event_focal_modelling_ready_final.csv")
}

v2_clean_records_path <- function() {
  reference_path <- who_don_v2_reference_dir("who_don_clean_records_seed.csv")
  if (file.exists(reference_path)) {
    return(reference_path)
  }
  who_don_clean_output_dir("records", "who_don_records_clean.csv")
}

v2_read_clean_final <- function() {
  v2_read_csv(v2_clean_final_path(), v2_required_final_cols)
}

v2_read_association_contract <- function() {
  v2_read_csv(who_don_v2_rules_dir("accepted_association_contract.csv"), v2_required_final_cols)
}

v2_read_records_source <- function() {
  v2_read_csv(
    who_don_v2_output_dir("records", "who_don_records_source.csv"),
    c("DonId", "record_id", "Title", "publication_datetime_utc", "article_url")
  )
}

v2_read_clean_records_seed <- function() {
  v2_read_csv(
    v2_clean_records_path(),
    c("DonId", "record_id", "Title", "publication_datetime_utc", "article_url")
  )
}

v2_read_records <- function() {
  records <- v2_read_csv(
    who_don_v2_output_dir("records", "who_don_records_clean.csv"),
    c("DonId", "record_id", "Title", "publication_datetime_utc", "article_url")
  )
  if (!"record_key" %in% names(records)) {
    records <- records %>%
      mutate(record_key = record_id)
  }
  records
}

v2_empty_review_decisions <- function() {
  tibble::tibble(
    review_id = character(),
    record_key = character(),
    decision_type = character(),
    decision_value = character(),
    confidence = character(),
    evidence_span = character(),
    reviewer = character(),
    review_source = character(),
    review_note = character(),
    review_date = as.Date(character())
  )
}

v2_first_present <- function(.data, cols) {
  cols <- intersect(cols, names(.data))
  if (length(cols) == 0) {
    rep(NA_character_, nrow(.data))
  } else {
    dplyr::coalesce(!!!lapply(cols, function(col) .data[[col]]))
  }
}
