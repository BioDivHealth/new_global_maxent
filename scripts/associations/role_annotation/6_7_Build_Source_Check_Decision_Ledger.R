#!/usr/bin/env Rscript
################################################################################
# 6_7_Build_Source_Check_Decision_Ledger.R
################################################################################
# Purpose: Build a review-only decision ledger for Deep Research candidate rows.
#          This script reads the consolidated candidate/source tables, preserves
#          user-added file_name metadata, and does not modify official role CSVs.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tidyr, tibble)

role_dir <- here::here("pathogen_association_data", "WHO", "role_annotation")
consolidated_dir <- file.path(role_dir, "deep_research_inputs", "consolidated")
papers_dir <- file.path(role_dir, "papers")
output_dir <- file.path(role_dir, "source_check")

candidate_queue_path <- file.path(consolidated_dir, "candidate_source_check_queue.csv")
source_request_path <- file.path(consolidated_dir, "candidate_source_request_list.csv")
unique_sources_path <- file.path(consolidated_dir, "candidate_unique_sources_to_fetch.csv")

required_paths <- c(candidate_queue_path, source_request_path, unique_sources_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required consolidated files: ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_stage_csv <- function(path) {
  read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

split_many <- function(x) {
  x <- coalesce(x, "")
  pieces <- unlist(str_split(x, "\\s*(?:[,;]|\\|)\\s*"), use.names = FALSE)
  pieces <- str_trim(pieces)
  pieces[pieces != ""]
}

candidate_queue <- read_stage_csv(candidate_queue_path) %>%
  mutate(candidate_row_id = paste0("candidate_", str_pad(row_number(), 3, pad = "0"))) %>%
  relocate(candidate_row_id, .before = 1)

source_request <- read_stage_csv(source_request_path)
unique_sources <- read_stage_csv(unique_sources_path)

if (!"file_name" %in% names(unique_sources)) {
  unique_sources$file_name <- NA_character_
}

source_file_lookup <- unique_sources %>%
  select(batch_id, source_id, source_title, source_url, file_name) %>%
  mutate(file_name = coalesce(file_name, ""))

if (!"file_name" %in% names(source_request)) {
  source_request$file_name <- NA_character_
}

source_request_with_files <- source_request %>%
  left_join(
    source_file_lookup,
    by = c("batch_id", "source_id", "source_title", "source_url"),
    suffix = c("", "_unique")
  ) %>%
  mutate(file_name = coalesce(.data$file_name_unique, .data$file_name)) %>%
  select(-any_of("file_name_unique"))

source_file_status <- source_request_with_files %>%
  select(
    batch_id,
    source_id,
    source_title,
    source_url,
    doi,
    pmid,
    pmcid,
    source_access,
    file_name
  ) %>%
  distinct() %>%
  mutate(file_name_piece = map(file_name, split_many)) %>%
  unnest_longer(file_name_piece, values_to = "file_name_piece", keep_empty = TRUE) %>%
  mutate(
    file_name_piece = if_else(is.na(file_name_piece), NA_character_, file_name_piece),
    local_pdf_path = if_else(
      is.na(file_name_piece) | file_name_piece == "",
      NA_character_,
      file.path(papers_dir, file_name_piece)
    ),
    local_pdf_exists = if_else(
      is.na(local_pdf_path),
      FALSE,
      file.exists(local_pdf_path)
    )
  ) %>%
  arrange(batch_id, source_id, file_name_piece)

candidate_sources_collapsed <- source_request_with_files %>%
  group_by(batch_id, disease_name, entity_type, entity_name, role_assignment, assignment_confidence) %>%
  summarise(
    source_ids = paste(unique(source_id), collapse = "; "),
    source_titles = paste(unique(source_title), collapse = " | "),
    source_urls = paste(unique(source_url), collapse = " | "),
    doi = paste(unique(na.omit(doi)), collapse = " | "),
    pmid = paste(unique(na.omit(pmid)), collapse = " | "),
    pmcid = paste(unique(na.omit(pmcid)), collapse = " | "),
    source_access = paste(unique(na.omit(source_access)), collapse = " | "),
    file_name = paste(unique(na.omit(file_name)), collapse = ", "),
    .groups = "drop"
  )

decision_ledger <- candidate_queue %>%
  left_join(
    candidate_sources_collapsed,
    by = c("batch_id", "disease_name", "entity_type", "entity_name", "role_assignment", "assignment_confidence")
  ) %>%
  mutate(
    local_pdf_paths = map_chr(file_name, function(value) {
      pieces <- split_many(value)
      if (length(pieces) == 0) return("")
      paste(file.path(papers_dir, pieces), collapse = " | ")
    }),
    local_pdf_status = map_chr(local_pdf_paths, function(value) {
      paths <- split_many(value)
      if (length(paths) == 0) return("no_local_pdf_expected_or_provided")
      exists <- file.exists(paths)
      if (all(exists)) return("all_local_pdfs_found")
      if (any(exists)) return("some_local_pdfs_found")
      "local_pdf_missing"
    }),
    source_checked = "no",
    source_check_method = "",
    evidence_found = "",
    checked_evidence_span = "",
    checked_evidence_location = "",
    checked_source_url = source_urls,
    checked_doi = doi,
    checked_pmid = pmid,
    checked_pmcid = pmcid,
    decision = "pending",
    accepted_role = "",
    accepted_confidence = "",
    accepted_evidence_scope = "",
    caveat = "",
    official_csv_target = "",
    decision_reason = "",
    reviewer = "Codex",
    review_date = as.character(Sys.Date()),
    import_ready = "no"
  ) %>%
  select(
    candidate_row_id,
    batch_order,
    phase,
    batch_id,
    disease_name,
    entity_type,
    entity_name,
    role_assignment,
    assignment_confidence,
    review_priority,
    evidence_source_ids,
    source_ids,
    source_titles,
    source_urls,
    doi,
    pmid,
    pmcid,
    source_access,
    file_name,
    local_pdf_paths,
    local_pdf_status,
    evidence_basis,
    review_reason,
    join_note,
    source_check_note,
    source_checked,
    source_check_method,
    evidence_found,
    checked_evidence_span,
    checked_evidence_location,
    checked_source_url,
    checked_doi,
    checked_pmid,
    checked_pmcid,
    decision,
    accepted_role,
    accepted_confidence,
    accepted_evidence_scope,
    caveat,
    official_csv_target,
    decision_reason,
    reviewer,
    review_date,
    import_ready
  )

write_csv(decision_ledger, file.path(output_dir, "candidate_source_check_decisions.csv"), na = "")
write_csv(source_request_with_files, file.path(output_dir, "candidate_source_request_list_with_files.csv"), na = "")
write_csv(source_file_status, file.path(output_dir, "source_file_status.csv"), na = "")

summary <- tibble(
  candidate_rows = nrow(decision_ledger),
  source_links = nrow(source_request_with_files),
  unique_source_file_rows = nrow(source_file_status),
  candidate_rows_with_any_file_name = sum(!is.na(decision_ledger$file_name) & decision_ledger$file_name != ""),
  candidate_rows_all_local_pdfs_found = sum(decision_ledger$local_pdf_status == "all_local_pdfs_found"),
  candidate_rows_some_local_pdfs_found = sum(decision_ledger$local_pdf_status == "some_local_pdfs_found"),
  candidate_rows_without_local_pdf = sum(decision_ledger$local_pdf_status == "no_local_pdf_expected_or_provided"),
  candidate_rows_missing_local_pdf = sum(decision_ledger$local_pdf_status == "local_pdf_missing")
)
write_csv(summary, file.path(output_dir, "source_check_summary.csv"), na = "")

readme_lines <- c(
  "# Source-Check Decision Ledger",
  "",
  "This folder contains review-only source-check artifacts for candidate Deep Research role rows.",
  "",
  "No official role evidence or assignment CSVs are modified by this workflow.",
  "",
  "Core files:",
  "",
  "- `candidate_source_check_decisions.csv`: one row per candidate role claim; fill decisions here.",
  "- `candidate_source_request_list_with_files.csv`: candidate-source links with the user-provided `file_name` metadata joined in.",
  "- `source_file_status.csv`: local PDF existence checks for each source/file pointer.",
  "- `source_check_summary.csv`: compact counts for the decision ledger.",
  "",
  "Decision vocabulary:",
  "",
  "- `accept`: source supports official evidence and assignment using existing vocabulary.",
  "- `accept_evidence_only`: useful source-backed evidence, but not an assignment-ready role row.",
  "- `defer`: taxonomy, role vocabulary, source-access, or interpretation issue remains.",
  "- `reject`: source does not support the proposed candidate role.",
  "",
  "Generated by:",
  "",
  "`Rscript scripts/associations/role_annotation/6_7_Build_Source_Check_Decision_Ledger.R`"
)
writeLines(readme_lines, file.path(output_dir, "README.md"), useBytes = TRUE)

message("Wrote source-check decision ledger.")
print(summary)
