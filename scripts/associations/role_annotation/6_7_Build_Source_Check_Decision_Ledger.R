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

source(here::here("scripts", "associations", "working_inputs.R"))

role_dir <- role_annotation_dir
consolidated_dir <- role_deep_research_consolidated_dir
papers_dir <- role_source_pdf_dir
output_dir <- role_source_check_dir

candidate_queue_path <- file.path(consolidated_dir, "candidate_source_check_queue.csv")
source_request_path <- file.path(consolidated_dir, "candidate_source_request_list.csv")
unique_sources_path <- file.path(consolidated_dir, "candidate_unique_sources_to_fetch.csv")
curated_decisions_path <- file.path(output_dir, "curated_source_check_decisions.csv")

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

repo_relative_path <- function(path) {
  repo_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  repo_prefix <- paste0(repo_root, "/")
  is_repo_path <- !is.na(normalized) & startsWith(normalized, repo_prefix)

  normalized[is_repo_path] <- substring(normalized[is_repo_path], nchar(repo_prefix) + 1L)
  normalized
}

curated_identity_columns <- c(
  "candidate_row_id",
  "batch_id",
  "disease_name",
  "entity_type",
  "entity_name",
  "role_assignment",
  "assignment_confidence"
)

curated_decision_columns <- c(
  "source_checked",
  "source_check_method",
  "evidence_found",
  "checked_evidence_span",
  "checked_evidence_location",
  "decision",
  "accepted_role",
  "accepted_confidence",
  "accepted_evidence_scope",
  "caveat",
  "official_csv_target",
  "decision_reason",
  "reviewer",
  "review_date",
  "import_ready"
)

apply_curated_source_check_decisions <- function(decision_ledger, curated_path) {
  if (!file.exists(curated_path)) {
    return(decision_ledger)
  }

  curated_decisions <- read_stage_csv(curated_path)
  required_columns <- c(curated_identity_columns, curated_decision_columns)
  missing_columns <- setdiff(required_columns, names(curated_decisions))
  if (length(missing_columns) > 0) {
    stop(
      "Curated source-check decisions are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  duplicated_ids <- unique(curated_decisions$candidate_row_id[duplicated(curated_decisions$candidate_row_id)])
  if (length(duplicated_ids) > 0) {
    stop(
      "Curated source-check decisions contain duplicate candidate_row_id values: ",
      paste(duplicated_ids, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_ids <- setdiff(curated_decisions$candidate_row_id, decision_ledger$candidate_row_id)
  if (length(unknown_ids) > 0) {
    stop(
      "Curated source-check decisions refer to candidate_row_id values absent from the regenerated ledger: ",
      paste(unknown_ids, collapse = ", "),
      call. = FALSE
    )
  }

  identity_comparison <- decision_ledger %>%
    select(all_of(curated_identity_columns)) %>%
    inner_join(
      curated_decisions %>% select(all_of(curated_identity_columns)),
      by = "candidate_row_id",
      suffix = c(".ledger", ".curated")
    )

  identity_mismatches <- map_dfr(setdiff(curated_identity_columns, "candidate_row_id"), function(column) {
    ledger_column <- paste0(column, ".ledger")
    curated_column <- paste0(column, ".curated")

    identity_comparison %>%
      filter(coalesce(.data[[ledger_column]], "") != coalesce(.data[[curated_column]], "")) %>%
      transmute(
        candidate_row_id,
        column = column,
        ledger_value = .data[[ledger_column]],
        curated_value = .data[[curated_column]]
      )
  })

  if (nrow(identity_mismatches) > 0) {
    mismatch_preview <- identity_mismatches %>%
      mutate(summary = paste0(candidate_row_id, ":", column)) %>%
      pull(summary) %>%
      head(10)

    stop(
      "Curated source-check decisions no longer match regenerated candidate identities: ",
      paste(mismatch_preview, collapse = ", "),
      call. = FALSE
    )
  }

  curated_values <- curated_decisions %>%
    select(candidate_row_id, all_of(curated_decision_columns)) %>%
    rename_with(~ paste0(.x, ".curated"), all_of(curated_decision_columns))

  filled_ledger <- decision_ledger %>%
    left_join(curated_values, by = "candidate_row_id")

  for (column in curated_decision_columns) {
    curated_column <- paste0(column, ".curated")
    filled_ledger[[column]] <- if_else(
      !is.na(filled_ledger[[curated_column]]),
      filled_ledger[[curated_column]],
      filled_ledger[[column]]
    )
    filled_ledger[[curated_column]] <- NULL
  }

  filled_ledger
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
    local_pdf_path_absolute = if_else(
      is.na(file_name_piece) | file_name_piece == "",
      NA_character_,
      file.path(papers_dir, file_name_piece)
    ),
    local_pdf_exists = if_else(
      is.na(local_pdf_path_absolute),
      FALSE,
      file.exists(local_pdf_path_absolute)
    ),
    local_pdf_path = if_else(
      is.na(local_pdf_path_absolute),
      NA_character_,
      repo_relative_path(local_pdf_path_absolute)
    )
  ) %>%
  select(-local_pdf_path_absolute) %>%
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
    file_name = paste(unique(file_name[!is.na(file_name) & file_name != ""]), collapse = ", "),
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
      paste(repo_relative_path(file.path(papers_dir, pieces)), collapse = " | ")
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

decision_ledger <- apply_curated_source_check_decisions(decision_ledger, curated_decisions_path)

write_csv(decision_ledger, file.path(output_dir, "candidate_source_check_decisions.csv"), na = "")
write_csv(source_request_with_files, file.path(output_dir, "candidate_source_request_list_with_files.csv"), na = "")
write_csv(source_file_status, file.path(output_dir, "source_file_status.csv"), na = "")

decision_summary <- decision_ledger %>%
  count(decision, import_ready, name = "n") %>%
  arrange(decision, import_ready)

write_csv(decision_summary, file.path(output_dir, "source_check_decision_summary.csv"), na = "")

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

progress_lines <- c(
  "# Source-Check Progress",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "Durable curated decisions:",
  "",
  "- `curated_source_check_decisions.csv` stores manual/source-checked curation fields.",
  "- `candidate_source_check_decisions.csv` is regenerated by merging fresh candidate/source metadata with those curated decisions.",
  "- Official role evidence and assignment CSVs are not modified by this ledger build.",
  "",
  "Decision summary:",
  "",
  paste(capture.output(print(decision_summary, n = Inf)), collapse = "\n")
)
writeLines(progress_lines, file.path(output_dir, "SOURCE_CHECK_PROGRESS.md"), useBytes = TRUE)

readme_lines <- c(
  "# Source-Check Decision Ledger",
  "",
  "This folder contains review-only source-check artifacts for candidate Deep Research role rows.",
  "",
  "No official role evidence or assignment CSVs are modified by this workflow.",
  "",
  "Core files:",
  "",
  "- `curated_source_check_decisions.csv`: durable manual/source-checked curation decisions keyed by candidate identity.",
  "- `candidate_source_check_decisions.csv`: regenerated one-row-per-candidate ledger with source metadata and curated decisions applied.",
  "- `candidate_source_request_list_with_files.csv`: candidate-source links with the user-provided `file_name` metadata joined in.",
  "- `source_file_status.csv`: local PDF existence checks for each source/file pointer.",
  "- `source_check_summary.csv`: compact counts for the decision ledger.",
  "- `source_check_decision_summary.csv`: counts by source-check decision and import-readiness.",
  "",
  "Decision vocabulary:",
  "",
  "- `accept`: source supports official evidence and assignment using existing vocabulary.",
  "- `accept_evidence_only`: useful source-backed evidence, but not an assignment-ready role row.",
  "- `defer`: taxonomy, role vocabulary, source-access, or interpretation issue remains.",
  "- `reject`: source does not support the proposed candidate role.",
  "",
  "Import status:",
  "",
  "- The source-check import script is idempotent. It skips rows whose",
  "  `source_check_candidate_id` is already present in the official role CSVs.",
  "- A rerun reporting `+0` official row deltas is expected after the accepted rows",
  "  have already been imported.",
  "- In the current package, 45 accepted/import-ready rows are already represented",
  "  in the official role evidence and assignment CSVs; 8 rows remain excluded as",
  "  evidence-only or deferred.",
  "",
  "Generated by:",
  "",
  "`Rscript scripts/associations/role_annotation/6_7_Build_Source_Check_Decision_Ledger.R`",
  "",
  "Import checked rows with:",
  "",
  "`Rscript scripts/associations/role_annotation/6_9_Import_Source_Checked_Role_Rows.R`"
)
writeLines(readme_lines, file.path(output_dir, "README.md"), useBytes = TRUE)

message("Wrote source-check decision ledger.")
print(summary)
print(decision_summary)
