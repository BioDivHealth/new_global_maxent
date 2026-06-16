#!/usr/bin/env Rscript
################################################################################
# 6_7_Build_Source_Check_Decision_Ledger.R
################################################################################
# Purpose: Build a review-only decision ledger for source-check candidate rows.
#          This script reads generated Deep Research candidates plus optional
#          manual candidate/source inputs, preserves user-added file_name
#          metadata, and does not modify official role CSVs.
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
source(here::here("scripts", "associations", "association_data_helpers.R"))

role_dir <- role_annotation_dir
consolidated_dir <- role_deep_research_consolidated_dir
papers_dir <- role_source_pdf_dir
input_dir <- role_source_check_input_dir
output_dir <- role_source_check_dir

candidate_queue_path <- file.path(consolidated_dir, "candidate_source_check_queue.csv")
source_request_path <- file.path(consolidated_dir, "candidate_source_request_list.csv")
unique_sources_path <- file.path(consolidated_dir, "candidate_unique_sources_to_fetch.csv")
manual_candidate_path <- file.path(input_dir, "manual_source_check_candidates.csv")
manual_source_path <- file.path(input_dir, "manual_source_check_sources.csv")
curated_decisions_path <- file.path(output_dir, "curated_source_check_decisions.csv")
role_gap_candidate_path <- role_gap_source_check_candidates_path()
candidate_id_overrides_path <- role_candidate_id_overrides_path()

required_paths <- c(candidate_queue_path, source_request_path, unique_sources_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required consolidated files: ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

read_stage_csv <- function(path) {
  read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

empty_stage_csv <- function(columns) {
  tibble(!!!setNames(rep(list(character()), length(columns)), columns))
}

read_optional_stage_csv <- function(path, columns) {
  if (!file.exists(path)) {
    return(empty_stage_csv(columns))
  }

  data <- read_stage_csv(path)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Manual source-check input is missing required columns in ", path, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data %>%
    select(all_of(columns))
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

candidate_queue_columns <- c(
  "candidate_row_id",
  "batch_order",
  "phase",
  "batch_id",
  "disease_name",
  "entity_type",
  "entity_name",
  "role_assignment",
  "assignment_confidence",
  "review_priority",
  "evidence_source_ids",
  "evidence_basis",
  "review_reason",
  "join_note",
  "source_check_note"
)

candidate_id_override_columns <- c(
  "candidate_id_key",
  "candidate_row_id",
  "batch_id",
  "disease_name",
  "entity_type",
  "entity_name",
  "role_assignment",
  "assignment_confidence"
)

source_request_columns <- c(
  "batch_order",
  "phase",
  "batch_id",
  "disease_name",
  "entity_type",
  "entity_name",
  "role_assignment",
  "assignment_confidence",
  "evidence_basis",
  "review_reason",
  "source_id",
  "source_lookup_status",
  "source_title",
  "authors_or_organization",
  "year",
  "source_type",
  "source_url",
  "doi",
  "pmid",
  "pmcid",
  "source_access",
  "rows_supported",
  "reliability_note",
  "source_check_note",
  "candidate_row",
  "file_name"
)

source_check_candidate_key <- function(batch_id, disease_name, entity_type,
                                       entity_name, role_assignment,
                                       assignment_confidence) {
  stable_candidate_id(
    "source_check",
    batch_id,
    disease_name,
    entity_type,
    entity_name,
    role_assignment,
    assignment_confidence
  )
}

read_candidate_id_overrides <- function(path) {
  overrides <- read_optional_stage_csv(path, candidate_id_override_columns)
  if (nrow(overrides) == 0) {
    return(overrides)
  }

  duplicated_keys <- overrides$candidate_id_key[duplicated(overrides$candidate_id_key)]
  if (length(duplicated_keys) > 0) {
    stop(
      "Candidate ID override keys are duplicated: ",
      paste(unique(duplicated_keys), collapse = ", "),
      call. = FALSE
    )
  }

  duplicated_ids <- overrides$candidate_row_id[duplicated(overrides$candidate_row_id)]
  if (length(duplicated_ids) > 0) {
    stop(
      "Candidate ID override values are duplicated: ",
      paste(unique(duplicated_ids), collapse = ", "),
      call. = FALSE
    )
  }

  overrides
}

validate_manual_candidate_ids <- function(manual_candidates) {
  if (nrow(manual_candidates) == 0) {
    return(invisible(NULL))
  }

  missing_ids <- manual_candidates %>%
    filter(is.na(candidate_row_id) | candidate_row_id == "") %>%
    pull(candidate_row_id)
  if (length(missing_ids) > 0) {
    stop("Manual source-check candidates must have stable candidate_row_id values.", call. = FALSE)
  }

  invalid_ids <- manual_candidates %>%
    filter(!str_detect(candidate_row_id, "^manual_[A-Za-z0-9_]+$")) %>%
    pull(candidate_row_id)
  if (length(invalid_ids) > 0) {
    stop(
      "Manual source-check candidate IDs must match ^manual_[A-Za-z0-9_]+$: ",
      paste(invalid_ids, collapse = ", "),
      call. = FALSE
    )
  }
}

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

candidate_id_overrides <- read_candidate_id_overrides(candidate_id_overrides_path)

deep_candidate_queue <- read_stage_csv(candidate_queue_path) %>%
  mutate(
    candidate_id_key = source_check_candidate_key(
      batch_id,
      disease_name,
      entity_type,
      entity_name,
      role_assignment,
      assignment_confidence
    )
  ) %>%
  left_join(
    candidate_id_overrides %>% select(candidate_id_key, candidate_row_id),
    by = "candidate_id_key"
  ) %>%
  mutate(
    candidate_row_id = coalesce(
      candidate_row_id,
      source_check_candidate_key(
        batch_id,
        disease_name,
        entity_type,
        entity_name,
        role_assignment,
        assignment_confidence
      )
    )
  ) %>%
  relocate(candidate_row_id, .before = 1) %>%
  select(-candidate_id_key) %>%
  select(all_of(candidate_queue_columns))

deep_source_request <- read_stage_csv(source_request_path)
unique_sources <- read_stage_csv(unique_sources_path)

manual_candidate_queue <- read_optional_stage_csv(manual_candidate_path, candidate_queue_columns)
manual_source_request <- read_optional_stage_csv(manual_source_path, source_request_columns)
role_gap_candidate_queue <- read_optional_stage_csv(role_gap_candidate_path, candidate_queue_columns)
validate_manual_candidate_ids(manual_candidate_queue)

candidate_queue <- bind_rows(deep_candidate_queue, manual_candidate_queue, role_gap_candidate_queue)

duplicated_candidate_ids <- unique(candidate_queue$candidate_row_id[duplicated(candidate_queue$candidate_row_id)])
if (length(duplicated_candidate_ids) > 0) {
  stop(
    "Source-check candidate inputs contain duplicate candidate_row_id values: ",
    paste(duplicated_candidate_ids, collapse = ", "),
    call. = FALSE
  )
}

if (!"file_name" %in% names(unique_sources)) {
  unique_sources$file_name <- NA_character_
}

source_file_lookup <- unique_sources %>%
  select(batch_id, source_id, source_title, source_url, file_name) %>%
  mutate(file_name = coalesce(file_name, ""))

if (!"file_name" %in% names(deep_source_request)) {
  deep_source_request$file_name <- NA_character_
}

deep_source_request_with_files <- deep_source_request %>%
  left_join(
    source_file_lookup,
    by = c("batch_id", "source_id", "source_title", "source_url"),
    suffix = c("", "_unique")
  ) %>%
  mutate(file_name = coalesce(.data$file_name_unique, .data$file_name)) %>%
  select(-any_of("file_name_unique")) %>%
  select(all_of(source_request_columns))

source_request_with_files <- bind_rows(
  deep_source_request_with_files,
  manual_source_request
)

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
  deep_research_candidate_rows = nrow(deep_candidate_queue),
  manual_candidate_rows = nrow(manual_candidate_queue),
  role_gap_candidate_rows = nrow(role_gap_candidate_queue),
  source_links = nrow(source_request_with_files),
  deep_research_source_links = nrow(deep_source_request_with_files),
  manual_source_links = nrow(manual_source_request),
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
  "- `candidate_source_check_decisions.csv` is regenerated by merging fresh Deep Research and manual candidate/source metadata with those curated decisions.",
  "- Official role evidence and assignment CSVs are not modified by this ledger build.",
  "",
  "Input summary:",
  "",
  paste0("- Deep Research candidates: ", nrow(deep_candidate_queue)),
  paste0("- Manual candidates: ", nrow(manual_candidate_queue)),
  paste0("- Role-gap candidates: ", nrow(role_gap_candidate_queue)),
  "",
  "Decision summary:",
  "",
  paste(capture.output(print(decision_summary, n = Inf)), collapse = "\n")
)
writeLines(progress_lines, file.path(output_dir, "SOURCE_CHECK_PROGRESS.md"), useBytes = TRUE)

readme_lines <- c(
  "# Source-Check Decision Ledger",
  "",
  "This folder contains review-only source-check artifacts for role rows from generated Deep Research outputs, durable manual candidate inputs, and generated role-gap candidates.",
  "Deep Research inputs are on-demand curation artifacts, not routine pipeline steps; the ledger consumes the current consolidated outputs when they exist.",
  "",
  "No official role evidence or assignment CSVs are modified by this workflow.",
  "",
  "Input files:",
  "",
  "- `input/manual_source_check_candidates.csv`: optional durable manual candidate rows, using stable `manual_*` candidate IDs.",
  "- `input/manual_source_check_sources.csv`: optional durable manual candidate-source links, including local PDF/text filename metadata when available.",
  "- `input/candidate_id_overrides.csv`: durable mapping that preserves historical generated candidate IDs when regenerated candidates have stable deterministic keys.",
  "- Deep Research candidate/source inputs are read from `pathogen_association_data/staged/role_annotation/deep_research_inputs/consolidated/`.",
  "- Role-gap candidates are read from `pathogen_association_data/staged/role_annotation/source_check_candidates/role_gap_source_check_candidates.csv` when present.",
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
  "- Current source-check decision counts are written to",
  "  `source_check_decision_summary.csv`.",
  "",
  "Generated by:",
  "",
  "`Rscript scripts/associations/role_annotation/source_check/01_build_source_check_decision_ledger.R`",
  "",
  "Import checked rows with:",
  "",
  "`Rscript scripts/associations/role_annotation/source_check/02_import_source_checked_role_rows.R`"
)
writeLines(readme_lines, file.path(output_dir, "README.md"), useBytes = TRUE)

message("Wrote source-check decision ledger.")
print(summary)
print(decision_summary)
