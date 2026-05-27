#!/usr/bin/env Rscript
################################################################################
# 6_6_Consolidate_Deep_Research_Reports.R
################################################################################
# Purpose: Consolidate the six reformatted Deep Research batch outputs into
#          reviewable all-batch staging artifacts without changing official role
#          evidence or assignment tables.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tibble, tidyr)

# ------------------------------------------------------------------------------|
#      Define paths ------------------------------------------------------------|
# ------------------------------------------------------------------------------|
role_dir <- here::here("pathogen_association_data", "WHO", "role_annotation")
input_root <- file.path(role_dir, "deep_research_inputs")
manifest_path <- file.path(input_root, "deep_research_reformat_manifest.csv")
output_dir <- file.path(input_root, "consolidated")
candidate_unique_sources_path <- file.path(output_dir, "candidate_unique_sources_to_fetch.csv")

if (!file.exists(manifest_path)) {
  stop(
    "Missing manifest. Run scripts/associations/role_annotation/",
    "6_5_Reformat_Deep_Research_Reports.R first.",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Load manifest and helpers ----------------------------------------------|
# ------------------------------------------------------------------------------|
manifest <- read_csv(manifest_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    batch_order = row_number(),
    phase = if_else(str_starts(batch_id, "v_"), "Phase V", "Phase N")
  )

read_optional_csv <- function(batch_id, filename) {
  path <- file.path(input_root, batch_id, "reformatted", filename)
  if (!file.exists(path)) {
    return(tibble())
  }

  read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(batch_id = batch_id, .before = 1)
}

read_all_batches <- function(filename) {
  combined <- map_dfr(manifest$batch_id, read_optional_csv, filename = filename)

  if (nrow(combined) == 0) {
    return(tibble())
  }

  if (!"source_report" %in% names(combined)) {
    combined$source_report <- NA_character_
  }
  if (!"phase" %in% names(combined)) {
    combined$phase <- NA_character_
  }

  combined %>%
    left_join(
      manifest %>%
        transmute(
          batch_id,
          batch_order,
          phase_manifest = phase,
          source_report_manifest = source_report
        ),
      by = "batch_id"
    ) %>%
    mutate(
      phase = coalesce(.data$phase, .data$phase_manifest),
      source_report = coalesce(.data$source_report, .data$source_report_manifest),
      .keep = "unused"
    ) %>%
    relocate(batch_order, phase, source_report, .after = batch_id)
}

normalize_confidence <- function(x) {
  case_when(
    is.na(x) | x == "" ~ NA_integer_,
    x == "high" ~ 1L,
    x == "medium" ~ 2L,
    x == "low-medium" ~ 2L,
    x == "low" ~ 3L,
    TRUE ~ 4L
  )
}

# Keep priority deterministic so the source-check queue is stable across runs.
assignment_priority <- function(action, confidence) {
  case_when(
    action == "candidate_add_after_source_check" & confidence %in% c("high", "medium", "low-medium") ~
      "P1_source_check_candidate",
    action == "candidate_add_after_source_check" ~
      "P2_source_check_low_or_unclear_candidate",
    action == "defer_vocabulary_or_taxonomy" ~
      "P3_taxonomy_or_vocabulary_decision",
    action == "evidence_only_group" ~
      "P4_evidence_only_group",
    action == "defer_insufficient_evidence" ~
      "P5_deferred_insufficient_evidence",
    action == "reject_or_negative_evidence_only" ~
      "P6_negative_or_unsupported_guardrail",
    action == "already_covered" ~
      "P7_already_covered_baseline",
    TRUE ~ "P8_other_review"
  )
}

# ------------------------------------------------------------------------------|
#      Load extracted batch tables --------------------------------------------|
# ------------------------------------------------------------------------------|
all_assignments <- read_all_batches("extracted_assignment_staging.csv") %>%
  mutate(
    confidence_rank = normalize_confidence(assignment_confidence),
    review_priority = assignment_priority(action, assignment_confidence),
    entity_type = if_else(is.na(entity_type) | entity_type == "", "host", entity_type)
  ) %>%
  arrange(
    match(
      review_priority,
      c(
        "P1_source_check_candidate",
        "P2_source_check_low_or_unclear_candidate",
        "P3_taxonomy_or_vocabulary_decision",
        "P4_evidence_only_group",
        "P5_deferred_insufficient_evidence",
        "P6_negative_or_unsupported_guardrail",
        "P7_already_covered_baseline",
        "P8_other_review"
      )
    ),
    confidence_rank,
    batch_order,
    disease_name,
    entity_type,
    entity_name
  )

all_host_evidence <- read_all_batches("extracted_host_role_evidence.csv")
all_vector_evidence <- read_all_batches("extracted_vector_role_evidence.csv")
all_sources <- read_all_batches("extracted_sources_used.csv")
all_quality_issues <- read_all_batches("extraction_quality_issues.csv")
all_candidate_presence <- read_all_batches("candidate_presence_check.csv")
all_vector_non_applicability <- read_all_batches("extracted_vector_non_applicability.csv")
all_cross_batch_summary <- read_all_batches("extracted_cross_batch_summary.csv")

# ------------------------------------------------------------------------------|
#      Build source-check queues ----------------------------------------------|
# ------------------------------------------------------------------------------|
candidate_source_check_queue <- all_assignments %>%
  filter(action == "candidate_add_after_source_check") %>%
  left_join(
    all_candidate_presence %>%
      select(batch_id, disease_name, entity_type, entity_name, join_note),
    by = c("batch_id", "disease_name", "entity_type", "entity_name")
  ) %>%
  mutate(
    join_note = coalesce(join_note, "not present in candidate_presence_check"),
    source_check_note = case_when(
      str_detect(evidence_source_ids, ";|,") ~ "multiple sources listed; verify each source supports the exact role grain",
      TRUE ~ "verify source URL/DOI/PMID and exact evidence span before import"
    )
  ) %>%
  arrange(confidence_rank, batch_order, disease_name, entity_type, entity_name)

candidate_source_request_list <- candidate_source_check_queue %>%
  mutate(
    source_id = str_split(coalesce(evidence_source_ids, ""), "\\s*[;,]\\s*")
  ) %>%
  unnest(source_id) %>%
  mutate(source_id = str_trim(source_id)) %>%
  filter(source_id != "") %>%
  left_join(
    all_sources %>%
      select(
        batch_id,
        source_id,
        source_disease_name = disease_name,
        source_title,
        authors_or_organization,
        year,
        source_type,
        source_url,
        doi,
        pmid,
        pmcid,
        source_access,
        rows_supported,
        reliability_note
      ),
    by = c("batch_id", "source_id")
  ) %>%
  mutate(
    source_lookup_status = if_else(
      is.na(source_title) | source_title == "",
      "missing_from_sources_used",
      "matched_sources_used"
    ),
    candidate_row = paste0(
      disease_name, " | ", entity_type, " | ", entity_name, " | ",
      role_assignment, " | ", assignment_confidence
    )
  ) %>%
  select(
    batch_order,
    phase,
    batch_id,
    disease_name,
    entity_type,
    entity_name,
    role_assignment,
    assignment_confidence,
    evidence_basis,
    review_reason,
    source_id,
    source_lookup_status,
    source_title,
    authors_or_organization,
    year,
    source_type,
    source_url,
    doi,
    pmid,
    pmcid,
    source_access,
    rows_supported,
    reliability_note,
    source_check_note,
    candidate_row
  ) %>%
  arrange(batch_order, disease_name, entity_type, entity_name, source_id)

candidate_unique_sources_to_fetch <- candidate_source_request_list %>%
  group_by(
    batch_order,
    phase,
    batch_id,
    source_id,
    source_lookup_status,
    source_title,
    authors_or_organization,
    year,
    source_type,
    source_url,
    doi,
    pmid,
    pmcid,
    source_access
  ) %>%
  summarise(
    candidate_rows_to_check = n(),
    candidate_rows = paste(unique(candidate_row), collapse = " | "),
  .groups = "drop"
) %>%
  arrange(batch_order, source_lookup_status, source_id)

# Carry forward manually entered file-name metadata when the source list is
# regenerated, so source-fetch progress is not lost between consolidation runs.
if (file.exists(candidate_unique_sources_path)) {
  existing_unique_sources <- read_csv(
    candidate_unique_sources_path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    na = c("", "NA")
  )

  if ("file_name" %in% names(existing_unique_sources)) {
    file_name_lookup <- existing_unique_sources %>%
      transmute(
        batch_id,
        source_title,
        source_url,
        file_name
      ) %>%
      filter(!is.na(file_name), file_name != "") %>%
      distinct()

    candidate_unique_sources_to_fetch <- candidate_unique_sources_to_fetch %>%
      left_join(
        file_name_lookup,
        by = c("batch_id", "source_title", "source_url")
      ) %>%
      relocate(file_name, .after = doi)
  }
}

# ------------------------------------------------------------------------------|
#      Build review summaries --------------------------------------------------|
# ------------------------------------------------------------------------------|
deferred_review_queue <- all_assignments %>%
  filter(action %in% c(
    "defer_vocabulary_or_taxonomy",
    "defer_insufficient_evidence",
    "reject_or_negative_evidence_only"
  )) %>%
  arrange(review_priority, confidence_rank, batch_order, disease_name, entity_type, entity_name)

baseline_or_group_rows <- all_assignments %>%
  filter(action %in% c("already_covered", "evidence_only_group")) %>%
  arrange(action, batch_order, disease_name, entity_type, entity_name)

action_summary <- all_assignments %>%
  count(batch_order, phase, batch_id, disease_name, action, name = "assignment_rows") %>%
  arrange(batch_order, disease_name, action)

priority_summary <- all_assignments %>%
  count(review_priority, entity_type, name = "assignment_rows") %>%
  arrange(review_priority, entity_type)

batch_summary <- manifest %>%
  select(batch_order, phase, batch_id, table_rows, quality_issue_rows) %>%
  left_join(
    all_assignments %>%
      count(batch_id, name = "assignment_rows"),
    by = "batch_id"
  ) %>%
  left_join(
    candidate_source_check_queue %>%
      count(batch_id, name = "candidate_source_check_rows"),
    by = "batch_id"
  ) %>%
  mutate(
    assignment_rows = coalesce(assignment_rows, 0L),
    candidate_source_check_rows = coalesce(candidate_source_check_rows, 0L)
  )

# ------------------------------------------------------------------------------|
#      Write consolidated CSV outputs -----------------------------------------|
# ------------------------------------------------------------------------------|
write_csv(all_assignments, file.path(output_dir, "all_assignment_staging.csv"), na = "")
write_csv(all_host_evidence, file.path(output_dir, "all_host_role_evidence.csv"), na = "")
write_csv(all_vector_evidence, file.path(output_dir, "all_vector_role_evidence.csv"), na = "")
write_csv(all_sources, file.path(output_dir, "all_sources_used.csv"), na = "")
write_csv(all_quality_issues, file.path(output_dir, "all_quality_issues.csv"), na = "")
write_csv(all_candidate_presence, file.path(output_dir, "all_candidate_presence_check.csv"), na = "")
write_csv(all_vector_non_applicability, file.path(output_dir, "all_vector_non_applicability.csv"), na = "")
write_csv(all_cross_batch_summary, file.path(output_dir, "all_cross_batch_summary.csv"), na = "")
write_csv(candidate_source_check_queue, file.path(output_dir, "candidate_source_check_queue.csv"), na = "")
write_csv(candidate_source_request_list, file.path(output_dir, "candidate_source_request_list.csv"), na = "")
write_csv(candidate_unique_sources_to_fetch, candidate_unique_sources_path, na = "")
write_csv(deferred_review_queue, file.path(output_dir, "deferred_review_queue.csv"), na = "")
write_csv(baseline_or_group_rows, file.path(output_dir, "baseline_or_group_rows.csv"), na = "")
write_csv(action_summary, file.path(output_dir, "assignment_action_summary.csv"), na = "")
write_csv(priority_summary, file.path(output_dir, "assignment_priority_summary.csv"), na = "")
write_csv(batch_summary, file.path(output_dir, "batch_summary.csv"), na = "")

# ------------------------------------------------------------------------------|
#      Write human-readable handoff files -------------------------------------|
# ------------------------------------------------------------------------------|
top_candidate_lines <- candidate_source_check_queue %>%
  mutate(
    line = paste0(
      "- ", disease_name, ": `", entity_name, "` (", entity_type, ", `",
      role_assignment, "`, ", assignment_confidence, ")"
    )
  ) %>%
  pull(line)

if (length(top_candidate_lines) > 40) {
  top_candidate_lines <- c(
    top_candidate_lines[seq_len(40)],
    paste0("- ...and ", length(top_candidate_lines) - 40, " more rows in `candidate_source_check_queue.csv`.")
  )
}

readme_lines <- c(
  "# Consolidated Deep Research Review",
  "",
  "This folder consolidates all six reformatted Deep Research batch reports into reviewable staging files.",
  "",
  "These are staging artifacts only. They do not update official role evidence or assignment CSVs.",
  "",
  "## Scope",
  "",
  paste0("- Reports consolidated: ", nrow(manifest)),
  paste0("- Extracted table rows across reports: ", sum(manifest$table_rows, na.rm = TRUE)),
  paste0("- Parser/QA issue rows across reports: ", sum(manifest$quality_issue_rows, na.rm = TRUE)),
  paste0("- Assignment staging rows: ", nrow(all_assignments)),
  paste0("- Candidate source-check rows: ", nrow(candidate_source_check_queue)),
  paste0("- Candidate source links to check: ", nrow(candidate_source_request_list)),
  paste0("- Unique candidate sources to fetch: ", nrow(candidate_unique_sources_to_fetch)),
  paste0("- Deferred/taxonomy/negative guardrail rows: ", nrow(deferred_review_queue)),
  "",
  "## Core Files",
  "",
  "- `candidate_source_check_queue.csv`: best next working file; source-check these before importing anything.",
  "- `candidate_source_request_list.csv`: candidate rows joined to their source IDs, URLs, DOIs, PMIDs, and PMCIDs.",
  "- `candidate_unique_sources_to_fetch.csv`: de-duplicated source list for fetching PDFs/pages once.",
  "- `all_assignment_staging.csv`: all assignment recommendations from all six reports, with `review_priority` added.",
  "- `deferred_review_queue.csv`: taxonomy/vocabulary, insufficient-evidence, and negative/unsupported rows.",
  "- `baseline_or_group_rows.csv`: already-covered and evidence-only group rows.",
  "- `all_host_role_evidence.csv` and `all_vector_role_evidence.csv`: extracted evidence tables from the reports.",
  "- `all_sources_used.csv`: source tables from all reports; source metadata still needs verification.",
  "- `all_quality_issues.csv`: parser/schema issues inherited from batch-level QA.",
  "- `batch_summary.csv`, `assignment_action_summary.csv`, and `assignment_priority_summary.csv`: compact overview tables.",
  "",
  "## Suggested Next Step",
  "",
  "Start by fetching sources from `candidate_unique_sources_to_fetch.csv`. Then work through `candidate_source_check_queue.csv`. For each row:",
  "",
  "1. Open the listed source or DOI/PMID/PMCID.",
  "2. Verify the source supports the exact disease, entity, role, taxonomic grain, geography, and caveat.",
  "3. Decide whether to accept as official evidence, keep as evidence-only, defer, or reject.",
  "4. Only after evidence acceptance, add any matching assignment row and run the role-readiness QA script.",
  "",
  "## Candidate Rows To Source-Check First",
  "",
  if (length(top_candidate_lines) == 0) "- No candidate rows found." else top_candidate_lines,
  "",
  "Generated by:",
  "",
  "`Rscript scripts/associations/role_annotation/6_6_Consolidate_Deep_Research_Reports.R`"
)

writeLines(readme_lines, file.path(output_dir, "CONSOLIDATED_DEEP_RESEARCH_REVIEW.md"), useBytes = TRUE)

source_lines <- candidate_unique_sources_to_fetch %>%
  mutate(
    line = paste0(
      "- `", source_id, "` | ", coalesce(source_title, "MISSING SOURCE METADATA"),
      " | URL: ", coalesce(source_url, ""),
      " | DOI: ", coalesce(doi, ""),
      " | PMID: ", coalesce(pmid, ""),
      " | PMCID: ", coalesce(pmcid, ""),
      " | candidate rows: ", candidate_rows_to_check
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# Candidate Source Request List",
    "",
    "Use this as the human-readable list of sources to fetch/check for candidate role rows.",
    "",
    "The machine-readable versions are:",
    "",
    "- `candidate_unique_sources_to_fetch.csv`",
    "- `candidate_source_request_list.csv`",
    "",
    "## Sources",
    "",
    source_lines
  ),
  file.path(output_dir, "SOURCE_REQUEST_LIST.md"),
  useBytes = TRUE
)

message("Wrote consolidated Deep Research staging artifacts.")
print(batch_summary)
print(priority_summary)
