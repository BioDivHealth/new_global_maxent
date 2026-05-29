#!/usr/bin/env Rscript
################################################################################
# 6_9_Import_Source_Checked_Role_Rows.R
################################################################################
# Purpose: Promote source-checked Deep Research role rows into the official
#          host/vector role evidence and assignment CSVs.
#
# Guardrails:
# - Import only source-check rows with decision == "accept" and import_ready.
# - Do not mutate the source-check decision ledger.
# - Keep official CSV schemas unchanged.
# - Keep every imported row reviewable and source-traceable.
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
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

is_true <- function(x) {
  as.character(x) %in% c("TRUE", "true", "True", "1", "yes", "Yes", "YES")
}

read_stage_csv <- function(path) {
  readr::read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(across(where(is.character), clean_text))
}

normalize_key <- function(x) {
  x %>%
    clean_text() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[._-]+", " ") %>%
    stringr::str_squish()
}

collapse_unique <- function(x, sep = "; ") {
  values <- clean_text(x)
  values <- unique(values[!is.na(values)])
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(values, collapse = sep)
}

first_non_empty <- function(...) {
  args <- list(...)
  if (length(args) == 0) {
    return(character())
  }

  out <- rep(NA_character_, length(args[[1]]))
  for (arg in args) {
    values <- clean_text(arg)
    fill <- is.na(out) & !is.na(values)
    out[fill] <- values[fill]
  }
  out
}

contains_target <- function(x, target) {
  stringr::str_detect(coalesce(x, ""), stringr::fixed(target))
}

make_note <- function(candidate_row_id, caveat, decision_reason, source_access,
                      file_name, source_check_method, checked_evidence_location) {
  fields <- c(
    paste0("source_check_candidate_id=", candidate_row_id),
    if (!is.na(caveat)) paste0("caveat=", caveat) else NA_character_,
    if (!is.na(decision_reason)) paste0("decision_reason=", decision_reason) else NA_character_,
    if (!is.na(source_access)) paste0("source_access=", source_access) else NA_character_,
    if (!is.na(file_name)) paste0("file_name=", file_name) else NA_character_,
    if (!is.na(source_check_method)) paste0("source_check_method=", source_check_method) else NA_character_,
    if (!is.na(checked_evidence_location)) {
      paste0("checked_evidence_location=", checked_evidence_location)
    } else {
      NA_character_
    }
  )
  paste(fields[!is.na(fields)], collapse = "; ")
}

extract_source_check_ids <- function(x) {
  ids <- stringr::str_extract_all(coalesce(x, ""), "source_check_candidate_id=candidate_[0-9]+")
  ids <- unlist(ids, use.names = FALSE)
  unique(stringr::str_remove(ids, "^source_check_candidate_id="))
}

classify_evidence_rows <- function(staged, existing, key_cols) {
  existing_ids <- extract_source_check_ids(existing$notes)
  existing_keys <- existing %>%
    distinct(across(all_of(key_cols))) %>%
    mutate(natural_duplicate = TRUE)

  staged %>%
    mutate(already_imported = candidate_row_id %in% existing_ids) %>%
    left_join(existing_keys, by = key_cols) %>%
    mutate(
      natural_duplicate = coalesce(natural_duplicate, FALSE),
      import_action = case_when(
        already_imported ~ "skipped_existing_source_check_id",
        natural_duplicate ~ "skipped_natural_duplicate",
        TRUE ~ "imported"
      )
    )
}

classify_assignment_rows <- function(staged, existing, key_cols) {
  existing_ids <- unique(clean_text(existing$evidence_record_ids))
  existing_keys <- existing %>%
    distinct(across(all_of(key_cols))) %>%
    mutate(natural_duplicate = TRUE)

  staged %>%
    mutate(already_imported = evidence_record_ids %in% existing_ids) %>%
    left_join(existing_keys, by = key_cols) %>%
    mutate(
      natural_duplicate = coalesce(natural_duplicate, FALSE),
      import_action = case_when(
        already_imported ~ "skipped_existing_source_check_id",
        natural_duplicate ~ "skipped_natural_duplicate",
        !evidence_available_for_assignment ~ "skipped_evidence_not_available",
        TRUE ~ "imported"
      )
    )
}

write_csv_strict <- function(data, path, schema) {
  if (!identical(names(data), schema)) {
    stop(
      "Refusing to write schema-changed CSV: ", path, "\nExpected: ",
      paste(schema, collapse = ", "), "\nActual: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  readr::write_csv(data, path, na = "")
}

stop_if_missing_required <- function(data, cols, label) {
  missing <- cols[!cols %in% names(data)]
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# ------------------------------------------------------------------------------|
#      Paths and preflight -----------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_root <- here::here()
source(here::here("scripts", "associations", "working_inputs.R"))

role_dir <- role_annotation_dir
source_check_dir <- role_source_check_dir
import_dir <- role_source_check_import_dir

paths <- list(
  data_decisions = file.path(repo_root, "docs", "DATA_DECISIONS.md"),
  role_plan = file.path(repo_root, "ROLE_EVIDENCE_FULL_CURATION_PLAN.md"),
  role_readme = file.path(role_dir, "README.md"),
  decisions = file.path(source_check_dir, "candidate_source_check_decisions.csv"),
  source_request = file.path(source_check_dir, "candidate_source_request_list_with_files.csv"),
  host_evidence = file.path(role_evidence_dir, "host_role_evidence.csv"),
  vector_evidence = file.path(role_evidence_dir, "vector_role_evidence.csv"),
  host_assignments = file.path(role_evidence_dir, "host_role_assignments.csv"),
  vector_assignments = file.path(role_evidence_dir, "vector_role_assignments.csv"),
  host_candidates = file.path(role_candidates_dir, "host_role_candidates.csv"),
  vector_candidates = role_vector_candidate_path("who"),
  host_taxonomy = file.path(who_virion_dir, "who_host_species_standardized.csv")
)

required_paths <- unlist(paths[c(
  "data_decisions", "role_plan", "role_readme", "decisions", "source_request",
  "host_evidence", "vector_evidence", "host_assignments", "vector_assignments"
)]);
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required import inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

dir.create(import_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Inputs ------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
decisions <- read_stage_csv(paths$decisions)
source_request <- read_stage_csv(paths$source_request)
host_evidence <- read_stage_csv(paths$host_evidence)
vector_evidence <- read_stage_csv(paths$vector_evidence)
host_assignments <- read_stage_csv(paths$host_assignments)
vector_assignments <- read_stage_csv(paths$vector_assignments)

host_evidence_schema <- names(host_evidence)
vector_evidence_schema <- names(vector_evidence)
host_assignment_schema <- names(host_assignments)
vector_assignment_schema <- names(vector_assignments)

stop_if_missing_required(
  decisions,
  c(
    "candidate_row_id", "disease_name", "entity_type", "entity_name",
    "role_assignment", "assignment_confidence", "source_titles", "source_urls",
    "source_access", "file_name", "source_check_method", "checked_evidence_span",
    "checked_evidence_location", "decision", "accepted_role", "accepted_confidence",
    "accepted_evidence_scope", "caveat", "official_csv_target", "decision_reason",
    "import_ready"
  ),
  "candidate_source_check_decisions.csv"
)

if (!"file_name" %in% names(decisions)) {
  stop("Expected file_name column is missing from candidate_source_check_decisions.csv", call. = FALSE)
}

accepted <- decisions %>%
  mutate(import_ready_flag = is_true(import_ready)) %>%
  filter(decision == "accept", import_ready_flag) %>%
  select(-import_ready_flag)

excluded <- decisions %>%
  mutate(import_ready_flag = is_true(import_ready)) %>%
  filter(!(decision == "accept" & import_ready_flag)) %>%
  select(-import_ready_flag)

if (nrow(accepted) != 45) {
  stop("Expected 45 accepted import-ready rows, found ", nrow(accepted), call. = FALSE)
}

if (nrow(excluded) != 8) {
  stop("Expected 8 excluded rows, found ", nrow(excluded), call. = FALSE)
}

missing_targets <- accepted %>%
  filter(is.na(official_csv_target)) %>%
  pull(candidate_row_id)
if (length(missing_targets) > 0) {
  stop("Accepted rows lack official_csv_target: ", paste(missing_targets, collapse = ", "), call. = FALSE)
}

source_type_lookup <- source_request %>%
  group_by(batch_id, disease_name, entity_type, entity_name, role_assignment, assignment_confidence) %>%
  summarise(
    source_type_collapsed = collapse_unique(source_type),
    source_request_file_name = collapse_unique(file_name),
    .groups = "drop"
  )

accepted <- accepted %>%
  left_join(
    source_type_lookup,
    by = c("batch_id", "disease_name", "entity_type", "entity_name", "role_assignment", "assignment_confidence")
  ) %>%
  mutate(
    source_type_import = first_non_empty(source_type_collapsed, source_check_method, source_access),
    file_name = first_non_empty(file_name, source_request_file_name)
  )

# ------------------------------------------------------------------------------|
#      Metadata fallbacks ------------------------------------------------------|
# ------------------------------------------------------------------------------|
disease_pathogen_lookup <- bind_rows(
  host_evidence %>% transmute(disease_name, source_pathogen, network_pathogen),
  host_assignments %>% transmute(disease_name, source_pathogen, network_pathogen),
  vector_evidence %>% transmute(disease_name, source_pathogen, network_pathogen = source_pathogen),
  vector_assignments %>% transmute(disease_name, source_pathogen, network_pathogen = source_pathogen)
) %>%
  filter(!is.na(disease_name), !is.na(source_pathogen)) %>%
  group_by(disease_name) %>%
  summarise(
    disease_source_pathogen = collapse_unique(source_pathogen),
    disease_network_pathogen = collapse_unique(network_pathogen),
    .groups = "drop"
  )

host_candidate_lookup <- if (file.exists(paths$host_candidates)) {
  read_stage_csv(paths$host_candidates) %>%
    transmute(
      disease_name,
      host_key = normalize_key(host),
      candidate_source_pathogen = active_source_pathogens,
      candidate_network_pathogen = network_pathogen,
      candidate_host_tax_id = host_tax_id
    ) %>%
    group_by(disease_name, host_key) %>%
    summarise(
      candidate_source_pathogen = collapse_unique(candidate_source_pathogen),
      candidate_network_pathogen = collapse_unique(candidate_network_pathogen),
      candidate_host_tax_id = collapse_unique(candidate_host_tax_id),
      .groups = "drop"
    )
} else {
  tibble(
    disease_name = character(),
    host_key = character(),
    candidate_source_pathogen = character(),
    candidate_network_pathogen = character(),
    candidate_host_tax_id = character()
  )
}

host_taxonomy_lookup <- if (file.exists(paths$host_taxonomy)) {
  host_taxonomy <- read_stage_csv(paths$host_taxonomy)
  bind_rows(
    host_taxonomy %>% transmute(host_key = normalize_key(Host), taxonomy_host_tax_id = HostTaxID),
    host_taxonomy %>% transmute(host_key = normalize_key(correct_name), taxonomy_host_tax_id = HostTaxID),
    host_taxonomy %>% transmute(host_key = normalize_key(host_species), taxonomy_host_tax_id = HostTaxID)
  ) %>%
    filter(!is.na(host_key), !is.na(taxonomy_host_tax_id)) %>%
    group_by(host_key) %>%
    summarise(taxonomy_host_tax_id = collapse_unique(taxonomy_host_tax_id), .groups = "drop")
} else {
  tibble(host_key = character(), taxonomy_host_tax_id = character())
}

vector_candidate_lookup <- if (file.exists(paths$vector_candidates)) {
  read_stage_csv(paths$vector_candidates) %>%
    transmute(
      disease_name,
      vector_key = normalize_key(vector_species),
      candidate_vector_join_key = vector_join_key
    ) %>%
    group_by(disease_name, vector_key) %>%
    summarise(candidate_vector_join_key = collapse_unique(candidate_vector_join_key), .groups = "drop")
} else {
  tibble(disease_name = character(), vector_key = character(), candidate_vector_join_key = character())
}

# ------------------------------------------------------------------------------|
#      Stage official rows -----------------------------------------------------|
# ------------------------------------------------------------------------------|
host_source <- accepted %>%
  filter(contains_target(official_csv_target, "host_role_evidence.csv")) %>%
  mutate(
    host = entity_name,
    host_key = normalize_key(entity_name)
  ) %>%
  left_join(host_candidate_lookup, by = c("disease_name", "host_key")) %>%
  left_join(host_taxonomy_lookup, by = "host_key") %>%
  left_join(disease_pathogen_lookup, by = "disease_name") %>%
  mutate(
    source_pathogen = first_non_empty(candidate_source_pathogen, disease_source_pathogen),
    network_pathogen = first_non_empty(candidate_network_pathogen, disease_network_pathogen, source_pathogen),
    host_tax_id = first_non_empty(candidate_host_tax_id, taxonomy_host_tax_id)
  )

vector_source <- accepted %>%
  filter(contains_target(official_csv_target, "vector_role_evidence.csv")) %>%
  mutate(
    vector_species = entity_name,
    vector_key = normalize_key(entity_name)
  ) %>%
  left_join(vector_candidate_lookup, by = c("disease_name", "vector_key")) %>%
  left_join(disease_pathogen_lookup, by = "disease_name") %>%
  mutate(
    source_pathogen = first_non_empty(disease_source_pathogen, disease_network_pathogen),
    vector_join_key = first_non_empty(candidate_vector_join_key, vector_key)
  )

required_host_values <- c("source_pathogen", "network_pathogen", "host", "accepted_role", "checked_evidence_span")
host_missing <- host_source %>%
  filter(if_any(all_of(required_host_values), is.na)) %>%
  pull(candidate_row_id)
if (length(host_missing) > 0) {
  stop("Host accepted rows still lack required import metadata: ", paste(host_missing, collapse = ", "), call. = FALSE)
}

required_vector_values <- c("source_pathogen", "vector_species", "vector_join_key", "accepted_role", "checked_evidence_span")
vector_missing <- vector_source %>%
  filter(if_any(all_of(required_vector_values), is.na)) %>%
  pull(candidate_row_id)
if (length(vector_missing) > 0) {
  stop("Vector accepted rows still lack required import metadata: ", paste(vector_missing, collapse = ", "), call. = FALSE)
}

host_evidence_staged <- host_source %>%
  transmute(
    candidate_row_id,
    disease_name,
    source_pathogen,
    network_pathogen,
    host,
    host_tax_id,
    role_claim = accepted_role,
    evidence_type = accepted_evidence_scope,
    evidence_direction = "supports",
    evidence_span = checked_evidence_span,
    source_citation = source_titles,
    source_url = source_urls,
    source_type = source_type_import,
    claim_confidence = accepted_confidence,
    needs_manual_review = "TRUE",
    notes = pmap_chr(
      list(candidate_row_id, caveat, decision_reason, source_access, file_name, source_check_method, checked_evidence_location),
      make_note
    )
  )

vector_evidence_staged <- vector_source %>%
  transmute(
    candidate_row_id,
    disease_name,
    source_pathogen,
    vector_species,
    vector_join_key,
    vector_role_claim = accepted_role,
    evidence_type = accepted_evidence_scope,
    evidence_direction = "supports",
    evidence_span = checked_evidence_span,
    source_citation = source_titles,
    source_url = source_urls,
    source_type = source_type_import,
    claim_confidence = accepted_confidence,
    needs_manual_review = "TRUE",
    notes = pmap_chr(
      list(candidate_row_id, caveat, decision_reason, source_access, file_name, source_check_method, checked_evidence_location),
      make_note
    )
  )

host_evidence_classified <- classify_evidence_rows(
  host_evidence_staged,
  host_evidence,
  c("disease_name", "host", "role_claim", "source_citation", "evidence_span")
)

vector_evidence_classified <- classify_evidence_rows(
  vector_evidence_staged,
  vector_evidence,
  c("disease_name", "vector_species", "vector_role_claim", "source_citation", "evidence_span")
)

host_imported_ids <- host_evidence_classified %>%
  filter(import_action %in% c("imported", "skipped_existing_source_check_id")) %>%
  pull(candidate_row_id)

vector_imported_ids <- vector_evidence_classified %>%
  filter(import_action %in% c("imported", "skipped_existing_source_check_id")) %>%
  pull(candidate_row_id)

host_assignments_staged <- host_source %>%
  transmute(
    candidate_row_id,
    disease_name,
    source_pathogen,
    network_pathogen,
    host,
    host_tax_id,
    host_role_assignment = accepted_role,
    assignment_status = "draft_source_backed",
    assignment_confidence = accepted_confidence,
    evidence_record_ids = paste0("host_role_evidence:source_check:", candidate_row_id),
    assignment_basis = decision_reason,
    needs_manual_review = "TRUE",
    review_notes = pmap_chr(
      list(candidate_row_id, caveat, decision_reason, source_access, file_name, source_check_method, checked_evidence_location),
      make_note
    ),
    evidence_available_for_assignment = candidate_row_id %in% host_imported_ids
  )

vector_assignments_staged <- vector_source %>%
  transmute(
    candidate_row_id,
    disease_name,
    source_pathogen,
    vector_species,
    vector_join_key,
    vector_role_assignment = accepted_role,
    assignment_status = "draft_source_backed",
    assignment_confidence = accepted_confidence,
    evidence_record_ids = paste0("vector_role_evidence:source_check:", candidate_row_id),
    assignment_basis = decision_reason,
    needs_manual_review = "TRUE",
    review_notes = pmap_chr(
      list(candidate_row_id, caveat, decision_reason, source_access, file_name, source_check_method, checked_evidence_location),
      make_note
    ),
    evidence_available_for_assignment = candidate_row_id %in% vector_imported_ids
  )

host_assignments_classified <- classify_assignment_rows(
  host_assignments_staged,
  host_assignments,
  c("disease_name", "host", "host_role_assignment", "evidence_record_ids")
)

vector_assignments_classified <- classify_assignment_rows(
  vector_assignments_staged,
  vector_assignments,
  c("disease_name", "vector_species", "vector_role_assignment", "evidence_record_ids")
)

# ------------------------------------------------------------------------------|
#      Write official CSVs -----------------------------------------------------|
# ------------------------------------------------------------------------------|
host_evidence_out <- bind_rows(
  host_evidence,
  host_evidence_classified %>%
    filter(import_action == "imported") %>%
    select(all_of(host_evidence_schema))
)

vector_evidence_out <- bind_rows(
  vector_evidence,
  vector_evidence_classified %>%
    filter(import_action == "imported") %>%
    select(all_of(vector_evidence_schema))
)

host_assignments_out <- bind_rows(
  host_assignments,
  host_assignments_classified %>%
    filter(import_action == "imported") %>%
    select(all_of(host_assignment_schema))
)

vector_assignments_out <- bind_rows(
  vector_assignments,
  vector_assignments_classified %>%
    filter(import_action == "imported") %>%
    select(all_of(vector_assignment_schema))
)

write_csv_strict(host_evidence_out, paths$host_evidence, host_evidence_schema)
write_csv_strict(vector_evidence_out, paths$vector_evidence, vector_evidence_schema)
write_csv_strict(host_assignments_out, paths$host_assignments, host_assignment_schema)
write_csv_strict(vector_assignments_out, paths$vector_assignments, vector_assignment_schema)

# ------------------------------------------------------------------------------|
#      Audit outputs -----------------------------------------------------------|
# ------------------------------------------------------------------------------|
role_vocab_summary <- bind_rows(
  tibble(
    table = "host_role_evidence",
    role = sort(unique(host_evidence_classified$role_claim)),
    was_present_before_import = role %in% unique(host_evidence$role_claim)
  ),
  tibble(
    table = "vector_role_evidence",
    role = sort(unique(vector_evidence_classified$vector_role_claim)),
    was_present_before_import = role %in% unique(vector_evidence$vector_role_claim)
  )
) %>%
  mutate(new_official_csv_role_value = !was_present_before_import)

new_role_lines <- role_vocab_summary %>%
  filter(new_official_csv_role_value) %>%
  transmute(line = paste0("- `", table, "`: `", role, "`")) %>%
  pull(line)

if (length(new_role_lines) == 0) {
  new_role_lines <- "- none"
}

row_actions <- bind_rows(
  host_evidence_classified %>%
    transmute(table = "host_role_evidence", candidate_row_id, disease_name, entity = host, role = role_claim, import_action),
  vector_evidence_classified %>%
    transmute(table = "vector_role_evidence", candidate_row_id, disease_name, entity = vector_species, role = vector_role_claim, import_action),
  host_assignments_classified %>%
    transmute(table = "host_role_assignments", candidate_row_id, disease_name, entity = host, role = host_role_assignment, import_action),
  vector_assignments_classified %>%
    transmute(table = "vector_role_assignments", candidate_row_id, disease_name, entity = vector_species, role = vector_role_assignment, import_action)
)

imported_rows <- row_actions %>% filter(import_action == "imported")

imported_row_summary <- row_actions %>%
  count(table, import_action, name = "rows") %>%
  arrange(table, import_action)

skipped_duplicate_summary <- row_actions %>%
  filter(import_action != "imported") %>%
  count(table, import_action, name = "rows") %>%
  arrange(table, import_action)

excluded_non_import_rows <- excluded %>%
  select(
    candidate_row_id,
    disease_name,
    entity_type,
    entity_name,
    role_assignment,
    decision,
    import_ready,
    official_csv_target,
    decision_reason,
    caveat
  )

readr::write_csv(imported_rows, file.path(import_dir, "imported_rows.csv"), na = "")
readr::write_csv(imported_row_summary, file.path(import_dir, "imported_row_summary.csv"), na = "")
readr::write_csv(skipped_duplicate_summary, file.path(import_dir, "skipped_duplicate_summary.csv"), na = "")
readr::write_csv(excluded_non_import_rows, file.path(import_dir, "excluded_non_import_rows.csv"), na = "")
readr::write_csv(role_vocab_summary, file.path(import_dir, "new_role_vocabulary_summary.csv"), na = "")

accepted_accounting <- row_actions %>%
  filter(table %in% c("host_role_evidence", "vector_role_evidence")) %>%
  count(import_action, name = "candidate_rows") %>%
  mutate(total_accepted_rows = sum(candidate_rows))

accepted_accounting_lines <- accepted_accounting %>%
  transmute(
    line = paste0(
      "- ", import_action, ": ", candidate_rows,
      " accepted evidence rows; total accepted rows: ", total_accepted_rows
    )
  ) %>%
  pull(line)

import_action_lines <- imported_row_summary %>%
  transmute(line = paste0("- `", table, "` ", import_action, ": ", rows)) %>%
  pull(line)

new_host_evidence <- nrow(host_evidence_out) - nrow(host_evidence)
new_vector_evidence <- nrow(vector_evidence_out) - nrow(vector_evidence)
new_host_assignments <- nrow(host_assignments_out) - nrow(host_assignments)
new_vector_assignments <- nrow(vector_assignments_out) - nrow(vector_assignments)

summary_lines <- c(
  "# Source-Checked Role Import Summary",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Guardrails",
  "",
  "- `docs/DATA_DECISIONS.md`, `ROLE_EVIDENCE_FULL_CURATION_PLAN.md`, and `role_annotation/README.md` were required preflight inputs.",
  "- Only `decision == \"accept\"` and `import_ready == TRUE` rows were considered for import.",
  "- Evidence-only and deferred rows were excluded.",
  "- Official role CSV schemas were preserved.",
  "- Imported rows are marked `needs_manual_review = TRUE`.",
  "",
  "## Official Row Deltas",
  "",
  paste0("- `host_role_evidence.csv`: +", new_host_evidence),
  paste0("- `vector_role_evidence.csv`: +", new_vector_evidence),
  paste0("- `host_role_assignments.csv`: +", new_host_assignments),
  paste0("- `vector_role_assignments.csv`: +", new_vector_assignments),
  "",
  "## Accepted Candidate Accounting",
  "",
  paste(accepted_accounting_lines, collapse = "\n"),
  "",
  "## Import Actions",
  "",
  paste(import_action_lines, collapse = "\n"),
  "",
  "## New Role Values",
  "",
  paste(new_role_lines, collapse = "\n"),
  "",
  "## Excluded Rows",
  "",
  paste0("- Excluded candidate rows: ", nrow(excluded_non_import_rows)),
  "- See `excluded_non_import_rows.csv` for evidence-only and deferred rows."
)

writeLines(summary_lines, file.path(import_dir, "SOURCE_CHECK_IMPORT_SUMMARY.md"), useBytes = TRUE)

# ------------------------------------------------------------------------------|
#      Validation --------------------------------------------------------------|
# ------------------------------------------------------------------------------|
if (sum(accepted_accounting$candidate_rows) != 45) {
  stop("Accepted candidate accounting did not total 45 rows.", call. = FALSE)
}

if (nrow(excluded_non_import_rows) != 8) {
  stop("Excluded-row accounting did not total 8 rows.", call. = FALSE)
}

host_assignment_join <- host_assignments_out %>%
  filter(stringr::str_detect(coalesce(evidence_record_ids, ""), "host_role_evidence:source_check:")) %>%
  left_join(
    host_evidence_out %>%
      transmute(
        disease_name,
        host,
        host_role_assignment = role_claim,
        evidence_present = TRUE
      ) %>%
      distinct(),
    by = c("disease_name", "host", "host_role_assignment")
  )

vector_assignment_join <- vector_assignments_out %>%
  filter(stringr::str_detect(coalesce(evidence_record_ids, ""), "vector_role_evidence:source_check:")) %>%
  left_join(
    vector_evidence_out %>%
      transmute(
        disease_name,
        vector_species,
        vector_role_assignment = vector_role_claim,
        evidence_present = TRUE
      ) %>%
      distinct(),
    by = c("disease_name", "vector_species", "vector_role_assignment")
  )

if (any(is.na(host_assignment_join$evidence_present))) {
  bad <- host_assignment_join %>%
    filter(is.na(evidence_present)) %>%
    pull(evidence_record_ids)
  stop("Host source-check assignments lack matching evidence: ", paste(bad, collapse = ", "), call. = FALSE)
}

if (any(is.na(vector_assignment_join$evidence_present))) {
  bad <- vector_assignment_join %>%
    filter(is.na(evidence_present)) %>%
    pull(evidence_record_ids)
  stop("Vector source-check assignments lack matching evidence: ", paste(bad, collapse = ", "), call. = FALSE)
}

message("Source-checked role import complete.")
message("Official row deltas:")
message("  host_role_evidence.csv: +", new_host_evidence)
message("  vector_role_evidence.csv: +", new_vector_evidence)
message("  host_role_assignments.csv: +", new_host_assignments)
message("  vector_role_assignments.csv: +", new_vector_assignments)
message("Import summary: ", file.path(import_dir, "SOURCE_CHECK_IMPORT_SUMMARY.md"))
