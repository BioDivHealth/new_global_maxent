#!/usr/bin/env Rscript
################################################################################
# 6_3_Build_Evidence_Readiness_QA.R
################################################################################
# Purpose: Build compact QA tables that make the current host, vector, GenBank,
#          WHO DON, and role-review evidence layers easier to inspect.
#
# Inputs : generated role_annotation, vector_screening, genbank_simple, and
#          disease_outbreak_news_v2 CSVs where present.
#
# Outputs: pathogen_association_data/evidence/role_annotation/qa/
#            evidence_layer_inventory.csv
#            disease_evidence_readiness.csv
#            vector_evidence_readiness_by_disease.csv
#            role_evidence_gap_summary.csv
#            role_evidence_claim_summary.csv
#            role_assignment_claim_summary.csv
#
# Notes  : These are QA/readiness summaries only. They do not assign biological
#          host or vector roles.
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

source(here::here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
}

read_csv_layer <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text))
}

prefer_existing_path <- function(primary, fallback) {
  if (file.exists(primary) || !file.exists(fallback)) {
    return(primary)
  }

  fallback
}

row_count <- function(x) {
  if (is.null(x)) {
    return(NA_integer_)
  }

  nrow(x)
}

col_count <- function(x) {
  if (is.null(x)) {
    return(NA_integer_)
  }

  ncol(x)
}

empty_disease_summary <- function() {
  tibble(disease_name = character())
}

full_join_disease_summaries <- function(tables) {
  tables <- purrr::compact(tables)

  if (length(tables) == 0) {
    return(empty_disease_summary())
  }

  purrr::reduce(tables, full_join, by = "disease_name")
}

count_non_missing_any <- function(data, cols) {
  present_cols <- intersect(cols, names(data))

  if (length(present_cols) == 0 || nrow(data) == 0) {
    return(rep(FALSE, nrow(data)))
  }

  rowSums(!is.na(data[present_cols])) > 0
}

sum_column <- function(data, col) {
  if (is.null(data) || !col %in% names(data)) {
    return(0)
  }

  sum(data[[col]], na.rm = TRUE)
}

summarise_role_evidence <- function(data, disease_col, signal_cols, prefix) {
  if (is.null(data) || !disease_col %in% names(data)) {
    return(empty_disease_summary())
  }

  data %>%
    mutate(
      disease_name = .data[[disease_col]],
      has_substantive_evidence = count_non_missing_any(pick(everything()), signal_cols)
    ) %>%
    filter(!is.na(disease_name)) %>%
    group_by(disease_name) %>%
    summarise(
      "{prefix}_rows" := dplyr::n(),
      "{prefix}_substantive_rows" := sum(has_substantive_evidence, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_assignments <- function(data, disease_col, assignment_col, prefix) {
  if (is.null(data) || !all(c(disease_col, assignment_col) %in% names(data))) {
    return(empty_disease_summary())
  }

  data %>%
    mutate(
      disease_name = .data[[disease_col]],
      has_assignment = !is.na(.data[[assignment_col]])
    ) %>%
    filter(!is.na(disease_name)) %>%
    group_by(disease_name) %>%
    summarise(
      "{prefix}_assignment_rows" := dplyr::n(),
      "{prefix}_substantive_assignment_rows" := sum(has_assignment, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
who_dir <- who_data_dir
role_dir <- role_annotation_dir
qa_dir <- role_qa_dir
vector_output_dir <- file.path(who_dir, "vector_screening", "outputs")
genbank_legacy_dir <- genbank_simple_legacy_dir
don_dir <- file.path(who_dir, "disease_outbreak_news_v2")

dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

genbank_summary_path <- prefer_existing_path(
  file.path(genbank_simple_evidence_dir, "genbank_readiness_disease_country_summary_standardized.csv"),
  prefer_existing_path(
    file.path(genbank_simple_evidence_dir, "genbank_disease_country_summary_standardized.csv"),
    prefer_existing_path(
      file.path(genbank_legacy_dir, "genbank_readiness_disease_country_summary_standardized.csv"),
      file.path(genbank_legacy_dir, "genbank_disease_country_summary_standardized.csv")
    )
  )
)
genbank_map_unmatched_path <- prefer_existing_path(
  file.path(genbank_simple_readiness_maps_dir, "genbank_disease_country_map_unmatched.csv"),
  prefer_existing_path(
    file.path(genbank_simple_standard_maps_dir, "genbank_disease_country_map_unmatched.csv"),
    prefer_existing_path(
      file.path(genbank_legacy_dir, "maps_readiness", "genbank_disease_country_map_unmatched.csv"),
      file.path(genbank_legacy_dir, "maps", "genbank_disease_country_map_unmatched.csv")
    )
  )
)

paths <- list(
  host_role_candidates = file.path(role_dir, "host_role_candidates.csv"),
  species_host_vector_roster = file.path(role_dir, "species_host_vector_roster.csv"),
  vector_competence_unmatched_review = file.path(
    vector_output_dir, "vector_competence_join_unmatched_review.csv"
  ),
  disease_vector_competence_annotated = file.path(
    vector_output_dir, "disease_vector_links_taxonomy_cleaned_competence_annotated.csv"
  ),
  genbank_disease_country_summary = genbank_summary_path,
  genbank_map_unmatched = genbank_map_unmatched_path,
  who_don_focal_modelling_ready = file.path(
    don_dir, "final", "who_don_modelling_ready.csv"
  ),
  host_role_evidence = file.path(role_dir, "host_role_evidence.csv"),
  vector_role_evidence = file.path(role_dir, "vector_role_evidence.csv"),
  host_role_assignments = file.path(role_dir, "host_role_assignments.csv"),
  vector_role_assignments = file.path(role_dir, "vector_role_assignments.csv")
)

paths_chr <- unlist(paths, use.names = FALSE)
layers <- purrr::map(paths, read_csv_layer)

repo_relative_path <- function(path) {
  repo_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  repo_prefix <- paste0(repo_root, "/")
  is_repo_path <- startsWith(normalized, repo_prefix)

  normalized[is_repo_path] <- substring(normalized[is_repo_path], nchar(repo_prefix) + 1L)
  normalized
}

# ------------------------------------------------------------------------------|
#      Inventory ---------------------------------------------------------------|
# ------------------------------------------------------------------------------|
inventory <- tibble(
  layer = names(paths),
  path = repo_relative_path(paths_chr),
  exists = file.exists(paths_chr),
  row_count = purrr::map_int(layers, row_count),
  column_count = purrr::map_int(layers, col_count)
) %>%
  mutate(
    status = case_when(
      !exists ~ "missing",
      row_count == 0 ~ "present_empty",
      TRUE ~ "present"
    ),
    notes = case_when(
      layer == "host_role_evidence" & row_count <= 1 ~
        "Currently only placeholder or near-empty evidence rows.",
      layer == "vector_role_evidence" & row_count == 0 ~
        "No vector role evidence rows yet.",
      layer %in% c("host_role_assignments", "vector_role_assignments") & row_count == 0 ~
        "No reviewed role assignments yet.",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------|
#      Disease-level readiness -------------------------------------------------|
# ------------------------------------------------------------------------------|
host_candidates <- layers$host_role_candidates
roster <- layers$species_host_vector_roster
vector_annotated <- layers$disease_vector_competence_annotated
genbank <- layers$genbank_disease_country_summary
genbank_unmatched <- layers$genbank_map_unmatched
who_don <- layers$who_don_focal_modelling_ready

host_summary <- if (is.null(host_candidates)) {
  empty_disease_summary()
} else {
  host_candidates %>%
    group_by(disease_name) %>%
    summarise(
      host_candidate_rows = dplyr::n(),
      human_host_candidate_rows = sum(is_true(is_human), na.rm = TRUE),
      livestock_host_candidate_rows = sum(is_true(is_livestock_like), na.rm = TRUE),
      bird_host_candidate_rows = sum(is_true(is_bird), na.rm = TRUE),
      rodent_host_candidate_rows = sum(is_true(is_rodent), na.rm = TRUE),
      bat_host_candidate_rows = sum(is_true(is_bat), na.rm = TRUE),
      .groups = "drop"
    )
}

roster_summary <- if (is.null(roster)) {
  empty_disease_summary()
} else {
  roster %>%
    group_by(disease_name) %>%
    summarise(
      roster_rows = dplyr::n(),
      roster_host_rows = sum(species_role == "host", na.rm = TRUE),
      roster_vector_rows = sum(species_role == "vector", na.rm = TRUE),
      roster_vectors_with_host_vector_evidence = sum(is_true(has_host_vector_evidence), na.rm = TRUE),
      roster_vectors_with_competence = sum(is_true(has_competence_evidence), na.rm = TRUE),
      roster_vectors_with_human_biting = sum(is_true(bites_humans), na.rm = TRUE),
      roster_taxonomy_caution_rows = sum(is_true(taxonomy_caution), na.rm = TRUE),
      .groups = "drop"
    )
}

vector_summary <- if (is.null(vector_annotated)) {
  empty_disease_summary()
} else {
  vector_annotated %>%
    group_by(disease_name) %>%
    summarise(
      disease_vector_rows = dplyr::n(),
      disease_vector_confirmed_rows = sum(best_evidence_level == "confirmed", na.rm = TRUE),
      disease_vector_probable_rows = sum(best_evidence_level == "probable", na.rm = TRUE),
      disease_vector_candidate_rows = sum(best_evidence_level == "candidate", na.rm = TRUE),
      vector_competence_joined_rows = sum(!is.na(vector_competence_status), na.rm = TRUE),
      vector_competence_transmission_yes_rows = sum(transmission_demonstrated == "yes", na.rm = TRUE),
      vector_competence_natural_infection_yes_rows = sum(natural_infection_reported == "yes", na.rm = TRUE),
      .groups = "drop"
    )
}

genbank_summary <- if (is.null(genbank)) {
  empty_disease_summary()
} else {
  genbank %>%
    group_by(disease_name = Disease_name) %>%
    summarise(
      genbank_country_rows = dplyr::n(),
      genbank_distinct_countries_or_territories = dplyr::n_distinct(country_standardized),
      genbank_records_with_country = sum(records_with_country, na.rm = TRUE),
      .groups = "drop"
    )
}

genbank_unmatched_summary <- if (is.null(genbank_unmatched)) {
  empty_disease_summary()
} else {
  genbank_unmatched %>%
    group_by(disease_name = Disease_name) %>%
    summarise(
      genbank_unmatched_map_rows = dplyr::n(),
      genbank_unmatched_records = sum(records_with_country, na.rm = TRUE),
      .groups = "drop"
    )
}

who_don_summary <- if (is.null(who_don)) {
  empty_disease_summary()
} else {
  disease_col <- dplyr::case_when(
    "disease_label_standard" %in% names(who_don) ~ "disease_label_standard",
    "disease_label_standard_refined" %in% names(who_don) ~ "disease_label_standard_refined",
    TRUE ~ NA_character_
  )
  confidence_col <- dplyr::case_when(
    "scope_confidence" %in% names(who_don) ~ "scope_confidence",
    "final_event_confidence" %in% names(who_don) ~ "final_event_confidence",
    TRUE ~ NA_character_
  )
  review_col <- dplyr::case_when(
    "needs_review" %in% names(who_don) ~ "needs_review",
    "final_needs_manual_review" %in% names(who_don) ~ "final_needs_manual_review",
    TRUE ~ NA_character_
  )

  if (is.na(disease_col)) {
    empty_disease_summary()
  } else {
    who_don %>%
      mutate(
        disease_name = .data[[disease_col]],
        confidence_value = if (!is.na(confidence_col)) .data[[confidence_col]] else NA_character_,
        needs_review_value = if (!is.na(review_col)) .data[[review_col]] else FALSE
      ) %>%
    filter(!is.na(disease_name)) %>%
    group_by(disease_name) %>%
    summarise(
      who_don_focal_country_rows = dplyr::n(),
      who_don_distinct_records = dplyr::n_distinct(record_key),
      who_don_distinct_countries = dplyr::n_distinct(country_standard),
        who_don_high_confidence_rows = sum(confidence_value == "high", na.rm = TRUE),
        who_don_manual_review_rows = sum(is_true(needs_review_value), na.rm = TRUE),
      .groups = "drop"
    )
  }
}

host_evidence_summary <- summarise_role_evidence(
  layers$host_role_evidence,
  "disease_name",
  c("role_claim", "evidence_span", "source_citation", "source_url"),
  "host_role_evidence"
)

vector_evidence_summary <- summarise_role_evidence(
  layers$vector_role_evidence,
  "disease_name",
  c("vector_role_claim", "evidence_span", "source_citation", "source_url"),
  "vector_role_evidence"
)

host_assignment_summary <- summarise_assignments(
  layers$host_role_assignments,
  "disease_name",
  "host_role_assignment",
  "host_role"
)

vector_assignment_summary <- summarise_assignments(
  layers$vector_role_assignments,
  "disease_name",
  "vector_role_assignment",
  "vector_role"
)

disease_readiness <- full_join_disease_summaries(
  list(
    host_summary,
    roster_summary,
    vector_summary,
    genbank_summary,
    genbank_unmatched_summary,
    who_don_summary,
    host_evidence_summary,
    vector_evidence_summary,
    host_assignment_summary,
    vector_assignment_summary
  )
)
readiness_numeric_cols <- c(
  "host_candidate_rows",
  "human_host_candidate_rows",
  "livestock_host_candidate_rows",
  "bird_host_candidate_rows",
  "rodent_host_candidate_rows",
  "bat_host_candidate_rows",
  "roster_rows",
  "roster_host_rows",
  "roster_vector_rows",
  "roster_vectors_with_host_vector_evidence",
  "roster_vectors_with_competence",
  "roster_vectors_with_human_biting",
  "roster_taxonomy_caution_rows",
  "disease_vector_rows",
  "disease_vector_confirmed_rows",
  "disease_vector_probable_rows",
  "disease_vector_candidate_rows",
  "vector_competence_joined_rows",
  "vector_competence_transmission_yes_rows",
  "vector_competence_natural_infection_yes_rows",
  "genbank_country_rows",
  "genbank_distinct_countries_or_territories",
  "genbank_records_with_country",
  "genbank_unmatched_map_rows",
  "genbank_unmatched_records",
  "who_don_focal_country_rows",
  "who_don_distinct_records",
  "who_don_distinct_countries",
  "who_don_high_confidence_rows",
  "who_don_manual_review_rows",
  "host_role_evidence_rows",
  "host_role_evidence_substantive_rows",
  "vector_role_evidence_rows",
  "vector_role_evidence_substantive_rows",
  "host_role_assignment_rows",
  "host_role_substantive_assignment_rows",
  "vector_role_assignment_rows",
  "vector_role_substantive_assignment_rows"
)

for (col in readiness_numeric_cols) {
  if (!col %in% names(disease_readiness)) {
    disease_readiness[[col]] <- 0
  }
}

disease_readiness <- disease_readiness %>%
  mutate(
    across(where(is.numeric), ~ tidyr::replace_na(.x, 0)),
    has_host_candidates = host_candidate_rows > 0,
    has_vector_rows = roster_vector_rows > 0 | disease_vector_rows > 0,
    has_competence_joined = vector_competence_joined_rows > 0 | roster_vectors_with_competence > 0,
    has_genbank_country_evidence = genbank_country_rows > 0,
    has_who_don_focal_evidence = who_don_focal_country_rows > 0,
    has_substantive_host_role_evidence = host_role_evidence_substantive_rows > 0,
    has_substantive_vector_role_evidence = vector_role_evidence_substantive_rows > 0,
    has_any_substantive_role_evidence =
      has_substantive_host_role_evidence | has_substantive_vector_role_evidence,
    role_assignment_status = case_when(
      host_role_substantive_assignment_rows > 0 | vector_role_substantive_assignment_rows > 0 ~
        "reviewed_assignment_present",
      has_any_substantive_role_evidence ~ "evidence_present_no_assignment",
      TRUE ~ "no_role_evidence_or_assignment"
    )
  ) %>%
  arrange(desc(has_vector_rows), disease_name)

# ------------------------------------------------------------------------------|
#      Vector readiness --------------------------------------------------------|
# ------------------------------------------------------------------------------|
roster_vector_support_summary <- if (all(
  c(
    "disease_name",
    "roster_vectors_with_host_vector_evidence",
    "roster_vectors_with_human_biting"
  ) %in% names(roster_summary)
)) {
  roster_summary %>%
    select(
      disease_name,
      roster_vectors_with_host_vector_evidence,
      roster_vectors_with_human_biting
    )
} else {
  tibble(
    disease_name = character(),
    roster_vectors_with_host_vector_evidence = numeric(),
    roster_vectors_with_human_biting = numeric()
  )
}

vector_readiness <- if (is.null(vector_annotated)) {
  empty_disease_summary()
} else {
  vector_annotated %>%
    group_by(disease_name) %>%
    summarise(
      vector_rows = dplyr::n(),
      vector_species_or_taxa = dplyr::n_distinct(vector_species_taxonomy_cleaned),
      confirmed_rows = sum(best_evidence_level == "confirmed", na.rm = TRUE),
      probable_rows = sum(best_evidence_level == "probable", na.rm = TRUE),
      candidate_rows = sum(best_evidence_level == "candidate", na.rm = TRUE),
      competence_joined_rows = sum(!is.na(vector_competence_status), na.rm = TRUE),
      competent_rows = sum(vector_competence_status == "competent", na.rm = TRUE),
      mixed_competence_rows = sum(vector_competence_status == "mixed", na.rm = TRUE),
      not_competent_rows = sum(vector_competence_status == "not_competent", na.rm = TRUE),
      unclear_competence_rows = sum(vector_competence_status == "unclear", na.rm = TRUE),
      transmission_yes_rows = sum(transmission_demonstrated == "yes", na.rm = TRUE),
      natural_infection_yes_rows = sum(natural_infection_reported == "yes", na.rm = TRUE),
      taxonomy_review_needed_rows = sum(is_true(review_needed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      roster_vector_support_summary,
      by = "disease_name"
    ) %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) %>%
    arrange(desc(vector_rows), disease_name)
}

# ------------------------------------------------------------------------------|
#      Evidence gap summary ----------------------------------------------------|
# ------------------------------------------------------------------------------|
role_gap_summary <- tibble(
  table = c(
    "host_role_evidence",
    "vector_role_evidence",
    "host_role_assignments",
    "vector_role_assignments"
  ),
  total_rows = c(
    row_count(layers$host_role_evidence),
    row_count(layers$vector_role_evidence),
    row_count(layers$host_role_assignments),
    row_count(layers$vector_role_assignments)
  ),
  substantive_rows = c(
    sum_column(host_evidence_summary, "host_role_evidence_substantive_rows"),
    sum_column(vector_evidence_summary, "vector_role_evidence_substantive_rows"),
    sum_column(host_assignment_summary, "host_role_substantive_assignment_rows"),
    sum_column(vector_assignment_summary, "vector_role_substantive_assignment_rows")
  )
) %>%
  mutate(
    status = case_when(
      is.na(total_rows) ~ "missing",
      substantive_rows > 0 ~ "substantive_rows_present",
      total_rows > 0 ~ "placeholder_or_empty_rows_only",
      TRUE ~ "empty"
    )
  )

# ------------------------------------------------------------------------------|
#      Role claim summaries ----------------------------------------------------|
# ------------------------------------------------------------------------------|
summarise_claims <- function(data, disease_col, entity_col, claim_col, confidence_col, review_col, prefix) {
  if (is.null(data) || !all(c(disease_col, entity_col, claim_col) %in% names(data))) {
    return(tibble())
  }

  confidence_col <- if (confidence_col %in% names(data)) confidence_col else NULL
  review_col <- if (review_col %in% names(data)) review_col else NULL

  data %>%
    transmute(
      table = prefix,
      disease_name = .data[[disease_col]],
      entity = .data[[entity_col]],
      claim = .data[[claim_col]],
      confidence = if (is.null(confidence_col)) NA_character_ else .data[[confidence_col]],
      needs_manual_review = if (is.null(review_col)) NA else is_true(.data[[review_col]])
    ) %>%
    filter(!is.na(disease_name), !is.na(claim)) %>%
    group_by(table, disease_name, claim, confidence) %>%
    summarise(
      rows = dplyr::n(),
      entities = dplyr::n_distinct(entity, na.rm = TRUE),
      manual_review_rows = sum(needs_manual_review, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(table, disease_name, claim, confidence)
}

role_evidence_claim_summary <- bind_rows(
  summarise_claims(
    layers$host_role_evidence,
    "disease_name",
    "host",
    "role_claim",
    "claim_confidence",
    "needs_manual_review",
    "host_role_evidence"
  ),
  summarise_claims(
    layers$vector_role_evidence,
    "disease_name",
    "vector_species",
    "vector_role_claim",
    "claim_confidence",
    "needs_manual_review",
    "vector_role_evidence"
  )
)

role_assignment_claim_summary <- bind_rows(
  summarise_claims(
    layers$host_role_assignments,
    "disease_name",
    "host",
    "host_role_assignment",
    "assignment_confidence",
    "needs_manual_review",
    "host_role_assignments"
  ),
  summarise_claims(
    layers$vector_role_assignments,
    "disease_name",
    "vector_species",
    "vector_role_assignment",
    "assignment_confidence",
    "needs_manual_review",
    "vector_role_assignments"
  )
)

# ------------------------------------------------------------------------------|
#      Write outputs -----------------------------------------------------------|
# ------------------------------------------------------------------------------|
write_csv(inventory, file.path(qa_dir, "evidence_layer_inventory.csv"), na = "")
write_csv(disease_readiness, file.path(qa_dir, "disease_evidence_readiness.csv"), na = "")
write_csv(vector_readiness, file.path(qa_dir, "vector_evidence_readiness_by_disease.csv"), na = "")
write_csv(role_gap_summary, file.path(qa_dir, "role_evidence_gap_summary.csv"), na = "")
write_csv(role_evidence_claim_summary, file.path(qa_dir, "role_evidence_claim_summary.csv"), na = "")
write_csv(role_assignment_claim_summary, file.path(qa_dir, "role_assignment_claim_summary.csv"), na = "")

message("Wrote QA inventory: ", file.path(qa_dir, "evidence_layer_inventory.csv"))
message("Wrote disease readiness QA: ", file.path(qa_dir, "disease_evidence_readiness.csv"))
message("Wrote vector readiness QA: ", file.path(qa_dir, "vector_evidence_readiness_by_disease.csv"))
message("Wrote role evidence gap summary: ", file.path(qa_dir, "role_evidence_gap_summary.csv"))
message("Wrote role evidence claim summary: ", file.path(qa_dir, "role_evidence_claim_summary.csv"))
message("Wrote role assignment claim summary: ", file.path(qa_dir, "role_assignment_claim_summary.csv"))
