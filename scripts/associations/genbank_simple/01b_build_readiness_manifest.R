# ------------------------------------------------------------------------------|
#      01b_build_readiness_manifest.R -----------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Build a reviewable expanded GenBank-simple manifest from the current
#          disease modelling readiness surface.
# Inputs : disease_modelling_readiness.csv
#          disease_modelling_readiness_full.csv
# Outputs: genbank_simple_readiness_manifest.csv
#          qa/genbank_simple_readiness_manifest_qa.csv
#
# Notes  : This script does not contact NCBI and does not overwrite the current
#          19-target GenBank-simple manifest. Set
#          GENBANK_SIMPLE_USE_LEGACY_19_MANIFEST=TRUE to add temporary
#          old-manifest provenance fields for comparison.
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))
source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Define input and output paths ------------------------------------------
# ------------------------------------------------------------------------------|
readiness_path <- file.path(readiness_dir, "disease_modelling_readiness.csv")
readiness_full_path <- file.path(readiness_dir, "disease_modelling_readiness_full.csv")
legacy_readiness_path <- here(
  "pathogen_association_data",
  "WHO",
  "role_annotation",
  "qa",
  "disease_modelling_readiness.csv"
)
legacy_readiness_full_path <- here(
  "pathogen_association_data",
  "WHO",
  "role_annotation",
  "qa",
  "disease_modelling_readiness_full.csv"
)
output_dir <- genbank_simple_dir
use_legacy_19_manifest <- parse_env_flag(
  "GENBANK_SIMPLE_USE_LEGACY_19_MANIFEST",
  default = FALSE
)
current_manifest_path <- genbank_simple_existing_file_path(
  output_dir,
  "genbank_simple_manifest.csv"
)
override_path <- genbank_simple_existing_file_path(
  output_dir,
  "genbank_readiness_query_overrides.csv"
)
readiness_manifest_path <- genbank_simple_file_path(
  output_dir,
  "genbank_simple_readiness_manifest.csv",
  create_parent = TRUE
)
readiness_qa_path <- genbank_simple_file_path(
  output_dir,
  "genbank_simple_readiness_manifest_qa.csv",
  create_parent = TRUE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Resolve legacy readiness fallbacks -------------------------------------
# ------------------------------------------------------------------------------|
if (!file.exists(readiness_path) && file.exists(legacy_readiness_path)) {
  readiness_path <- legacy_readiness_path
}

if (!file.exists(readiness_full_path) && file.exists(legacy_readiness_full_path)) {
  readiness_full_path <- legacy_readiness_full_path
}

# ------------------------------------------------------------------------------|
#      Define required schemas -------------------------------------------------
# ------------------------------------------------------------------------------|
required_readiness_cols <- c(
  "analysis_unit_id",
  "readiness_disease_name",
  "analysis_unit_label",
  "modelling_scope_status",
  "recommended_next_action",
  "pathogen_species_name",
  "pathogen_taxid",
  "genbank_distinct_countries_or_territories",
  "who_don_distinct_countries"
)
required_full_cols <- c(
  "analysis_unit_id",
  "source_pathogen",
  "source_msl39_viral_name",
  "matched_pathogen_names",
  "matched_taxids",
  "match_review_flag",
  "shared_species_proxy_flag",
  "match_review_notes",
  "pathogen_species_name_source",
  "pathogen_taxid_source"
)
required_override_cols <- c(
  "analysis_unit_id",
  "review_decision",
  "override_query_pathogen_label",
  "override_pathogen_taxid",
  "override_reason"
)

# ------------------------------------------------------------------------------|
#      Helper functions --------------------------------------------------------
# ------------------------------------------------------------------------------|
stop_if_missing <- function(data, cols, label) {
  missing_cols <- setdiff(cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(data)
}

has_concrete_influenza_subtype <- function(x) {
  x <- clean_text(x)
  !is.na(x) & stringr::str_detect(x, "\\(H[0-9]+N[0-9]+\\)")
}

collapse_flags <- function(x) {
  x <- clean_text(x)
  x <- unique(x[!is.na(x)])

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

collapse_semicolon_flags <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  collapse_flags(unlist(strsplit(x, ";\\s*")))
}

target_status <- function(status) {
  status <- clean_text(status)

  dplyr::case_when(
    any(status == "defer_scope_guardrail", na.rm = TRUE) ~ "defer_scope_guardrail",
    any(status == "defer_no_species_name", na.rm = TRUE) ~ "defer_no_species_name",
    any(status == "review_before_future_retrieval", na.rm = TRUE) ~
      "review_before_future_retrieval",
    TRUE ~ "ready_for_future_retrieval"
  )
}

# ------------------------------------------------------------------------------|
#      Load readiness tables ---------------------------------------------------
# ------------------------------------------------------------------------------|
readiness <- read_csv(readiness_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))
readiness_full <- read_csv(readiness_full_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

stop_if_missing(readiness, required_readiness_cols, "disease_modelling_readiness.csv")
stop_if_missing(readiness_full, required_full_cols, "disease_modelling_readiness_full.csv")

if (anyDuplicated(readiness$analysis_unit_id) > 0) {
  stop("disease_modelling_readiness.csv has duplicated analysis_unit_id values.", call. = FALSE)
}

if (anyDuplicated(readiness_full$analysis_unit_id) > 0) {
  stop("disease_modelling_readiness_full.csv has duplicated analysis_unit_id values.", call. = FALSE)
}

# ------------------------------------------------------------------------------|
#      Load optional legacy manifest and manual query overrides ---------------
# ------------------------------------------------------------------------------|
current_manifest <- if (use_legacy_19_manifest && file.exists(current_manifest_path)) {
  read_csv(current_manifest_path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text)) %>%
    transmute(
      current_target_id = target_id,
      current_manifest_pathogen = Pathogens,
      current_manifest_disease = Disease_name,
      current_manifest_query_used = query_used,
      current_manifest_query_strategy = query_strategy
    ) %>%
    distinct(
      current_manifest_pathogen,
      current_manifest_disease,
      .keep_all = TRUE
    )
} else {
  tibble(
    current_target_id = character(),
    current_manifest_pathogen = character(),
    current_manifest_disease = character(),
    current_manifest_query_used = character(),
    current_manifest_query_strategy = character()
  )
}

query_overrides <- if (file.exists(override_path)) {
  read_csv(
    override_path,
    col_types = cols(.default = col_character()),
    na = c("", "NA")
  ) %>%
    mutate(across(where(is.character), clean_text)) %>%
    stop_if_missing(required_override_cols, "genbank_readiness_query_overrides.csv") %>%
    select(any_of(c(required_override_cols, "review_notes", "source_url"))) %>%
    filter(!is.na(review_decision)) %>%
    mutate(
      review_decision = stringr::str_to_lower(review_decision)
    )
} else {
  tibble(
    analysis_unit_id = character(),
    review_decision = character(),
    override_query_pathogen_label = character(),
    override_pathogen_taxid = character(),
    override_reason = character()
  )
}

unexpected_review_decisions <- setdiff(
  unique(query_overrides$review_decision),
  c("approve", "override", "defer")
)
if (length(unexpected_review_decisions) > 0) {
  stop(
    "Unexpected review_decision values in genbank_readiness_query_overrides.csv: ",
    paste(unexpected_review_decisions, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(query_overrides$analysis_unit_id) > 0) {
  stop("genbank_readiness_query_overrides.csv has duplicated analysis_unit_id values.", call. = FALSE)
}

# ------------------------------------------------------------------------------|
#      Build one QA row per readiness analysis unit ---------------------------
# ------------------------------------------------------------------------------|
readiness_qa <- readiness %>%
  left_join(
    readiness_full %>%
      select(all_of(required_full_cols)),
    by = "analysis_unit_id"
  ) %>%
  left_join(query_overrides, by = "analysis_unit_id") %>%
  mutate(
    match_review_flag = as_logical_flag(match_review_flag),
    shared_species_proxy_flag = as_logical_flag(shared_species_proxy_flag),
    has_manual_genbank_override = review_decision %in% c("approve", "override"),
    manual_genbank_defer = review_decision == "defer",
    manual_genbank_defer = dplyr::coalesce(manual_genbank_defer, FALSE),
    genbank_distinct_countries_or_territories = suppressWarnings(
      as.numeric(genbank_distinct_countries_or_territories)
    ),
    genbank_distinct_countries_or_territories = dplyr::coalesce(
      genbank_distinct_countries_or_territories,
      0
    ),
    missing_species_name = is.na(pathogen_species_name),
    source_pathogen_concrete_influenza = has_concrete_influenza_subtype(source_pathogen),
    query_pathogen_label = case_when(
      has_manual_genbank_override & !is.na(override_query_pathogen_label) ~
        override_query_pathogen_label,
      source_pathogen_concrete_influenza ~ source_pathogen,
      !missing_species_name ~ pathogen_species_name,
      TRUE ~ NA_character_
    ),
    query_pathogen_source = case_when(
      has_manual_genbank_override & !is.na(override_query_pathogen_label) ~
        "manual_genbank_readiness_override",
      source_pathogen_concrete_influenza ~ "source_pathogen_concrete_influenza_subtype",
      !missing_species_name ~ "pathogen_species_name",
      TRUE ~ NA_character_
    ),
    pathogen_taxid = case_when(
      has_manual_genbank_override & !is.na(override_pathogen_taxid) ~ override_pathogen_taxid,
      TRUE ~ pathogen_taxid
    ),
    pathogen_taxid_source = case_when(
      has_manual_genbank_override & !is.na(override_pathogen_taxid) ~
        "manual_genbank_readiness_override",
      TRUE ~ pathogen_taxid_source
    ),
    target_id = if_else(
      !is.na(query_pathogen_label),
      sanitize_filename(query_pathogen_label),
      NA_character_
    ),
    coronavirus_scope_deferred = is_coronavirus_excluded(
      query_pathogen_label,
      readiness_disease_name
    ),
    coronavirus_scope_deferred = dplyr::coalesce(coronavirus_scope_deferred, FALSE),
    broad_influenza_without_subtype =
      is_broad_or_unwanted_influenza(query_pathogen_label) &
      !source_pathogen_concrete_influenza,
    broad_influenza_without_subtype = dplyr::coalesce(
      broad_influenza_without_subtype,
      FALSE
    ),
    non_include_scope = modelling_scope_status != "include",
    non_include_scope = dplyr::coalesce(non_include_scope, FALSE),
    taxid_missing = is.na(pathogen_taxid),
    already_has_genbank_country_evidence =
      genbank_distinct_countries_or_territories > 0
  )

duplicate_species_groups <- readiness_qa %>%
  filter(!missing_species_name) %>%
  count(pathogen_species_name, name = "species_group_rows") %>%
  filter(species_group_rows > 1)

# ------------------------------------------------------------------------------|
#      Add status, QA flags, and optional legacy-manifest matches -------------
# ------------------------------------------------------------------------------|
readiness_qa <- readiness_qa %>%
  left_join(duplicate_species_groups, by = "pathogen_species_name") %>%
  mutate(
    species_group_rows = dplyr::coalesce(species_group_rows, 1L),
    duplicate_species_group = species_group_rows > 1
  ) %>%
  left_join(
    current_manifest,
    by = c(
      "query_pathogen_label" = "current_manifest_pathogen",
      "readiness_disease_name" = "current_manifest_disease"
    )
  ) %>%
  mutate(
    existing_19_target = !is.na(current_target_id),
    row_manifest_status = case_when(
      manual_genbank_defer ~ "review_before_future_retrieval",
      missing_species_name ~ "defer_no_species_name",
      coronavirus_scope_deferred | broad_influenza_without_subtype ~
        "defer_scope_guardrail",
      has_manual_genbank_override & !taxid_missing ~ "ready_for_future_retrieval",
      taxid_missing | shared_species_proxy_flag | match_review_flag ~
        "review_before_future_retrieval",
      TRUE ~ "ready_for_future_retrieval"
    )
  ) %>%
  rowwise() %>%
  mutate(
    qa_flags = collapse_flags(c(
      if (missing_species_name) "missing_species_name" else NA_character_,
      if (non_include_scope) "non_include_scope" else NA_character_,
      if (coronavirus_scope_deferred) "coronavirus_scope_deferred" else NA_character_,
      if (broad_influenza_without_subtype) "broad_influenza_without_subtype" else NA_character_,
      if (duplicate_species_group) "duplicate_species_group" else NA_character_,
      if (shared_species_proxy_flag) "shared_species_proxy" else NA_character_,
      if (match_review_flag) "match_review_flag" else NA_character_,
      if (has_manual_genbank_override) "manual_genbank_override" else NA_character_,
      if (manual_genbank_defer) "manual_genbank_defer" else NA_character_,
      if (taxid_missing) "taxid_missing" else NA_character_,
      if (already_has_genbank_country_evidence) "already_has_genbank_country_evidence" else NA_character_,
      if (use_legacy_19_manifest && existing_19_target) "existing_19_target" else NA_character_
    )),
    row_manifest_status_reason = collapse_flags(c(
      if (has_manual_genbank_override & !taxid_missing) "manual GenBank query override applied" else NA_character_,
      if (manual_genbank_defer) "manual GenBank query review deferred" else NA_character_,
      if (missing_species_name) "pathogen_species_name missing in slim readiness" else NA_character_,
      if (coronavirus_scope_deferred) "coronavirus scope deferred" else NA_character_,
      if (broad_influenza_without_subtype) "broad influenza label lacks concrete subtype" else NA_character_,
      if (taxid_missing) "pathogen_taxid missing" else NA_character_,
      if (shared_species_proxy_flag & !has_manual_genbank_override) {
        "shared species proxy requires review"
      } else {
        NA_character_
      },
      if (match_review_flag & !has_manual_genbank_override) {
        "pathogen match review flag set"
      } else {
        NA_character_
      },
      if (!missing_species_name &&
          !coronavirus_scope_deferred &&
          !broad_influenza_without_subtype &&
          !taxid_missing &&
          !shared_species_proxy_flag &&
          !match_review_flag &&
          !has_manual_genbank_override) {
        "species-level GenBank query ready; modelling_scope_status retained as QA only"
      } else {
        NA_character_
      }
    ))
  ) %>%
  ungroup() %>%
  select(
    analysis_unit_id,
    readiness_disease_name,
    analysis_unit_label,
    modelling_scope_status,
    recommended_next_action,
    pathogen_species_name,
    pathogen_taxid,
    source_pathogen,
    source_msl39_viral_name,
    matched_pathogen_names,
    matched_taxids,
    match_review_flag,
    shared_species_proxy_flag,
    match_review_notes,
    pathogen_species_name_source,
    pathogen_taxid_source,
    review_decision,
    override_query_pathogen_label,
    override_pathogen_taxid,
    override_reason,
    query_pathogen_label,
    query_pathogen_source,
    target_id,
    row_manifest_status,
    row_manifest_status_reason,
    qa_flags,
    genbank_distinct_countries_or_territories,
    who_don_distinct_countries,
    current_target_id,
    current_manifest_query_used,
    current_manifest_query_strategy
  ) %>%
  arrange(
    is.na(query_pathogen_label),
    query_pathogen_label,
    readiness_disease_name,
    analysis_unit_id
  )

# ------------------------------------------------------------------------------|
#      Collapse QA rows into unique future retrieval targets ------------------
# ------------------------------------------------------------------------------|
readiness_manifest <- readiness_qa %>%
  filter(!is.na(query_pathogen_label)) %>%
  group_by(target_id, query_pathogen_label, query_pathogen_source) %>%
  summarise(
    pathogen_species_name = collapse_unique(pathogen_species_name),
    pathogen_taxid = collapse_unique(pathogen_taxid),
    current_target_id = collapse_unique(current_target_id),
    current_manifest_query_used = collapse_unique(current_manifest_query_used),
    current_manifest_query_strategy = collapse_unique(current_manifest_query_strategy),
    analysis_unit_ids = collapse_unique(analysis_unit_id),
    readiness_disease_names = collapse_unique(readiness_disease_name),
    readiness_row_count = dplyr::n(),
    manifest_status = target_status(row_manifest_status),
    manifest_status_reason = collapse_semicolon_flags(row_manifest_status_reason),
    qa_flags = collapse_semicolon_flags(qa_flags),
    existing_genbank_country_evidence =
      any(genbank_distinct_countries_or_territories > 0, na.rm = TRUE),
    existing_genbank_distinct_countries_or_territories =
      sum(genbank_distinct_countries_or_territories, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    source_db = "nuccore",
    query_used = build_simple_query(
      query_pathogen_label,
      unlist(strsplit(dplyr::coalesce(pathogen_taxid, ""), ";\\s*"))
    ),
    query_strategy = case_when(
      is_allowed_influenza_target(query_pathogen_label) ~
        "influenza_subtype_constrained_full_retrieval",
      !is.na(pathogen_taxid) ~ "taxid_full_retrieval",
      TRUE ~ "organism_name_full_retrieval"
    ),
    query_source = case_when(
      is_allowed_influenza_target(query_pathogen_label) ~ "simple_subtype_guardrail",
      TRUE ~ "readiness_generated_query"
    ),
    retrieval_policy = "future_full_deterministic_pagination",
    allow_sampling = FALSE
  ) %>%
  ungroup() %>%
  select(
    target_id,
    query_pathogen_label,
    query_pathogen_source,
    pathogen_species_name,
    pathogen_taxid,
    source_db,
    query_used,
    query_strategy,
    query_source,
    retrieval_policy,
    allow_sampling,
    analysis_unit_ids,
    readiness_disease_names,
    readiness_row_count,
    manifest_status,
    manifest_status_reason,
    qa_flags,
    existing_genbank_country_evidence,
    existing_genbank_distinct_countries_or_territories,
    current_target_id,
    current_manifest_query_strategy
  ) %>%
  arrange(query_pathogen_label, target_id)

duplicate_target_ids <- readiness_manifest %>%
  count(target_id, name = "target_id_rows") %>%
  filter(target_id_rows > 1)

# ------------------------------------------------------------------------------|
#      Validate and write outputs ---------------------------------------------
# ------------------------------------------------------------------------------|
if (nrow(duplicate_target_ids) > 0) {
  stop(
    "Readiness manifest target_id values are not unique: ",
    paste(duplicate_target_ids$target_id, collapse = ", "),
    call. = FALSE
  )
}

write_csv(readiness_manifest, readiness_manifest_path)
write_csv(readiness_qa, readiness_qa_path)

message("Wrote readiness manifest rows: ", nrow(readiness_manifest))
message("Wrote readiness manifest QA rows: ", nrow(readiness_qa))
message("Legacy 19-target manifest provenance enabled: ", use_legacy_19_manifest)
message(
  "QA species-level rows: ",
  sum(!is.na(readiness_qa$pathogen_species_name))
)
message(
  "Manifest statuses: ",
  paste(
    readiness_manifest %>%
      count(manifest_status) %>%
      transmute(label = paste0(manifest_status, "=", n)) %>%
      pull(label),
    collapse = "; "
  )
)
