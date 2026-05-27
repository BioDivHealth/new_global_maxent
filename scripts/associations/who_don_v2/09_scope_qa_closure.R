library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

who_don_v2_ensure_dirs()

completed_review_dir <- function(...) {
  who_don_v2_qa_archive_dir("completed_review_surfaces", ...)
}

input_path <- completed_review_dir("v2_scope_adjudication_candidates_enriched.csv")
workpack_path <- completed_review_dir("v2_scope_qa_closure_workpack.csv")
summary_path <- completed_review_dir("v2_scope_qa_closure_summary.csv")
manifest_path <- completed_review_dir("v2_scope_qa_closure_manifest.csv")
review_path <- who_don_v2_output_dir("review", "v2_scope_adjudication_review_decisions.csv")

required_cols <- c(
  "review_id",
  "record_key",
  "country_standard",
  "disease_standard",
  "scope_review_bucket",
  "review_priority",
  "evidence_text",
  "claim_type",
  "association_scope",
  "scope_rule_id",
  "llm_use_policy"
)

review_decision_cols <- c(
  "decision_id",
  "scope_review_key",
  "review_id",
  "record_key",
  "country_standard",
  "disease_standard",
  "scope_review_bucket",
  "scope_workpack",
  "pattern_group_key",
  "evidence_text",
  "review_decision",
  "decision_confidence",
  "decision_reason",
  "reviewer",
  "review_date",
  "rule_candidate",
  "proposed_rule_id",
  "review_note",
  "llm_use_policy"
)

text_norm <- function(x) {
  x %>%
    coalesce("") %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[0-9]+", "<num>") %>%
    str_squish()
}

stable_hash <- function(x) {
  vapply(
    x,
    function(value) {
      if (is.na(value)) {
        value <- ""
      }
      bytes <- as.integer(charToRaw(enc2utf8(value)))
      hash <- 2166136261
      for (byte in bytes) {
        hash <- (hash * 16777619 + byte) %% 4294967296
      }
      sprintf("%010.0f", hash)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

add_missing_review_cols <- function(x) {
  missing_cols <- setdiff(review_decision_cols, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA_character_
  }
  x %>%
    mutate(across(all_of(review_decision_cols), as.character)) %>%
    select(all_of(review_decision_cols))
}

scope_workpack_for <- function(scope_review_bucket, disease_standard) {
  disease <- text_norm(disease_standard)
  case_when(
    scope_review_bucket == "uncertain_language_needs_sampling" &
      str_detect(disease, regex("influenza|respiratory|mers|sars", ignore_case = TRUE)) ~
      "A_respiratory_influenza",
    scope_review_bucket == "uncertain_language_needs_sampling" &
      str_detect(disease, regex("cholera|salmonellosis|listeriosis|hepatitis e|gastro|food|shigellosis|melamine", ignore_case = TRUE)) ~
      "B_gastro_foodborne",
    scope_review_bucket == "uncertain_language_needs_sampling" &
      str_detect(disease, regex("ebola|haemorrhagic|hemorrhagic|lassa|marburg|rift valley|plague|anthrax|crimean|cchf", ignore_case = TRUE)) ~
      "C_haemorrhagic_zoonotic",
    scope_review_bucket == "uncertain_language_needs_sampling" ~
      "D_vector_vpd_other",
    scope_review_bucket == "import_or_exposure_context_language" ~
      "E_import_exposure_context",
    scope_review_bucket == "surveillance_lab_or_partner_language" ~
      "F_surveillance_lab_partner_context",
    scope_review_bucket == "background_or_historical_language" ~
      "G_background_historical_context",
    TRUE ~ "H_other_scope_review"
  )
}

scope_context_pattern <- paste(
  c(
    "recommendations by",
    "for more information",
    "specific recommendations",
    "product distribution",
    "trace forward",
    "hecolin",
    "licensed .*china",
    "notification of the cases",
    "notified of the cases",
    "no special restrictions on travel or trade",
    "travel or trade",
    "travel advice",
    "decision affects imports",
    "\\bimports from\\b",
    "\\bexported from\\b",
    "history of travel",
    "travelled to",
    "returned from",
    "exposure in",
    "exposures in",
    "source of infection",
    "surrounding countries",
    "neighbou?r",
    "\\bborder",
    "at-risk",
    "previous outbreaks?",
    "previously",
    "historical",
    "endemic",
    "in the past",
    "since <num>.*sporadic",
    "laboratory analyses",
    "reference laboratory",
    "laboratory network",
    "samples? (was|were )?sent",
    "sequencing",
    "sequence analysis",
    "genomic",
    "supporting the response",
    "support the response",
    "technical support",
    "who office",
    "field assessments",
    "task force consists",
    "public health professionals from",
    "team members were drawn",
    "conference in",
    "study presented",
    "transferred to",
    "vaccine doses sent",
    "free of cholera",
    "no cholera cases",
    "vector .*reported from",
    "aedes albopictus.*reported from",
    "removed from areas with recent local",
    "government.*continuing.*source",
    "french society of paediatrics",
    "available at:",
    "\\bhttp",
    "factsheet"
  ),
  collapse = "|"
)

scope_direct_event_pattern <- paste(
  c(
    "of the <num> cases confirmed to date in .* <num> have been fatal",
    "altogether, <num> cases, <num> of them fatal, have been reported in",
    "confirmed cases? (had|have|has) been reported in",
    "laboratory-confirmed cases? (had|have|has) been reported in",
    "new cases? (were|was) reported in",
    "has reported a total of <num> cases",
    "reported a total of <num> cases",
    "has confirmed the country.?s .* case",
    "reported a confirmed case of",
    "human cases? of .* continue to occur in",
    "has officially declared the epidemic",
    "first case of .* diagnosed in",
    "countries recently reporting new or increased .* activity are",
    "has been detected in .* sewage samples",
    "new emergence of .* in ",
    "virus circulation .* reported in",
    "outbreaks? (was|were|has been|have been) reported in",
    "outbreaks? (is|are) occurring in",
    "declared an outbreak"
  ),
  collapse = "|"
)

title_like_event_hit <- function(evidence_norm, country_standard, disease_standard) {
  country_norm <- text_norm(country_standard)
  country_pattern <- str_replace_all(country_norm, "([\\W])", "\\\\\\1")
  country_hit <- mapply(
    function(text, pattern) {
      pattern != "" && str_detect(text, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE))
    },
    evidence_norm,
    country_pattern,
    USE.NAMES = FALSE
  )
  disease_hit <- str_detect(
    evidence_norm,
    regex(
      paste(
        c(
          "cholera",
          "ebola",
          "marburg",
          "lassa",
          "plague",
          "yellow fever",
          "dengue",
          "meningococcal",
          "west nile",
          "japanese encephalitis",
          "hepatitis e",
          "shigellosis",
          "rift valley",
          "crimean-congo",
          "anthrax",
          "mayaro",
          "poliovirus"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  country_hit &
    disease_hit &
    str_count(evidence_norm, "\\S+") <= 10L &
    str_detect(evidence_norm, regex("\\bin\\s+", ignore_case = TRUE)) &
    !str_detect(text_norm(disease_standard), regex("influenza|respiratory|mers|sars", ignore_case = TRUE))
}

scope_pattern_decisions <- function(scope_review_bucket, evidence_norm, country_standard, disease_standard) {
  direct_event_hit <- str_detect(evidence_norm, regex(scope_direct_event_pattern, ignore_case = TRUE))
  context_hit <- scope_review_bucket %in% c(
    "import_or_exposure_context_language",
    "surveillance_lab_or_partner_language",
    "background_or_historical_language"
  ) | str_detect(evidence_norm, regex(scope_context_pattern, ignore_case = TRUE))
  direct_event_hit <- direct_event_hit |
    title_like_event_hit(evidence_norm, country_standard, disease_standard)

  direct_event_hit <- direct_event_hit & !context_hit

  tibble::tibble(
    review_decision = case_when(
      direct_event_hit ~ "accept_event_pattern",
      context_hit ~ "reject_context_pattern",
      TRUE ~ "defer_insufficient_evidence_closed"
    ),
    decision_confidence = case_when(
      direct_event_hit ~ "high",
      context_hit ~ "high",
      TRUE ~ "low"
    ),
    decision_reason = case_when(
      direct_event_hit ~
        "Repeated direct event-country wording: the evidence span itself reports confirmed cases, deaths, or an outbreak in the candidate country.",
      context_hit ~
        "Repeated non-event scope context: the evidence span is import/travel/exposure, surveillance/lab/partner, historical/background, neighbouring-country, trade, or support wording.",
      TRUE ~
        "Scope QA review found no conservative repeated event or context pattern in the available evidence span; close as insufficient evidence with no production effect."
    ),
    reviewer = "codex_scope_qa_closure",
    review_date = format(Sys.Date(), "%Y-%m-%d"),
    rule_candidate = case_when(
      direct_event_hit ~ "yes",
      context_hit ~ "yes",
      TRUE ~ "no"
    ),
    proposed_rule_id = case_when(
      direct_event_hit ~ "scope_direct_event_country_pattern",
      context_hit ~ "scope_non_event_context_pattern",
      TRUE ~ NA_character_
    ),
    review_note = case_when(
      direct_event_hit ~
        "Use only through the deterministic claim/scope policy layer; do not apply as a row-level modelling override.",
      context_hit ~
        "Use only through the deterministic claim/scope policy layer; these rows should not become focal-event modelling evidence.",
      TRUE ~
        "Closed rather than inferred from disease, country, section, or record title alone."
    )
  )
}

allowed_review_decisions <- c(
  "accept_event_pattern",
  "reject_context_pattern",
  "needs_full_article_review",
  "needs_rule_change",
  "defer_ambiguous",
  "defer_insufficient_evidence_closed"
)

scope_candidates <- v2_read_csv(input_path, required_cols)

if (anyDuplicated(scope_candidates$review_id) > 0) {
  stop("Expected unique review_id values in ", input_path, call. = FALSE)
}

scope_workpack <- scope_candidates %>%
  mutate(
    evidence_norm = text_norm(evidence_text),
    scope_workpack = scope_workpack_for(scope_review_bucket, disease_standard),
    scope_review_key = review_id,
    pattern_group_key = paste(scope_workpack, disease_standard, evidence_norm, sep = "::"),
    stable_order_key = stable_hash(scope_review_key)
  ) %>%
  add_count(pattern_group_key, name = "pattern_group_rows") %>%
  arrange(scope_workpack, stable_order_key, scope_review_key)

expected_workpack_counts <- tibble::tribble(
  ~scope_workpack, ~expected_rows,
  "A_respiratory_influenza", 532L,
  "B_gastro_foodborne", 229L,
  "C_haemorrhagic_zoonotic", 309L,
  "D_vector_vpd_other", 353L,
  "E_import_exposure_context", 107L,
  "F_surveillance_lab_partner_context", 61L,
  "G_background_historical_context", 36L
)

observed_workpack_counts <- scope_workpack %>%
  count(scope_workpack, name = "rows")

unexpected_counts <- expected_workpack_counts %>%
  left_join(observed_workpack_counts, by = "scope_workpack") %>%
  mutate(rows = coalesce(rows, 0L)) %>%
  filter(rows != expected_rows)

if (nrow(unexpected_counts) > 0) {
  message(
    "Scope closure workpack counts differ from the accepted 1627-row baseline; ",
    "continuing because deterministic scope rules may have changed the current QA surface."
  )
}

new_review_rows <- scope_workpack %>%
  bind_cols(scope_pattern_decisions(.$scope_review_bucket, .$evidence_norm, .$country_standard, .$disease_standard)) %>%
  transmute(
    decision_id = paste("scope_qa_closure", stable_order_key, sep = "::"),
    scope_review_key,
    review_id,
    record_key,
    country_standard,
    disease_standard,
    scope_review_bucket,
    scope_workpack,
    pattern_group_key,
    evidence_text,
    review_decision,
    decision_confidence,
    decision_reason,
    reviewer,
    review_date,
    rule_candidate,
    proposed_rule_id,
    review_note,
    llm_use_policy = "not_llm_input"
  )

existing_review_rows <- if (file.exists(review_path)) {
  existing <- add_missing_review_cols(v2_read_csv(review_path))
  bad_decisions <- existing %>%
    filter(
      !is.na(review_decision),
      review_decision != "",
      !review_decision %in% allowed_review_decisions
    ) %>%
    distinct(review_decision)
  if (nrow(bad_decisions) > 0) {
    stop(
      "Unexpected review_decision values in ", review_path, ": ",
      paste(bad_decisions$review_decision, collapse = ", "),
      call. = FALSE
    )
  }
  existing
} else {
  tibble::tibble(!!!stats::setNames(rep(list(character()), length(review_decision_cols)), review_decision_cols))
}

review_rows <- existing_review_rows %>%
  semi_join(scope_workpack %>% distinct(scope_review_key), by = "scope_review_key") %>%
  bind_rows(
    new_review_rows %>%
      anti_join(existing_review_rows, by = "scope_review_key")
  ) %>%
  left_join(
    new_review_rows %>%
      transmute(
        scope_review_key,
        default_review_decision = review_decision,
        default_decision_confidence = decision_confidence,
        default_decision_reason = decision_reason,
        default_reviewer = reviewer,
        default_review_date = review_date,
        default_rule_candidate = rule_candidate,
        default_proposed_rule_id = proposed_rule_id,
        default_review_note = review_note,
        default_llm_use_policy = llm_use_policy
      ),
    by = "scope_review_key"
  ) %>%
  mutate(
    generated_to_default = reviewer == "codex_scope_qa_closure" &
      !is.na(default_review_decision) &
      default_review_decision != "",
    review_decision = if_else(
      is.na(review_decision) | review_decision == "" | generated_to_default,
      default_review_decision,
      review_decision
    ),
    decision_confidence = if_else(
      is.na(decision_confidence) | decision_confidence == "" | generated_to_default,
      default_decision_confidence,
      decision_confidence
    ),
    decision_reason = if_else(
      is.na(decision_reason) | decision_reason == "" | generated_to_default,
      default_decision_reason,
      decision_reason
    ),
    reviewer = if_else(is.na(reviewer) | reviewer == "" | generated_to_default, default_reviewer, reviewer),
    review_date = if_else(is.na(review_date) | review_date == "" | generated_to_default, default_review_date, review_date),
    rule_candidate = if_else(
      is.na(rule_candidate) | rule_candidate == "" | generated_to_default,
      default_rule_candidate,
      rule_candidate
    ),
    proposed_rule_id = if_else(
      is.na(proposed_rule_id) | proposed_rule_id == "" | generated_to_default,
      default_proposed_rule_id,
      proposed_rule_id
    ),
    review_note = if_else(
      is.na(review_note) | review_note == "" | generated_to_default,
      default_review_note,
      review_note
    ),
    llm_use_policy = if_else(
      is.na(llm_use_policy) | llm_use_policy == "" | generated_to_default,
      default_llm_use_policy,
      llm_use_policy
    )
  ) %>%
  select(all_of(review_decision_cols)) %>%
  arrange(scope_workpack, scope_review_key)

if (nrow(review_rows) != nrow(scope_workpack)) {
  stop("Scope decision rows do not align one-to-one with the closure workpack.", call. = FALSE)
}
if (anyDuplicated(review_rows$scope_review_key) > 0) {
  stop("Duplicate scope_review_key values found in review decisions.", call. = FALSE)
}

workpack_summary <- scope_workpack %>%
  count(scope_workpack, scope_review_bucket, name = "rows") %>%
  arrange(scope_workpack, scope_review_bucket)

decision_summary <- review_rows %>%
  count(scope_workpack, review_decision, decision_confidence, proposed_rule_id, name = "rows") %>%
  arrange(scope_workpack, review_decision, desc(rows))

summary <- bind_rows(
  workpack_summary %>%
    transmute(
      summary_type = "workpack_rows",
      scope_workpack,
      scope_review_bucket,
      review_decision = NA_character_,
      decision_confidence = NA_character_,
      proposed_rule_id = NA_character_,
      rows
    ),
  decision_summary %>%
    transmute(
      summary_type = "decision_rows",
      scope_workpack,
      scope_review_bucket = NA_character_,
      review_decision,
      decision_confidence,
      proposed_rule_id,
      rows
    )
)

manifest <- tibble::tibble(
  metric = c(
    "generated_at_utc",
    "input_file",
    "input_file_md5",
    "total_scope_rows",
    "uncertain_workpack_rows",
    "context_workpack_rows",
    "review_decision_file",
    "llm_policy",
    "baseline_count_status"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    input_path,
    unname(tools::md5sum(input_path)),
    as.character(nrow(scope_workpack)),
    as.character(sum(scope_workpack$scope_review_bucket == "uncertain_language_needs_sampling")),
    as.character(sum(scope_workpack$scope_review_bucket != "uncertain_language_needs_sampling")),
    review_path,
    "not_llm_input",
    if_else(nrow(unexpected_counts) == 0, "matches_accepted_1627_baseline", "changed_after_deterministic_scope_rules")
  ),
  note = c(
    "Timestamp for this regeneration; row identity is driven by review_id, not this value.",
    "Optional scope QA surface from 07_quality_tightening_review_surfaces.R.",
    "Checksum of the input scope QA surface used for this closure pass.",
    "All current scope adjudication rows assigned exactly one closure workpack.",
    "Rows split into four disease-oriented workpacks for read-only subagent review.",
    "Rows from import/exposure, surveillance/lab/partner, and background/historical context buckets.",
    "Durable optional scope review decisions; production does not read this file.",
    "This closure surface is not a broad LLM queue.",
    "Records whether the current QA surface still matches the pre-closure 1627-row baseline."
  )
)

v2_write_csv(scope_workpack, workpack_path)
v2_write_csv(summary, summary_path)
v2_write_csv(manifest, manifest_path)
v2_write_csv(review_rows, review_path)

message(
  "Wrote scope QA closure workpack: ",
  nrow(scope_workpack),
  " rows; decisions: ",
  nrow(review_rows)
)
