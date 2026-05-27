library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

who_don_v2_ensure_dirs()

completed_review_dir <- function(...) {
  who_don_v2_qa_archive_dir("completed_review_surfaces", ...)
}

input_path <- completed_review_dir("v2_native_new_country_priority_review.csv")
sample_path <- completed_review_dir("v2_medium_native_new_country_sample.csv")
manifest_path <- completed_review_dir("v2_medium_native_new_country_sample_manifest.csv")
review_path <- who_don_v2_output_dir("review", "v2_medium_native_new_country_review_decisions.csv")

required_cols <- c(
  "record_key",
  "country_standard",
  "country_claim_type",
  "country_evidence_location",
  "country_confidence",
  "disease_standard",
  "Title",
  "local_evidence_text",
  "country_evidence_text",
  "alias_type",
  "native_new_value_score",
  "priority_bucket",
  "llm_use_policy"
)

review_decision_cols <- c(
  "decision_id",
  "candidate_key",
  "record_key",
  "country_standard",
  "disease_standard",
  "pattern_group_key",
  "local_evidence_text",
  "country_evidence_location",
  "native_new_value_score",
  "review_decision",
  "decision_confidence",
  "decision_reason",
  "reviewer",
  "review_date",
  "rule_candidate",
  "proposed_rule_id",
  "review_note"
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

reported_case_pattern <- paste(
  c(
    "number of laboratory confirmed cases and deaths",
    "countries officially reporting cases",
    "following countries have reported laboratory confirmed cases",
    "reported laboratory confirmed cases with no deaths",
          "reported laboratory-confirmed cases originating in the following countries",
          "also reported laboratory-confirmed cases",
          "cases have now also been notified from",
          "has now been reported in .* in addition to",
          "laboratory[- ]confirmed cases of new influenza .* officially reported to who",
          "update on countries and cases .* cumulative total",
          "new cases were reported (in|from)",
          "new probable sars cases were also reported",
          "cases are now being reported in",
          "countries .* newly reported",
          "countries reported their first suspected cases",
          "reported their first (suspected |probable |confirmed )?cases",
          "reported their first probable cases"
        ),
        collapse = "|"
      )

non_event_context_pattern <- paste(
  c(
    "specific recommendations by national regulatory authorities",
    "listed by country",
    "exported to",
    "contributing personnel and materials",
    "experts come from",
    "worked closely with",
    "collaborating centre",
    "collaborating center",
    "reference laboratory",
    "but no illness or death in humans",
    "reported no influenza activity",
    "no influenza activity was reported",
    "low influenza activity was reported",
    "overall levels of ili remain low",
    "activity remains low",
    "team members were drawn",
    "direct financial support",
    "médecins sans frontières",
    "medecins sans frontieres",
    "are assisting",
    "laboratory network",
    "originating from",
    "spread from",
    "migration of people",
    "neighbou?r",
    "preparedness",
    "readiness"
  ),
  collapse = "|"
)

medium_country_pattern_decisions <- function(evidence_norm) {
  reported_case_hit <- str_detect(evidence_norm, regex(reported_case_pattern, ignore_case = TRUE))
  negative_report_hit <- str_detect(
    evidence_norm,
    regex(
      paste(
        c(
          "reported no influenza activity",
          "no influenza activity was reported",
          "no similar pattern",
          "no cases",
          "without reporting cases",
          "list of definitions of qualitative indicators"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  reported_case_hit <- reported_case_hit & !negative_report_hit
  non_event_context_hit <- str_detect(evidence_norm, regex(non_event_context_pattern, ignore_case = TRUE))

  tibble::tibble(
    review_decision = case_when(
      reported_case_hit ~ "accept_pattern",
      non_event_context_hit ~ "reject_pattern",
      TRUE ~ "defer_insufficient_evidence_closed"
    ),
    decision_confidence = case_when(
      reported_case_hit ~ "high",
      non_event_context_hit ~ "high",
      TRUE ~ "low"
    ),
    decision_reason = case_when(
      reported_case_hit ~
        "Repeated explicit country case-report pattern: the evidence span itself reports laboratory-confirmed cases, notified cases, or reported cases for the listed countries.",
      non_event_context_hit ~
        "Repeated non-event context pattern: the country appears in partner, regulator, export, laboratory, expert-origin, or background wording rather than as an outbreak country.",
      TRUE ~
        "Full medium-country review found no conservative repeated accept/reject pattern in the sampled evidence; close as insufficient evidence with no production effect."
    ),
    reviewer = "codex_full_medium_country_review",
    review_date = format(Sys.Date(), "%Y-%m-%d"),
    rule_candidate = case_when(
      reported_case_hit ~ "yes",
      non_event_context_hit ~ "no",
      TRUE ~ "no"
    ),
    proposed_rule_id = case_when(
      reported_case_hit ~ "medium_native_new_reported_cases_policy",
      TRUE ~ NA_character_
    ),
    review_note = case_when(
      reported_case_hit ~
        "Promote only through a narrow country-adoption policy requiring local-event medium-confidence native-new evidence plus one of these explicit reported-case phrases.",
      non_event_context_hit ~
        "Reviewed as a sample-level rejection pattern; not promoted as a production extraction rule in this pass.",
      TRUE ~
        "Closed rather than inferred from section, disease, or country alone."
    )
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

allowed_review_decisions <- c(
  "accept_pattern",
  "reject_pattern",
  "needs_full_article_review",
  "needs_rule_change",
  "defer_ambiguous",
  "defer_insufficient_evidence_closed"
)

large_pattern_threshold <- 10L
rule_candidate_pattern_threshold <- 5L
per_disease_n <- 5L
high_density_record_threshold <- 5L
high_density_record_n <- 5L

native_new <- v2_read_csv(input_path, required_cols)

medium <- native_new %>%
  filter(priority_bucket == "medium_value_local_event")

if (nrow(medium) == 0) {
  stop("No medium_value_local_event rows found in ", input_path, call. = FALSE)
}

medium_enriched <- medium %>%
  mutate(
    evidence_norm = text_norm(coalesce(local_evidence_text, country_evidence_text)),
    title_norm = text_norm(Title),
    candidate_key = paste(record_key, country_standard, country_evidence_location, evidence_norm, sep = "::"),
    pattern_group_key = paste(disease_standard, country_evidence_location, alias_type, evidence_norm, sep = "::"),
    stable_order_key = stable_hash(candidate_key)
  ) %>%
  add_count(pattern_group_key, name = "pattern_group_rows") %>%
  add_count(record_key, name = "native_new_countries_in_record") %>%
  mutate(
    deterministic_rule_candidate = pattern_group_rows >= rule_candidate_pattern_threshold &
      !str_detect(
        evidence_norm,
        regex(
          paste(
            c(
              "travel", "import", "history", "previous", "surveillance",
              "laborator", "sequence", "sequencing", "partner", "preparedness",
              "readiness", "border", "neighbou?r", "region", "regional",
              "globally", "worldwide"
            ),
            collapse = "|"
          )
        )
      )
  ) %>%
  arrange(stable_order_key, candidate_key)

sample_membership <- bind_rows(
  medium_enriched %>%
    filter(pattern_group_rows >= large_pattern_threshold) %>%
    transmute(candidate_key, sample_rule = "large_repeated_pattern_group"),
  medium_enriched %>%
    filter(deterministic_rule_candidate) %>%
    transmute(candidate_key, sample_rule = "potential_rule_candidate_group"),
  medium_enriched %>%
    group_by(disease_standard) %>%
    arrange(stable_order_key, candidate_key, .by_group = TRUE) %>%
    slice_head(n = per_disease_n) %>%
    ungroup() %>%
    transmute(candidate_key, sample_rule = "per_disease_stable_sample"),
  medium_enriched %>%
    filter(native_new_countries_in_record >= high_density_record_threshold) %>%
    group_by(record_key) %>%
    arrange(stable_order_key, candidate_key, .by_group = TRUE) %>%
    slice_head(n = high_density_record_n) %>%
    ungroup() %>%
    transmute(candidate_key, sample_rule = "high_density_record_sample")
) %>%
  distinct() %>%
  group_by(candidate_key) %>%
  summarise(sample_reason = str_c(sort(unique(sample_rule)), collapse = "; "), .groups = "drop")

medium_sample <- medium_enriched %>%
  inner_join(sample_membership, by = "candidate_key") %>%
  arrange(sample_reason, stable_order_key, candidate_key)

new_review_rows <- medium_sample %>%
  bind_cols(medium_country_pattern_decisions(.$evidence_norm)) %>%
  transmute(
    decision_id = paste("medium_native_new_country", stable_order_key, sep = "::"),
    candidate_key,
    record_key,
    country_standard,
    disease_standard,
    pattern_group_key,
    local_evidence_text,
    country_evidence_location,
    native_new_value_score = as.character(native_new_value_score),
    review_decision,
    decision_confidence,
    decision_reason,
    reviewer,
    review_date,
    rule_candidate = coalesce(rule_candidate, if_else(deterministic_rule_candidate, "yes", "no")),
    proposed_rule_id,
    review_note
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
  semi_join(medium_sample %>% distinct(candidate_key), by = "candidate_key") %>%
  bind_rows(
    new_review_rows %>%
      anti_join(existing_review_rows, by = "candidate_key")
  ) %>%
  left_join(
    new_review_rows %>%
      transmute(
        candidate_key,
        default_review_decision = review_decision,
        default_decision_confidence = decision_confidence,
        default_decision_reason = decision_reason,
        default_reviewer = reviewer,
        default_review_date = review_date,
        default_rule_candidate = rule_candidate,
        default_proposed_rule_id = proposed_rule_id,
        default_review_note = review_note
      ),
    by = "candidate_key"
  ) %>%
  mutate(
    generated_defer_to_default = review_decision == "defer_ambiguous" &
      reviewer %in% c("codex_pattern_review", "codex_full_medium_country_review") &
      !is.na(default_review_decision) &
      default_review_decision != "",
    review_decision = if_else(
      is.na(review_decision) | review_decision == "" | generated_defer_to_default,
      default_review_decision,
      review_decision
    ),
    decision_confidence = if_else(
      is.na(decision_confidence) | decision_confidence == "" | generated_defer_to_default,
      default_decision_confidence,
      decision_confidence
    ),
    decision_reason = if_else(
      is.na(decision_reason) | decision_reason == "" | generated_defer_to_default,
      default_decision_reason,
      decision_reason
    ),
    reviewer = if_else(is.na(reviewer) | reviewer == "" | generated_defer_to_default, default_reviewer, reviewer),
    review_date = if_else(is.na(review_date) | review_date == "" | generated_defer_to_default, default_review_date, review_date),
    rule_candidate = if_else(
      is.na(rule_candidate) | rule_candidate == "" | generated_defer_to_default,
      default_rule_candidate,
      rule_candidate
    ),
    proposed_rule_id = if_else(
      is.na(proposed_rule_id) | proposed_rule_id == "" | generated_defer_to_default,
      default_proposed_rule_id,
      proposed_rule_id
    ),
    review_note = if_else(
      is.na(review_note) | review_note == "" | generated_defer_to_default,
      default_review_note,
      review_note
    )
  ) %>%
  select(all_of(review_decision_cols)) %>%
  arrange(candidate_key, decision_id)

sample_rule_summary <- sample_membership %>%
  tidyr::separate_rows(sample_reason, sep = "; ") %>%
  count(sample_reason, name = "rows") %>%
  arrange(sample_reason)

manifest <- bind_rows(
  tibble::tibble(
    metric = c(
      "generated_at_utc",
      "input_file",
      "input_file_md5",
      "total_medium_rows",
      "total_sampled_rows",
      "diseases_represented",
      "countries_represented",
      "records_represented",
      "pattern_groups_represented",
      "large_pattern_threshold",
      "rule_candidate_pattern_threshold",
      "per_disease_sample_n",
      "high_density_record_threshold",
      "high_density_record_sample_n",
      "review_decision_file"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      input_path,
      unname(tools::md5sum(input_path)),
      as.character(nrow(medium_enriched)),
      as.character(nrow(medium_sample)),
      as.character(n_distinct(medium_sample$disease_standard)),
      as.character(n_distinct(medium_sample$country_standard)),
      as.character(n_distinct(medium_sample$record_key)),
      as.character(n_distinct(medium_sample$pattern_group_key)),
      as.character(large_pattern_threshold),
      as.character(rule_candidate_pattern_threshold),
      as.character(per_disease_n),
      as.character(high_density_record_threshold),
      as.character(high_density_record_n),
      review_path
    ),
    note = c(
      "Timestamp for this regeneration; sample row ordering is driven by stable_order_key, not this value.",
      "Source review surface filtered to priority_bucket == medium_value_local_event.",
      "Checksum of the input review surface used for this sample.",
      "All medium native-new country candidates available for sampling.",
      "Distinct candidate_key rows selected by deterministic sample rules.",
      "Distinct disease_standard values in the sample.",
      "Distinct country_standard values in the sample.",
      "Distinct DON record keys in the sample.",
      "Distinct pattern_group_key values in the sample.",
      "All rows from pattern groups at or above this row count are sampled.",
      "Repeated groups at or above this count are flagged as potential deterministic rule candidates when context blockers are absent.",
      "Stable-hash rows sampled per disease after other grouping fields are materialized.",
      "Records at or above this native-new row count are high-density records.",
      "Stable-hash rows sampled per high-density record.",
      "Durable manual-review template; blank review_decision values have no production effect."
    )
  ),
  sample_rule_summary %>%
    transmute(
      metric = paste0("sample_rule_rows_", sample_reason),
      value = as.character(rows),
      note = "Rows selected by this deterministic sample rule before overlap collapse."
    )
)

v2_write_csv(medium_sample, sample_path)
v2_write_csv(manifest, manifest_path)
v2_write_csv(review_rows, review_path)

message(
  "Wrote medium native-new country sample: ",
  nrow(medium_sample),
  " rows from ",
  nrow(medium_enriched),
  " medium candidates"
)
