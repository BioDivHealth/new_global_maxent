# ------------------------------------------------------------------------------
# 9_WHO_DON_OpenAI_Cleanup.R
# ------------------------------------------------------------------------------
# Purpose: Clean, rank, and deduplicate OpenAI WHO DON country candidate rows
#          into a best-evidence table while preserving full-span provenance and
#          review flags for weaker support.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
batch_dir <- who_don_output_path("openai_batch")

candidates_path <- file.path(batch_dir, "who_don_openai_country_candidates_full.csv")
ranked_path <- file.path(batch_dir, "who_don_openai_country_candidates_ranked_full.csv")
best_path <- file.path(batch_dir, "who_don_openai_country_best_full.csv")
record_summary_path <- file.path(batch_dir, "who_don_openai_country_record_summary_full.csv")
summary_path <- file.path(batch_dir, "who_don_openai_cleanup_summary_full.csv")

records <- read_who_don_records("who_don_records_clean.csv") %>%
  who_don_with_record_key()

candidates <- readr::read_csv(
  candidates_path,
  show_col_types = FALSE,
  na = c("", "NA")
)

normalize_country_key <- function(x) {
  normalize_match_text(x)
}

contains_country_literal <- function(span, country) {
  span_key <- normalize_match_text(span)
  country_key <- normalize_country_key(country)

  ifelse(
    is.na(span_key) | is.na(country_key),
    FALSE,
    stringr::str_detect(
      span_key,
      stringr::regex(
        paste0("(?<![a-z])", escape_regex(country_key), "(?![a-z])"),
        ignore_case = TRUE
      )
    )
  )
}

span_quality_label <- function(span, country) {
  span_key <- normalize_match_text(span)
  country_key <- normalize_country_key(country)

  dplyr::case_when(
    is.na(country) | country == "" ~ "no_country",
    is.na(span) | span == "" ~ "missing_span",
    !is.na(span_key) & !is.na(country_key) & span_key == country_key ~ "country_name_only",
    contains_country_literal(span, country) ~ "country_literal",
    stringr::str_detect(
      dplyr::coalesce(span_key, ""),
      stringr::regex(
        "reported in|reported from|outbreak in|confirmed in|identified in|detected in|notified who|notified who of|ministry of health|national ihr focal point",
        ignore_case = TRUE
      )
    ) ~ "event_phrase_without_country_literal",
    stringr::str_detect(
      dplyr::coalesce(span_key, ""),
      stringr::regex(
        "province|city|region|district|state|mainland|travel|visited|returned|origin|tour",
        ignore_case = TRUE
      )
    ) ~ "indirect_location_context",
    TRUE ~ "other"
  )
}

reasoning_rank <- function(reasoning_label) {
  dplyr::case_when(
    reasoning_label == "explicit_event_country" ~ 1L,
    reasoning_label == "background_only" ~ 2L,
    reasoning_label == "no_country" ~ 9L,
    TRUE ~ 5L
  )
}

span_quality_rank <- function(span_quality) {
  dplyr::case_when(
    span_quality == "country_literal" ~ 1L,
    span_quality == "event_phrase_without_country_literal" ~ 2L,
    span_quality == "indirect_location_context" ~ 3L,
    span_quality == "other" ~ 4L,
    span_quality == "country_name_only" ~ 5L,
    span_quality == "missing_span" ~ 6L,
    span_quality == "no_country" ~ 9L,
    TRUE ~ 8L
  )
}

needs_manual_review_flag <- function(span_quality, reasoning_label, confidence, span) {
  dplyr::case_when(
    reasoning_label == "no_country" ~ FALSE,
    confidence %in% c("medium", "low") ~ TRUE,
    span_quality %in% c("country_name_only", "indirect_location_context", "other", "missing_span") ~ TRUE,
    is.na(span) | nchar(span) < 10 ~ TRUE,
    TRUE ~ FALSE
  )
}

ranked_candidates <- candidates %>%
  left_join(
    records %>%
      select(record_key, DonId, record_id, publication_datetime_utc, PublicationDateAndTime),
    by = "record_key"
  ) %>%
  mutate(
    span_quality = span_quality_label(evidence_span, country_standard),
    span_contains_country_literal = contains_country_literal(evidence_span, country_standard),
    confidence_rank = who_don_confidence_rank(confidence),
    reasoning_rank = reasoning_rank(reasoning_label),
    span_quality_rank = span_quality_rank(span_quality),
    span_nchar = nchar(dplyr::coalesce(evidence_span, "")),
    needs_manual_review = needs_manual_review_flag(
      span_quality = span_quality,
      reasoning_label = reasoning_label,
      confidence = confidence,
      span = evidence_span
    ),
    ranking_score = dplyr::case_when(
      reasoning_label == "explicit_event_country" & confidence == "high" & span_quality == "country_literal" ~ 1L,
      reasoning_label == "explicit_event_country" & confidence == "high" & span_quality == "event_phrase_without_country_literal" ~ 2L,
      reasoning_label == "explicit_event_country" & confidence == "high" & span_quality == "indirect_location_context" ~ 3L,
      reasoning_label == "explicit_event_country" & confidence == "high" ~ 4L,
      reasoning_label == "explicit_event_country" & confidence == "medium" ~ 5L,
      reasoning_label == "background_only" & confidence == "high" ~ 6L,
      reasoning_label == "background_only" & confidence == "medium" ~ 7L,
      reasoning_label == "explicit_event_country" & confidence == "low" ~ 8L,
      reasoning_label == "no_country" ~ 99L,
      TRUE ~ 50L
    )
  ) %>%
  arrange(record_key, country_standard, ranking_score, span_quality_rank, confidence_rank)

best_country_rows <- ranked_candidates %>%
  filter(!is.na(country_standard), country_standard != "") %>%
  group_by(record_key, country_standard) %>%
  summarise(
    DonId = first(stats::na.omit(DonId)),
    record_id = first(stats::na.omit(record_id)),
    Title = first(stats::na.omit(Title)),
    article_url = first(stats::na.omit(article_url)),
    publication_datetime_utc = first(stats::na.omit(publication_datetime_utc)),
    best_evidence_span = first(stats::na.omit(evidence_span)),
    best_reasoning_label = first(stats::na.omit(reasoning_label)),
    best_confidence = first(stats::na.omit(confidence)),
    best_span_quality = first(stats::na.omit(span_quality)),
    span_contains_country_literal = first(span_contains_country_literal),
    needs_manual_review = any(needs_manual_review, na.rm = TRUE),
    n_candidate_spans = dplyr::n(),
    all_evidence_spans = safe_collapse_unique(evidence_span),
    all_reasoning_labels = safe_collapse_unique(reasoning_label),
    all_confidences = safe_collapse_unique(confidence),
    all_span_qualities = safe_collapse_unique(span_quality),
    .groups = "drop"
  )

record_summary <- ranked_candidates %>%
  group_by(record_key) %>%
  summarise(
    DonId = first(stats::na.omit(DonId)),
    record_id = first(stats::na.omit(record_id)),
    Title = first(stats::na.omit(Title)),
    article_url = first(stats::na.omit(article_url)),
    publication_datetime_utc = first(stats::na.omit(publication_datetime_utc)),
    llm_country_count = dplyr::n_distinct(country_standard[!is.na(country_standard) & country_standard != ""]),
    llm_country_list = safe_collapse_unique(country_standard[country_standard != ""]),
    any_manual_review = any(needs_manual_review, na.rm = TRUE),
    any_no_country = any(reasoning_label == "no_country", na.rm = TRUE),
    .groups = "drop"
  )

summary_tbl <- tibble::tribble(
  ~metric, ~value,
  "candidate_rows", nrow(ranked_candidates),
  "candidate_records", dplyr::n_distinct(ranked_candidates$record_key),
  "candidate_records_with_country", dplyr::n_distinct(ranked_candidates$record_key[!is.na(ranked_candidates$country_standard) & ranked_candidates$country_standard != ""]),
  "best_country_rows", nrow(best_country_rows),
  "best_records_with_country", dplyr::n_distinct(best_country_rows$record_key),
  "records_no_country", dplyr::n_distinct(record_summary$record_key[record_summary$llm_country_count == 0]),
  "best_rows_needing_manual_review", sum(best_country_rows$needs_manual_review, na.rm = TRUE)
)

readr::write_csv(ranked_candidates, ranked_path, na = "")
readr::write_csv(best_country_rows, best_path, na = "")
readr::write_csv(record_summary, record_summary_path, na = "")
readr::write_csv(summary_tbl, summary_path, na = "")

message_if(
  "Saved ranked OpenAI WHO DON candidates to: ",
  ranked_path,
  " | rows=",
  nrow(ranked_candidates),
  verbose = verbose
)
message_if(
  "Saved best OpenAI WHO DON country table to: ",
  best_path,
  " | rows=",
  nrow(best_country_rows),
  verbose = verbose
)
message_if(
  "Saved OpenAI WHO DON record summary to: ",
  record_summary_path,
  verbose = verbose
)
message_if(
  "Saved OpenAI WHO DON cleanup summary to: ",
  summary_path,
  verbose = verbose
)
