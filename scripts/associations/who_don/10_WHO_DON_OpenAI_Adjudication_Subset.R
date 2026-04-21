# ------------------------------------------------------------------------------
# 10_WHO_DON_OpenAI_Adjudication_Subset.R
# ------------------------------------------------------------------------------
# Purpose: Build a stricter second-pass OpenAI adjudication subset from the
#          uncertain WHO DON LLM outputs.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, jsonlite, purrr, readr, stringr, tibble)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
batch_dir <- who_don_output_path("openai_batch")

records_path <- who_don_output_path("who_don_records_clean.csv")
best_path <- file.path(batch_dir, "who_don_openai_country_best_full.csv")
record_summary_path <- file.path(batch_dir, "who_don_openai_country_record_summary_full.csv")

adjudication_csv_path <- file.path(batch_dir, "who_don_openai_adjudication_input.csv")
adjudication_jsonl_path <- file.path(batch_dir, "who_don_openai_adjudication_input.jsonl")
adjudication_summary_path <- file.path(batch_dir, "who_don_openai_adjudication_summary.csv")

include_no_country <- tolower(Sys.getenv("WHO_DON_INCLUDE_NO_COUNTRY", "true")) %in% c(
  "true", "1", "yes", "y"
)
max_text_chars <- suppressWarnings(
  as.integer(Sys.getenv("WHO_DON_ADJUDICATION_TEXT_CHARS", "12000"))
)
if (is.na(max_text_chars) || max_text_chars <= 0) {
  max_text_chars <- 12000L
}

records <- readr::read_csv(records_path, show_col_types = FALSE, na = c("", "NA")) %>%
  who_don_with_record_key() %>%
  mutate(
    publication_datetime_utc = dplyr::coalesce(publication_datetime_utc, PublicationDateAndTime),
    summary_text = truncate_text(summary_text, max_chars = max_text_chars),
    overview_text = truncate_text(overview_text, max_chars = max_text_chars),
    response_text = truncate_text(response_text, max_chars = max_text_chars)
  )

best_rows <- readr::read_csv(best_path, show_col_types = FALSE, na = c("", "NA"))
record_summary <- readr::read_csv(record_summary_path, show_col_types = FALSE, na = c("", "NA"))

review_records <- record_summary %>%
  mutate(
    include_for_adjudication = any_manual_review | (include_no_country & llm_country_count == 0),
    adjudication_reason = dplyr::case_when(
      any_manual_review & llm_country_count == 0 ~ "manual_review_and_no_country",
      any_manual_review ~ "manual_review",
      llm_country_count == 0 ~ "no_country",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(include_for_adjudication) %>%
  select(
    record_key,
    DonId,
    record_id,
    Title,
    article_url,
    publication_datetime_utc,
    llm_country_count,
    llm_country_list,
    any_manual_review,
    any_no_country,
    adjudication_reason
  )

best_rows_nested <- best_rows %>%
  mutate(
    best_evidence_span = dplyr::coalesce(best_evidence_span, ""),
    best_reasoning_label = dplyr::coalesce(best_reasoning_label, ""),
    best_confidence = dplyr::coalesce(best_confidence, ""),
    best_span_quality = dplyr::coalesce(best_span_quality, ""),
    needs_manual_review = dplyr::coalesce(needs_manual_review, FALSE)
  ) %>%
  group_by(record_key) %>%
  summarise(
    current_best_country_rows = list(
      purrr::transpose(list(
        country_standard = country_standard,
        best_evidence_span = best_evidence_span,
        best_reasoning_label = best_reasoning_label,
        best_confidence = best_confidence,
        best_span_quality = best_span_quality,
        needs_manual_review = needs_manual_review
      ))
    ),
    manual_review_country_rows = list(
      purrr::transpose(list(
        country_standard = country_standard[needs_manual_review],
        best_evidence_span = best_evidence_span[needs_manual_review],
        best_reasoning_label = best_reasoning_label[needs_manual_review],
        best_confidence = best_confidence[needs_manual_review],
        best_span_quality = best_span_quality[needs_manual_review],
        needs_manual_review = needs_manual_review[needs_manual_review]
      ))
    ),
    current_best_country_list = safe_collapse_unique(country_standard),
    manual_review_country_count = sum(needs_manual_review, na.rm = TRUE),
    .groups = "drop"
  )

adjudication_records <- review_records %>%
  left_join(
    records %>%
      select(
        record_key,
        Title,
        article_url,
        publication_datetime_utc,
        summary_text,
        overview_text,
        response_text
      ),
    by = "record_key",
    suffix = c("_summary", "")
  ) %>%
  left_join(best_rows_nested, by = "record_key") %>%
  mutate(
    Title = dplyr::coalesce(Title, Title_summary),
    article_url = dplyr::coalesce(article_url, article_url_summary),
    publication_datetime_utc = dplyr::coalesce(
      publication_datetime_utc,
      publication_datetime_utc_summary
    ),
    current_best_country_rows = purrr::map(
      current_best_country_rows,
      ~ if (is.null(.x)) list() else .x
    ),
    manual_review_country_rows = purrr::map(
      manual_review_country_rows,
      ~ if (is.null(.x)) list() else .x
    ),
    current_best_country_list = dplyr::coalesce(current_best_country_list, ""),
    manual_review_country_count = dplyr::coalesce(manual_review_country_count, 0L),
    adjudication_notes = dplyr::case_when(
      adjudication_reason == "manual_review_and_no_country" ~
        "Strict adjudication: validate the uncertain country evidence and decide whether any country should be retained.",
      adjudication_reason == "manual_review" ~
        "Strict adjudication: keep only countries with strong explicit support and return one best span per retained country.",
      adjudication_reason == "no_country" ~
        "Strict adjudication: confirm whether the record truly has no explicit country evidence, or recover explicit countries if present.",
      TRUE ~
        "Strict adjudication: keep only strongly supported explicit countries."
    )
  ) %>%
  select(
    record_key,
    DonId,
    record_id,
    Title,
    article_url,
    publication_datetime_utc,
    llm_country_count,
    llm_country_list,
    any_manual_review,
    any_no_country,
    adjudication_reason,
    current_best_country_list,
    manual_review_country_count,
    summary_text,
    overview_text,
    response_text,
    current_best_country_rows,
    manual_review_country_rows,
    adjudication_notes
  ) %>%
  arrange(desc(any_manual_review), desc(any_no_country), record_key)

jsonl_rows <- purrr::pmap(
  adjudication_records,
  function(
    record_key,
    DonId,
    record_id,
    Title,
    article_url,
    publication_datetime_utc,
    llm_country_count,
    llm_country_list,
    any_manual_review,
    any_no_country,
    adjudication_reason,
    current_best_country_list,
    manual_review_country_count,
    summary_text,
    overview_text,
    response_text,
    current_best_country_rows,
    manual_review_country_rows,
    adjudication_notes
  ) {
    list(
      record_key = record_key,
      llm_status = "pending",
      extraction_contract_version = "who_don_country_v2_adjudication",
      model_provider = "openai",
      input = list(
        Title = Title,
        article_url = article_url,
        summary_text = summary_text,
        overview_text = overview_text,
        response_text = response_text,
        llm_country_count = llm_country_count,
        current_best_country_list = current_best_country_list,
        current_best_country_rows = current_best_country_rows,
        manual_review_country_rows = manual_review_country_rows,
        adjudication_reason = adjudication_reason,
        adjudication_notes = adjudication_notes
      )
    )
  }
)

summary_tbl <- tibble::tribble(
  ~metric, ~value,
  "records_for_adjudication", nrow(adjudication_records),
  "manual_review_records", sum(adjudication_records$any_manual_review, na.rm = TRUE),
  "no_country_records", sum(adjudication_records$llm_country_count == 0, na.rm = TRUE),
  "include_no_country", include_no_country,
  "max_text_chars", max_text_chars
)

readr::write_csv(adjudication_records, adjudication_csv_path, na = "")
writeLines(
  purrr::map_chr(
    jsonl_rows,
    ~ jsonlite::toJSON(.x, auto_unbox = TRUE, null = "null", na = "null")
  ),
  adjudication_jsonl_path
)
readr::write_csv(summary_tbl, adjudication_summary_path, na = "")

message_if(
  "Saved WHO DON adjudication input CSV to: ",
  adjudication_csv_path,
  " | rows=",
  nrow(adjudication_records),
  verbose = verbose
)
message_if(
  "Saved WHO DON adjudication input JSONL to: ",
  adjudication_jsonl_path,
  " | rows=",
  length(jsonl_rows),
  verbose = verbose
)
message_if(
  "Saved WHO DON adjudication summary to: ",
  adjudication_summary_path,
  verbose = verbose
)
