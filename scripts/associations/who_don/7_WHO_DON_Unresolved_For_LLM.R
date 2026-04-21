# ------------------------------------------------------------------------------
# 7_WHO_DON_Unresolved_For_LLM.R
# ------------------------------------------------------------------------------
# Purpose: Export a compact unresolved WHO DON subset for later OpenAI batch
#          extraction with structured outputs, while excluding poor LLM targets.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, jsonlite, purrr, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
records <- read_who_don_records("who_don_records_clean.csv") %>%
  who_don_with_record_key()

best_evidence <- readr::read_csv(
  who_don_output_path("who_don_country_evidence_best.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  who_don_with_record_key()

gazetteer <- load_who_don_geography_gazetteer()
country_pattern <- build_country_match_pattern(gazetteer)
ambiguous_terms <- gazetteer %>%
  filter(is_ambiguous) %>%
  pull(alias_key) %>%
  unique()

country_allowlist <- gazetteer %>%
  filter(geography_type == "country", !is_ambiguous, !is.na(country_standard)) %>%
  pull(country_standard) %>%
  unique() %>%
  sort()

matched_ids <- unique(stats::na.omit(best_evidence$record_key[!is.na(best_evidence$country_standard)]))
unresolved_records <- records %>%
  filter(!record_key %in% matched_ids)

text_has_country <- function(text) {
  text <- clean_scalar_text(text)

  vapply(
    text,
    FUN.VALUE = logical(1),
    function(one_text) {
      if (is.na(one_text)) {
        return(FALSE)
      }

      normalized_text <- normalize_match_text(one_text)

      if (is.na(normalized_text) || normalized_text == "") {
        return(FALSE)
      }

      stringr::str_detect(
        normalized_text,
        stringr::regex(country_pattern, ignore_case = TRUE)
      )
    }
  )
}

text_has_ambiguous_only <- function(text) {
  text <- clean_scalar_text(text)

  vapply(
    text,
    FUN.VALUE = logical(1),
    function(one_text) {
      if (is.na(one_text)) {
        return(FALSE)
      }

      normalized_text <- normalize_match_text(one_text)

      if (is.na(normalized_text) || normalized_text == "" || text_has_country(one_text)) {
        return(FALSE)
      }

      any(stringr::str_detect(normalized_text, stringr::fixed(ambiguous_terms)))
    }
  )
}

llm_candidates <- unresolved_records %>%
  mutate(
    title_suffix = who_don_title_suffix(Title),
    region_hint = who_don_region_hint(Title, title_suffix),
    slug_hint = who_don_article_slug(article_url),
    summary_first = who_don_first_sentences(summary_text, n_sentences = 2L),
    overview_first = who_don_first_sentences(overview_text, n_sentences = 2L),
    response_first = who_don_first_sentences(response_text, n_sentences = 2L),
    summary_text_short = truncate_text(summary_text, max_chars = 1500L),
    overview_text_short = truncate_text(overview_text, max_chars = 1500L),
    response_text_short = truncate_text(response_text, max_chars = 1000L),
    combined_text_nchar = nchar(dplyr::coalesce(summary_text_short, "")) +
      nchar(dplyr::coalesce(overview_text_short, "")),
    has_country_in_first_sentences = text_has_country(summary_first) |
      text_has_country(overview_first),
    has_country_anywhere = text_has_country(summary_text_short) |
      text_has_country(overview_text_short) |
      text_has_country(response_text_short),
    has_ambiguous_only = text_has_ambiguous_only(summary_first) |
      text_has_ambiguous_only(overview_first),
    is_region_or_global = !is.na(region_hint),
    is_multi_country_only = stringr::str_detect(
      normalize_match_text(dplyr::coalesce(title_suffix, Title)),
      stringr::regex("multi country|multi-country", ignore_case = TRUE)
    ),
    exclusion_bucket = case_when(
      combined_text_nchar < 80 ~ "too_little_text",
      is_region_or_global & !has_country_in_first_sentences ~ "global_or_region_only",
      is_multi_country_only & !has_country_in_first_sentences ~ "multi_country_only",
      has_ambiguous_only & !has_country_anywhere ~ "ambiguous_only",
      TRUE ~ NA_character_
    ),
    llm_status = ifelse(is.na(exclusion_bucket), "pending", NA_character_)
  )

llm_input <- llm_candidates %>%
  filter(is.na(exclusion_bucket)) %>%
  transmute(
    record_key,
    DonId,
    record_id,
    Title,
    publication_datetime_utc,
    article_url,
    summary_text_short,
    overview_text_short,
    response_text_short,
    slug_hint,
    region_hint,
    llm_status = "pending"
  )

llm_jsonl_records <- llm_candidates %>%
  filter(is.na(exclusion_bucket)) %>%
  mutate(
    summary_text = truncate_text(summary_text, max_chars = 2500L),
    overview_text = truncate_text(overview_text, max_chars = 2500L),
    response_text = truncate_text(response_text, max_chars = 1500L)
  ) %>%
  transmute(
    payload = purrr::pmap(
      list(
        record_key,
        Title,
        article_url,
        summary_text,
        overview_text,
        response_text,
        slug_hint,
        region_hint
      ),
      function(
        record_key,
        Title,
        article_url,
        summary_text,
        overview_text,
        response_text,
        slug_hint,
        region_hint
      ) {
        list(
          record_key = record_key,
          llm_status = "pending",
          extraction_contract_version = "who_don_country_v1",
          model_provider = "openai",
          input = list(
            Title = Title,
            article_url = article_url,
            summary_text = summary_text,
            overview_text = overview_text,
            response_text = response_text,
            slug_hint = slug_hint,
            region_hint = region_hint
          ),
          expected_output_schema = list(
            record_key = "string",
            has_country_evidence = "boolean",
            countries = "array<string>",
            evidence_spans = "array<string>",
            reasoning_label = c(
              "explicit_event_country",
              "background_only",
              "regional_only",
              "no_country"
            ),
            confidence = c("high", "medium", "low")
          ),
          constraints = list(
            country_allowlist = country_allowlist,
            allow_no_country = TRUE,
            require_verbatim_evidence_spans = TRUE
          )
        )
      }
    )
  )

jsonl_lines <- purrr::map_chr(
  llm_jsonl_records$payload,
  ~ jsonlite::toJSON(.x, auto_unbox = TRUE, null = "null")
)

exclusion_summary <- tibble(
  exclusion_bucket = c(
    "global_or_region_only",
    "multi_country_only",
    "too_little_text",
    "ambiguous_only",
    "still_candidate_for_manual_review"
  )
) %>%
  left_join(
    llm_candidates %>%
      filter(!is.na(exclusion_bucket)) %>%
      count(exclusion_bucket, name = "n_records"),
    by = "exclusion_bucket"
  ) %>%
  mutate(n_records = dplyr::coalesce(n_records, 0L))

csv_path <- who_don_output_path("who_don_unresolved_llm_input.csv")
jsonl_path <- who_don_output_path("who_don_unresolved_llm_input.jsonl")
summary_path <- who_don_output_path("who_don_unresolved_llm_exclusion_summary.csv")

readr::write_csv(llm_input, csv_path, na = "")
writeLines(jsonl_lines, jsonl_path)
readr::write_csv(exclusion_summary, summary_path, na = "")

message_if(
  "Saved WHO DON unresolved LLM CSV input to: ",
  csv_path,
  " | rows=",
  nrow(llm_input),
  verbose = verbose
)
message_if(
  "Saved WHO DON unresolved LLM JSONL input to: ",
  jsonl_path,
  verbose = verbose
)
message_if(
  "Saved WHO DON unresolved LLM exclusion summary to: ",
  summary_path,
  verbose = verbose
)
