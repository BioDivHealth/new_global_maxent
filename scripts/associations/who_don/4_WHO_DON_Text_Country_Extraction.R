# ------------------------------------------------------------------------------
# 4_WHO_DON_Text_Country_Extraction.R
# ------------------------------------------------------------------------------
# Purpose: Extract country mentions from cleaned WHO DON narrative text, score
#          them by outbreak-context phrases, and preserve section-level
#          provenance for later review.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
records <- read_who_don_records("who_don_records_clean.csv")
gazetteer <- load_who_don_geography_gazetteer()
country_pattern <- build_country_match_pattern(gazetteer)
section_cols <- who_don_section_columns()
event_patterns <- who_don_event_patterns()

country_lookup <- gazetteer %>%
  filter(geography_type == "country", !is_ambiguous) %>%
  select(alias, alias_key, country_standard) %>%
  distinct()

ambiguous_terms <- gazetteer %>%
  filter(is_ambiguous) %>%
  pull(alias_key) %>%
  unique()

extract_section_country_hits <- function(record_row, section_name, section_text) {
  if (is.na(section_text)) {
    return(tibble())
  }

  normalized_text <- normalize_match_text(section_text)

  if (is.na(normalized_text) || normalized_text == "") {
    return(tibble())
  }

  locs <- stringr::str_locate_all(
    normalized_text,
    stringr::regex(country_pattern, ignore_case = TRUE)
  )[[1]]

  if (nrow(locs) == 0) {
    return(tibble())
  }

  matches <- stringr::str_sub(normalized_text, locs[, 1], locs[, 2])

  hits <- tibble(
    match_start = locs[, 1],
    match_end = locs[, 2],
    matched_text = matches
  ) %>%
    mutate(
      alias_key = normalize_match_text(matched_text)
    ) %>%
    left_join(country_lookup, by = "alias_key") %>%
    mutate(
      context_window = purrr::map2_chr(
        match_start,
        match_end,
        ~ extract_context_window(normalized_text, .x, .y)
      ),
      extraction_basis = purrr::map_chr(
        context_window,
        function(one_context) {
          hit <- event_patterns[stringr::str_detect(one_context, stringr::fixed(event_patterns))]

          if (length(hit) == 0) {
            return("country_mention")
          }

          hit[[1]]
        }
      ),
      confidence = case_when(
        extraction_basis != "country_mention" ~ "high",
        section_name %in% c("summary_text", "overview_text", "epidemiology_text", "response_text") ~ "medium",
        TRUE ~ "low"
      ),
      needs_manual_review = confidence == "low",
      DonId = record_row$DonId,
      record_id = dplyr::coalesce(record_row$record_id, record_row$DonId, record_row$Id),
      Title = record_row$Title,
      publication_datetime_utc = record_row$publication_datetime_utc,
      article_url = record_row$article_url,
      section_name = section_name,
      extraction_source = "text"
    ) %>%
    filter(!is.na(country_standard)) %>%
    transmute(
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      article_url,
      country = matched_text,
      country_standard,
      section_name,
      matched_text,
      context_window,
      extraction_source,
      extraction_basis,
      confidence,
      needs_manual_review
    ) %>%
    distinct()

  hits
}

message_if("Extracting text-based WHO DON country evidence...", verbose = verbose)

text_hits <- purrr::map_dfr(
  seq_len(nrow(records)),
  function(i) {
    one_record <- records[i, ]

    purrr::map_dfr(
      section_cols,
      function(section_name) {
        extract_section_country_hits(
          record_row = one_record,
          section_name = section_name,
          section_text = one_record[[section_name]][[1]]
        )
      }
    )
  }
)

text_hits <- text_hits %>%
  distinct(
    DonId,
    country_standard,
    section_name,
    extraction_basis,
    matched_text,
    .keep_all = TRUE
  )

ambiguous_mentions <- purrr::map_dfr(
  seq_len(nrow(records)),
  function(i) {
    one_record <- records[i, ]

    purrr::map_dfr(
      section_cols,
      function(section_name) {
        section_text <- one_record[[section_name]][[1]]

        if (is.na(section_text)) {
          return(tibble())
        }

        normalized_text <- normalize_match_text(section_text)

        if (is.na(normalized_text) || normalized_text == "") {
          return(tibble())
        }

        purrr::map_dfr(
          ambiguous_terms,
          function(term_key) {
            locs <- stringr::str_locate_all(
              normalized_text,
              stringr::regex(
                paste0("(?<![a-z])", escape_regex(term_key), "(?![a-z])"),
                ignore_case = TRUE
              )
            )[[1]]

            if (nrow(locs) == 0) {
              return(tibble())
            }

            tibble(
              DonId = one_record$DonId,
              record_id = dplyr::coalesce(one_record$record_id, one_record$DonId, one_record$Id),
              Title = one_record$Title,
              section_name = section_name,
              ambiguous_text = term_key,
              context_window = purrr::map2_chr(
                locs[, 1],
                locs[, 2],
                ~ extract_context_window(normalized_text, .x, .y)
              )
            )
          }
        )
      }
    )
  }
)

text_path <- who_don_output_path("who_don_text_country_evidence.csv")
ambiguous_path <- who_don_output_path("who_don_country_qa_ambiguous_mentions.csv")

readr::write_csv(text_hits, text_path, na = "")
readr::write_csv(ambiguous_mentions, ambiguous_path, na = "")

message_if(
  "Saved text-based WHO DON evidence to: ",
  text_path,
  " | rows=",
  nrow(text_hits),
  verbose = verbose
)
