# ------------------------------------------------------------------------------
# 4_1_WHO_DON_Deterministic_Pass2.R
# ------------------------------------------------------------------------------
# Purpose: Run a conservative second deterministic recovery pass on unresolved
#          WHO DON records using old title patterns, first-sentence event
#          phrases, and corroborated article slug hints.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, stringr, tibble, tidyr)

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

unmatched_records <- readr::read_csv(
  who_don_output_path("who_don_country_qa_unmatched_records.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)

gazetteer <- load_who_don_geography_gazetteer()
country_pattern <- build_country_match_pattern(gazetteer)

strict_event_patterns <- c(
  "reported in",
  "reported from",
  "detected in",
  "confirmed in",
  "identified in",
  "notified by",
  "national ihr focal point for",
  "ministry of health of"
)

country_lookup <- gazetteer %>%
  filter(geography_type == "country", !is_ambiguous) %>%
  select(alias, alias_key, country_standard) %>%
  distinct()

region_lookup <- gazetteer %>%
  filter(geography_type == "who_region") %>%
  select(alias, alias_key, geography_type) %>%
  distinct()

current_best_ids <- unique(stats::na.omit(best_evidence$record_key[!is.na(best_evidence$country_standard)]))
qa_unmatched_ids <- unique(stats::na.omit(unmatched_records$record_key))
unresolved_ids <- unique(c(qa_unmatched_ids, setdiff(records$record_key, current_best_ids)))

records_unresolved <- records %>%
  filter(record_key %in% unresolved_ids)

match_exact_country <- function(raw_text) {
  who_don_match_country_alias(raw_text, gazetteer)
}

split_geo_candidates <- function(raw_text) {
  raw_text <- clean_scalar_text(raw_text)

  if (is.na(raw_text)) {
    return(character(0))
  }

  if (nrow(match_exact_country(raw_text)) > 0) {
    return(raw_text)
  }

  try_split <- function(pattern) {
    pieces <- stringr::str_split(raw_text, pattern, simplify = FALSE)[[1]]
    pieces <- clean_scalar_text(pieces)
    pieces <- pieces[!is.na(pieces)]

    if (length(pieces) <= 1) {
      return(NULL)
    }

    ok <- purrr::map_lgl(pieces, ~ nrow(match_exact_country(.x)) > 0)

    if (all(ok)) {
      return(pieces)
    }

    NULL
  }

  for (pattern in c(
    stringr::regex("\\s*;\\s*"),
    stringr::regex("\\s*/\\s*"),
    stringr::regex("\\s*,\\s*"),
    stringr::regex("\\s+and\\s+", ignore_case = TRUE)
  )) {
    split_hit <- try_split(pattern)

    if (!is.null(split_hit)) {
      return(split_hit)
    }
  }

  raw_text
}

extract_title_country_candidates <- function(title) {
  title <- clean_scalar_text(title)
  empty_tbl <- tibble(
    country = character(),
    country_standard = character(),
    alias = character(),
    alias_key = character(),
    extraction_basis = character(),
    raw_match = character()
  )

  if (is.na(title)) {
    return(empty_tbl)
  }

  patterns <- tribble(
    ~pattern, ~basis,
    "^(?:.+?[-–—:]\\s*)?(?:global\\s+)?update\\s*\\((.+)\\)$", "title_update_parentheses",
    "^(?:.+?[-–—:]\\s*)?situation\\s+in\\s+(.+)$", "title_situation_in",
    "^(?:.+?[-–—:]\\s*)?confirmed\\s+in\\s+(.+)$", "title_confirmed_in",
    "^(?:.+?[-–—:]\\s*)?detected\\s+in\\s+(.+)$", "title_detected_in",
    "^(?:.+?[-–—:]\\s*)?reported\\s+in\\s+(.+)$", "title_reported_in"
  )

  purrr::map_dfr(
    seq_len(nrow(patterns)),
    function(i) {
      hit <- stringr::str_match(title, stringr::regex(patterns$pattern[[i]], ignore_case = TRUE))
      captured <- clean_scalar_text(hit[, 2])

      if (is.na(captured)) {
        return(empty_tbl)
      }

      pieces <- split_geo_candidates(captured)

      purrr::map_dfr(
        pieces,
        function(one_piece) {
          matched <- match_exact_country(one_piece)

          if (nrow(matched) == 0) {
            return(empty_tbl)
          }

          matched %>%
            mutate(
              extraction_basis = patterns$basis[[i]],
              raw_match = captured
            )
        }
      )
    }
  ) %>%
    distinct(country_standard, extraction_basis, .keep_all = TRUE)
}

extract_country_mentions <- function(section_text, require_event_phrase = FALSE) {
  section_text <- clean_scalar_text(section_text)
  empty_tbl <- tibble(
    match_start = integer(),
    match_end = integer(),
    matched_text = character(),
    alias_key = character(),
    alias = character(),
    country_standard = character(),
    context_window = character(),
    extraction_basis = character()
  )

  if (is.na(section_text)) {
    return(empty_tbl)
  }

  normalized_text <- normalize_match_text(section_text)

  if (is.na(normalized_text) || normalized_text == "") {
    return(empty_tbl)
  }

  locs <- stringr::str_locate_all(
    normalized_text,
    stringr::regex(country_pattern, ignore_case = TRUE)
  )[[1]]

  if (nrow(locs) == 0) {
    return(empty_tbl)
  }

  tibble(
    match_start = locs[, 1],
    match_end = locs[, 2],
    matched_text = stringr::str_sub(normalized_text, locs[, 1], locs[, 2])
  ) %>%
    mutate(alias_key = normalize_match_text(matched_text)) %>%
    left_join(country_lookup, by = "alias_key") %>%
    filter(!is.na(country_standard)) %>%
    mutate(
      context_window = purrr::map2_chr(
        match_start,
        match_end,
        ~ extract_context_window(normalized_text, .x, .y)
      ),
      extraction_basis = purrr::map_chr(
        context_window,
        function(one_context) {
          hit <- strict_event_patterns[
            stringr::str_detect(one_context, stringr::fixed(strict_event_patterns))
          ]

          if (length(hit) == 0) {
            return("country_mention")
          }

          hit[[1]]
        }
      )
    ) %>%
    filter(!require_event_phrase | extraction_basis != "country_mention") %>%
    distinct(country_standard, extraction_basis, matched_text, .keep_all = TRUE)
}

extract_slug_country_mentions <- function(article_url) {
  slug <- who_don_article_slug(article_url)
  slug_text <- clean_scalar_text(slug)
  empty_tbl <- tibble(
    matched_text = character(),
    alias_key = character(),
    alias = character(),
    country_standard = character()
  )

  if (is.na(slug_text)) {
    return(empty_tbl)
  }

  normalized_slug <- normalize_match_text(stringr::str_replace_all(slug_text, "[-_]", " "))

  if (is.na(normalized_slug) || normalized_slug == "") {
    return(empty_tbl)
  }

  locs <- stringr::str_locate_all(
    normalized_slug,
    stringr::regex(country_pattern, ignore_case = TRUE)
  )[[1]]

  if (nrow(locs) == 0) {
    return(empty_tbl)
  }

  tibble(
    matched_text = stringr::str_sub(normalized_slug, locs[, 1], locs[, 2])
  ) %>%
    mutate(alias_key = normalize_match_text(matched_text)) %>%
    left_join(country_lookup, by = "alias_key") %>%
    filter(!is.na(country_standard)) %>%
    distinct(country_standard, matched_text)
}

extract_region_hits <- function(section_text) {
  section_text <- clean_scalar_text(section_text)

  if (is.na(section_text)) {
    return(character(0))
  }

  normalized_text <- normalize_match_text(section_text)

  region_lookup %>%
    filter(!is.na(alias_key)) %>%
    filter(stringr::str_detect(normalized_text, stringr::fixed(alias_key))) %>%
    pull(alias) %>%
    unique()
}

build_review_row <- function(record_row, title_hits, first_summary_hits, first_overview_hits, slug_hits) {
  title_suffix <- who_don_title_suffix(record_row$Title)
  summary_first <- who_don_first_sentences(record_row$summary_text, n_sentences = 2L)
  overview_first <- who_don_first_sentences(record_row$overview_text, n_sentences = 2L)
  region_hint <- who_don_region_hint(record_row$Title, title_suffix)

  tibble(
    record_key = record_row$record_key,
    DonId = record_row$DonId,
    record_id = record_row$record_id,
    Title = record_row$Title,
    publication_datetime_utc = record_row$publication_datetime_utc,
    article_url = record_row$article_url,
    title_suffix = title_suffix,
    region_hint = region_hint,
    slug_hint = who_don_article_slug(record_row$article_url),
    summary_first = summary_first,
    overview_first = overview_first,
    title_country_candidates = safe_collapse_unique(title_hits$country_standard),
    summary_first_countries = safe_collapse_unique(first_summary_hits$country_standard),
    overview_first_countries = safe_collapse_unique(first_overview_hits$country_standard),
    slug_country_candidates = safe_collapse_unique(slug_hits$country_standard),
    summary_region_hits = safe_collapse_unique(extract_region_hits(summary_first)),
    overview_region_hits = safe_collapse_unique(extract_region_hits(overview_first))
  )
}

message_if(
  "Running WHO DON deterministic pass 2 on unresolved records: ",
  nrow(records_unresolved),
  verbose = verbose
)

pass2_rows <- purrr::map_dfr(
  seq_len(nrow(records_unresolved)),
  function(i) {
    one_record <- records_unresolved[i, ]
    title_suffix <- who_don_title_suffix(one_record$Title)
    region_hint <- who_don_region_hint(one_record$Title, title_suffix)
    summary_first <- who_don_first_sentences(one_record$summary_text, n_sentences = 2L)
    overview_first <- who_don_first_sentences(one_record$overview_text, n_sentences = 2L)

    title_hits <- extract_title_country_candidates(one_record$Title) %>%
      transmute(
        record_key = one_record$record_key,
        DonId = one_record$DonId,
        record_id = one_record$record_id,
        Title = one_record$Title,
        publication_datetime_utc = one_record$publication_datetime_utc,
        article_url = one_record$article_url,
        country = country_standard,
        country_standard,
        section_name = "title",
        matched_text = raw_match,
        context_window = one_record$Title,
        extraction_source = "deterministic_pass2",
        extraction_basis,
        confidence = "high",
        needs_manual_review = FALSE
      )

    summary_hits <- extract_country_mentions(summary_first, require_event_phrase = TRUE) %>%
      transmute(
        record_key = one_record$record_key,
        DonId = one_record$DonId,
        record_id = one_record$record_id,
        Title = one_record$Title,
        publication_datetime_utc = one_record$publication_datetime_utc,
        article_url = one_record$article_url,
        country = matched_text,
        country_standard,
        section_name = "summary_first_sentences",
        matched_text,
        context_window,
        extraction_source = "deterministic_pass2",
        extraction_basis = paste0("summary_first_sentences:", extraction_basis),
        confidence = "high",
        needs_manual_review = FALSE
      )

    overview_hits <- extract_country_mentions(overview_first, require_event_phrase = TRUE) %>%
      transmute(
        record_key = one_record$record_key,
        DonId = one_record$DonId,
        record_id = one_record$record_id,
        Title = one_record$Title,
        publication_datetime_utc = one_record$publication_datetime_utc,
        article_url = one_record$article_url,
        country = matched_text,
        country_standard,
        section_name = "overview_first_sentences",
        matched_text,
        context_window,
        extraction_source = "deterministic_pass2",
        extraction_basis = paste0("overview_first_sentences:", extraction_basis),
        confidence = "high",
        needs_manual_review = FALSE
      )

    slug_hits <- extract_slug_country_mentions(one_record$article_url)
    corroborated_countries <- unique(c(
      title_hits$country_standard,
      extract_country_mentions(summary_first, require_event_phrase = FALSE)$country_standard,
      extract_country_mentions(overview_first, require_event_phrase = FALSE)$country_standard
    ))

    slug_rows <- slug_hits %>%
      filter(country_standard %in% corroborated_countries) %>%
      transmute(
        record_key = one_record$record_key,
        DonId = one_record$DonId,
        record_id = one_record$record_id,
        Title = one_record$Title,
        publication_datetime_utc = one_record$publication_datetime_utc,
        article_url = one_record$article_url,
        country = country_standard,
        country_standard,
        section_name = "article_url_slug",
        matched_text = matched_text,
        context_window = paste(
          "slug:",
          who_don_article_slug(one_record$article_url),
          "| title:",
          clean_scalar_text(one_record$Title),
          "| summary_first:",
          clean_scalar_text(summary_first)
        ),
        extraction_source = "deterministic_pass2",
        extraction_basis = "slug_plus_title_or_first_sentence",
        confidence = "medium",
        needs_manual_review = FALSE
      )

    all_rows <- bind_rows(title_hits, summary_hits, overview_hits, slug_rows) %>%
      distinct(record_key, country_standard, section_name, extraction_basis, confidence, .keep_all = TRUE)

    if (!is.na(region_hint) && nrow(all_rows) == 0) {
      return(tibble())
    }

    all_rows
  }
)

pass2_rows <- pass2_rows %>%
  filter(!is.na(country_standard)) %>%
  distinct()

review_tbl <- purrr::map_dfr(
  seq_len(nrow(records_unresolved)),
  function(i) {
    one_record <- records_unresolved[i, ]
    title_hits <- extract_title_country_candidates(one_record$Title)
    summary_hits <- extract_country_mentions(
      who_don_first_sentences(one_record$summary_text, n_sentences = 2L),
      require_event_phrase = FALSE
    )
    overview_hits <- extract_country_mentions(
      who_don_first_sentences(one_record$overview_text, n_sentences = 2L),
      require_event_phrase = FALSE
    )
    slug_hits <- extract_slug_country_mentions(one_record$article_url)

    build_review_row(one_record, title_hits, summary_hits, overview_hits, slug_hits)
  }
)

evidence_path <- who_don_output_path("who_don_deterministic_pass2_evidence.csv")
review_path <- who_don_output_path("who_don_deterministic_pass2_review.csv")

readr::write_csv(pass2_rows, evidence_path, na = "")
readr::write_csv(review_tbl, review_path, na = "")

message_if(
  "Saved WHO DON deterministic pass 2 evidence to: ",
  evidence_path,
  " | rows=",
  nrow(pass2_rows),
  verbose = verbose
)
message_if(
  "Saved WHO DON deterministic pass 2 review table to: ",
  review_path,
  verbose = verbose
)
