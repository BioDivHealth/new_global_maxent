# ------------------------------------------------------------------------------
# 5_WHO_DON_Country_Evidence_Union.R
# ------------------------------------------------------------------------------
# Purpose: Combine title-based and text-based WHO DON country evidence into a
#          full audit table and a best-available country summary.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE

title_evidence <- readr::read_csv(
  who_don_output_path("who_don_title_country_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)

pass2_evidence <- readr::read_csv(
  who_don_output_path("who_don_deterministic_pass2_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)

text_evidence <- readr::read_csv(
  who_don_output_path("who_don_text_country_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)

title_union <- title_evidence %>%
  transmute(
    record_key = who_don_record_key(record_id, DonId),
    DonId,
    record_id,
    Title,
    publication_datetime_utc,
    article_url,
    country,
    country_standard,
    section_name = NA_character_,
    matched_text = raw_match,
    context_window = raw_match,
    extraction_source,
    extraction_basis,
    confidence,
    needs_manual_review = FALSE
  )

pass2_union <- pass2_evidence %>%
  transmute(
    record_key,
    DonId,
    record_id,
    Title,
    publication_datetime_utc,
    article_url,
    country,
    country_standard,
    section_name,
    matched_text,
    context_window,
    extraction_source,
    extraction_basis,
    confidence,
    needs_manual_review
  )

text_union <- text_evidence %>%
  transmute(
    record_key = who_don_record_key(record_id, DonId),
    DonId,
    record_id,
    Title,
    publication_datetime_utc,
    article_url,
    country,
    country_standard,
    section_name,
    matched_text,
    context_window,
    extraction_source,
    extraction_basis,
    confidence,
    needs_manual_review
  )

union_tbl <- bind_rows(title_union, pass2_union, text_union) %>%
  filter(!is.na(country_standard)) %>%
  distinct()

best_tbl <- union_tbl %>%
  mutate(
    ranking_score = case_when(
      extraction_source == "title" & confidence == "high" ~ 1L,
      extraction_source == "deterministic_pass2" & confidence == "high" ~ 2L,
      extraction_source == "text" & confidence == "high" ~ 3L,
      extraction_source == "deterministic_pass2" & confidence == "medium" ~ 4L,
      extraction_source == "text" & confidence == "medium" ~ 5L,
      extraction_source == "text" & confidence == "low" ~ 6L,
      TRUE ~ 9L
    )
  ) %>%
  group_by(record_key, DonId, record_id, Title, publication_datetime_utc, article_url, country_standard) %>%
  arrange(ranking_score, who_don_confidence_rank(confidence), .by_group = TRUE) %>%
  summarise(
    country = first(stats::na.omit(country)),
    best_confidence = first(confidence),
    best_extraction_source = first(extraction_source),
    all_extraction_sources = safe_collapse_unique(extraction_source),
    all_sections = safe_collapse_unique(section_name),
    needs_manual_review = any(needs_manual_review, na.rm = TRUE),
    .groups = "drop"
  )

union_path <- who_don_output_path("who_don_country_evidence_union.csv")
best_path <- who_don_output_path("who_don_country_evidence_best.csv")

readr::write_csv(union_tbl, union_path, na = "")
readr::write_csv(best_tbl, best_path, na = "")

message_if(
  "Saved WHO DON country evidence union to: ",
  union_path,
  " | rows=",
  nrow(union_tbl),
  verbose = verbose
)
message_if(
  "Saved WHO DON best country evidence to: ",
  best_path,
  " | rows=",
  nrow(best_tbl),
  verbose = verbose
)
