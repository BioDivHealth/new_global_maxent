# ------------------------------------------------------------------------------
# 6_WHO_DON_QA_Summary.R
# ------------------------------------------------------------------------------
# Purpose: Summarize WHO DON country-recovery coverage and write review-focused
#          QA tables for unmatched and ambiguous records.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE

records <- read_who_don_records("who_don_records_clean.csv")
records <- records %>%
  who_don_with_record_key()
title_evidence <- readr::read_csv(
  who_don_output_path("who_don_title_country_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)
title_evidence <- title_evidence %>%
  who_don_with_record_key()
pass2_evidence <- readr::read_csv(
  who_don_output_path("who_don_deterministic_pass2_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)
pass2_evidence <- pass2_evidence %>%
  who_don_with_record_key()
text_evidence <- readr::read_csv(
  who_don_output_path("who_don_text_country_evidence.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)
text_evidence <- text_evidence %>%
  who_don_with_record_key()
best_evidence <- readr::read_csv(
  who_don_output_path("who_don_country_evidence_best.csv"),
  show_col_types = FALSE,
  na = c("", "NA")
)
best_evidence <- best_evidence %>%
  who_don_with_record_key()

title_country_ids <- unique(stats::na.omit(title_evidence$record_key[!is.na(title_evidence$country_standard)]))
pass2_country_ids <- unique(stats::na.omit(pass2_evidence$record_key[!is.na(pass2_evidence$country_standard)]))
text_country_ids <- unique(stats::na.omit(text_evidence$record_key[!is.na(text_evidence$country_standard)]))
final_country_ids <- unique(stats::na.omit(best_evidence$record_key[!is.na(best_evidence$country_standard)]))

summary_tbl <- tibble::tribble(
  ~metric, ~value,
  "total_records", dplyr::n_distinct(records$record_key),
  "records_with_title_country_evidence", length(title_country_ids),
  "records_with_deterministic_pass2_country_evidence", length(pass2_country_ids),
  "records_with_any_text_country_evidence", length(text_country_ids),
  "records_with_any_final_country_evidence", length(final_country_ids),
  "records_still_unmatched", dplyr::n_distinct(records$record_key) - length(final_country_ids)
)

unmatched_records <- records %>%
  filter(!record_key %in% final_country_ids) %>%
  transmute(
    record_key,
    DonId,
    record_id,
    Title,
    PublicationDateAndTime,
    article_url,
    title_suffix = who_don_title_suffix(Title),
    region_hint = who_don_region_hint(Title, who_don_title_suffix(Title)),
    slug_hint = who_don_article_slug(article_url)
  )

unmatched_titles <- unmatched_records %>%
  mutate(title_suffix_key = normalize_match_text(title_suffix)) %>%
  count(title_suffix, title_suffix_key, sort = TRUE, name = "n_records")

summary_path <- who_don_output_path("who_don_country_qa_summary.csv")
unmatched_titles_path <- who_don_output_path("who_don_country_qa_unmatched_titles.csv")
unmatched_records_path <- who_don_output_path("who_don_country_qa_unmatched_records.csv")

readr::write_csv(summary_tbl, summary_path, na = "")
readr::write_csv(unmatched_titles, unmatched_titles_path, na = "")
readr::write_csv(unmatched_records, unmatched_records_path, na = "")

message_if(
  "Saved WHO DON QA summary to: ",
  summary_path,
  " | unmatched records=",
  nrow(unmatched_records),
  verbose = verbose
)
