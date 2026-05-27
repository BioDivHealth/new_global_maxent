source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_extraction.R"))

who_don_v2_ensure_dirs()

native <- v2_extract_native_country_candidates()

v2_write_csv(
  native$candidates,
  who_don_v2_output_dir("candidates", "who_don_country_candidates_native.csv")
)
v2_write_csv(
  native$review_queue,
  who_don_v2_output_dir("review", "v2_country_candidate_review_queue.csv")
)
v2_write_stage_diagnostic(
  native$summary,
  "v2_native_country_extraction_summary.csv"
)

message(
  "Wrote native country candidates: ",
  nrow(native$candidates),
  " rows; ",
  nrow(native$review_queue),
  " review rows"
)
