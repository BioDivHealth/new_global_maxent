source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_extraction.R"))

who_don_v2_ensure_dirs()

native <- v2_extract_native_disease_candidates()

v2_write_csv(
  native$candidates,
  who_don_v2_output_dir("candidates", "who_don_disease_candidates_native.csv")
)
v2_write_csv(
  native$review_queue,
  who_don_v2_output_dir("review", "influenza_refinement_review_queue.csv")
)
v2_write_stage_diagnostic(
  native$summary,
  "v2_native_disease_extraction_summary.csv"
)

message(
  "Wrote native disease candidates: ",
  nrow(native$candidates),
  " rows; ",
  nrow(native$review_queue),
  " review rows"
)
