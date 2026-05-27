source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_compare.R"))

who_don_v2_ensure_dirs()

comparison <- v2_compare_disease_candidates()

v2_write_stage_diagnostic(
  comparison$summary,
  "v2_native_disease_vs_seeded_summary.csv"
)
v2_write_stage_diagnostic(
  comparison$summary_by_disease,
  "v2_native_disease_vs_seeded_summary_by_disease.csv"
)
v2_write_stage_diagnostic(
  comparison$summary_by_anchor,
  "v2_native_disease_vs_seeded_summary_by_anchor.csv"
)
v2_write_stage_diagnostic(
  comparison$diff,
  "v2_native_disease_vs_seeded_diff.csv"
)
v2_write_stage_diagnostic(
  comparison$unmatched_seeded,
  "v2_native_disease_unmatched_seeded.csv"
)
v2_write_stage_diagnostic(
  comparison$new_native,
  "v2_native_disease_new_candidates.csv"
)
v2_write_stage_diagnostic(
  comparison$native_new_candidate_review,
  "v2_native_new_candidate_noise_review.csv"
)
v2_write_stage_diagnostic(
  comparison$influenza_change_review,
  "v2_native_influenza_change_review.csv"
)
v2_write_stage_diagnostic(
  comparison$subtype_changed_classification,
  "v2_influenza_subtype_changed_classification.csv"
)
v2_write_csv(
  comparison$influenza_review_queue,
  who_don_v2_output_dir("review", "influenza_refinement_review_queue.csv")
)
v2_write_csv(
  comparison$adoption_decisions,
  who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv")
)

message(
  "Wrote native-vs-seeded disease comparison: ",
  nrow(comparison$diff),
  " comparison rows; ",
  nrow(comparison$influenza_change_review),
  " influenza review rows; ",
  nrow(comparison$adoption_decisions),
  " adoption decision rows"
)
