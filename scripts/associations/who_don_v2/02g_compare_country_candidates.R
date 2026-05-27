source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_extraction.R"))

who_don_v2_ensure_dirs()

comparison <- v2_compare_country_candidates()
adoption_decisions <- v2_country_adoption_decisions(comparison$diff)

v2_write_csv(
  comparison$summary,
  who_don_v2_output_dir("qa", "v2_native_country_vs_accepted_summary.csv")
)
v2_write_stage_diagnostic(
  comparison$summary_by_claim,
  "v2_native_country_vs_accepted_summary_by_claim.csv"
)
v2_write_stage_diagnostic(
  comparison$diff,
  "v2_native_country_vs_accepted_diff.csv"
)
v2_write_stage_diagnostic(
  comparison$native_new,
  "v2_native_country_new_candidates.csv"
)
v2_write_stage_diagnostic(
  comparison$unmatched_accepted,
  "v2_native_country_unmatched_accepted.csv"
)
v2_write_csv(
  adoption_decisions,
  who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv")
)

message(
  "Wrote native-vs-accepted country comparison: ",
  nrow(comparison$diff),
  " comparison rows; ",
  nrow(adoption_decisions),
  " adoption decision rows"
)
