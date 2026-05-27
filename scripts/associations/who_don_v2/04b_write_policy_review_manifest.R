source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_policy.R"))

who_don_v2_ensure_dirs()

manifest <- v2_policy_review_manifest()
v2_write_csv(manifest, who_don_v2_output_dir("qa", "v2_policy_review_decision_manifest.csv"))

message("Wrote policy/review decision manifest: ", nrow(manifest), " rows")
