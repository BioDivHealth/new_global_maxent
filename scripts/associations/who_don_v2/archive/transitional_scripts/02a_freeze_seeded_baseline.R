source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_baseline.R"))

baseline <- v2_write_seeded_baseline()

message(
  "Wrote seeded v2 baseline: ",
  nrow(baseline$manifest),
  " files; ",
  nrow(baseline$counts),
  " CSV count rows"
)

