library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_count_rows_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  nrow(v2_read_csv(path))
}

v2_policy_review_manifest <- function() {
  tibble::tribble(
    ~decision_layer, ~decision_kind, ~durable_input_path, ~generated_output_path, ~notes,
    "country", "candidate_adoption", who_don_v2_scripts_dir("rules", "country_candidate_policy_decisions.csv"), who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv"), "Country adoption decisions are generated from native-vs-accepted diff categories and durable policy rows.",
    "disease", "candidate_adoption", who_don_v2_scripts_dir("rules", "disease_candidate_policy_decisions.csv"), who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv"), "Disease adoption decisions are generated from native-vs-seeded diffs plus durable disease policy rows.",
    "influenza", "subtype_policy", who_don_v2_scripts_dir("rules", "influenza_subtype_policy_decisions.csv"), who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv"), "Influenza subtype policy decisions feed the generated disease adoption table.",
    "scope", "manual_adjudication", who_don_v2_output_dir("review", "who_don_review_decisions_accepted.csv"), who_don_v2_output_dir("review", "who_don_review_decisions_applied.csv"), "Accepted scope sidecar decisions are applied to generated evidence without overwriting the durable input file."
  ) %>%
    mutate(
      durable_input_exists = file.exists(durable_input_path),
      generated_output_exists = file.exists(generated_output_path),
      durable_input_rows = vapply(durable_input_path, v2_count_rows_if_exists, integer(1)),
      generated_output_rows = vapply(generated_output_path, v2_count_rows_if_exists, integer(1))
    )
}
