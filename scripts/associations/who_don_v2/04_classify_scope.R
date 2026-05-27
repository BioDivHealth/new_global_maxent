library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_claims.R"))

who_don_v2_ensure_dirs()

evidence <- v2_read_csv(
  who_don_v2_output_dir("evidence", "who_don_association_evidence.csv"),
  c("record_key", "country_standard", "disease_standard", "association_scope", "scope_rule_id")
)
evidence <- evidence %>% mutate(evidence_row_id = row_number())
claims <- v2_build_claims_from_evidence(evidence)
evidence <- v2_apply_claim_scope(evidence, claims)
v2_write_csv(claims, who_don_v2_output_dir("evidence", "who_don_claims.csv"))

seeded_decisions <- evidence %>%
  filter(grepl("openai|manual|override|review", paste(country_source_method, disease_source_method, review_status), ignore.case = TRUE)) %>%
  transmute(
    review_id = paste(record_key, country_standard, disease_standard, "scope", sep = "::"),
    record_key,
    decision_type = "scope_classification",
    decision_value = association_scope,
    confidence = scope_confidence,
    evidence_span = scope_evidence_text,
    reviewer = "clean_pipeline",
    review_source = "clean_final_seed",
    review_note = paste("Accepted in clean final via", country_source_method, disease_source_method),
    review_date = as.Date(NA)
  ) %>%
  distinct()

review_queue <- v2_review_queue_from_evidence(evidence)
adjudication_candidates <- v2_scope_adjudication_candidates(review_queue)
decisions_path <- who_don_v2_output_dir("review", "who_don_review_decisions_accepted.csv")
decisions <- if (file.exists(decisions_path)) {
  v2_read_csv(decisions_path)
} else {
  v2_empty_review_decisions()
}
if (!file.exists(decisions_path)) {
  v2_write_csv(decisions, decisions_path)
}

applied <- v2_apply_review_decisions(evidence, decisions)

v2_write_csv(review_queue, who_don_v2_output_dir("review", "who_don_review_queue.csv"))
v2_write_csv(
  adjudication_candidates,
  who_don_v2_output_dir("review", "who_don_scope_adjudication_candidates.csv")
)
v2_write_csv(seeded_decisions, who_don_v2_output_dir("review", "who_don_review_decisions_seeded_from_clean.csv"))
v2_write_csv(applied$evidence, who_don_v2_output_dir("review", "who_don_review_decisions_applied.csv"))
v2_write_stage_diagnostic(applied$change_log, "v2_review_decision_change_log.csv")
v2_write_stage_diagnostic(applied$unmatched, "v2_review_decision_unmatched.csv")

rule_summary <- applied$evidence %>%
  count(scope_rule_id, final_association_scope, final_scope_confidence, name = "rows") %>%
  arrange(scope_rule_id, final_association_scope)
v2_write_stage_diagnostic(rule_summary, "v2_rule_hit_summary.csv")

queue_summary <- review_queue %>%
  count(review_surface, llm_use_policy, review_task, current_decision, reason_for_review, name = "rows")
v2_write_csv(queue_summary, who_don_v2_output_dir("qa", "v2_review_queue_summary.csv"))

claim_summary <- claims %>%
  count(claim_type, claim_scope, claim_confidence, name = "rows") %>%
  arrange(claim_type, claim_scope, claim_confidence)
v2_write_stage_diagnostic(claim_summary, "v2_claim_type_summary.csv")

message(
  "Applied v2 scope classification and review sidecar: ",
  nrow(applied$evidence),
  " rows; ",
  nrow(adjudication_candidates),
  " scope adjudication candidates; ",
  nrow(claims),
  " claim rows"
)
