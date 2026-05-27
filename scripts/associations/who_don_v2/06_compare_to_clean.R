library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))

who_don_v2_ensure_dirs()

clean_audit_dir <- function(...) {
  who_don_v2_qa_archive_dir("clean_audit", ...)
}

clean_final <- v2_read_clean_final()
v2_audit <- v2_read_csv(
  who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv"),
  c("record_key", "country_standard", "disease_standard", "association_scope")
)

clean_keyed <- clean_final %>%
  transmute(
    record_key,
    country_standard = coalesce(country_standard_final, country_standard),
    disease_standard = coalesce(disease_label_standard_refined, disease_label_standard, disease_label_clean),
    clean_scope = v2_scope_from_clean(pick(everything()))
  ) %>%
  distinct()

v2_keyed <- v2_audit %>%
  transmute(
    record_key,
    country_standard,
    disease_standard,
    v2_scope = association_scope,
    v2_scope_rule_id = scope_rule_id,
    v2_source_method = source_method,
    v2_country_adoption_decision = country_adoption_decision,
    v2_country_adoption_note = country_adoption_note
  ) %>%
  distinct()

adoption_path <- who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv")
adoption_explanation <- if (file.exists(adoption_path)) {
  v2_read_csv(
    adoption_path,
    c("record_key", "disease_standard", "diff_category", "adoption_decision", "adoption_source", "decision_note")
  ) %>%
    group_by(record_key, disease_standard) %>%
    summarise(
      disease_adoption_diff_category = paste(sort(unique(diff_category)), collapse = " | "),
      disease_adoption_decision = paste(sort(unique(adoption_decision)), collapse = " | "),
      disease_adoption_source = paste(sort(unique(adoption_source)), collapse = " | "),
      disease_adoption_note = paste(sort(unique(decision_note)), collapse = " | "),
      .groups = "drop"
    )
} else {
  tibble::tibble(
    record_key = character(),
    disease_standard = character(),
    disease_adoption_diff_category = character(),
    disease_adoption_decision = character(),
    disease_adoption_source = character(),
    disease_adoption_note = character()
  )
}

row_diff <- full_join(
  clean_keyed,
  v2_keyed,
  by = c("record_key", "country_standard", "disease_standard")
) %>%
  left_join(adoption_explanation, by = c("record_key", "disease_standard")) %>%
  mutate(
    difference_category = case_when(
      !is.na(clean_scope) & !is.na(v2_scope) & clean_scope == v2_scope ~ "exact_match",
      !is.na(clean_scope) & !is.na(v2_scope) & clean_scope != v2_scope &
        grepl("^claim_type:", coalesce(v2_scope_rule_id, "")) &
        v2_scope != "uncertain_focality" ~ "scope_changed_by_claim_policy",
      !is.na(clean_scope) & !is.na(v2_scope) & clean_scope != v2_scope ~ "scope_changed_needs_review",
      !is.na(clean_scope) & is.na(v2_scope) &
        disease_adoption_decision == "reject_seeded_weak" ~ "clean_removed_by_policy",
      !is.na(clean_scope) & is.na(v2_scope) ~ "clean_removed_unexplained",
      is.na(clean_scope) & !is.na(v2_scope) &
        disease_adoption_decision == "accept_native" ~ "v2_added_by_policy",
      is.na(clean_scope) & !is.na(v2_scope) &
        v2_country_adoption_decision == "accept_native_reviewed" ~ "v2_added_by_country_policy",
      is.na(clean_scope) & !is.na(v2_scope) ~ "v2_added_unexplained",
      TRUE ~ "needs_review"
    )
  )

summary <- row_diff %>%
  count(difference_category, name = "rows") %>%
  arrange(difference_category)

scope_change_review <- row_diff %>%
  filter(!is.na(clean_scope), !is.na(v2_scope), clean_scope != v2_scope)

missing_clean_rows <- row_diff %>% filter(!is.na(clean_scope), is.na(v2_scope))
new_rows <- row_diff %>% filter(is.na(clean_scope), !is.na(v2_scope))

set.seed(20260501)
sample_n_safe <- function(x, n) {
  if (nrow(x) <= n) x else dplyr::slice_sample(x, n = n)
}
review_sample <- bind_rows(
  row_diff %>% filter(difference_category == "exact_match") %>% sample_n_safe(50),
  missing_clean_rows %>% sample_n_safe(50),
  new_rows %>% sample_n_safe(50),
  scope_change_review %>% sample_n_safe(50)
)

v2_write_csv(row_diff, clean_audit_dir("v2_vs_clean_row_diff.csv"))
v2_write_csv(summary, clean_audit_dir("v2_vs_clean_summary.csv"))
v2_write_csv(scope_change_review, clean_audit_dir("v2_scope_change_review.csv"))
v2_write_csv(missing_clean_rows, clean_audit_dir("v2_missing_clean_rows.csv"))
v2_write_csv(new_rows, clean_audit_dir("v2_new_rows.csv"))
v2_write_csv(review_sample, clean_audit_dir("v2_deterministic_review_sample.csv"))

if (file.exists(adoption_path)) {
  adoption_decisions <- v2_read_csv(
    adoption_path,
    c("diff_category", "adoption_decision", "adoption_source", "review_priority")
  )

  adoption_summary <- adoption_decisions %>%
    count(diff_category, adoption_decision, adoption_source, review_priority, name = "rows") %>%
    arrange(adoption_decision, desc(rows), diff_category)

  manual_review_rows <- adoption_decisions %>%
    filter(adoption_decision == "needs_manual_review")

  adoption_gate <- tibble::tibble(
    metric = c(
      "clean_vs_v2_exact_rows",
      "clean_vs_v2_non_exact_rows",
      "native_adoption_decision_rows",
      "native_accept_rows",
      "native_reject_context_rows",
      "seeded_reject_weak_rows",
      "native_keep_seeded_rows",
      "native_manual_review_rows",
      "native_adoption_ready"
    ),
    value = c(
      sum(row_diff$difference_category == "exact_match"),
      sum(row_diff$difference_category != "exact_match"),
      nrow(adoption_decisions),
      sum(adoption_decisions$adoption_decision == "accept_native"),
      sum(adoption_decisions$adoption_decision == "reject_native_context"),
      sum(adoption_decisions$adoption_decision == "reject_seeded_weak"),
      sum(adoption_decisions$adoption_decision == "keep_seeded"),
      nrow(manual_review_rows),
      as.integer(nrow(manual_review_rows) == 0)
    ),
    note = c(
      "Current v2 final exports are built from reviewed disease-adoption evidence.",
      "Non-exact rows should be explained by disease adoption policy after native adoption.",
      "Rows in native-vs-seeded comparison requiring a non-exact decision.",
      "Native disease candidates accepted by deterministic policy.",
      "Native disease candidates rejected as context/non-event.",
      "Seeded disease candidates rejected because native evidence supports a more specific or cleaner event label.",
      "Seeded disease candidates kept by policy.",
      "Rows still blocking full native disease adoption.",
      "1 means native disease candidates can replace seeded candidates without unresolved review rows."
    )
  )

  manual_review_summary <- manual_review_rows %>%
    count(diff_category, review_priority, name = "rows") %>%
    arrange(desc(rows), diff_category)

  v2_write_csv(adoption_summary, clean_audit_dir("v2_native_adoption_decision_summary.csv"))
  v2_write_csv(adoption_gate, clean_audit_dir("v2_native_adoption_gate.csv"))
  v2_write_csv(manual_review_summary, clean_audit_dir("v2_native_adoption_manual_review_summary.csv"))
}

strict <- identical(Sys.getenv("WHO_DON_V2_STRICT_COMPARE"), "1")
if (strict && any(summary$difference_category != "exact_match")) {
  stop("Strict v2 comparison found non-exact rows. See qa/archive/clean_audit/v2_vs_clean_row_diff.csv", call. = FALSE)
}

message("Wrote v2 comparison: ", paste(summary$difference_category, summary$rows, collapse = "; "))
