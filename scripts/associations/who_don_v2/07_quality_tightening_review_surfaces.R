library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_extraction.R"))

who_don_v2_ensure_dirs()

completed_review_dir <- function(...) {
  who_don_v2_qa_archive_dir("completed_review_surfaces", ...)
}

country_comparison <- v2_compare_country_candidates()

v2_quality_text <- function(x) {
  str_to_lower(str_squish(coalesce(as.character(x), "")))
}

v2_quality_gap_category <- function(country_standard, evidence_text, accepted_scope) {
  evidence <- v2_quality_text(evidence_text)
  scope <- v2_quality_text(accepted_scope)

  case_when(
    evidence == "" ~ "missing_accepted_evidence_text",
    str_detect(evidence, regex("map of|number of laboratory|laboratory-confirmed cases|alerts have been reported from|\\b[0-9]+\\b.*\\b[0-9]+\\b", ignore_case = TRUE)) ~
      "table_or_case_count_listing",
    str_detect(scope, regex("historical|background|lab|surveillance", ignore_case = TRUE)) ~
      "non_focal_context_legacy_exception",
    country_standard == "Ireland" &
      str_detect(evidence, regex("\\bnorthern ireland\\b", ignore_case = TRUE)) ~
      "non_focal_context_legacy_exception",
    country_standard == "Mongolia" &
      str_detect(evidence, regex("\\binner mongolia\\b", ignore_case = TRUE)) ~
      "non_focal_context_legacy_exception",
    str_detect(evidence, regex("history of .*transmission and outbreaks|previous|historical|endemic|region|regional|neighbouring|neighboring", ignore_case = TRUE)) ~
      "non_focal_context_legacy_exception",
    str_count(evidence, "\\S+") <= 6L &
      (
        str_detect(evidence, fixed(str_to_lower(country_standard))) |
          str_detect(evidence, regex(paste0("\\bin\\s+", stringr::str_replace_all(str_to_lower(country_standard), "([\\W])", "\\\\\\1"), "$"), ignore_case = TRUE))
      ) ~
      "legacy_country_or_title_only_span",
    str_detect(evidence, fixed(str_to_lower(country_standard))) &
      str_detect(evidence, regex("authorities have confirmed|first .*confirmed cases|ministry of health( has)? reported|national ihr focal point|ha?emorrhagic fever in", ignore_case = TRUE)) ~
      "legacy_event_span_not_recovered",
    str_detect(evidence, regex("holland|zaire|libyan arab jamahiriya|burma", ignore_case = TRUE)) ~
      "historical_country_name",
    str_detect(evidence, regex("west bank and gaza strip|gaza strip|hong kong", ignore_case = TRUE)) ~
      "subnational_or_territory_mapped_to_country",
    str_detect(evidence, regex("^\\s*(u\\s*s|viet\\s*nam|republic of korea)\\s*$", ignore_case = TRUE)) ~
      "legacy_country_or_title_only_span",
    str_detect(evidence, regex("canadian authorities", ignore_case = TRUE)) ~
      "legacy_event_span_not_recovered",
    str_detect(evidence, regex("korea, republic of|islamic republic of iran|uk", ignore_case = TRUE)) ~
      "who_or_common_name_variant",
    str_detect(evidence, regex("luxemburg|hezegovina|côte d'ivoir|cote d`ivoire|cap verde", ignore_case = TRUE)) ~
      "spelling_or_typo_variant",
    str_detect(evidence, regex("cote d.?ivoire", ignore_case = TRUE)) ~
      "spelling_or_typo_variant",
    str_detect(evidence, fixed(str_to_lower(country_standard))) ~
      "native_rule_or_standardization_gap",
    str_count(evidence, "\\S+") <= 2L ~
      "short_accepted_span_needs_context",
    TRUE ~ "manual_review_needed"
  )
}

v2_quality_gap_action <- function(gap_category) {
  case_when(
    gap_category %in% c(
      "historical_country_name",
      "subnational_or_territory_mapped_to_country",
      "who_or_common_name_variant",
      "spelling_or_typo_variant",
      "native_rule_or_standardization_gap"
    ) ~ "review_for_country_alias_or_standardization_rule",
    gap_category %in% c(
      "legacy_country_or_title_only_span",
      "legacy_event_span_not_recovered",
      "table_or_case_count_listing",
      "non_focal_context_legacy_exception"
    ) ~
      "keep_legacy_exception_unless_event_support_is_confirmed",
    gap_category == "short_accepted_span_needs_context" ~
      "inspect_record_context_before_rule_change",
    TRUE ~ "manual_review"
  )
}

country_gaps <- country_comparison$unmatched_accepted %>%
  mutate(
    gap_category = v2_quality_gap_category(
      country_standard,
      accepted_country_evidence_text,
      accepted_country_scope
    ),
    recommended_action = v2_quality_gap_action(gap_category),
    review_decision_class = case_when(
      gap_category %in% c(
        "historical_country_name",
        "who_or_common_name_variant",
        "spelling_or_typo_variant"
      ) ~ "needs_alias_rule",
      gap_category %in% c("subnational_or_territory_mapped_to_country") ~
        "needs_standardization_rule",
      gap_category %in% c("native_rule_or_standardization_gap", "short_accepted_span_needs_context") ~
        "needs_context_review",
      gap_category %in% c(
        "legacy_country_or_title_only_span",
        "legacy_event_span_not_recovered",
        "table_or_case_count_listing",
        "non_focal_context_legacy_exception"
      ) ~
        "keep_legacy_exception",
      TRUE ~ "needs_context_review"
    ),
    review_priority = case_when(
      gap_category %in% c(
        "historical_country_name",
        "subnational_or_territory_mapped_to_country",
        "who_or_common_name_variant",
        "spelling_or_typo_variant",
        "native_rule_or_standardization_gap"
      ) ~ "high_rule_review",
      accepted_country_scope == "focal_event_country" ~ "high_focal_legacy_exception_review",
      gap_category == "manual_review_needed" ~ "medium_manual_review",
      TRUE ~ "low_context_or_listing_review"
    )
  ) %>%
  arrange(review_priority, gap_category, country_standard, record_key)

country_gap_summary <- country_gaps %>%
  count(gap_category, recommended_action, review_priority, name = "rows") %>%
  arrange(review_priority, desc(rows), gap_category)

scope_candidates <- v2_read_csv(
  who_don_v2_output_dir("review", "who_don_scope_adjudication_candidates.csv")
)
final_audit <- v2_read_csv(
  who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv")
)

record_disease_lookup <- final_audit %>%
  distinct(record_key, disease_standard) %>%
  filter(!is.na(disease_standard), disease_standard != "") %>%
  arrange(record_key, disease_standard) %>%
  group_by(record_key) %>%
  summarise(
    disease_standard = str_c(unique(disease_standard), collapse = "; "),
    .groups = "drop"
  )

scope_enriched <- scope_candidates %>%
  left_join(
    final_audit %>%
      select(
        record_key,
        country_standard,
        disease_standard,
        claim_type,
        scope_rule_id,
        association_scope,
        scope_confidence,
        scope_evidence_text,
        country_evidence_text,
        disease_evidence_text,
        source_method,
        country_adoption_decision
      ),
    by = c("record_key", "country_standard", "disease_standard")
  ) %>%
  mutate(
    evidence_lower = v2_quality_text(evidence_text),
    scope_review_bucket = case_when(
      str_detect(evidence_lower, regex("confirmed cases|laboratory-confirmed|outbreak|has reported|reported .*cases|deaths?", ignore_case = TRUE)) &
        !str_detect(evidence_lower, regex("travel|import|history|previous|surveillance|sequence|laboratory in|reference laboratory", ignore_case = TRUE)) ~
        "possible_focal_event_language",
      str_detect(evidence_lower, regex("travel|import|returned from|history of travel|exposure", ignore_case = TRUE)) ~
        "import_or_exposure_context_language",
      str_detect(evidence_lower, regex("previous|historical|endemic|globally|worldwide|neighbouring|neighboring", ignore_case = TRUE)) ~
        "background_or_historical_language",
      str_detect(evidence_lower, regex("surveillance|sequence|sequencing|genomic|laboratory|reference laboratory|partner", ignore_case = TRUE)) ~
        "surveillance_lab_or_partner_language",
      TRUE ~ "uncertain_language_needs_sampling"
    ),
    review_priority = case_when(
      scope_review_bucket == "possible_focal_event_language" ~ "high_scope_review",
      country_adoption_decision == "accept_legacy_exception" ~ "medium_legacy_country_scope_review",
      scope_review_bucket == "uncertain_language_needs_sampling" ~ "medium_scope_sampling",
      TRUE ~ "low_context_scope_review"
    ),
    scope_policy_decision_class = case_when(
      scope_review_bucket == "possible_focal_event_language" ~
        "candidate_for_deterministic_event_claim_rule",
      scope_review_bucket == "import_or_exposure_context_language" ~
        "candidate_for_import_or_exposure_rule",
      scope_review_bucket == "background_or_historical_language" ~
        "candidate_for_background_rule",
      scope_review_bucket == "surveillance_lab_or_partner_language" ~
        "candidate_for_surveillance_or_lab_rule",
      TRUE ~ "needs_sampling_or_manual_review"
    )
  ) %>%
  select(-evidence_lower) %>%
  arrange(review_priority, scope_review_bucket, disease_standard, country_standard, record_key)

scope_summary <- bind_rows(
  scope_enriched %>%
    count(claim_type, scope_rule_id, association_scope, review_priority, scope_review_bucket, name = "rows") %>%
    mutate(summary_type = "claim_rule_scope_bucket"),
  scope_enriched %>%
    count(disease_standard, review_priority, scope_review_bucket, name = "rows") %>%
    mutate(summary_type = "disease_bucket"),
  scope_enriched %>%
    count(country_standard, review_priority, scope_review_bucket, name = "rows") %>%
    mutate(summary_type = "country_bucket"),
  scope_enriched %>%
    count(reason_for_review, review_priority, scope_review_bucket, name = "rows") %>%
    mutate(summary_type = "reason_bucket")
) %>%
  arrange(summary_type, review_priority, desc(rows))

scope_sample <- scope_enriched %>%
  group_by(review_priority, scope_review_bucket, disease_standard) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(review_priority, scope_review_bucket, disease_standard, record_key)

native_new <- country_comparison$native_new %>%
  left_join(record_disease_lookup, by = "record_key") %>%
  mutate(disease_standard = coalesce(disease_standard, "record_disease_unresolved")) %>%
  mutate(
    native_new_value_score = case_when(country_claim_type == "local_event" ~ 4L, TRUE ~ 0L) +
      case_when(country_evidence_location == "title" ~ 4L, country_evidence_location %in% c("summary", "overview", "epidemiology") ~ 2L, TRUE ~ 0L) +
      case_when(country_confidence == "high" ~ 3L, country_confidence == "medium" ~ 2L, country_confidence == "review" ~ -2L, TRUE ~ 0L) +
      case_when(country_claim_type %in% c("background_context", "lab_or_partner_context", "surveillance_or_sequence_context") ~ -3L, TRUE ~ 0L) +
      case_when(alias_type %in% c("abbreviation", "countrycode_un.name.en") ~ -1L, TRUE ~ 0L),
    priority_bucket = case_when(
      native_new_value_score >= 10L ~ "high_value_title_or_strong_local_event",
      native_new_value_score >= 6L ~ "medium_value_local_event",
      country_claim_type == "local_event" ~ "low_value_local_event_needs_context",
      TRUE ~ "context_or_uncertain_keep_rejected"
    ),
    initial_review_class = case_when(
      priority_bucket == "high_value_title_or_strong_local_event" ~
        "manual_review_adopt_or_rule_update_candidate",
      priority_bucket == "medium_value_local_event" ~
        "manual_review_sample_before_adoption",
      priority_bucket == "low_value_local_event_needs_context" ~
        "inspect_context_before_any_adoption",
      TRUE ~ "keep_rejected_unreviewed"
    ),
    review_decision_class = case_when(
      priority_bucket == "high_value_title_or_strong_local_event" ~ "adopt",
      priority_bucket == "medium_value_local_event" ~ "needs_further_evidence",
      priority_bucket == "low_value_local_event_needs_context" ~ "needs_context_review",
      TRUE ~ "reject"
    ),
    llm_use_policy = case_when(
      priority_bucket %in% c("high_value_title_or_strong_local_event", "medium_value_local_event") ~
        "eligible_only_after_manual_subset_selection",
      TRUE ~ "not_llm_input"
    )
  ) %>%
  arrange(desc(native_new_value_score), priority_bucket, country_standard, record_key)

native_new_summary <- native_new %>%
  count(priority_bucket, initial_review_class, llm_use_policy, country_claim_type, country_evidence_location, name = "rows") %>%
  arrange(priority_bucket, desc(rows))

adjudication_subset_candidates <- bind_rows(
  scope_enriched %>%
    filter(review_priority %in% c("high_scope_review", "medium_scope_sampling")) %>%
    transmute(
      queue_type = "scope_adjudication_candidate",
      review_id,
      record_key,
      DonId,
      Title,
      article_url,
      country_standard,
      disease_standard,
      priority_bucket = review_priority,
      review_reason = scope_review_bucket,
      evidence_text,
      suggested_decision_surface = "manual_scope_adjudication",
      llm_use_policy = "eligible_only_after_manual_subset_selection"
    ) %>%
    slice_head(n = 100),
  native_new %>%
    filter(priority_bucket %in% c("high_value_title_or_strong_local_event", "medium_value_local_event")) %>%
    transmute(
      queue_type = "native_new_country_candidate",
      review_id = paste(record_key, country_standard, "native_new_country", sep = "::"),
      record_key,
      DonId,
      Title,
      article_url,
      country_standard,
      disease_standard,
      priority_bucket,
      review_reason = initial_review_class,
      evidence_text = local_evidence_text,
      suggested_decision_surface = "manual_country_candidate_adoption",
      llm_use_policy
    ) %>%
    slice_head(n = 100)
) %>%
  arrange(queue_type, priority_bucket, record_key, country_standard)

quality_manifest <- tibble::tibble(
  artifact = c(
    "qa/archive/completed_review_surfaces/v2_country_recovery_gap_review.csv",
    "qa/archive/completed_review_surfaces/v2_country_recovery_gap_summary.csv",
    "qa/archive/completed_review_surfaces/v2_scope_adjudication_candidates_enriched.csv",
    "qa/archive/completed_review_surfaces/v2_scope_adjudication_summary.csv",
    "qa/archive/completed_review_surfaces/v2_scope_adjudication_review_sample.csv",
    "qa/archive/completed_review_surfaces/v2_native_new_country_priority_review.csv",
    "qa/archive/completed_review_surfaces/v2_native_new_country_priority_summary.csv",
    "v2_targeted_adjudication_subset_candidates.csv"
  ),
  rows = c(
    nrow(country_gaps),
    nrow(country_gap_summary),
    nrow(scope_enriched),
    nrow(scope_summary),
    nrow(scope_sample),
    nrow(native_new),
    nrow(native_new_summary),
    nrow(adjudication_subset_candidates)
  ),
  note = c(
    "Categorized accepted record-country rows not recovered natively.",
    "Summary of native-country recovery gap categories.",
    "Scope adjudication candidates joined to final claim/scope metadata.",
    "Scope adjudication summaries by claim/rule, disease, country, and reason.",
    "Deterministic sample from major scope adjudication buckets.",
    "Ranked native-new country candidates for selective review.",
    "Summary of ranked native-new country candidate buckets.",
    "Small candidate surface for future manual/LLM adjudication; not routine input."
  )
)

v2_write_csv(country_gaps, completed_review_dir("v2_country_recovery_gap_review.csv"))
v2_write_csv(country_gap_summary, completed_review_dir("v2_country_recovery_gap_summary.csv"))
v2_write_csv(scope_enriched, completed_review_dir("v2_scope_adjudication_candidates_enriched.csv"))
v2_write_csv(scope_summary, completed_review_dir("v2_scope_adjudication_summary.csv"))
v2_write_csv(scope_sample, completed_review_dir("v2_scope_adjudication_review_sample.csv"))
v2_write_csv(native_new, completed_review_dir("v2_native_new_country_priority_review.csv"))
v2_write_csv(native_new_summary, completed_review_dir("v2_native_new_country_priority_summary.csv"))
v2_write_csv(
  adjudication_subset_candidates,
  who_don_v2_output_dir("review", "v2_targeted_adjudication_subset_candidates.csv")
)
v2_write_csv(quality_manifest, completed_review_dir("v2_quality_tightening_manifest.csv"))

message(
  "Wrote v2 quality tightening review surfaces: ",
  nrow(country_gaps),
  " country gaps; ",
  nrow(scope_enriched),
  " scope candidates; ",
  nrow(native_new),
  " native-new country candidates"
)
