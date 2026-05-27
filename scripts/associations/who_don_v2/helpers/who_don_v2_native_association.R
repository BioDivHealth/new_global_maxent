suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_extraction.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_compare.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_claims.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_final_shaping.R"))

v2_native_empty_to_na <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::na_if(x, "")
}

v2_native_norm_key <- function(x) {
  dplyr::coalesce(v2_native_empty_to_na(x), "")
}

v2_native_logical_from_text <- function(x) {
  x <- stringr::str_to_lower(as.character(x))
  dplyr::coalesce(x %in% c("true", "t", "1", "yes"), FALSE)
}

v2_native_scope_from_country_claim <- function(country_claim_type) {
  dplyr::case_when(
    country_claim_type == "imported_case" ~ "imported_case_country",
    country_claim_type == "exposure_origin" ~ "travel_or_import_context_country",
    country_claim_type == "background_context" ~ "historical_or_background_context_country",
    country_claim_type == "lab_or_partner_context" ~ "lab_or_partner_context_country",
    country_claim_type == "surveillance_or_sequence_context" ~ "surveillance_or_sequence_context_country",
    TRUE ~ "uncertain_focality"
  )
}

v2_accepted_native_countries <- function() {
  native_accept_decisions <- c("accept_native", "accept_native_reviewed")

  v2_read_csv(
    who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv"),
    c(
      "record_key", "country_standard", "present_native", "adoption_decision",
      "native_country_raw", "native_country_evidence_text",
      "native_country_evidence_full_text", "native_country_evidence_location",
      "native_country_source_method", "native_country_rule_id",
      "native_country_confidence", "native_country_claim_type",
      "native_country_claim_reason"
    )
  ) %>%
    filter(present_native, adoption_decision %in% native_accept_decisions) %>%
    transmute(
      record_key,
      country_standard,
      country_raw = native_country_raw,
      country_evidence_text = coalesce(
        native_country_evidence_text,
        native_country_evidence_full_text
      ),
      country_evidence_location = native_country_evidence_location,
      country_source_method = native_country_source_method,
      country_rule_id = native_country_rule_id,
      country_confidence = native_country_confidence,
      country_needs_review = FALSE,
      country_claim_type = native_country_claim_type,
      country_claim_reason = native_country_claim_reason,
      country_adoption_decision = adoption_decision,
      country_adoption_decision_id = decision_id,
      country_adoption_note = decision_note
    ) %>%
    arrange(record_key, country_standard, country_adoption_decision) %>%
    distinct(record_key, country_standard, .keep_all = TRUE)
}

v2_accepted_native_diseases <- function() {
  comparison <- v2_compare_disease_candidates()

  native_accept_decisions <- comparison$adoption_decisions %>%
    filter(adoption_decision == "accept_native", present_native) %>%
    mutate(
      disease_key = v2_native_norm_key(disease_standard),
      influenza_type_key = v2_native_norm_key(influenza_type),
      influenza_subtype_key = v2_native_norm_key(influenza_subtype)
    ) %>%
    transmute(
      record_key,
      disease_key,
      influenza_type_key,
      influenza_subtype_key,
      native_disease_adoption_decision = adoption_decision,
      native_disease_adoption_decision_id = decision_id,
      native_disease_adoption_note = decision_note
    ) %>%
    distinct()

  comparison$diff %>%
    filter(present_native) %>%
    mutate(
      disease_key = v2_native_norm_key(disease_standard),
      influenza_type_key = v2_native_norm_key(influenza_type),
      influenza_subtype_key = v2_native_norm_key(influenza_subtype)
    ) %>%
    left_join(
      native_accept_decisions,
      by = c(
        "record_key",
        "disease_key",
        "influenza_type_key",
        "influenza_subtype_key"
      )
    ) %>%
    filter(
      diff_category == "exact_match" |
        !is.na(native_disease_adoption_decision)
    ) %>%
    transmute(
      record_key,
      DonId = coalesce(DonId_native, DonId_seeded),
      record_id = coalesce(record_id_native, record_id_seeded),
      disease_raw = disease_raw_native,
      disease_standard = v2_native_empty_to_na(disease_standard),
      disease_evidence_text = disease_evidence_text_native,
      disease_evidence_location = disease_evidence_location_native,
      disease_source_method = disease_source_method_native,
      disease_rule_id = disease_rule_id_native,
      disease_confidence = disease_confidence_native,
      disease_needs_review = v2_native_logical_from_text(disease_needs_review_native),
      influenza_type = v2_native_empty_to_na(influenza_type),
      influenza_subtype = v2_native_empty_to_na(influenza_subtype),
      native_disease_diff_category = diff_category,
      native_disease_adoption_decision = coalesce(
        native_disease_adoption_decision,
        "native_exact_match"
      ),
      native_disease_adoption_decision_id,
      native_disease_adoption_note = coalesce(
        native_disease_adoption_note,
        "Native disease candidate exactly matches the seeded comparison key."
      )
    ) %>%
    arrange(record_key, disease_standard, influenza_type, influenza_subtype) %>%
    distinct(record_key, disease_standard, influenza_type, influenza_subtype, .keep_all = TRUE)
}

v2_build_native_association_from_layers <- function(records, countries, diseases) {
  countries %>%
    inner_join(
      diseases,
      by = "record_key",
      relationship = "many-to-many"
    ) %>%
    left_join(
      records %>%
        select(
          record_key,
          record_DonId = DonId,
          record_record_id = record_id,
          Title,
          publication_datetime_utc,
          article_url
        ),
      by = "record_key"
    ) %>%
    mutate(
      DonId = coalesce(DonId, record_DonId),
      association_scope = v2_native_scope_from_country_claim(country_claim_type),
      scope_confidence = case_when(
        association_scope == "uncertain_focality" ~ "review",
        TRUE ~ coalesce(country_confidence, disease_confidence, "medium")
      ),
      scope_rule_id = paste0("native_country_claim:", coalesce(country_claim_type, "unknown")),
      scope_reason = country_claim_reason,
      scope_evidence_text = coalesce(
        country_evidence_text,
        disease_evidence_text
      ),
      source_method = case_when(
        country_adoption_decision %in% c("accept_native", "accept_native_reviewed") &
          native_disease_adoption_decision == "accept_native" ~
          "native_only_country_disease_reviewed_adoption",
        country_adoption_decision %in% c("accept_native", "accept_native_reviewed") ~
          "native_only_country_seeded_exact_disease",
        TRUE ~ "native_only_association"
      ),
      needs_review = disease_needs_review | association_scope == "uncertain_focality",
      review_status = if_else(needs_review, "queued", "not_needed"),
      review_decision_id = native_disease_adoption_decision_id,
      review_note = native_disease_adoption_note,
      clean_association_scope = NA_character_,
      clean_don_country_report_scope = NA_character_,
      clean_final_country_role = NA_character_,
      clean_final_event_country_flag = NA,
      clean_strict_focal_event_country_flag = NA,
      clean_disease_label_clean = NA_character_,
      clean_disease_label_standard = NA_character_,
      clean_disease_label_standard_refined = NA_character_
    ) %>%
    transmute(
      record_key,
      DonId,
      record_id = coalesce(record_id, record_record_id),
      Title,
      publication_datetime_utc,
      article_url,
      country_raw,
      country_standard,
      country_evidence_text,
      country_evidence_location,
      country_source_method,
      country_rule_id,
      country_confidence,
      disease_raw,
      disease_standard,
      disease_evidence_text,
      disease_evidence_location,
      disease_source_method,
      disease_rule_id,
      disease_confidence,
      influenza_type,
      influenza_subtype,
      association_scope,
      scope_confidence,
      scope_rule_id,
      scope_reason,
      scope_evidence_text,
      source_method,
      needs_review,
      review_status,
      review_decision_id,
      review_note,
      clean_association_scope,
      clean_don_country_report_scope,
      clean_final_country_role,
      clean_final_event_country_flag,
      clean_strict_focal_event_country_flag,
      clean_disease_label_clean,
      clean_disease_label_standard,
      clean_disease_label_standard_refined,
      country_claim_type,
      country_claim_reason,
      country_adoption_decision,
      country_adoption_decision_id,
      country_adoption_note,
      native_disease_diff_category,
      native_disease_adoption_decision,
      native_disease_adoption_decision_id,
      native_disease_adoption_note
    ) %>%
    distinct()
}

v2_build_native_association_evidence <- function(records = v2_read_records()) {
  native_countries <- v2_accepted_native_countries()
  native_diseases <- v2_accepted_native_diseases()

  v2_build_native_association_from_layers(
    records = records,
    countries = native_countries,
    diseases = native_diseases
  )
}

v2_apply_native_association_scope <- function(evidence) {
  evidence <- evidence %>%
    mutate(evidence_row_id = row_number())

  claims <- v2_build_claims_from_evidence(evidence)
  scoped <- v2_apply_claim_scope(evidence, claims)
  review_queue <- v2_review_queue_from_evidence(scoped)
  adjudication_candidates <- v2_scope_adjudication_candidates(review_queue)

  decisions_path <- who_don_v2_output_dir("review", "who_don_review_decisions_accepted.csv")
  decisions <- if (file.exists(decisions_path)) {
    v2_read_csv(decisions_path)
  } else {
    v2_empty_review_decisions()
  }

  applied <- v2_apply_review_decisions(scoped, decisions)

  list(
    evidence = applied$evidence,
    claims = claims,
    review_queue = review_queue,
    adjudication_candidates = adjudication_candidates,
    review_change_log = applied$change_log,
    review_unmatched = applied$unmatched
  )
}

v2_native_association_audit_shape <- function(evidence) {
  v2_shape_final_scope_audit(evidence)
}

v2_build_native_association_audit <- function(records = v2_read_records()) {
  native_countries <- v2_accepted_native_countries()
  native_diseases <- v2_accepted_native_diseases()

  evidence <- v2_build_native_association_from_layers(
    records = records,
    countries = native_countries,
    diseases = native_diseases
  )
  scoped <- v2_apply_native_association_scope(evidence)
  association_audit <- v2_native_association_audit_shape(scoped$evidence)

  list(
    evidence = evidence,
    scoped_evidence = scoped$evidence,
    association_audit = association_audit,
    claims = scoped$claims,
    review_queue = scoped$review_queue,
    adjudication_candidates = scoped$adjudication_candidates,
    review_change_log = scoped$review_change_log,
    review_unmatched = scoped$review_unmatched,
    accepted_native_countries = native_countries,
    accepted_native_diseases = native_diseases
  )
}
