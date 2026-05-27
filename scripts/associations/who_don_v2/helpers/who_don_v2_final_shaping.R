suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

v2_final_scope_audit_input_cols <- c(
  "record_key",
  "DonId",
  "record_id",
  "Title",
  "publication_datetime_utc",
  "article_url",
  "country_standard",
  "disease_standard",
  "disease_raw",
  "influenza_type",
  "influenza_subtype",
  "final_association_scope",
  "final_scope_confidence",
  "scope_rule_id",
  "scope_reason",
  "scope_evidence_text",
  "claim_type",
  "claim_id",
  "claim_rule_id",
  "claim_evidence_text",
  "claim_provenance",
  "country_evidence_text",
  "disease_evidence_text",
  "source_method",
  "country_source_method",
  "country_rule_id",
  "country_claim_type",
  "country_claim_reason",
  "country_adoption_decision",
  "country_adoption_decision_id",
  "country_adoption_note",
  "needs_review",
  "review_status",
  "review_decision_id",
  "final_review_source",
  "final_review_note",
  "clean_association_scope",
  "clean_don_country_report_scope",
  "clean_final_country_role",
  "clean_final_event_country_flag",
  "clean_strict_focal_event_country_flag"
)

v2_check_final_shaping_cols <- function(evidence) {
  missing_cols <- setdiff(v2_final_scope_audit_input_cols, names(evidence))
  if (length(missing_cols) > 0) {
    stop(
      "Final shaping evidence missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

v2_shape_final_scope_audit <- function(evidence) {
  v2_check_final_shaping_cols(evidence)

  evidence %>%
    transmute(
      record_key,
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      article_url,
      country_standard,
      disease_standard,
      disease_raw,
      influenza_type,
      influenza_subtype,
      association_scope = final_association_scope,
      scope_confidence = final_scope_confidence,
      scope_rule_id,
      scope_reason,
      scope_evidence_text,
      claim_type,
      claim_id,
      claim_rule_id,
      claim_evidence_text,
      claim_provenance,
      country_evidence_text,
      disease_evidence_text,
      source_method,
      country_source_method,
      country_rule_id,
      country_claim_type,
      country_claim_reason,
      country_adoption_decision,
      country_adoption_decision_id,
      country_adoption_note,
      needs_review,
      review_status,
      review_decision_id,
      final_review_source,
      final_review_note,
      clean_association_scope,
      clean_don_country_report_scope,
      clean_final_country_role,
      clean_final_event_country_flag,
      clean_strict_focal_event_country_flag
    )
}

v2_shape_modelling_ready <- function(audit) {
  required_cols <- c(
    "record_key",
    "DonId",
    "record_id",
    "Title",
    "publication_datetime_utc",
    "article_url",
    "country_standard",
    "disease_standard",
    "association_scope",
    "scope_confidence",
    "scope_rule_id",
    "scope_reason",
    "claim_type",
    "influenza_type",
    "influenza_subtype",
    "source_method",
    "needs_review"
  )
  missing_cols <- setdiff(required_cols, names(audit))
  if (length(missing_cols) > 0) {
    stop(
      "Final audit missing modelling-ready column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  audit %>%
    filter(association_scope == "focal_event_country") %>%
    transmute(
      record_key,
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      article_url,
      country_standard,
      disease_label_standard = disease_standard,
      association_scope,
      scope_confidence,
      scope_rule_id,
      scope_reason,
      claim_type,
      influenza_type,
      influenza_subtype,
      source_method,
      needs_review
    )
}

v2_shape_final_outputs <- function(evidence,
                                   compatibility_exports_enabled = FALSE,
                                   compatibility_final_rows = NA_integer_,
                                   compatibility_modelling_rows = NA_integer_,
                                   compatibility_note = NULL) {
  audit <- v2_shape_final_scope_audit(evidence)
  modelling <- v2_shape_modelling_ready(audit)

  if (is.null(compatibility_note)) {
    compatibility_note <- "skipped; set WHO_DON_V2_WRITE_COMPAT=1 to refresh clean-shaped compatibility exports"
  }

  final_qa <- tibble(
    metric = c(
      "audit_rows",
      "modelling_rows",
      "compatibility_exports_enabled",
      "compatibility_final_rows",
      "compatibility_modelling_rows"
    ),
    value = c(
      nrow(audit),
      nrow(modelling),
      as.integer(isTRUE(compatibility_exports_enabled)),
      compatibility_final_rows,
      compatibility_modelling_rows
    ),
    note = c(
      "canonical v2 final audit rows",
      "canonical v2 focal event-country rows",
      compatibility_note,
      compatibility_note,
      compatibility_note
    )
  )

  list(
    audit = audit,
    modelling = modelling,
    final_qa = final_qa
  )
}
