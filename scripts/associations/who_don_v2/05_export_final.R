library(dplyr)
library(jsonlite)

source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_rules.R"))

who_don_v2_ensure_dirs()

evidence <- v2_read_csv(
  who_don_v2_output_dir("review", "who_don_review_decisions_applied.csv"),
  c("record_key", "country_standard", "disease_standard", "final_association_scope")
)

write_compatibility_exports <- identical(Sys.getenv("WHO_DON_V2_WRITE_COMPAT"), "1")

audit <- evidence %>%
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

modelling <- audit %>%
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

v2_write_csv(audit, who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv"))
v2_write_csv(modelling, who_don_v2_output_dir("final", "who_don_modelling_ready.csv"))

compatibility_final_rows <- NA_integer_
compatibility_modelling_rows <- NA_integer_
compatibility_note <- "skipped; set WHO_DON_V2_WRITE_COMPAT=1 to refresh clean-shaped compatibility exports"

if (write_compatibility_exports) {
  clean_final <- v2_read_clean_final()
  clean_modelling <- v2_read_csv(v2_clean_modelling_path())

  # Compatibility exports intentionally keep the current clean schemas and values.
  v2_write_csv(
    clean_final,
    who_don_v2_output_dir("final", "who_don_country_disease_event_focal_scope_evidence_final.csv")
  )
  v2_write_csv(
    clean_modelling,
    who_don_v2_output_dir("final", "who_don_country_disease_event_focal_modelling_ready_final.csv")
  )

  compatibility_final_rows <- nrow(clean_final)
  compatibility_modelling_rows <- nrow(clean_modelling)
  compatibility_note <- "written from v2-local clean reference seeds"
}

final_qa <- tibble::tibble(
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
    as.integer(write_compatibility_exports),
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
v2_write_csv(final_qa, who_don_v2_output_dir("qa", "v2_final_export_summary.csv"))

message(
  "Wrote v2 final exports: ",
  nrow(audit),
  " audit rows, ",
  nrow(modelling),
  " modelling rows; compatibility exports ",
  if (write_compatibility_exports) "written" else "skipped"
)
