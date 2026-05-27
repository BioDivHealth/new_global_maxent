library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_is_blank <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

v2_check_required_values <- function(x, cols, dataset) {
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0) {
    stop(dataset, " missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tibble::tibble(
    check = paste0(dataset, "_no_blank_", cols),
    failed_rows = vapply(cols, function(col) sum(v2_is_blank(x[[col]])), integer(1)),
    status = if_else(failed_rows == 0L, "pass", "fail"),
    note = paste("Required non-blank column:", cols)
  )
}

v2_production_output_specs <- function() {
  tibble::tibble(
    output_role = c(
      "association_evidence",
      "accepted_association_contract_reference",
      "claim_evidence",
      "records_source_fixture",
      "country_candidates_native",
      "country_candidate_review_queue",
      "country_adoption_decisions",
      "country_native_vs_accepted_summary",
      "disease_adoption_decisions",
      "disease_rule_model",
      "scope_review_surface",
      "scope_adjudication_candidates",
      "scope_review_decisions_applied",
      "policy_review_decision_manifest",
      "final_scope_audit",
      "modelling_ready",
      "web_rows_json",
      "web_meta_json",
      "final_export_summary",
      "scope_review_queue_summary",
      "production_checks"
    ),
    path = c(
      who_don_v2_output_dir("evidence", "who_don_association_evidence.csv"),
      who_don_v2_rules_dir("accepted_association_contract.csv"),
      who_don_v2_output_dir("evidence", "who_don_claims.csv"),
      who_don_v2_output_dir("records", "who_don_records_source.csv"),
      who_don_v2_output_dir("candidates", "who_don_country_candidates_native.csv"),
      who_don_v2_output_dir("review", "v2_country_candidate_review_queue.csv"),
      who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv"),
      who_don_v2_output_dir("qa", "v2_native_country_vs_accepted_summary.csv"),
      who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv"),
      who_don_v2_rules_dir("disease_rule_model.csv"),
      who_don_v2_output_dir("review", "who_don_review_queue.csv"),
      who_don_v2_output_dir("review", "who_don_scope_adjudication_candidates.csv"),
      who_don_v2_output_dir("review", "who_don_review_decisions_applied.csv"),
      who_don_v2_output_dir("qa", "v2_policy_review_decision_manifest.csv"),
      who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv"),
      who_don_v2_output_dir("final", "who_don_modelling_ready.csv"),
      who_don_v2_output_dir("web", "who_don_web.json"),
      who_don_v2_output_dir("web", "who_don_meta.json"),
      who_don_v2_output_dir("qa", "v2_final_export_summary.csv"),
      who_don_v2_output_dir("qa", "v2_review_queue_summary.csv"),
      who_don_v2_output_dir("qa", "v2_production_checks.csv")
    )
  )
}

v2_write_output_manifest <- function(specs, manifest_path = who_don_v2_output_dir("final", "who_don_v2_output_manifest.csv")) {
  existing <- file.exists(specs$path)
  if (!all(existing)) {
    missing <- specs$path[!existing]
    stop("Cannot write manifest; missing production outputs: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  manifest_root <- normalizePath(who_don_v2_output_dir(), mustWork = TRUE)
  output_manifest <- specs %>%
    mutate(
      relative_path = sub(paste0("^", manifest_root, "/?"), "", normalizePath(path)),
      rows = vapply(path, function(one_path) {
        if (grepl("[.]csv$", one_path)) {
          nrow(v2_read_csv(one_path))
        } else {
          NA_integer_
        }
      }, integer(1)),
      bytes = file.info(path)$size,
      md5 = unname(tools::md5sum(path)),
      frozen_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )

  v2_write_csv(output_manifest, manifest_path)
  invisible(output_manifest)
}

v2_validate_production_outputs <- function() {
  audit <- v2_read_csv(
    who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv"),
    c("record_key", "country_standard", "disease_standard", "association_scope")
  )
  modelling <- v2_read_csv(
    who_don_v2_output_dir("final", "who_don_modelling_ready.csv"),
    c("association_scope")
  )
  adoption <- v2_read_csv(
    who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv"),
    c("adoption_decision")
  )
  country_adoption <- v2_read_csv(
    who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv"),
    c("adoption_decision")
  )
  review_queue <- v2_read_csv(
    who_don_v2_output_dir("review", "who_don_review_queue.csv"),
    c("review_surface", "llm_use_policy")
  )
  country_candidates <- v2_read_csv(
    who_don_v2_output_dir("candidates", "who_don_country_candidates_native.csv"),
    c("record_key", "country_standard", "country_claim_type", "local_evidence_text")
  )

  required_value_checks <- v2_check_required_values(
    audit,
    c("record_key", "country_standard", "disease_standard", "association_scope"),
    "final_scope_audit"
  )

  provenance_cols <- c(
    "scope_evidence_text",
    "country_evidence_text",
    "disease_evidence_text",
    "source_method",
    "scope_rule_id",
    "country_adoption_decision",
    "claim_type"
  )
  missing_provenance_cols <- setdiff(provenance_cols, names(audit))
  if (length(missing_provenance_cols) > 0) {
    stop(
      "final_scope_audit missing expected evidence/provenance columns: ",
      paste(missing_provenance_cols, collapse = ", "),
      call. = FALSE
    )
  }

  checks <- bind_rows(
    required_value_checks,
    v2_check_required_values(
      country_candidates,
      c("record_key", "country_standard", "country_claim_type", "local_evidence_text"),
      "native_country_candidates"
    ),
    tibble::tibble(
      check = c(
        "modelling_rows_are_focal_event_country",
        "final_scope_audit_no_blank_scope_evidence_text",
        "final_scope_audit_no_blank_country_evidence_text",
        "final_scope_audit_no_blank_source_method",
        "final_scope_audit_no_blank_scope_rule_id",
        "final_scope_audit_no_blank_country_adoption_decision",
        "final_scope_audit_no_blank_claim_type",
        "native_reviewed_rows_no_blank_disease_evidence_text",
        "disease_adoption_no_needs_manual_review",
        "country_adoption_no_needs_manual_review",
        "review_queue_has_no_routine_llm_policy"
      ),
      failed_rows = c(
        sum(modelling$association_scope != "focal_event_country" | v2_is_blank(modelling$association_scope)),
        sum(v2_is_blank(audit$scope_evidence_text)),
        sum(v2_is_blank(audit$country_evidence_text)),
        sum(v2_is_blank(audit$source_method)),
        sum(v2_is_blank(audit$scope_rule_id)),
        sum(v2_is_blank(audit$country_adoption_decision)),
        sum(v2_is_blank(audit$claim_type)),
        audit %>%
          filter(source_method %in% c(
            "v2_native_reviewed_adoption",
            "v2_native_country_disease_reviewed_adoption"
          )) %>%
          summarise(failed_rows = sum(v2_is_blank(disease_evidence_text)), .groups = "drop") %>%
          pull(failed_rows),
        sum(adoption$adoption_decision == "needs_manual_review", na.rm = TRUE),
        sum(country_adoption$adoption_decision == "needs_manual_review", na.rm = TRUE),
        sum(!review_queue$llm_use_policy %in% c("not_llm_input", "eligible_only_after_manual_subset_selection"))
      ),
      note = c(
        "Modelling output is restricted to focal event-country rows.",
        "Scope evidence text is required for final audit rows.",
        "Country evidence text is required for final audit rows.",
        "Source method is required for final audit rows.",
        "Scope rule ID is required for final audit rows.",
        "Country adoption decision is required for final audit rows.",
        "Claim type is required for final audit rows.",
        "Native-reviewed disease additions require disease evidence text.",
        "Disease adoption decisions must have no unresolved manual-review rows.",
        "Country adoption decisions must have no unresolved manual-review rows.",
        "Review queue may expose adjudication candidates, but not routine LLM submission rows."
      )
    )
  ) %>%
    mutate(status = if_else(failed_rows == 0L, "pass", "fail"))

  v2_write_csv(checks, who_don_v2_output_dir("qa", "v2_production_checks.csv"))

  failed <- checks %>% filter(status == "fail")
  if (nrow(failed) > 0) {
    stop(
      "WHO DON v2 production checks failed: ",
      paste(failed$check, failed$failed_rows, sep = "=", collapse = "; "),
      call. = FALSE
    )
  }

  invisible(checks)
}
