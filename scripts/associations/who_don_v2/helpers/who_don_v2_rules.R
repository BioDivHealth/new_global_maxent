library(dplyr)
library(stringr)
library(tidyr)
library(readr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_write_seed_csv_if_missing <- function(x, path) {
  if (!file.exists(path)) {
    v2_write_csv(x, path)
  }
  invisible(path)
}

v2_seed_rule_tables <- function() {
  who_don_v2_ensure_dirs()

  country_aliases <- tibble::tibble(
    alias = character(),
    country_standard = character(),
    alias_type = character(),
    is_ambiguous = logical(),
    priority = integer(),
    notes = character()
  )

  disease_aliases <- tibble::tibble(
    alias = character(),
    disease_standard = character(),
    disease_group = character(),
    alias_type = character(),
    priority = integer(),
    is_generic_label = logical(),
    notes = character()
  )

  scope_rules <- tibble::tribble(
    ~scope_rule_id, ~scope, ~pattern, ~pattern_location, ~country_must_be_mentioned, ~disease_must_be_mentioned, ~priority, ~confidence, ~reason_label, ~notes,
    "manual_review_scope", "manual_review_applied", "", "review", FALSE, FALSE, 1L, "high", "accepted review decision", "Reserved for sidecar review decisions.",
    "strict_focal_flag", "focal_event_country", "", "clean_final_field", FALSE, FALSE, 10L, "high", "strict focal flag from clean baseline", "Compatibility seed from accepted clean final field.",
    "probable_focal_flag", "focal_event_country", "", "clean_final_field", FALSE, FALSE, 20L, "medium", "probable focal flag from clean baseline", "Compatibility seed from accepted clean final field.",
    "imported_scope", "imported_case_country", "imported|travel", "clean_scope_or_reason", FALSE, FALSE, 30L, "medium", "import or travel context", "Initial table-driven pattern for later native v2 classification.",
    "context_scope", "historical_or_background_context_country", "background|context|neighbour|neighbor|historic|surveillance", "clean_scope_or_reason", FALSE, FALSE, 40L, "medium", "background or context country", "Initial table-driven pattern for later native v2 classification.",
    "not_event_country", "not_final_event_country", "not.event|non.event|not final", "clean_scope_or_reason", FALSE, FALSE, 50L, "medium", "not a final event country", "Initial table-driven pattern for later native v2 classification.",
    "uncertain_scope", "uncertain_focality", "", "fallback", FALSE, FALSE, 999L, "review", "uncertain scope", "Fallback when no stronger rule is available."
  )

  influenza_rules <- tibble::tribble(
    ~rule_id, ~pattern, ~influenza_type, ~influenza_subtype, ~specificity, ~priority, ~notes,
    "flu_h5n1", "\\bA?\\(?H5N1\\)?\\b", "A", "H5N1", "subtype", 1L, "Subtype evidence in label or evidence span.",
    "flu_h5n2", "\\bA?\\(?H5N2\\)?\\b", "A", "H5N2", "subtype", 2L, "Subtype evidence in label or evidence span.",
    "flu_h5n5", "\\bA?\\(?H5N5\\)?\\b", "A", "H5N5", "subtype", 3L, "Subtype evidence in label or evidence span.",
    "flu_h5n6", "\\bA?\\(?H5N6\\)?\\b", "A", "H5N6", "subtype", 4L, "Subtype evidence in label or evidence span.",
    "flu_h5n8", "\\bA?\\(?H5N8\\)?\\b", "A", "H5N8", "subtype", 5L, "Subtype evidence in label or evidence span.",
    "flu_h7n2", "\\bA?\\(?H7N2\\)?\\b", "A", "H7N2", "subtype", 6L, "Subtype evidence in label or evidence span.",
    "flu_h7n4", "\\bA?\\(?H7N4\\)?\\b", "A", "H7N4", "subtype", 7L, "Subtype evidence in label or evidence span.",
    "flu_h7n7", "\\bA?\\(?H7N7\\)?\\b", "A", "H7N7", "subtype", 8L, "Subtype evidence in label or evidence span.",
    "flu_h7n9", "\\bA?\\(?H7N9\\)?\\b", "A", "H7N9", "subtype", 9L, "Subtype evidence in label or evidence span.",
    "flu_h9n2", "\\bA?\\(?H9N2\\)?\\b", "A", "H9N2", "subtype", 10L, "Subtype evidence in label or evidence span.",
    "flu_h1n1", "\\bA?\\(?H1N1\\)?\\b", "influenza", "H1N1", "subtype", 11L, "Subtype evidence in label or evidence span.",
    "flu_h1n2", "\\bA?\\(?H1N2\\)?\\b", "influenza", "H1N2", "subtype", 12L, "Subtype evidence in label or evidence span.",
    "flu_h3n2", "\\bA?\\(?H3N2\\)?\\b", "influenza", "H3N2", "subtype", 13L, "Subtype evidence in label or evidence span.",
    "flu_h3n8", "\\bA?\\(?H3N8\\)?\\b", "A", "H3N8", "subtype", 14L, "Subtype evidence in label or evidence span.",
    "flu_h10n3", "\\bA?\\(?H10N3\\)?\\b", "A", "H10N3", "subtype", 15L, "Subtype evidence in label or evidence span.",
    "flu_h10n5", "\\bA?\\(?H10N5\\)?\\b", "A", "H10N5", "subtype", 16L, "Subtype evidence in label or evidence span.",
    "flu_h1", "\\bH1\\b", "A", "H1", "subtype_family", 48L, "H1 family evidence without explicit N subtype.",
    "flu_h5", "\\bH5\\b", "A", "H5", "subtype_family", 50L, "H5 family evidence without explicit N subtype.",
    "flu_generic", "influenza", "influenza", NA_character_, "generic", 100L, "Generic influenza mention."
  )

  source_ranks <- tibble::tribble(
    ~source_method, ~rank, ~confidence, ~notes,
    "accepted_review", 1L, "high", "Manual or LLM decision accepted as sidecar input.",
    "clean_final", 2L, "high", "Accepted clean final evidence used to seed v2.",
    "title_rule", 3L, "high", "Reserved for native v2 country/disease extraction.",
    "body_rule", 4L, "medium", "Reserved for native v2 country/disease extraction.",
    "record_level", 5L, "low", "Reserved for weak record-level association evidence."
  )

  influenza_label_standardization <- tibble::tribble(
    ~influenza_subtype, ~canonical_disease_standard, ~canonical_influenza_type, ~canonical_influenza_subtype, ~notes,
    "H1N1", "Influenza (H1N1)", "influenza", "H1N1", "Keep established seasonal/pandemic H1N1 label convention.",
    "H1N2", "Influenza (H1N2)", "influenza", "H1N2", "Keep established H1N2 label convention.",
    "H1", "Influenza A(H1)", "A", "H1", "Keep H1 family label when N subtype is not available.",
    "H3N2", "Influenza (H3N2)", "influenza", "H3N2", "Keep established H3N2 label convention.",
    "H5", "Influenza (H5 subtype)", "A", "H5", "Keep H5 family label when N subtype is not available.",
    "H5N1", "Influenza A(H5N1)", "A", "H5N1", "Use Influenza A(...) for avian/specific A subtypes.",
    "H5N2", "Influenza A(H5N2)", "A", "H5N2", "Use Influenza A(...) for avian/specific A subtypes.",
    "H5N5", "Influenza A(H5N5)", "A", "H5N5", "Use Influenza A(...) for avian/specific A subtypes.",
    "H5N6", "Influenza A(H5N6)", "A", "H5N6", "Use Influenza A(...) for avian/specific A subtypes.",
    "H5N8", "Influenza A(H5N8)", "A", "H5N8", "Use Influenza A(...) for avian/specific A subtypes.",
    "H7", "Influenza (H7 subtype)", "A", "H7", "Keep H7 family label when N subtype is not available.",
    "H7N2", "Influenza A(H7N2)", "A", "H7N2", "Use Influenza A(...) for avian/specific A subtypes.",
    "H7N4", "Influenza A(H7N4)", "A", "H7N4", "Use Influenza A(...) for avian/specific A subtypes.",
    "H7N7", "Influenza A(H7N7)", "A", "H7N7", "Use Influenza A(...) for avian/specific A subtypes.",
    "H7N9", "Influenza A(H7N9)", "A", "H7N9", "Use Influenza A(...) for avian/specific A subtypes.",
    "H9N2", "Influenza A(H9N2)", "A", "H9N2", "Use Influenza A(...) for avian/specific A subtypes.",
    "H10N3", "Influenza A(H10N3)", "A", "H10N3", "Use Influenza A(...) for avian/specific A subtypes.",
    "H10N5", "Influenza A(H10N5)", "A", "H10N5", "Use Influenza A(...) for avian/specific A subtypes."
  )

  v2_write_seed_csv_if_missing(country_aliases, who_don_v2_rules_dir("country_aliases.csv"))
  v2_write_seed_csv_if_missing(disease_aliases, who_don_v2_rules_dir("disease_aliases.csv"))
  v2_write_seed_csv_if_missing(scope_rules, who_don_v2_rules_dir("scope_rules.csv"))
  v2_write_seed_csv_if_missing(influenza_rules, who_don_v2_rules_dir("influenza_subtype_rules.csv"))
  v2_write_seed_csv_if_missing(
    influenza_label_standardization,
    who_don_v2_rules_dir("influenza_label_standardization.csv")
  )
  v2_write_seed_csv_if_missing(source_ranks, who_don_v2_rules_dir("evidence_source_ranks.csv"))
}

v2_fallback_record_key_from_url <- function(article_url) {
  slug <- stringr::str_replace(
    as.character(article_url),
    "^https?://[^/]+/emergencies/disease-outbreak-news/item/",
    ""
  )
  slug <- stringr::str_replace_all(slug, "[^A-Za-z0-9]+", "-")
  slug <- stringr::str_replace_all(slug, "^-+|-+$", "")
  slug <- stringr::str_to_lower(slug)
  dplyr::if_else(
    !is.na(slug) & slug != "",
    paste0("url-", slug),
    NA_character_
  )
}

v2_materialize_reference_inputs <- function() {
  who_don_v2_ensure_dirs()

  copy_if_missing <- function(source_path, reference_path, required_cols = character()) {
    if (file.exists(reference_path)) {
      return(invisible(FALSE))
    }
    seed <- v2_read_csv(source_path, required_cols)
    v2_write_csv(seed, reference_path)
    invisible(TRUE)
  }

  copy_if_missing(
    who_don_clean_output_dir("final", "who_don_country_disease_event_focal_scope_evidence_final.csv"),
    who_don_v2_reference_dir("who_don_clean_final_seed.csv"),
    v2_required_final_cols
  )
  copy_if_missing(
    who_don_clean_output_dir("final", "who_don_country_disease_event_focal_modelling_ready_final.csv"),
    who_don_v2_reference_dir("who_don_clean_modelling_seed.csv")
  )
  copy_if_missing(
    who_don_clean_output_dir("records", "who_don_records_clean.csv"),
    who_don_v2_reference_dir("who_don_clean_records_seed.csv"),
    c("DonId", "record_id", "Title", "publication_datetime_utc", "article_url")
  )
}

v2_clean_record_key_map <- function(association_contract = v2_read_association_contract()) {
  association_contract %>%
    filter(
      !is.na(article_url),
      article_url != "",
      !is.na(record_key),
      record_key != ""
    ) %>%
    distinct(article_url, clean_record_key = record_key)
}

v2_materialize_records <- function(
  source_records,
  association_contract = v2_read_association_contract()
) {
  key_map <- v2_clean_record_key_map(association_contract)

  records <- source_records %>%
    mutate(
      source_record_id = record_id,
      fallback_record_key = v2_fallback_record_key_from_url(article_url)
    ) %>%
    left_join(key_map, by = "article_url") %>%
    mutate(
      record_key = coalesce(clean_record_key, fallback_record_key, record_id),
      record_key_source = case_when(
        !is.na(clean_record_key) & clean_record_key != "" ~ "clean_final_article_url",
        !is.na(fallback_record_key) & fallback_record_key != "" ~ "article_url_fallback",
        TRUE ~ "record_id_fallback"
      )
    ) %>%
    select(
      record_key,
      record_key_source,
      source_record_id,
      everything(),
      -clean_record_key,
      -fallback_record_key
    )

  duplicated_keys <- records %>%
    count(record_key, name = "record_rows") %>%
    filter(record_rows > 1)

  if (nrow(duplicated_keys) > 0) {
    stop(
      "v2 record_key repair still produced duplicated record keys: ",
      paste(head(duplicated_keys$record_key, 20), collapse = ", "),
      call. = FALSE
    )
  }

  records
}

v2_scope_from_clean <- function(x) {
  scope_text <- str_to_lower(paste(
    x$association_scope,
    x$don_country_report_scope,
    x$final_country_role,
    x$final_reasoning_label,
    sep = " | "
  ))

  dplyr::case_when(
    !is.na(x$strict_focal_event_country_flag) & x$strict_focal_event_country_flag ~ "focal_event_country",
    "probable_focal_event_country_flag" %in% names(x) &
      !is.na(x$probable_focal_event_country_flag) &
      x$probable_focal_event_country_flag ~ "focal_event_country",
    str_detect(scope_text, "import|travel") ~ "imported_case_country",
    str_detect(scope_text, "lab|partner") ~ "lab_or_partner_context_country",
    str_detect(scope_text, "surveillance|sequence") ~ "surveillance_or_sequence_context_country",
    str_detect(scope_text, "background|context|neighbo|historic") ~ "historical_or_background_context_country",
    str_detect(scope_text, "not.event|non.event|not final") ~ "not_final_event_country",
    TRUE ~ "uncertain_focality"
  )
}

v2_scope_rule_id <- function(scope, x) {
  dplyr::case_when(
    scope == "focal_event_country" &
      !is.na(x$strict_focal_event_country_flag) &
      x$strict_focal_event_country_flag ~ "strict_focal_flag",
    scope == "focal_event_country" ~ "probable_focal_flag",
    scope == "imported_case_country" ~ "imported_scope",
    scope == "lab_or_partner_context_country" ~ "context_scope",
    scope == "surveillance_or_sequence_context_country" ~ "context_scope",
    scope == "historical_or_background_context_country" ~ "context_scope",
    scope == "not_final_event_country" ~ "not_event_country",
    TRUE ~ "uncertain_scope"
  )
}

v2_scope_confidence <- function(scope, x) {
  dplyr::case_when(
    scope == "focal_event_country" &
      !is.na(x$strict_focal_event_country_flag) &
      x$strict_focal_event_country_flag ~ "high",
    scope == "uncertain_focality" ~ "review",
    TRUE ~ dplyr::coalesce(x$focal_scope_confidence, x$final_event_confidence, "medium")
  )
}

v2_country_candidates_from_clean <- function(clean_final) {
  clean_final %>%
    transmute(
      record_key,
      DonId,
      record_id,
      country_raw = country_standard,
      country_standard = coalesce(country_standard_final, country_standard),
      country_evidence_text = v2_first_present(
        pick(everything()),
        c("final_evidence_span", "best_evidence_span", "evidence_sentence")
      ),
      country_evidence_location = coalesce(evidence_section, "clean_final"),
      country_source_method = coalesce(final_source, country_source_stage, "clean_final"),
      country_rule_id = coalesce(final_reasoning_label, best_reasoning_label, "clean_final_country"),
      country_confidence = coalesce(final_event_confidence, event_confidence, "medium"),
      country_needs_review = coalesce(final_needs_manual_review, needs_manual_review, FALSE)
    ) %>%
    distinct()
}

v2_disease_candidates_from_clean <- function(clean_final) {
  clean_final %>%
    transmute(
      record_key,
      DonId,
      record_id,
      disease_raw = disease_label_raw,
      disease_standard = coalesce(
        disease_label_standard_refined,
        disease_label_standard,
        disease_label_clean
      ),
      disease_evidence_text = coalesce(
        influenza_subtype_evidence_span,
        final_evidence_span,
        evidence_sentence
      ),
      disease_evidence_location = coalesce(evidence_section, "clean_final"),
      disease_source_method = coalesce(disease_refinement_source, disease_match_basis, "clean_final"),
      disease_rule_id = coalesce(disease_match_basis, "clean_final_disease"),
      disease_confidence = coalesce(disease_match_confidence, "medium"),
      disease_needs_review = coalesce(disease_refinement_needs_review, broad_needs_manual_review, FALSE),
      influenza_type,
      influenza_subtype
    ) %>%
    distinct()
}

v2_association_evidence_from_clean <- function(clean_final) {
  scope <- v2_scope_from_clean(clean_final)

  clean_final %>%
    mutate(
      v2_scope = scope,
      scope_rule_id = v2_scope_rule_id(scope, pick(everything())),
      scope_confidence = v2_scope_confidence(scope, pick(everything()))
    ) %>%
    transmute(
      record_key,
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      article_url,
      country_raw = country_standard,
      country_standard = coalesce(country_standard_final, country_standard),
      country_evidence_text = coalesce(final_evidence_span, best_evidence_span, evidence_sentence),
      country_evidence_location = coalesce(evidence_section, "clean_final"),
      country_source_method = coalesce(final_source, country_source_stage, "clean_final"),
      country_rule_id = coalesce(final_reasoning_label, best_reasoning_label, "clean_final_country"),
      country_confidence = coalesce(final_event_confidence, event_confidence, "medium"),
      disease_raw = disease_label_raw,
      disease_standard = coalesce(disease_label_standard_refined, disease_label_standard, disease_label_clean),
      disease_evidence_text = coalesce(influenza_subtype_evidence_span, final_evidence_span, evidence_sentence),
      disease_evidence_location = coalesce(evidence_section, "clean_final"),
      disease_source_method = coalesce(disease_refinement_source, disease_match_basis, "clean_final"),
      disease_rule_id = coalesce(disease_match_basis, "clean_final_disease"),
      disease_confidence = coalesce(disease_match_confidence, "medium"),
      influenza_type,
      influenza_subtype,
      association_scope = v2_scope,
      scope_confidence,
      scope_rule_id,
      scope_reason = coalesce(focal_scope_reasoning_label, final_reasoning_label, best_reasoning_label),
      scope_evidence_text = coalesce(final_evidence_span, best_evidence_span, evidence_sentence),
      source_method = "clean_final_seed",
      needs_review = coalesce(final_needs_manual_review, focal_scope_needs_review, disease_refinement_needs_review, FALSE),
      review_status = if_else(needs_review, "queued", "not_needed"),
      review_decision_id = NA_character_,
      review_note = NA_character_,
      clean_association_scope = association_scope,
      clean_don_country_report_scope = don_country_report_scope,
      clean_final_country_role = final_country_role,
      clean_final_event_country_flag = final_event_country_flag,
      clean_strict_focal_event_country_flag = strict_focal_event_country_flag,
      clean_disease_label_clean = disease_label_clean,
      clean_disease_label_standard = disease_label_standard,
      clean_disease_label_standard_refined = disease_label_standard_refined
    )
}

v2_adoption_key <- function(x) {
  x %>%
    mutate(
      influenza_type_key = coalesce(influenza_type, ""),
      influenza_subtype_key = coalesce(influenza_subtype, "")
    )
}

v2_association_evidence_from_reviewed_candidates <- function(clean_final, adoption_decisions) {
  base <- v2_association_evidence_from_clean(clean_final)

  required <- c(
    "decision_id", "record_key", "disease_standard", "influenza_type",
    "influenza_subtype", "present_native", "adoption_decision",
    "diff_category", "decision_note", "native_evidence_text",
    "native_source_method"
  )
  missing_cols <- setdiff(required, names(adoption_decisions))
  if (length(missing_cols) > 0) {
    stop(
      "Disease adoption decisions missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  rejected_seeded <- adoption_decisions %>%
    filter(adoption_decision == "reject_seeded_weak") %>%
    v2_adoption_key() %>%
    distinct(record_key, disease_standard, influenza_type_key, influenza_subtype_key)

  base_reviewed <- base %>%
    v2_adoption_key() %>%
    anti_join(
      rejected_seeded,
      by = c("record_key", "disease_standard", "influenza_type_key", "influenza_subtype_key")
    ) %>%
    select(-influenza_type_key, -influenza_subtype_key) %>%
    mutate(
      source_method = if_else(
        paste(record_key, disease_standard, coalesce(influenza_type, ""), coalesce(influenza_subtype, ""), sep = "::") %in%
          paste(
            adoption_decisions$record_key,
            adoption_decisions$disease_standard,
            coalesce(adoption_decisions$influenza_type, ""),
            coalesce(adoption_decisions$influenza_subtype, ""),
            sep = "::"
          ),
        "v2_seeded_reviewed_adoption",
        source_method
      )
    )

  country_rows <- base %>%
    group_by(record_key, country_standard) %>%
    slice(1) %>%
    ungroup()

  accepted_native <- adoption_decisions %>%
    filter(adoption_decision == "accept_native", present_native) %>%
    transmute(
      record_key,
      adopted_disease_standard = disease_standard,
      adopted_influenza_type = influenza_type,
      adopted_influenza_subtype = influenza_subtype,
      adoption_decision_id = decision_id,
      adoption_diff_category = diff_category,
      adoption_note = decision_note,
      adopted_disease_evidence_text = native_evidence_text,
      adopted_disease_source_method = coalesce(native_source_method, "native_adoption")
    ) %>%
    distinct()

  native_added <- country_rows %>%
    inner_join(accepted_native, by = "record_key", relationship = "many-to-many") %>%
    mutate(
      disease_raw = adopted_disease_standard,
      disease_standard = adopted_disease_standard,
      disease_evidence_text = adopted_disease_evidence_text,
      disease_evidence_location = "native_reviewed_candidate",
      disease_source_method = adopted_disease_source_method,
      disease_rule_id = adoption_diff_category,
      disease_confidence = "high",
      influenza_type = adopted_influenza_type,
      influenza_subtype = adopted_influenza_subtype,
      source_method = "v2_native_reviewed_adoption",
      needs_review = FALSE,
      review_status = "not_needed",
      review_decision_id = adoption_decision_id,
      review_note = adoption_note,
      clean_disease_label_clean = NA_character_,
      clean_disease_label_standard = NA_character_,
      clean_disease_label_standard_refined = NA_character_
    ) %>%
    select(all_of(names(base_reviewed)))

  bind_rows(base_reviewed, native_added) %>%
    distinct()
}

v2_review_queue_from_evidence <- function(evidence) {
  evidence %>%
    filter(needs_review | association_scope == "uncertain_focality") %>%
    transmute(
      review_id = paste(record_key, country_standard, disease_standard, "scope", sep = "::"),
      record_key,
      DonId,
      Title,
      article_url,
      country_standard,
      disease_standard,
      review_task = "scope_classification",
      evidence_text = scope_evidence_text,
      current_decision = association_scope,
      allowed_decisions = paste(
        c(
          "focal_event_country",
          "secondary_local_transmission_country",
          "imported_case_country",
          "travel_or_import_context_country",
          "historical_or_background_context_country",
          "lab_or_partner_context_country",
          "surveillance_or_sequence_context_country",
          "reported_from_other_countries_background",
          "not_final_event_country",
          "uncertain_focality"
        ),
        collapse = "|"
      ),
      reason_for_review = if_else(
        association_scope == "uncertain_focality",
        "uncertain deterministic scope",
        "clean evidence carried review flag"
      ),
      review_surface = if_else(
        association_scope == "uncertain_focality",
        "adjudication_candidate",
        "audit_only"
      ),
      llm_use_policy = if_else(
        association_scope == "uncertain_focality",
        "eligible_only_after_manual_subset_selection",
        "not_llm_input"
      )
    ) %>%
    distinct()
}

v2_scope_adjudication_candidates <- function(review_queue) {
  review_queue %>%
    filter(review_surface == "adjudication_candidate") %>%
    arrange(record_key, country_standard, disease_standard)
}

v2_apply_review_decisions <- function(evidence, decisions) {
  if (nrow(decisions) == 0) {
    return(list(
      evidence = evidence %>%
        mutate(
          final_association_scope = association_scope,
          final_scope_confidence = scope_confidence,
          final_review_source = NA_character_,
          final_review_note = NA_character_
        ),
      change_log = tibble::tibble(),
      unmatched = decisions
    ))
  }

  required <- c("review_id", "record_key", "decision_type", "decision_value")
  missing_cols <- setdiff(required, names(decisions))
  if (length(missing_cols) > 0) {
    stop("Review decisions missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  deduped <- decisions %>%
    filter(decision_type == "scope_classification") %>%
    distinct(review_id, .keep_all = TRUE)

  keyed <- evidence %>%
    mutate(review_id_join = paste(record_key, country_standard, disease_standard, "scope", sep = "::"))

  joined <- keyed %>%
    left_join(
      deduped %>%
        transmute(
          review_id_join = review_id,
          review_decision_id = review_id,
          decision_value,
          review_confidence = confidence,
          final_review_source = review_source,
          final_review_note = review_note
        ),
      by = "review_id_join"
    )

  changed <- joined %>%
    filter(!is.na(decision_value), decision_value != association_scope) %>%
    transmute(
      review_decision_id,
      record_key,
      country_standard,
      disease_standard,
      before_scope = association_scope,
      after_scope = decision_value,
      review_source = final_review_source,
      review_note = final_review_note
    )

  unmatched <- deduped %>%
    anti_join(keyed, by = c("review_id" = "review_id_join"))

  list(
    evidence = joined %>%
      transmute(
        !!!syms(names(evidence)),
        final_association_scope = coalesce(decision_value, association_scope),
        final_scope_confidence = coalesce(review_confidence, scope_confidence),
        final_review_source,
        final_review_note
      ),
    change_log = changed,
    unmatched = unmatched
  )
}
