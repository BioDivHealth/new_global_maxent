library(dplyr)
library(stringr)
library(tidyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_disease_compare_key <- c(
  "record_key",
  "disease_standard",
  "influenza_type",
  "influenza_subtype"
)

v2_influenza_label_standardization_cols <- c(
  "influenza_subtype",
  "canonical_disease_standard",
  "canonical_influenza_type",
  "canonical_influenza_subtype"
)

v2_influenza_policy_decision_cols <- c(
  "record_key",
  "disease_standard",
  "influenza_subtype",
  "policy_decision",
  "policy_diff_category",
  "policy_review_priority",
  "policy_note"
)

v2_disease_candidate_policy_decision_cols <- v2_influenza_policy_decision_cols

v2_normalize_compare_value <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::if_else(is.na(x), "", x)
}

v2_is_influenza_label <- function(x) {
  stringr::str_detect(x, stringr::regex("influenza", ignore_case = TRUE))
}

v2_is_specific_influenza <- function(disease_standard, influenza_subtype) {
  v2_is_influenza_label(disease_standard) & !is.na(influenza_subtype) & influenza_subtype != ""
}

v2_is_generic_influenza <- function(disease_standard, influenza_subtype) {
  v2_is_influenza_label(disease_standard) & (is.na(influenza_subtype) | influenza_subtype == "")
}

v2_is_wild_polio_context <- function(local_text) {
  str_detect(
    str_to_lower(coalesce(local_text, "")),
    regex(
      paste(
        c(
          "no wild poliovirus",
          "last (indigenous )?(case of )?wild",
          "prior to this outbreak",
          "has not had a wild",
          "one of only",
          "classified as",
          "epi-centre",
          "epicentre",
          "no case of wild"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

v2_is_wild_polio_event <- function(local_text) {
  str_detect(
    str_to_lower(coalesce(local_text, "")),
    regex(
      paste(
        c(
          "cases? (of |due to |caused by )?wild poliovirus",
          "wild poliovirus type [123].*(case|outbreak|isolat|import|confirm)",
          "confirmed wild poliovirus",
          "poliomyelitis due to wild poliovirus",
          "wild poliovirus.*isolat",
          "wild poliovirus.*outbreak"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

v2_is_dengue_context <- function(local_text) {
  str_detect(
    str_to_lower(coalesce(local_text, "")),
    regex(
      paste(
        c(
          "similar to dengue",
          "dengue-like",
          "test(ed)? (for )?dengue.*negative",
          "dengue.*negative",
          "suspicion of dengue",
          "suspected dengue",
          "compatible with.*dengue",
          "including dengue",
          "dengue and leptospirosis.*present",
          "circulation of zika.*dengue",
          "other pathogens.*dengue",
          "differential"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

v2_is_dengue_event <- function(local_text) {
  str_detect(
    str_to_lower(coalesce(local_text, "")),
    regex(
      paste(
        c(
          "dengue (fever )?(outbreak|epidemic)",
          "cases? of dengue",
          "dengue cases?",
          "reported .*dengue",
          "increase .*dengue",
          "number of cases of dengue",
          "confirmed .*dengue",
          "positive for dengue"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

v2_native_new_syndrome_labels <- c(
  "Acute diarrhoea",
  "Acute febrile illness",
  "Acute haemorrhagic fever syndrome",
  "Acute respiratory syndrome",
  "Acute watery diarrhoea",
  "Bloody diarrhoea",
  "Diarrhoeal disease",
  "Haemorrhagic fever syndrome"
)

v2_is_native_context_pattern <- function(local_text) {
  str_detect(
    str_to_lower(coalesce(local_text, "")),
    regex(
      paste(
        c(
          "needs? to be excluded",
          "ruled out",
          "excluded",
          "has excluded",
          "negative results",
          "tested negative",
          "were negative",
          "not positive",
          "none .* positive",
          "differential testing",
          "similar to",
          "same family as",
          "family includes",
          "family .* include",
          "vector of",
          "vectors? of",
          "initially suspected",
          "initially thought",
          "presumptive diagnosis",
          "other causes under investigation",
          "other pathogens",
          "multiple disease outbreaks",
          "multiple ongoing emergencies",
          "concurrent(ly)? .*outbreak",
          "coincides with",
          "same area as",
          "previous reports?",
          "previous outbreaks?",
          "site of an epidemic",
          "for comparison",
          "retrospectively listed",
          "also affected by",
          "in the context of",
          "last case of",
          "last reported",
          "co-circulat",
          "circulation of",
          "coinfection with",
          "also transmit",
          "living with hiv",
          "hiv/aids control",
          "present in central america",
          "co-circulation",
          "distinguish .* from",
          "shares? .*clinical signs",
          "misdiagnosed",
          "no cases? of",
          "vaccination campaign",
          "immuni[sz]ation",
          "for all arrivals",
          "visitors? arriving",
          "pilgrimage",
          "umra",
          "does not recommend",
          "clinically diagnosed as",
          "not limited to",
          "severe .* is more likely",
          "complications",
          "control program"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

v2_is_seeded_covid_coronavirus_context <- function(text) {
  str_detect(
    str_to_lower(coalesce(text, "")),
    regex("2019-ncov|sars-cov-2|covid-19|wuhan|hubei province", ignore_case = TRUE)
  )
}

v2_is_seeded_meningitis_weak <- function(text) {
  text <- str_to_lower(coalesce(text, ""))
  str_detect(
    text,
    regex(
      paste(
        c(
          "viral meningitis",
          "fungal meningitis",
          "suspected fungal meningitis",
          "meningoencephalitis",
          "aseptic meningitis"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  ) |
    (
      str_detect(text, regex("meningitis", ignore_case = TRUE)) &
        !str_detect(text, regex("meningococ|meningitidis|epidemic meningitis", ignore_case = TRUE))
    )
}

v2_seeded_text_has_subtype <- function(text, subtype) {
  text <- str_to_lower(coalesce(text, ""))
  subtype <- str_to_lower(coalesce(subtype, ""))
  out <- rep(FALSE, length(text))
  valid <- !is.na(subtype) & subtype != "" & !is.na(text) & text != ""
  out[valid] <- str_detect(text[valid], fixed(subtype[valid]))
  out
}

v2_apply_influenza_compare_standardization <- function(x, influenza_standardization) {
  missing_cols <- setdiff(v2_influenza_label_standardization_cols, names(influenza_standardization))
  if (length(missing_cols) > 0) {
    stop(
      "Missing influenza standardization columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  influenza_standardization <- influenza_standardization %>%
    transmute(
      influenza_subtype_lookup = v2_normalize_compare_value(influenza_subtype),
      canonical_disease_standard = v2_normalize_compare_value(canonical_disease_standard),
      canonical_influenza_type = v2_normalize_compare_value(canonical_influenza_type),
      canonical_influenza_subtype = v2_normalize_compare_value(canonical_influenza_subtype)
    ) %>%
    distinct(influenza_subtype_lookup, .keep_all = TRUE)

  x %>%
    mutate(
      disease_standard_original = disease_standard,
      influenza_type_original = influenza_type,
      influenza_subtype_original = influenza_subtype,
      influenza_subtype_lookup = v2_normalize_compare_value(influenza_subtype),
      is_influenza_before_standardization = v2_is_influenza_label(disease_standard)
    ) %>%
    left_join(influenza_standardization, by = "influenza_subtype_lookup") %>%
    mutate(
      disease_standard = case_when(
        is_influenza_before_standardization & influenza_subtype_lookup != "" &
          canonical_disease_standard != "" ~ canonical_disease_standard,
        TRUE ~ disease_standard
      ),
      influenza_type = case_when(
        is_influenza_before_standardization & influenza_subtype_lookup != "" &
          canonical_influenza_type != "" ~ canonical_influenza_type,
        is_influenza_before_standardization & influenza_subtype_lookup == "" ~ "influenza",
        TRUE ~ influenza_type
      ),
      influenza_subtype = case_when(
        is_influenza_before_standardization & influenza_subtype_lookup != "" &
          canonical_influenza_subtype != "" ~ canonical_influenza_subtype,
        TRUE ~ influenza_subtype
      )
    ) %>%
    select(
      -influenza_subtype_lookup,
      -canonical_disease_standard,
      -canonical_influenza_type,
      -canonical_influenza_subtype,
      -is_influenza_before_standardization
    )
}

v2_extract_influenza_subtype_from_label <- function(x) {
  label <- str_to_upper(coalesce(as.character(x), ""))
  specific <- str_match(label, "\\bH[0-9]+N[0-9]+\\b")[, 1]
  family <- str_match(label, "\\bH[0-9]+\\b")[, 1]
  na_if(coalesce(specific, family), "")
}

v2_add_influenza_compare_keys <- function(
  x,
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  missing_cols <- setdiff(v2_disease_compare_key, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Missing disease comparison columns for influenza compare keys: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  missing_standardization_cols <- setdiff(
    v2_influenza_label_standardization_cols,
    names(influenza_standardization)
  )
  if (length(missing_standardization_cols) > 0) {
    stop(
      "Missing influenza standardization columns: ",
      paste(missing_standardization_cols, collapse = ", "),
      call. = FALSE
    )
  }

  influenza_standardization <- influenza_standardization %>%
    transmute(
      influenza_subtype_lookup = v2_normalize_compare_value(influenza_subtype),
      canonical_disease_standard = v2_normalize_compare_value(canonical_disease_standard),
      canonical_influenza_type = v2_normalize_compare_value(canonical_influenza_type),
      canonical_influenza_subtype = v2_normalize_compare_value(canonical_influenza_subtype)
    ) %>%
    distinct(influenza_subtype_lookup, .keep_all = TRUE)

  x %>%
    mutate(
      disease_original_key = v2_normalize_compare_value(disease_standard),
      influenza_type_original_key = v2_normalize_compare_value(influenza_type),
      influenza_subtype_original_key = v2_normalize_compare_value(influenza_subtype),
      influenza_subtype_label_key = v2_normalize_compare_value(
        v2_extract_influenza_subtype_from_label(disease_standard)
      ),
      influenza_subtype_lookup = coalesce(
        na_if(influenza_subtype_original_key, ""),
        na_if(influenza_subtype_label_key, ""),
        ""
      ),
      is_influenza_compare_label = v2_is_influenza_label(disease_standard)
    ) %>%
    left_join(influenza_standardization, by = "influenza_subtype_lookup") %>%
    mutate(
      disease_compare_key = case_when(
        is_influenza_compare_label & influenza_subtype_lookup != "" &
          canonical_disease_standard != "" ~ canonical_disease_standard,
        is_influenza_compare_label & influenza_subtype_lookup == "" ~ "Influenza",
        TRUE ~ disease_original_key
      ),
      influenza_type_compare_key = case_when(
        is_influenza_compare_label & influenza_subtype_lookup != "" &
          canonical_influenza_type != "" ~ canonical_influenza_type,
        is_influenza_compare_label & influenza_subtype_lookup == "" ~ "influenza",
        TRUE ~ influenza_type_original_key
      ),
      influenza_subtype_compare_key = case_when(
        is_influenza_compare_label & influenza_subtype_lookup != "" &
          canonical_influenza_subtype != "" ~ canonical_influenza_subtype,
        is_influenza_compare_label ~ influenza_subtype_lookup,
        TRUE ~ influenza_subtype_original_key
      ),
      influenza_compare_changed =
        disease_compare_key != disease_original_key |
          influenza_type_compare_key != influenza_type_original_key |
          influenza_subtype_compare_key != influenza_subtype_original_key,
      influenza_compare_note = case_when(
        !is_influenza_compare_label ~ NA_character_,
        influenza_compare_changed & influenza_subtype_lookup != "" ~ paste0(
          "Influenza comparison key standardized from subtype ",
          influenza_subtype_lookup,
          "."
        ),
        influenza_compare_changed ~ "Generic influenza comparison key standardized.",
        TRUE ~ "Influenza comparison key already canonical."
      )
    ) %>%
    select(
      -disease_original_key,
      -influenza_type_original_key,
      -influenza_subtype_original_key,
      -influenza_subtype_label_key,
      -influenza_subtype_lookup,
      -canonical_disease_standard,
      -canonical_influenza_type,
      -canonical_influenza_subtype,
      -is_influenza_compare_label
    )
}

v2_collapse_disease_candidates <- function(x, source_label, influenza_standardization) {
  missing_cols <- setdiff(v2_disease_compare_key, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Missing disease comparison columns in ", source_label, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  optional_cols <- c(
    "local_evidence_text",
    "event_anchor_class",
    "event_anchor_reason",
    "candidate_promotion_status",
    "candidate_promotion_reason",
    "influenza_subtype_evidence_scope"
  )
  for (col in setdiff(optional_cols, names(x))) {
    x[[col]] <- NA_character_
  }

  x %>%
    v2_apply_influenza_compare_standardization(influenza_standardization) %>%
    mutate(
      across(all_of(v2_disease_compare_key), v2_normalize_compare_value),
      disease_needs_review = as.character(disease_needs_review),
      disease_evidence_text = stringr::str_squish(as.character(disease_evidence_text)),
      local_evidence_text = stringr::str_squish(as.character(local_evidence_text))
    ) %>%
    group_by(across(all_of(v2_disease_compare_key))) %>%
    summarise(
      candidate_rows = n(),
      DonId = dplyr::first(na.omit(DonId)),
      record_id = dplyr::first(na.omit(record_id)),
      disease_standard_original = paste(sort(unique(na.omit(disease_standard_original))), collapse = " | "),
      influenza_type_original = paste(sort(unique(na.omit(influenza_type_original))), collapse = " | "),
      influenza_subtype_original = paste(sort(unique(na.omit(influenza_subtype_original))), collapse = " | "),
      disease_raw = paste(sort(unique(na.omit(disease_raw))), collapse = " | "),
      disease_evidence_location = paste(sort(unique(na.omit(disease_evidence_location))), collapse = " | "),
      disease_source_method = paste(sort(unique(na.omit(disease_source_method))), collapse = " | "),
      disease_rule_id = paste(sort(unique(na.omit(disease_rule_id))), collapse = " | "),
      disease_confidence = paste(sort(unique(na.omit(disease_confidence))), collapse = " | "),
      disease_needs_review = paste(sort(unique(na.omit(disease_needs_review))), collapse = " | "),
      disease_evidence_text = paste(head(unique(na.omit(disease_evidence_text)), 3), collapse = " || "),
      local_evidence_text = paste(head(unique(na.omit(local_evidence_text)), 3), collapse = " || "),
      event_anchor_class = paste(sort(unique(na.omit(event_anchor_class))), collapse = " | "),
      event_anchor_reason = paste(head(unique(na.omit(event_anchor_reason)), 3), collapse = " || "),
      candidate_promotion_status = paste(sort(unique(na.omit(candidate_promotion_status))), collapse = " | "),
      candidate_promotion_reason = paste(head(unique(na.omit(candidate_promotion_reason)), 3), collapse = " || "),
      influenza_subtype_evidence_scope = paste(sort(unique(na.omit(influenza_subtype_evidence_scope))), collapse = " | "),
      .groups = "drop"
    ) %>%
    mutate(
      candidate_source = source_label,
      is_influenza = v2_is_influenza_label(disease_standard),
      is_specific_influenza = v2_is_specific_influenza(disease_standard, influenza_subtype),
      is_generic_influenza = v2_is_generic_influenza(disease_standard, influenza_subtype)
    )
}

v2_clean_resolution_flags <- function(resolution_seed) {
  resolution_seed %>%
    transmute(
      record_key = v2_normalize_compare_value(record_key),
      seeded_disease_standard = v2_normalize_compare_value(resolved_disease_label),
      seeded_influenza_type = v2_normalize_compare_value(influenza_type),
      seeded_influenza_subtype = v2_normalize_compare_value(influenza_subtype),
      raw_was_generic_influenza = v2_normalize_compare_value(disease_label_raw) == "Influenza",
      subtype_evidence_present = v2_normalize_compare_value(influenza_subtype_evidence_span) != "",
      disease_refinement_source = v2_normalize_compare_value(disease_refinement_source)
    ) %>%
    group_by(
      record_key,
      seeded_disease_standard,
      seeded_influenza_type,
      seeded_influenza_subtype
    ) %>%
    summarise(
      seeded_from_generic_influenza = any(raw_was_generic_influenza),
      seeded_has_subtype_evidence = any(subtype_evidence_present),
      seeded_refinement_sources = paste(sort(unique(na.omit(disease_refinement_source))), collapse = " | "),
      .groups = "drop"
    )
}

v2_record_level_influenza_context <- function(seed_collapsed, native_collapsed) {
  seeded <- seed_collapsed %>%
    filter(is_influenza) %>%
    group_by(record_key) %>%
    summarise(
      seeded_has_generic_influenza = any(is_generic_influenza),
      seeded_has_specific_influenza = any(is_specific_influenza),
      seeded_influenza_subtypes = paste(sort(unique(influenza_subtype[influenza_subtype != ""])), collapse = "|"),
      .groups = "drop"
    )

  native <- native_collapsed %>%
    filter(is_influenza) %>%
    group_by(record_key) %>%
    summarise(
      native_has_generic_influenza = any(is_generic_influenza),
      native_has_specific_influenza = any(is_specific_influenza),
      native_influenza_subtypes = paste(sort(unique(influenza_subtype[influenza_subtype != ""])), collapse = "|"),
      .groups = "drop"
    )

  full_join(seeded, native, by = "record_key") %>%
    mutate(
      across(
        c(
          seeded_has_generic_influenza,
          seeded_has_specific_influenza,
          native_has_generic_influenza,
          native_has_specific_influenza
        ),
        ~ dplyr::coalesce(.x, FALSE)
      ),
      across(
        c(seeded_influenza_subtypes, native_influenza_subtypes),
        ~ dplyr::coalesce(.x, "")
      )
    )
}

v2_record_level_native_disease_context <- function(native_collapsed) {
  native_collapsed %>%
    group_by(record_key) %>%
    summarise(
      native_has_specific_polio = any(
        disease_standard %in% c(
          "Wild poliovirus",
          "Vaccine-derived poliovirus",
          "Circulating vaccine-derived poliovirus type 1",
          "Circulating vaccine-derived poliovirus type 2"
        )
      ),
      native_polio_labels = paste(
        sort(unique(disease_standard[
          disease_standard %in% c(
            "Wild poliovirus",
            "Vaccine-derived poliovirus",
            "Circulating vaccine-derived poliovirus type 1",
            "Circulating vaccine-derived poliovirus type 2"
          )
        ])),
        collapse = "|"
      ),
      .groups = "drop"
    )
}

v2_read_influenza_policy_decisions <- function(path = who_don_v2_rules_dir("influenza_subtype_policy_decisions.csv")) {
  if (!file.exists(path)) {
    return(tibble::tibble(
      record_key = character(),
      disease_standard = character(),
      influenza_subtype = character(),
      policy_decision = character(),
      policy_diff_category = character(),
      policy_review_priority = character(),
      policy_note = character()
    ))
  }

  v2_read_csv(path, v2_influenza_policy_decision_cols) %>%
    mutate(
      record_key = v2_normalize_compare_value(record_key),
      disease_standard = v2_normalize_compare_value(disease_standard),
      influenza_subtype = v2_normalize_compare_value(influenza_subtype)
    )
}

v2_read_disease_candidate_policy_decisions <- function(
  path = who_don_v2_rules_dir("disease_candidate_policy_decisions.csv")
) {
  if (!file.exists(path)) {
    return(tibble::tibble(
      record_key = character(),
      disease_standard = character(),
      influenza_subtype = character(),
      policy_decision = character(),
      policy_diff_category = character(),
      policy_review_priority = character(),
      policy_note = character()
    ))
  }

  v2_read_csv(path, v2_disease_candidate_policy_decision_cols) %>%
    mutate(
      record_key = v2_normalize_compare_value(record_key),
      disease_standard = v2_normalize_compare_value(disease_standard),
      influenza_subtype = v2_normalize_compare_value(influenza_subtype)
    )
}

v2_combine_candidate_policy_decisions <- function(...) {
  policies <- bind_rows(...) %>%
    mutate(
      record_key = v2_normalize_compare_value(record_key),
      disease_standard = v2_normalize_compare_value(disease_standard),
      influenza_subtype = v2_normalize_compare_value(influenza_subtype)
    )

  duplicated_keys <- policies %>%
    count(record_key, disease_standard, influenza_subtype, name = "policy_rows") %>%
    filter(policy_rows > 1)

  if (nrow(duplicated_keys) > 0) {
    stop(
      "Duplicate disease candidate policy decisions for: ",
      paste(
        paste(
          duplicated_keys$record_key,
          duplicated_keys$disease_standard,
          duplicated_keys$influenza_subtype,
          sep = " / "
        ),
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  policies
}

v2_candidate_adoption_decisions <- function(diff) {
  diff %>%
    filter(diff_category != "exact_match") %>%
    transmute(
      decision_id = paste(record_key, disease_standard, influenza_type, influenza_subtype, diff_category, sep = "::"),
      record_key,
      DonId = coalesce(DonId_native, DonId_seeded),
      record_id = coalesce(record_id_native, record_id_seeded),
      disease_standard,
      influenza_type,
      influenza_subtype,
      diff_category,
      review_priority,
      present_seeded,
      present_native,
      candidate_rows_seeded,
      candidate_rows_native,
      adoption_decision = case_when(
        str_detect(diff_category, "^accepted_native") ~ "accept_native",
        str_detect(diff_category, "^accepted_seeded") ~ "keep_seeded",
        diff_category == "acceptable_alias_cleanup" ~ "keep_seeded",
        str_detect(diff_category, "^rejected_native") ~ "reject_native_context",
        str_detect(diff_category, "^rejected_seeded") ~ "reject_seeded_weak",
        str_detect(diff_category, "^deferred_native") ~ "reject_native_context",
        str_detect(diff_category, "^manual_review") ~ "needs_manual_review",
        diff_category == "comparison_key_mismatch" ~ "standardize",
        TRUE ~ "needs_manual_review"
      ),
      adoption_source = case_when(
        str_detect(diff_category, "^accepted_native") ~ "comparison_policy",
        str_detect(diff_category, "^accepted_seeded") ~ "comparison_policy",
        str_detect(diff_category, "^rejected_native") ~ "comparison_policy",
        str_detect(diff_category, "^rejected_seeded") ~ "comparison_policy",
        str_detect(diff_category, "^deferred_native") ~ "row_policy",
        str_detect(diff_category, "^manual_review") ~ "review_queue",
        diff_category == "acceptable_alias_cleanup" ~ "comparison_policy",
        diff_category == "comparison_key_mismatch" ~ "standardization_review",
        TRUE ~ "review_queue"
      ),
      decision_note = case_when(
        !is.na(policy_note) & policy_note != "" ~ policy_note,
        diff_category == "accepted_seeded_legacy_exception_policy" ~
          "Conservative adoption policy keeps this seeded clean candidate as an explicit legacy exception because native extraction did not recover a replacement.",
        diff_category == "rejected_native_unadopted_candidate_policy" ~
          "Conservative adoption policy does not adopt native-new candidates unless an explicit deterministic acceptance policy supports them.",
        adoption_decision == "accept_native" ~ "Native candidate is accepted by deterministic comparison policy.",
        adoption_decision == "keep_seeded" ~ "Seeded candidate is kept by deterministic comparison policy.",
        adoption_decision == "reject_native_context" ~ "Native candidate is rejected because evidence is contextual or non-event.",
        adoption_decision == "reject_seeded_weak" ~ "Seeded candidate is rejected because native evidence supports a more specific or cleaner event label.",
        adoption_decision == "standardize" ~ "Native and seeded candidates require key/label standardization review.",
        TRUE ~ "Candidate difference requires manual review before native adoption."
      ),
      seeded_evidence_text = disease_evidence_text_seeded,
      native_evidence_text = disease_evidence_text_native,
      seeded_source_method = disease_source_method_seeded,
      native_source_method = disease_source_method_native
    ) %>%
    arrange(
      match(adoption_decision, c("needs_manual_review", "standardize", "accept_native", "keep_seeded", "reject_native_context", "reject_seeded_weak")),
      record_key,
      disease_standard
    )
}

v2_compare_disease_candidates <- function(
  seeded_candidates = v2_read_csv(who_don_v2_output_dir("candidates", "who_don_disease_candidates.csv")),
  native_candidates = v2_read_csv(who_don_v2_output_dir("candidates", "who_don_disease_candidates_native.csv")),
  resolution_seed = v2_read_csv(who_don_v2_rules_dir("disease_resolution_seed_from_clean.csv")),
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv")),
  influenza_policy_decisions = v2_read_influenza_policy_decisions(),
  disease_candidate_policy_decisions = v2_read_disease_candidate_policy_decisions()
) {
  seeded <- v2_collapse_disease_candidates(seeded_candidates, "seeded", influenza_standardization)
  native <- v2_collapse_disease_candidates(native_candidates, "native", influenza_standardization)
  clean_flags <- v2_clean_resolution_flags(resolution_seed)
  influenza_context <- v2_record_level_influenza_context(seeded, native)
  native_disease_context <- v2_record_level_native_disease_context(native)
  candidate_policy_decisions <- v2_combine_candidate_policy_decisions(
    influenza_policy_decisions,
    disease_candidate_policy_decisions
  )

  diff <- full_join(
    seeded %>% rename_with(~ paste0(.x, "_seeded"), -all_of(v2_disease_compare_key)),
    native %>% rename_with(~ paste0(.x, "_native"), -all_of(v2_disease_compare_key)),
    by = v2_disease_compare_key
  ) %>%
    left_join(
      clean_flags,
      by = c(
        "record_key" = "record_key",
        "disease_standard" = "seeded_disease_standard",
        "influenza_type" = "seeded_influenza_type",
        "influenza_subtype" = "seeded_influenza_subtype"
      )
    ) %>%
    left_join(
      candidate_policy_decisions,
      by = c("record_key", "disease_standard", "influenza_subtype")
    ) %>%
    left_join(influenza_context, by = "record_key") %>%
    left_join(native_disease_context, by = "record_key") %>%
    mutate(
      present_seeded = !is.na(candidate_rows_seeded),
      present_native = !is.na(candidate_rows_native),
      is_influenza = v2_is_influenza_label(disease_standard),
      is_specific_influenza = v2_is_specific_influenza(disease_standard, influenza_subtype),
      is_generic_influenza = v2_is_generic_influenza(disease_standard, influenza_subtype),
      seeded_from_generic_influenza = coalesce(seeded_from_generic_influenza, FALSE),
      seeded_has_subtype_evidence = coalesce(seeded_has_subtype_evidence, FALSE),
      seeded_has_specific_influenza = coalesce(seeded_has_specific_influenza, FALSE),
      seeded_has_generic_influenza = coalesce(seeded_has_generic_influenza, FALSE),
      native_has_specific_influenza = coalesce(native_has_specific_influenza, FALSE),
      native_has_generic_influenza = coalesce(native_has_generic_influenza, FALSE),
      native_has_specific_polio = coalesce(native_has_specific_polio, FALSE),
      native_polio_labels = coalesce(native_polio_labels, ""),
      seeded_influenza_subtypes = coalesce(seeded_influenza_subtypes, ""),
      native_influenza_subtypes = coalesce(native_influenza_subtypes, ""),
      native_wild_polio_context = disease_standard == "Wild poliovirus" &
        v2_is_wild_polio_context(local_evidence_text_native),
      native_wild_polio_event = disease_standard == "Wild poliovirus" &
        v2_is_wild_polio_event(local_evidence_text_native),
      native_dengue_context = disease_standard == "Dengue" &
        v2_is_dengue_context(local_evidence_text_native),
      native_dengue_event = disease_standard == "Dengue" &
        v2_is_dengue_event(local_evidence_text_native),
      native_context_pattern = v2_is_native_context_pattern(local_evidence_text_native),
      same_record_influenza_subtype_in_other_source = case_when(
        present_seeded & !present_native & is_influenza & influenza_subtype != "" ~
          str_detect(paste0("|", native_influenza_subtypes, "|"), fixed(paste0("|", influenza_subtype, "|"))),
        !present_seeded & present_native & is_influenza & influenza_subtype != "" ~
          str_detect(paste0("|", seeded_influenza_subtypes, "|"), fixed(paste0("|", influenza_subtype, "|"))),
        TRUE ~ FALSE
      ),
      diff_category = case_when(
        present_seeded & present_native ~ "exact_match",
        !is.na(policy_diff_category) & policy_diff_category != "" ~ policy_diff_category,
        same_record_influenza_subtype_in_other_source ~ "comparison_key_mismatch",
        !present_seeded & present_native & is_specific_influenza &
          !seeded_has_specific_influenza &
          str_detect(coalesce(influenza_subtype_evidence_scope_native, ""), "title|local_event") ~
          "accepted_native_specific_influenza_policy",
        !present_seeded & present_native & disease_standard == "Encephalitis" ~
          "rejected_native_clinical_manifestation_policy",
        !present_seeded & present_native & disease_standard == "Wild poliovirus" &
          native_wild_polio_context ~
          "rejected_native_wild_polio_context_policy",
        !present_seeded & present_native & disease_standard == "Wild poliovirus" &
          native_wild_polio_event & !native_wild_polio_context ~
          "accepted_native_specific_polio_policy",
        !present_seeded & present_native & disease_standard == "Dengue" &
          native_dengue_context ~
          "rejected_native_dengue_context_policy",
        !present_seeded & present_native & disease_standard == "Dengue" &
          native_dengue_event & !native_dengue_context ~
          "accepted_native_dengue_event_policy",
        !present_seeded & present_native & disease_standard == "Influenza" &
          !str_detect(coalesce(event_anchor_class_native, ""), "event_title") ~
          "rejected_native_generic_influenza_context_policy",
        !present_seeded & present_native & native_context_pattern &
          !str_detect(coalesce(event_anchor_class_native, ""), "event_title") ~
          "rejected_native_context_pattern_policy",
        !present_seeded & present_native &
          disease_standard %in% v2_native_new_syndrome_labels &
          !str_detect(coalesce(event_anchor_class_native, ""), "event_title") ~
          "rejected_native_syndrome_context_policy",
        !present_seeded & present_native &
          str_detect(coalesce(event_anchor_class_native, ""), "event_title") &
          disease_standard != "Encephalitis" ~
          "accepted_native_title_disease_policy",
        present_seeded & !present_native & disease_standard == "Influenza" &
          native_has_specific_influenza ~
          "rejected_seeded_generic_influenza_superseded_policy",
        present_seeded & !present_native & disease_standard == "Poliomyelitis" &
          native_has_specific_polio ~
          "rejected_seeded_generic_polio_superseded_policy",
        present_seeded & !present_native & disease_standard == "Meningococcal disease" &
          v2_is_seeded_meningitis_weak(disease_evidence_text_seeded) ~
          "rejected_seeded_generic_meningitis_policy",
        present_seeded & !present_native &
          disease_standard %in% c("Middle East Respiratory Syndrome (MERS)", "Severe Acute Respiratory Syndrome (SARS)") &
          v2_is_seeded_covid_coronavirus_context(disease_evidence_text_seeded) ~
          "rejected_seeded_covid_coronavirus_context_policy",
        present_seeded & !present_native & is_specific_influenza &
          native_has_generic_influenza &
          v2_seeded_text_has_subtype(disease_evidence_text_seeded, influenza_subtype) ~
          "accepted_seeded_specific_influenza_policy",
        present_seeded & !present_native & is_specific_influenza &
          native_has_generic_influenza & seeded_from_generic_influenza &
          !seeded_has_subtype_evidence ~ "acceptable_alias_cleanup",
        present_seeded & !present_native & is_specific_influenza &
          native_has_generic_influenza ~ "manual_review_native_more_generic",
        !present_seeded & present_native & is_generic_influenza &
          seeded_has_specific_influenza ~ "manual_review_native_more_generic",
        !present_seeded & present_native & is_specific_influenza &
          seeded_has_generic_influenza ~ "manual_review_native_more_specific",
        present_seeded & !present_native & is_specific_influenza &
          native_has_specific_influenza & seeded_influenza_subtypes != native_influenza_subtypes ~
          "influenza_subtype_changed",
        !present_seeded & present_native & is_specific_influenza &
          seeded_has_specific_influenza & seeded_influenza_subtypes != native_influenza_subtypes ~
          "influenza_subtype_changed",
        present_seeded & !present_native ~ "accepted_seeded_legacy_exception_policy",
        !present_seeded & present_native ~ "rejected_native_unadopted_candidate_policy",
        TRUE ~ "needs_review"
      ),
      review_priority = case_when(
        diff_category == "exact_match" ~ "none",
        !is.na(policy_review_priority) & policy_review_priority != "" ~ policy_review_priority,
        diff_category == "comparison_key_mismatch" ~ "low",
        diff_category == "acceptable_alias_cleanup" ~ "low",
        diff_category == "accepted_native_specific_influenza_policy" ~ "low",
        diff_category == "accepted_native_dengue_event_policy" ~ "low",
        diff_category == "accepted_native_specific_polio_policy" ~ "low",
        diff_category == "accepted_native_title_disease_policy" ~ "low",
        diff_category == "accepted_seeded_legacy_exception_policy" ~ "low",
        diff_category == "accepted_seeded_specific_influenza_policy" ~ "low",
        diff_category == "rejected_seeded_generic_influenza_superseded_policy" ~ "low",
        diff_category == "rejected_seeded_generic_polio_superseded_policy" ~ "low",
        diff_category == "rejected_seeded_generic_meningitis_policy" ~ "low",
        diff_category == "rejected_seeded_covid_coronavirus_context_policy" ~ "low",
        diff_category == "rejected_native_clinical_manifestation_policy" ~ "low",
        diff_category == "rejected_native_dengue_context_policy" ~ "low",
        diff_category == "rejected_native_generic_influenza_context_policy" ~ "low",
        diff_category == "rejected_native_context_pattern_policy" ~ "low",
        diff_category == "rejected_native_syndrome_context_policy" ~ "low",
        diff_category == "rejected_native_wild_polio_context_policy" ~ "low",
        diff_category == "rejected_native_unadopted_candidate_policy" ~ "low",
        diff_category %in% c(
          "manual_review_native_more_generic",
          "manual_review_native_more_specific",
          "manual_review_native_missing_seeded",
          "influenza_subtype_changed"
        ) ~ "medium",
        diff_category == "manual_review_native_new_candidate" ~ "medium",
        TRUE ~ "medium"
      )
    ) %>%
    arrange(diff_category, record_key, disease_standard, influenza_type, influenza_subtype)

  unmatched_seeded <- diff %>%
    filter(present_seeded, !present_native) %>%
    arrange(review_priority, record_key, disease_standard)

  new_native <- diff %>%
    filter(!present_seeded, present_native) %>%
    arrange(review_priority, record_key, disease_standard)

  native_new_candidate_review <- new_native %>%
    mutate(
      native_new_review_class = case_when(
        !is.na(policy_decision) & policy_decision == "reject_native_context" ~
          "reviewed_rejected_native_context",
        !is.na(policy_decision) & policy_decision == "defer_native_differential_diagnosis" ~
          "reviewed_deferred_native_differential",
        diff_category == "accepted_native_specific_influenza_policy" ~
          "reviewed_accepted_native_specific_influenza",
        diff_category %in% c("accepted_native_h3n2_policy", "accepted_native_h5_family_policy") ~
          "reviewed_accepted_native_specific_influenza",
        diff_category == "accepted_native_title_disease_policy" ~
          "reviewed_accepted_native_title_disease",
        diff_category == "rejected_native_clinical_manifestation_policy" ~
          "reviewed_rejected_native_clinical_manifestation",
        diff_category == "accepted_native_specific_polio_policy" ~
          "reviewed_accepted_native_specific_polio",
        diff_category == "rejected_native_wild_polio_context_policy" ~
          "reviewed_rejected_native_wild_polio_context",
        diff_category == "accepted_native_dengue_event_policy" ~
          "reviewed_accepted_native_dengue_event",
        diff_category == "rejected_native_dengue_context_policy" ~
          "reviewed_rejected_native_dengue_context",
        diff_category == "rejected_native_generic_influenza_context_policy" ~
          "reviewed_rejected_native_generic_influenza_context",
        diff_category == "rejected_native_context_pattern_policy" ~
          "reviewed_rejected_native_context_pattern",
        diff_category == "rejected_native_syndrome_context_policy" ~
          "reviewed_rejected_native_syndrome_context",
        diff_category == "comparison_key_mismatch" ~
          "comparison_key_mismatch",
        diff_category == "manual_review_native_new_candidate" ~
          "manual_review_native_new_candidate",
        is_specific_influenza &
          str_detect(coalesce(influenza_subtype_evidence_scope_native, ""), "title|local_event") ~
          "specific_influenza_event_candidate",
        disease_standard %in% c("Acute respiratory infection", "COVID-19", "Malaria", "Respiratory illness") &
          str_detect(coalesce(event_anchor_class_native, ""), "event_title|event_summary") ~
          "background_prone_label_with_strong_anchor",
        disease_standard %in% c("Acute respiratory infection", "COVID-19", "Malaria", "Respiratory illness") ~
          "background_prone_label_needs_review",
        str_detect(coalesce(event_anchor_class_native, ""), "event_overview") ~
          "event_overview_needs_review",
        TRUE ~ "native_new_candidate_needs_review"
      ),
      native_new_review_note = case_when(
        native_new_review_class == "comparison_key_mismatch" ~
          "Native and seeded share the same record-level influenza subtype but differ in label/type standardization.",
        native_new_review_class == "manual_review_native_new_candidate" ~
          "Native-new disease candidate remains plausible but needs disease-specific manual review before adoption.",
        native_new_review_class == "reviewed_rejected_native_context" ~
          policy_note,
        native_new_review_class == "reviewed_deferred_native_differential" ~
          policy_note,
        native_new_review_class == "reviewed_accepted_native_specific_influenza" ~
          "Specific native influenza subtype has title/local-event evidence and no unresolved competing seeded subtype.",
        native_new_review_class == "reviewed_accepted_native_title_disease" ~
          "Disease alias appears in the DON title and is accepted as a native title disease candidate.",
        native_new_review_class == "reviewed_rejected_native_clinical_manifestation" ~
          "Generic encephalitis is treated as a clinical manifestation/context label unless it was already seeded or title-supported as the event disease.",
        native_new_review_class == "reviewed_accepted_native_specific_polio" ~
          "Wild poliovirus has explicit case, outbreak, confirmation, importation, or isolation wording and no negative/historical wild-polio context trigger.",
        native_new_review_class == "reviewed_rejected_native_wild_polio_context" ~
          "Wild poliovirus appears in negative, historical, endemicity, or prior-outbreak context rather than as the DON event label.",
        native_new_review_class == "reviewed_accepted_native_dengue_event" ~
          "Dengue has explicit case, outbreak, epidemic, increase, confirmation, or positive-test wording and no dengue-context trigger.",
        native_new_review_class == "reviewed_rejected_native_dengue_context" ~
          "Dengue appears in differential diagnosis, negative testing, similarity, co-circulation, or regional context rather than as the DON event label.",
        native_new_review_class == "reviewed_rejected_native_generic_influenza_context" ~
          "Generic influenza appears only in non-title overview/summary context and is not adopted as a native-new event disease.",
        native_new_review_class == "reviewed_rejected_native_context_pattern" ~
          "Disease alias appears in explicit non-event context such as ruled-out, negative, differential, similar-to, same-family, vector, or concurrent-outbreak wording.",
        native_new_review_class == "reviewed_rejected_native_syndrome_context" ~
          "Broad clinical syndrome or symptom label appears only outside the DON title and is not adopted as a native-new event disease.",
        native_new_review_class == "specific_influenza_event_candidate" ~
          "Specific influenza subtype has native title/local evidence but no seeded candidate key.",
        native_new_review_class == "background_prone_label_with_strong_anchor" ~
          "Background-prone label survived the promotion gate because title or summary supports it.",
        native_new_review_class == "background_prone_label_needs_review" ~
          "Background-prone label remains as native-new and should be checked manually.",
        native_new_review_class == "event_overview_needs_review" ~
          "Native-new label has overview evidence but was absent from seeded v2.",
        TRUE ~ "Native-new label needs manual review before adoption."
      )
    ) %>%
    select(
      native_new_review_class,
      native_new_review_note,
      everything()
    )

  influenza_change_review <- diff %>%
    filter(is_influenza, diff_category != "exact_match") %>%
    arrange(review_priority, record_key, disease_standard)

  subtype_changed_classification <- influenza_change_review %>%
    filter(diff_category == "influenza_subtype_changed") %>%
    mutate(
      subtype_change_review_class = case_when(
        !is.na(policy_decision) & policy_decision == "accept_native_h5_family_event_label" ~
          "accepted_native_h5_family_policy",
        record_key %in% c("2021-DON354", "2023-DON434", "2023-DON453") ~
          "native_h5_family_preferred_over_seeded_h5n1_background",
        record_key == "2024-DON504" ~
          "native_multi_subtype_title_gap",
        record_key == "6243fca4-7a8c-450c-8f52-260d88c7a422" ~
          "manual_review_seasonal_h1_h3n2_record",
        record_key == "9c1da48b-b69d-4843-83b8-cc5bafdd5854" ~
          "manual_review_h5_h5n1_poultry_context",
        TRUE ~ "manual_review_influenza_subtype_change"
      ),
      subtype_change_recommended_action = case_when(
        subtype_change_review_class == "accepted_native_h5_family_policy" ~
          policy_note,
        subtype_change_review_class == "native_h5_family_preferred_over_seeded_h5n1_background" ~
          "Prefer native H5 family event label; keep seeded H5N1 as background/context evidence unless manual review says otherwise.",
        subtype_change_review_class == "native_multi_subtype_title_gap" ~
          "Update native extraction to emit one candidate per explicit title subtype for coinfection records.",
        subtype_change_review_class == "manual_review_seasonal_h1_h3n2_record" ~
          "Manual review needed; record includes broad seasonal influenza text with H1 and H3N2 evidence.",
        subtype_change_review_class == "manual_review_h5_h5n1_poultry_context" ~
          "Manual review needed; distinguish event subtype from poultry/background context.",
        TRUE ~ "Manual review needed before adoption."
      )
    ) %>%
    select(
      subtype_change_review_class,
      subtype_change_recommended_action,
      everything()
    )

  influenza_review_queue <- influenza_change_review %>%
    transmute(
      review_id = paste(record_key, disease_standard, influenza_type, influenza_subtype, diff_category, sep = "::"),
      record_key,
      DonId = coalesce(DonId_native, DonId_seeded),
      record_id = coalesce(record_id_native, record_id_seeded),
      disease_standard,
      influenza_type,
      influenza_subtype,
      diff_category,
      review_priority,
      seeded_influenza_subtypes,
      native_influenza_subtypes,
      seeded_from_generic_influenza,
      seeded_has_subtype_evidence,
      seeded_evidence_text = disease_evidence_text_seeded,
      native_evidence_text = disease_evidence_text_native,
      current_decision = "review_before_adoption",
      allowed_decisions = "accept_native|keep_seeded|standardize|reject_native|needs_llm_or_manual",
      reason_for_review = case_when(
        diff_category == "acceptable_alias_cleanup" ~
          "Seeded subtype came from generic influenza clean-final provenance without explicit subtype evidence.",
        diff_category == "manual_review_native_more_generic" ~
          "Native extraction kept influenza broader than seeded v2.",
        diff_category == "manual_review_native_more_specific" ~
          "Native extraction found more specific influenza subtype evidence than seeded v2.",
        diff_category == "influenza_subtype_changed" ~
          "Seeded and native extraction disagree on influenza subtype.",
        diff_category == "comparison_key_mismatch" ~
          "Seeded and native have the same record-level influenza subtype but different label/type keys.",
        diff_category == "accepted_native_h5_family_policy" ~
          "Reviewed policy decision accepts native H5 family event label and treats seeded H5N1 as background/context evidence.",
        diff_category == "accepted_native_specific_influenza_policy" ~
          "Reviewed rule accepts native specific influenza subtype because title/local-event evidence supports it and seeded has no competing subtype.",
        TRUE ~ "Influenza candidate differs between seeded and native disease extraction."
      )
    ) %>%
    distinct()

  summary <- diff %>%
    count(diff_category, review_priority, name = "comparison_rows") %>%
    arrange(match(review_priority, c("high", "medium", "low", "none")), diff_category)

  summary_by_disease <- diff %>%
    count(diff_category, disease_standard, name = "comparison_rows") %>%
    arrange(diff_category, desc(comparison_rows), disease_standard)

  summary_by_anchor <- diff %>%
    count(
      diff_category,
      event_anchor_class_native,
      influenza_subtype_evidence_scope_native,
      name = "comparison_rows"
    ) %>%
    arrange(diff_category, desc(comparison_rows))

  adoption_decisions <- v2_candidate_adoption_decisions(diff)

  list(
    diff = diff,
    summary = summary,
    summary_by_disease = summary_by_disease,
    summary_by_anchor = summary_by_anchor,
    unmatched_seeded = unmatched_seeded,
    new_native = new_native,
    native_new_candidate_review = native_new_candidate_review,
    influenza_change_review = influenza_change_review,
    subtype_changed_classification = subtype_changed_classification,
    influenza_review_queue = influenza_review_queue,
    adoption_decisions = adoption_decisions
  )
}
