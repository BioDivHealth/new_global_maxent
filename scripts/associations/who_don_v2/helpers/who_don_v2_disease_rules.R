library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))

v2_influenza_specific_pattern <- "(h[0-9]+n[0-9]+|a\\(h[0-9]+n[0-9]+\\))"
v2_influenza_h_only_pattern <- "(^|[^a-z0-9])h[0-9]+([^a-z0-9]|$)"

v2_resolved_disease_label <- function(x) {
  dplyr::coalesce(
    x$disease_label_standard_refined,
    x$disease_label_standard,
    x$disease_label_clean
  )
}

v2_disease_resolution_seed_from_clean <- function(clean_final) {
  clean_final %>%
    mutate(resolved_disease_label = v2_resolved_disease_label(pick(everything()))) %>%
    transmute(
      record_key,
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      country_standard,
      disease_label_raw,
      disease_label_clean,
      disease_label_standard,
      disease_label_standard_refined,
      resolved_disease_label,
      influenza_type,
      influenza_subtype,
      influenza_subtype_candidates,
      influenza_subtype_evidence_span,
      disease_refinement_source,
      final_source,
      source_note = case_when(
        disease_label_raw != resolved_disease_label ~
          "Clean final raw disease label was resolved or refined before v2 seeding.",
        TRUE ~ "Clean final raw disease label matched resolved label."
      )
    ) %>%
    distinct()
}

v2_is_safe_disease_alias <- function(alias, disease_standard, influenza_type, influenza_subtype) {
  alias_norm <- str_to_lower(str_squish(alias))
  standard_norm <- str_to_lower(str_squish(disease_standard))
  is_influenza <- str_detect(alias_norm, "influenza") | str_detect(standard_norm, "influenza") |
    !is.na(influenza_type) | !is.na(influenza_subtype)
  alias_has_specific <- str_detect(alias_norm, v2_influenza_specific_pattern)
  alias_has_h_only <- str_detect(alias_norm, v2_influenza_h_only_pattern)
  standard_has_specific <- str_detect(standard_norm, v2_influenza_specific_pattern)
  standard_has_h_only <- str_detect(standard_norm, v2_influenza_h_only_pattern)

  case_when(
    !is_influenza ~ TRUE,
    alias_norm == "influenza" & (standard_has_specific | standard_has_h_only) ~ FALSE,
    alias_has_h_only & standard_has_specific & !alias_has_specific ~ FALSE,
    standard_has_specific & !alias_has_specific ~ FALSE,
    TRUE ~ TRUE
  )
}

v2_disease_alias_type <- function(alias, disease_standard, is_generic_label) {
  alias_norm <- str_to_lower(str_squish(alias))
  standard_norm <- str_to_lower(str_squish(disease_standard))

  case_when(
    alias_norm == "influenza" ~ "influenza_generic",
    str_detect(alias_norm, v2_influenza_specific_pattern) |
      str_detect(alias_norm, v2_influenza_h_only_pattern) ~ "influenza_subtype",
    str_detect(alias_norm, "^[A-Z0-9() -]{2,}$") & alias != disease_standard ~ "abbreviation",
    is_generic_label ~ "syndrome",
    alias_norm == standard_norm ~ "exact_label",
    TRUE ~ "synonym"
  )
}

v2_safe_disease_aliases <- function(disease_candidates) {
  seeded_aliases <- disease_candidates %>%
    mutate(
      safe_alias = v2_is_safe_disease_alias(disease_raw, disease_standard, influenza_type, influenza_subtype),
      alias_type = v2_disease_alias_type(
        disease_raw,
        disease_standard,
        is.na(influenza_subtype) & str_detect(disease_standard, regex("influenza", ignore_case = TRUE))
      ),
      is_generic_label = is.na(influenza_subtype) &
        str_detect(disease_standard, regex("influenza|syndrome|unknown|respiratory|diarrh", ignore_case = TRUE)),
      requires_subtype_evidence = alias_type == "influenza_subtype",
      allowed_without_subtype = alias_type != "influenza_subtype",
      notes = case_when(
        alias_type == "influenza_generic" ~ "Safe generic influenza alias; subtype must be resolved by explicit subtype evidence.",
        alias_type == "influenza_subtype" ~ "Safe only because the alias itself contains subtype evidence.",
        TRUE ~ "Safe native disease alias derived from accepted clean disease candidates."
      )
    ) %>%
    filter(
      safe_alias,
      !is.na(disease_raw),
      disease_raw != "",
      !is.na(disease_standard),
      disease_standard != ""
    ) %>%
    transmute(
      alias = disease_raw,
      disease_standard,
      disease_group = if_else(
        str_detect(disease_standard, regex("influenza", ignore_case = TRUE)),
        "influenza",
        NA_character_
      ),
      alias_type,
      priority = case_when(
        alias_type == "exact_label" ~ 10L,
        alias_type == "influenza_subtype" ~ 20L,
        alias_type == "synonym" ~ 30L,
        alias_type == "abbreviation" ~ 40L,
        alias_type == "influenza_generic" ~ 90L,
        TRUE ~ 100L
      ),
      is_generic_label,
      requires_subtype_evidence,
      allowed_without_subtype,
      notes
    ) %>%
    distinct() %>%
    arrange(priority, alias, disease_standard)

  manual_aliases <- tibble::tribble(
    ~alias, ~disease_standard, ~disease_group, ~alias_type, ~priority, ~is_generic_label, ~requires_subtype_evidence, ~allowed_without_subtype, ~notes,
    "MERS-CoV", "Middle East Respiratory Syndrome (MERS)", NA_character_, "abbreviation", 25L, TRUE, FALSE, TRUE, "Safe DON wording variant for MERS; keep syndrome-level and reviewable.",
    "Middle East respiratory syndrome coronavirus", "Middle East Respiratory Syndrome (MERS)", NA_character_, "synonym", 25L, TRUE, FALSE, TRUE, "Safe DON wording variant for MERS.",
    "Middle East respiratory syndrome coronavirus (MERS-CoV)", "Middle East Respiratory Syndrome (MERS)", NA_character_, "synonym", 25L, TRUE, FALSE, TRUE, "Safe DON wording variant for MERS.",
    "Novel coronavirus infection", "Middle East Respiratory Syndrome (MERS)", NA_character_, "synonym", 80L, TRUE, FALSE, TRUE, "Review-only MERS-era wording; event anchoring prevents COVID-era promotion.",
    "Hemorrhagic fever syndrome", "Haemorrhagic fever syndrome", NA_character_, "synonym", 25L, TRUE, FALSE, TRUE, "US spelling variant for broad haemorrhagic fever syndrome.",
    "Acute hemorrhagic fever syndrome", "Acute haemorrhagic fever syndrome", NA_character_, "synonym", 25L, TRUE, FALSE, TRUE, "US spelling variant for broad haemorrhagic fever syndrome.",
    "Ebola haemorrhagic fever", "Ebola virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Older DON wording for Ebola virus disease.",
    "Ebola hemorrhagic fever", "Ebola virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "US spelling variant for older DON Ebola wording.",
    "Ebola haemorragic fever", "Ebola virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON typo for older Ebola haemorrhagic fever wording.",
    "Ebola disease", "Ebola virus disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "DON title wording for Ebola virus disease, including Sudan ebolavirus-era reports.",
    "Ebola infection", "Ebola virus disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Older DON title wording for Ebola virus disease.",
    "Ebola fever", "Ebola virus disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Older DON evidence wording for Ebola virus disease.",
    "Ebola outbreak", "Ebola virus disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "DON title wording for Ebola virus disease outbreaks.",
    "Sudan virus disease", "Sudan virus disease (Ebola virus disease)", NA_character_, "synonym", 20L, FALSE, FALSE, TRUE, "Short DON title wording for Sudan virus disease.",
    "Sudan ebolavirus", "Sudan virus disease (Ebola virus disease)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen wording for Sudan virus disease.",
    "Sudan Ebola virus", "Sudan virus disease (Ebola virus disease)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Older DON pathogen wording for Sudan virus disease.",
    "Marburg haemorrhagic fever", "Marburg virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Older DON wording for Marburg virus disease.",
    "Marburg hemorrhagic fever", "Marburg virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "US spelling variant for older DON Marburg wording.",
    "Meningooccal disease", "Meningococcal disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON typo for meningococcal disease.",
    "Meningococcal meningitidis", "Meningococcal disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Older DON title wording for meningococcal disease.",
    "meningococcal meningitis", "Meningococcal disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON wording for meningococcal disease.",
    "cerebrospinal meningitis", "Meningococcal disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Older DON wording for epidemic meningococcal meningitis; generic/viral meningitis remains intentionally unmatched.",
    "epidemic meningitis", "Meningococcal disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Older DON wording for epidemic meningococcal meningitis; generic/viral meningitis remains intentionally unmatched.",
    "viral haemorrhagic fever", "Haemorrhagic fever syndrome", NA_character_, "syndrome", 35L, TRUE, FALSE, TRUE, "Broad DON syndrome wording; keep reviewable.",
    "viral hemorrhagic fever", "Haemorrhagic fever syndrome", NA_character_, "syndrome", 35L, TRUE, FALSE, TRUE, "US spelling variant for broad DON syndrome wording.",
    "haemorrhagic fever", "Haemorrhagic fever syndrome", NA_character_, "syndrome", 80L, TRUE, FALSE, TRUE, "Generic broad DON syndrome wording; keep reviewable.",
    "hemorrhagic fever", "Haemorrhagic fever syndrome", NA_character_, "syndrome", 80L, TRUE, FALSE, TRUE, "US spelling variant for generic broad DON syndrome wording.",
    "severe acute respiratory syndrome", "Severe Acute Respiratory Syndrome (SARS)", NA_character_, "synonym", 25L, TRUE, FALSE, TRUE, "Expanded SARS wording; abbreviation-only matching is avoided to prevent SARS-CoV-2 false positives.",
    "SARS", "Severe Acute Respiratory Syndrome (SARS)", NA_character_, "abbreviation", 40L, TRUE, FALSE, TRUE, "Historical SARS abbreviation; alias regex blocks SARS-CoV hyphenated matches.",
    "acute respiratory syndrome", "Acute respiratory syndrome", NA_character_, "syndrome", 35L, TRUE, FALSE, TRUE, "DON syndrome wording; promotion gate keeps this title/event anchored.",
    "respiratory illnesses", "Respiratory illness", NA_character_, "syndrome", 35L, TRUE, FALSE, TRUE, "Plural DON title wording for respiratory illness; promotion gate keeps this title/event anchored.",
    "monkeypox", "Mpox (Monkeypox)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Older DON wording for mpox.",
    "mpox", "Mpox (Monkeypox)", NA_character_, "synonym", 20L, FALSE, FALSE, TRUE, "Current DON wording for mpox.",
    "MPXV", "Mpox (Monkeypox)", NA_character_, "abbreviation", 35L, FALSE, FALSE, TRUE, "Common abbreviation for monkeypox virus in DON mpox reports.",
    "West Nile virus", "West Nile fever", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON wording for West Nile fever.",
    "West Nile virus infection", "West Nile fever", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON wording for West Nile fever.",
    "Nipah virus", "Nipah virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON wording for Nipah virus disease.",
    "Crimean-Congo haemorrhagic fever", "Crimean-Congo hemorrhagic fever", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "British spelling variant for CCHF.",
    "CCHF", "Crimean-Congo hemorrhagic fever", NA_character_, "abbreviation", 40L, FALSE, FALSE, TRUE, "Standard abbreviation for Crimean-Congo hemorrhagic fever.",
    "E.coli O157:H7", "E. coli O157 infection", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "E. coli O157:H7", "E. coli O157 infection", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "E.coli serotype O157:H7", "E. coli O157 infection", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "E. coli serotype O157:H7", "E. coli O157 infection", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "E.coli O157", "E. coli O157 infection", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "E. coli O157", "E. coli O157 infection", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Common DON pathogen wording for E. coli O157 infection.",
    "acute infective myocarditis", "Myocarditis associated with enterovirus infection", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "DON title wording for myocarditis associated with enterovirus infection.",
    "Enterovirus D68", "Enterovirus D68 respiratory disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON pathogen wording for enterovirus D68 respiratory disease.",
    "EV-D68", "Enterovirus D68 respiratory disease", NA_character_, "abbreviation", 35L, FALSE, FALSE, TRUE, "Standard abbreviation for enterovirus D68.",
    "hand, foot and mouth disease", "Hand foot and mouth disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Punctuated DON wording for hand foot and mouth disease.",
    "HFMD", "Hand foot and mouth disease", NA_character_, "abbreviation", 35L, FALSE, FALSE, TRUE, "Standard abbreviation for hand foot and mouth disease.",
    "OROV", "Oropouche fever", NA_character_, "abbreviation", 35L, FALSE, FALSE, TRUE, "Standard abbreviation for Oropouche virus/febrile disease in DON reports.",
    "Oropouche virus", "Oropouche fever", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen wording for Oropouche fever.",
    "Shiga bacillus", "Shigellosis (bacillary dysentery)", NA_character_, "synonym", 35L, FALSE, FALSE, TRUE, "DON pathogen wording for bacillary dysentery.",
    "Shigella dysenteriae", "Shigellosis (bacillary dysentery)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen wording for bacillary dysentery.",
    "Shigella sonnei", "Shigellosis (bacillary dysentery)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen wording for Shigella sonnei dysentery reports.",
    "XDR Shigella", "Shigellosis (bacillary dysentery)", NA_character_, "synonym", 35L, FALSE, FALSE, TRUE, "DON shorthand for extensively drug-resistant Shigella reports.",
    "Zika virus infection", "Zika virus disease", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Common DON wording for Zika virus disease.",
    "Zika virus", "Zika virus disease", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "Common DON wording for Zika virus disease.",
    "Acute watery diarrhoeal syndrome", "Acute watery diarrhoea", NA_character_, "syndrome", 25L, TRUE, FALSE, TRUE, "Older DON wording for acute watery diarrhoea.",
    "Diarrhoeal diseases", "Diarrhoeal disease", NA_character_, "syndrome", 30L, TRUE, FALSE, TRUE, "Older DON plural title wording for diarrhoeal disease.",
    "Salmonella Agona", "Salmonellosis", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen/serovar wording for Salmonella Agona infection reports.",
    "Salmonella infections", "Salmonellosis", NA_character_, "synonym", 30L, FALSE, FALSE, TRUE, "DON title wording for Salmonella infection reports.",
    "Klebsiella pneumoniae", "Klebsiella infection (pneumonia sepsis)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Pathogen wording for hypervirulent Klebsiella pneumoniae reports.",
    "hypervirulent Klebsiella", "Klebsiella infection (pneumonia sepsis)", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON title wording for hypervirulent Klebsiella reports.",
    "acute hepatitis of unknown aetiology", "Severe acute hepatitis of unknown aetiology", NA_character_, "syndrome", 25L, TRUE, FALSE, TRUE, "DON wording variant for severe acute hepatitis of unknown aetiology.",
    "acute hepatitis of unknown etiology", "Severe acute hepatitis of unknown aetiology", NA_character_, "syndrome", 25L, TRUE, FALSE, TRUE, "US spelling variant for severe acute hepatitis of unknown aetiology.",
    "hepatitis of unknown origin", "Severe acute hepatitis of unknown aetiology", NA_character_, "syndrome", 35L, TRUE, FALSE, TRUE, "DON wording variant for severe acute hepatitis of unknown aetiology.",
    "Polio", "Poliomyelitis", NA_character_, "synonym", 35L, FALSE, FALSE, TRUE, "Common DON title wording for poliomyelitis.",
    "vaccine-derived poliovirus", "Vaccine-derived poliovirus", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON wording for vaccine-derived poliovirus.",
    "vaccine derived poliovirus", "Vaccine-derived poliovirus", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "Hyphenless DON wording for vaccine-derived poliovirus.",
    "circulating vaccine derived poliovirus", "Vaccine-derived poliovirus", NA_character_, "synonym", 25L, FALSE, FALSE, TRUE, "DON title wording for circulating vaccine-derived poliovirus."
  )

  bind_rows(seeded_aliases, manual_aliases) %>%
    distinct(alias, disease_standard, alias_type, .keep_all = TRUE) %>%
    arrange(priority, alias, disease_standard)
}

v2_validate_disease_aliases <- function(disease_aliases) {
  required_cols <- c(
    "alias", "disease_standard", "disease_group", "alias_type", "priority",
    "is_generic_label", "requires_subtype_evidence", "allowed_without_subtype", "notes"
  )

  issues <- list()
  missing_cols <- setdiff(required_cols, names(disease_aliases))
  if (length(missing_cols) > 0) {
    issues[["missing_columns"]] <- tibble::tibble(
      severity = "blocking",
      issue_type = "missing_required_columns",
      alias = NA_character_,
      disease_standard = NA_character_,
      detail = paste(missing_cols, collapse = ", ")
    )
    return(bind_rows(issues))
  }

  add_issue <- function(.data, issue_type, detail_text, severity = "blocking") {
    .data %>%
      transmute(
        severity = severity,
        issue_type = issue_type,
        alias = as.character(alias),
        disease_standard = as.character(disease_standard),
        detail = detail_text
      )
  }

  issues[["empty_alias"]] <- disease_aliases %>%
    filter(is.na(alias) | alias == "") %>%
    add_issue("empty_alias", "Alias is missing.")

  issues[["empty_standard"]] <- disease_aliases %>%
    filter(is.na(disease_standard) | disease_standard == "") %>%
    add_issue("empty_disease_standard", "Disease standard is missing.")

  issues[["bad_priority"]] <- disease_aliases %>%
    filter(is.na(priority) | suppressWarnings(is.na(as.integer(priority)))) %>%
    add_issue("invalid_priority", "Priority must be numeric/integer.")

  issues[["bad_logical"]] <- disease_aliases %>%
    filter(
      is.na(is_generic_label) |
        is.na(requires_subtype_evidence) |
        is.na(allowed_without_subtype)
    ) %>%
    add_issue("invalid_logical_flags", "Logical flags must be non-missing.")

  issues[["duplicates"]] <- disease_aliases %>%
    count(alias, disease_standard, alias_type, name = "n") %>%
    filter(n > 1) %>%
    transmute(
      severity = "blocking",
      issue_type = "duplicate_alias_rule",
      alias = as.character(alias),
      disease_standard = as.character(disease_standard),
      detail = paste0("Duplicate rule count: ", n)
    )

  generic_multi <- disease_aliases %>%
    mutate(alias_norm = str_to_lower(str_squish(alias))) %>%
    filter(alias_type == "influenza_generic" | alias_norm == "influenza") %>%
    group_by(alias) %>%
    summarise(
      standards = paste(sort(unique(disease_standard)), collapse = " | "),
      n_standards = n_distinct(disease_standard),
      .groups = "drop"
    ) %>%
    filter(n_standards > 1) %>%
    transmute(
      severity = "blocking",
      issue_type = "generic_alias_multiple_targets",
      alias,
      disease_standard = standards,
      detail = "Generic alias maps to multiple disease standards."
    )
  issues[["generic_multi"]] <- generic_multi

  issues[["influenza_to_h5n1"]] <- disease_aliases %>%
    mutate(alias_norm = str_to_lower(str_squish(alias))) %>%
    filter(alias_norm == "influenza", disease_standard %in% c("Influenza A(H5N1)", "Influenza (H5 subtype)")) %>%
    add_issue("unsafe_influenza_generic_to_subtype", "Plain Influenza must not map directly to H5/H5N1.")

  issues[["subtype_requires_evidence"]] <- disease_aliases %>%
    filter(alias_type == "influenza_subtype", !requires_subtype_evidence) %>%
    add_issue("influenza_subtype_without_evidence_requirement", "Subtype aliases must require subtype evidence.")

  bind_rows(issues) %>%
    arrange(severity, issue_type, alias, disease_standard)
}

v2_build_disease_rule_model <- function(
  disease_aliases,
  influenza_rules = v2_read_csv(who_don_v2_rules_dir("influenza_subtype_rules.csv")),
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  alias_model <- disease_aliases %>%
    transmute(
      rule_id = paste0("alias:", row_number()),
      rule_type = "alias",
      disease_group = coalesce(
        disease_group,
        if_else(str_detect(disease_standard, regex("influenza", ignore_case = TRUE)), "influenza", "other")
      ),
      standard_label = disease_standard,
      alias,
      subtype_pattern = NA_character_,
      influenza_type = NA_character_,
      influenza_subtype = NA_character_,
      specificity_rank = case_when(
        alias_type == "influenza_subtype" ~ 1L,
        alias_type == "exact_label" ~ 2L,
        alias_type %in% c("synonym", "abbreviation") ~ 3L,
        alias_type == "influenza_generic" ~ 8L,
        is_generic_label ~ 9L,
        TRUE ~ 5L
      ),
      requires_title_or_event_anchor = is_generic_label |
        alias_type == "influenza_generic" |
        disease_standard %in% c("Acute respiratory infection", "COVID-19", "Malaria", "Respiratory illness"),
      requires_subtype_evidence,
      background_exclusion_pattern = if_else(
        is_generic_label |
          disease_standard %in% c("Acute respiratory infection", "COVID-19", "Malaria", "Respiratory illness"),
        "history of|historical|previously|differential diagnosis|advice|preparedness|surveillance|background",
        NA_character_
      ),
      priority,
      notes
    )

  subtype_model <- influenza_rules %>%
    left_join(
      influenza_standardization %>%
        transmute(
          influenza_subtype,
          standard_label = canonical_disease_standard,
          canonical_influenza_type,
          canonical_influenza_subtype
        ),
      by = "influenza_subtype"
    ) %>%
    transmute(
      rule_id,
      rule_type = "influenza_subtype",
      disease_group = "influenza",
      standard_label = if_else(
        specificity == "generic",
        "Influenza",
        coalesce(standard_label, paste0("Influenza A(", influenza_subtype, ")"))
      ),
      alias = NA_character_,
      subtype_pattern = pattern,
      influenza_type = coalesce(canonical_influenza_type, influenza_type),
      influenza_subtype = canonical_influenza_subtype,
      specificity_rank = case_when(
        specificity == "subtype" ~ 1L,
        specificity == "subtype_family" ~ 2L,
        specificity == "generic" ~ 9L,
        TRUE ~ 5L
      ),
      requires_title_or_event_anchor = specificity != "generic",
      requires_subtype_evidence = specificity != "generic",
      background_exclusion_pattern = "history of|historical|previously|differential diagnosis|surveillance|background",
      priority,
      notes
    )

  bind_rows(alias_model, subtype_model) %>%
    arrange(disease_group, specificity_rank, priority, rule_type, standard_label)
}

v2_prepare_disease_rules <- function(association_contract = v2_read_association_contract()) {
  who_don_v2_ensure_dirs()

  disease_candidates <- v2_disease_candidates_from_clean(association_contract)
  resolution_seed <- v2_disease_resolution_seed_from_clean(association_contract)
  disease_aliases <- v2_safe_disease_aliases(disease_candidates)
  validation <- v2_validate_disease_aliases(disease_aliases)
  disease_rule_model <- v2_build_disease_rule_model(disease_aliases)

  v2_write_csv(disease_aliases, who_don_v2_rules_dir("disease_aliases.csv"))
  v2_write_csv(resolution_seed, who_don_v2_rules_dir("disease_resolution_seed_from_clean.csv"))
  v2_write_csv(disease_rule_model, who_don_v2_rules_dir("disease_rule_model.csv"))
  v2_write_csv(validation, who_don_v2_output_dir("qa", "v2_disease_rule_validation.csv"))

  if (any(validation$severity == "blocking")) {
    stop(
      "Blocking disease rule validation issues found. See qa/v2_disease_rule_validation.csv",
      call. = FALSE
    )
  }

  invisible(list(
    disease_aliases = disease_aliases,
    disease_rule_model = disease_rule_model,
    resolution_seed = resolution_seed,
    validation = validation
  ))
}
