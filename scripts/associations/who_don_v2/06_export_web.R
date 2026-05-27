library(dplyr)
library(stringr)
library(jsonlite)
library(countrycode)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

who_don_v2_ensure_dirs()

audit_path <- who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv")
modelling_path <- who_don_v2_output_dir("final", "who_don_modelling_ready.csv")
web_data_path <- who_don_v2_output_dir("web", "who_don_web.json")
web_meta_path <- who_don_v2_output_dir("web", "who_don_meta.json")

generic_disease_labels <- c(
  "Influenza",
  "Acute respiratory syndrome",
  "Acute respiratory infection",
  "Haemorrhagic fever syndrome",
  "Acute haemorrhagic fever syndrome",
  "Vaccine-derived poliovirus",
  "Acute febrile illness",
  "Acute fever and rash syndrome",
  "Acute watery diarrhoea",
  "Bloody diarrhoea",
  "Cluster of unexplained community deaths",
  "Acute neurological syndrome",
  "Acute encephalitis syndrome"
)

iso3_fallback <- c(
  "Kosovo" = "XKX",
  "Micronesia" = "FSM"
)

blank_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

coalesce_chr <- function(...) {
  dplyr::coalesce(!!!lapply(list(...), blank_to_na))
}

normalise_bool <- function(x, default = FALSE) {
  out <- case_when(
    is.na(x) ~ default,
    is.logical(x) ~ x,
    str_to_lower(as.character(x)) %in% c("true", "t", "1", "yes") ~ TRUE,
    str_to_lower(as.character(x)) %in% c("false", "f", "0", "no") ~ FALSE,
    TRUE ~ default
  )
  as.logical(out)
}

normalise_scope_confidence <- function(x) {
  case_when(
    x %in% c("high", "medium", "manual_override_high", "manual_override_keep_non_strict") ~ x,
    TRUE ~ "not_applicable"
  )
}

derive_event_layer_confidence <- function(scope_confidence) {
  case_when(
    scope_confidence %in% c("high", "manual_override_high") ~ "high",
    scope_confidence == "medium" ~ "medium",
    TRUE ~ "low"
  )
}

simplify_source <- function(x) {
  case_when(
    str_detect(x, regex("openai|llm", ignore_case = TRUE)) ~ "llm",
    TRUE ~ "deterministic"
  )
}

derive_disease_group <- function(disease_label) {
  case_when(
    str_detect(disease_label, "^Influenza(\\b|\\s|$)") ~ "Influenza",
    disease_label %in% c(
      "Poliomyelitis",
      "Wild poliovirus",
      "Vaccine-derived poliovirus",
      "Circulating vaccine-derived poliovirus type 1",
      "Circulating vaccine-derived poliovirus type 2"
    ) ~ "Poliovirus / poliomyelitis",
    disease_label %in% c(
      "Ebola virus disease",
      "Sudan virus disease (Ebola virus disease)"
    ) ~ "Ebola virus disease",
    disease_label %in% c(
      "Hantavirus disease",
      "Hantavirus pulmonary syndrome",
      "Seoul virus disease"
    ) ~ "Hantavirus",
    disease_label %in% c(
      "Enterovirus infection",
      "Enterovirus-Echovirus 11 infection",
      "Myocarditis associated with enterovirus infection",
      "Hand foot and mouth disease"
    ) ~ "Enterovirus",
    disease_label %in% c(
      "Hepatitis A",
      "Hepatitis E",
      "Severe acute hepatitis of unknown aetiology"
    ) ~ "Hepatitis",
    disease_label %in% c(
      "Undiagnosed illness",
      "Undiagnosed disease",
      "Undiagnosed febrile illness",
      "Unknown illness",
      "Pneumonia of unknown cause",
      "Cluster of unexplained community deaths"
    ) ~ "Undiagnosed / unknown illness",
    TRUE ~ NA_character_
  )
}

derive_disease_group_display <- function(disease_group) {
  case_when(
    disease_group == "Influenza" ~ "Influenza [all]",
    disease_group == "Poliovirus / poliomyelitis" ~ "Poliovirus / poliomyelitis [all]",
    disease_group == "Ebola virus disease" ~ "Ebola virus disease [all]",
    disease_group == "Hantavirus" ~ "Hantavirus [all]",
    disease_group == "Enterovirus" ~ "Enterovirus [all]",
    disease_group == "Hepatitis" ~ "Hepatitis [all]",
    disease_group == "Undiagnosed / unknown illness" ~ "Undiagnosed / unknown illness [all]",
    TRUE ~ NA_character_
  )
}

derive_scope_group <- function(scope) {
  case_when(
    scope == "focal_event_country" ~ "focal",
    scope %in% c(
      "imported_case_country",
      "secondary_local_transmission_country",
      "travel_or_import_context_country"
    ) ~ "imported_or_secondary",
    scope %in% c(
      "historical_or_background_context_country",
      "lab_or_partner_context_country",
      "surveillance_or_sequence_context_country"
    ) ~ "background_or_context",
    scope == "not_final_event_country" ~ "not_event_country",
    scope == "uncertain_focality" ~ "review",
    TRUE ~ "other"
  )
}

derive_scope_group_display <- function(scope_group) {
  case_when(
    scope_group == "focal" ~ "Focal outbreak/reporting country",
    scope_group == "imported_or_secondary" ~ "Imported/secondary transmission country",
    scope_group == "background_or_context" ~ "Background/context mention",
    scope_group == "not_event_country" ~ "Not event-country layer",
    scope_group == "review" ~ "Needs review",
    TRUE ~ "Other"
  )
}

derive_scope_display <- function(scope) {
  case_when(
    scope == "focal_event_country" ~ "Focal event country",
    scope == "imported_case_country" ~ "Imported case country",
    scope == "secondary_local_transmission_country" ~ "Secondary local transmission country",
    scope == "travel_or_import_context_country" ~ "Travel/import context country",
    scope == "historical_or_background_context_country" ~ "Historical/background context country",
    scope == "lab_or_partner_context_country" ~ "Lab/partner context country",
    scope == "surveillance_or_sequence_context_country" ~ "Surveillance/sequence context country",
    scope == "not_final_event_country" ~ "Not in final event-country layer",
    scope == "uncertain_focality" ~ "Uncertain focality",
    TRUE ~ scope
  )
}

normalise_influenza_type_display <- function(influenza_type) {
  case_when(
    is.na(influenza_type) | influenza_type == "" ~ "A",
    str_to_lower(influenza_type) == "influenza" ~ "A",
    TRUE ~ influenza_type
  )
}

canonical_disease_display <- function(disease_label, influenza_type, influenza_subtype) {
  influenza_type_display <- normalise_influenza_type_display(influenza_type)

  case_when(
    !is.na(influenza_type) & !is.na(influenza_subtype) &
      str_detect(disease_label, "^Influenza") &
      str_detect(influenza_subtype, "^H[0-9]+N[0-9]+$") ~
      paste0("Influenza ", influenza_type_display, "(", influenza_subtype, ")"),
    !is.na(influenza_type) & !is.na(influenza_subtype) &
      str_detect(disease_label, "^Influenza") &
      str_detect(influenza_subtype, "^H[0-9]+$") ~
      paste0("Influenza ", influenza_type_display, "(", influenza_subtype, " subtype)"),
    str_detect(disease_label, "^Influenza \\(H[0-9]+N[0-9]+\\)$") ~
      str_replace(disease_label, "^Influenza \\((H[0-9]+N[0-9]+)\\)$", "Influenza A(\\1)"),
    str_detect(disease_label, "^Influenza \\(H[0-9]+ subtype\\)$") ~
      str_replace(disease_label, "^Influenza \\((H[0-9]+ subtype)\\)$", "Influenza A(\\1)"),
    TRUE ~ disease_label
  )
}

audit <- v2_read_csv(
  audit_path,
  c(
    "record_key",
    "country_standard",
    "disease_standard",
    "association_scope",
    "scope_evidence_text",
    "country_evidence_text",
    "disease_evidence_text",
    "source_method",
    "scope_rule_id"
  )
)

modelling <- v2_read_csv(
  modelling_path,
  c("record_key", "country_standard", "disease_label_standard", "association_scope")
)

strict_keys <- modelling %>%
  transmute(
    strict_key = paste(record_key, country_standard, disease_label_standard, sep = "||")
  ) %>%
  distinct()

web_rows <- audit %>%
  mutate(
    strict_key = paste(record_key, country_standard, disease_standard, sep = "||"),
    include_in_strict_modelling = strict_key %in% strict_keys$strict_key,
    date = as.character(as.Date(str_sub(publication_datetime_utc, 1, 10))),
    don_id = coalesce_chr(DonId, record_id, record_key),
    iso3 = countrycode::countrycode(
      country_standard,
      origin = "country.name",
      destination = "iso3c",
      warn = FALSE
    ),
    iso3 = coalesce_chr(iso3, unname(iso3_fallback[country_standard])),
    disease_label_web = canonical_disease_display(
      disease_standard,
      blank_to_na(influenza_type),
      blank_to_na(influenza_subtype)
    ),
    disease_group = derive_disease_group(disease_label_web),
    disease_group_display = derive_disease_group_display(disease_group),
    scope_group = derive_scope_group(association_scope),
    scope_group_display = derive_scope_group_display(scope_group),
    scope_display = derive_scope_display(association_scope),
    scope_confidence_web = normalise_scope_confidence(scope_confidence),
    event_layer_confidence = derive_event_layer_confidence(scope_confidence_web),
    event_layer_needs_review = normalise_bool(needs_review),
    final_don_ready_for_downstream = include_in_strict_modelling,
    evidence = coalesce_chr(
      scope_evidence_text,
      claim_evidence_text,
      country_evidence_text,
      disease_evidence_text
    ),
    evidence_sentence = coalesce_chr(scope_evidence_text, claim_evidence_text),
    trigger_phrase = coalesce_chr(claim_rule_id, scope_rule_id),
    source_simple = simplify_source(source_method),
    was_adjudicated = str_detect(
      coalesce_chr(source_method, final_review_source, review_status, ""),
      regex("review|manual|adjudicat|override|openai|llm", ignore_case = TRUE)
    )
  ) %>%
  transmute(
    record_key,
    don_id,
    record_id,
    title = Title,
    date,
    article_url,
    country = country_standard,
    country_original = country_standard,
    iso3,
    disease = disease_label_web,
    disease_display = disease_label_web,
    disease_group,
    disease_group_display,
    disease_standard,
    disease_refined = disease_standard,
    disease_downstream = disease_standard,
    disease_raw,
    influenza_type = blank_to_na(influenza_type),
    influenza_subtype = blank_to_na(influenza_subtype),
    is_generic_label = disease_label_web %in% generic_disease_labels,
    country_role = case_when(
      include_in_strict_modelling ~ "event_country",
      scope_group == "review" ~ "needs_review",
      TRUE ~ "context_country"
    ),
    event_country = include_in_strict_modelling,
    include_in_strict_modelling,
    probable_focal_event_country = include_in_strict_modelling,
    scope = association_scope,
    don_country_report_scope = association_scope,
    scope_display,
    scope_group,
    scope_group_display,
    is_context_only = scope_group %in% c("background_or_context", "not_event_country", "imported_or_secondary"),
    confidence = scope_confidence_web,
    scope_confidence = scope_confidence_web,
    event_layer_confidence,
    event_layer_reasoning = scope_rule_id,
    event_layer_needs_review,
    evidence,
    reasoning = scope_reason,
    event_reasoning = scope_rule_id,
    source = source_simple,
    needs_review = event_layer_needs_review,
    legacy_event_review_flag = event_layer_needs_review,
    final_don_ready_for_downstream,
    disease_refinement_review_status = "resolved_or_not_needed",
    evidence_sentence,
    trigger_phrase,
    exclusion_phrase = NA_character_,
    focal_scope_needs_review = event_layer_needs_review,
    manual_focal_review_decision = review_decision_id,
    manual_focal_review_class = review_status,
    manual_focal_review_note = final_review_note,
    manual_focal_review_override_applied = !is.na(review_decision_id) & review_decision_id != "",
    disease_refinement_source = source_method,
    disease_refinement_needs_review = FALSE,
    influenza_subtype_candidates = blank_to_na(influenza_subtype),
    influenza_subtype_evidence_span = if_else(
      !is.na(blank_to_na(influenza_subtype)),
      blank_to_na(disease_evidence_text),
      NA_character_
    ),
    was_adjudicated,
    adj_changed_role = FALSE,
    adj_changed_flag = FALSE,
    adj_changed_conf = FALSE
  ) %>%
  arrange(date, record_key, country, disease)

if (any(is.na(web_rows$record_key) | web_rows$record_key == "")) {
  stop("Web export has blank record_key values.", call. = FALSE)
}
if (any(is.na(web_rows$country) | web_rows$country == "")) {
  stop("Web export has blank country values.", call. = FALSE)
}
if (any(is.na(web_rows$disease) | web_rows$disease == "")) {
  stop("Web export has blank disease values.", call. = FALSE)
}

strict_n <- sum(web_rows$include_in_strict_modelling %in% TRUE, na.rm = TRUE)
if (strict_n != nrow(modelling)) {
  stop(
    "Strict web-row count does not match modelling rows: ",
    strict_n,
    " vs ",
    nrow(modelling),
    call. = FALSE
  )
}

count_rows <- function(x, ...) {
  x %>%
    count(..., name = "n") %>%
    arrange(desc(n), ...)
}

disease_counts <- count_rows(web_rows, disease)
strict_disease_counts <- web_rows %>%
  filter(include_in_strict_modelling) %>%
  count(disease, name = "n") %>%
  arrange(desc(n), disease)

disease_group_counts <- web_rows %>%
  filter(!is.na(disease_group)) %>%
  count(disease_group, disease_group_display, name = "n") %>%
  arrange(desc(n), disease_group)

strict_disease_group_counts <- web_rows %>%
  filter(include_in_strict_modelling, !is.na(disease_group)) %>%
  count(disease_group, disease_group_display, name = "n") %>%
  arrange(desc(n), disease_group)

country_counts <- web_rows %>%
  distinct(record_key, country, iso3) %>%
  count(country, iso3, name = "n") %>%
  arrange(desc(n), country)

strict_country_counts <- web_rows %>%
  filter(include_in_strict_modelling) %>%
  distinct(record_key, country, iso3) %>%
  count(country, iso3, name = "n") %>%
  arrange(desc(n), country)

scope_counts <- web_rows %>%
  count(don_country_report_scope, scope_display, scope_group, scope_group_display, name = "n") %>%
  arrange(scope_group, desc(n), scope_display)

scope_group_counts <- web_rows %>%
  count(scope_group, scope_group_display, name = "n") %>%
  arrange(desc(n), scope_group)

valid_dates <- blank_to_na(web_rows$date)
meta <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  n_rows_total = nrow(web_rows),
  n_rows_strict = strict_n,
  n_rows_event = sum(web_rows$event_country %in% TRUE, na.rm = TRUE),
  n_don_records = n_distinct(web_rows$don_id),
  n_don_records_strict = n_distinct(web_rows$don_id[web_rows$include_in_strict_modelling %in% TRUE]),
  n_diseases = n_distinct(web_rows$disease),
  n_diseases_strict = n_distinct(web_rows$disease[web_rows$include_in_strict_modelling %in% TRUE]),
  n_countries = n_distinct(web_rows$country),
  n_countries_strict = n_distinct(web_rows$country[web_rows$include_in_strict_modelling %in% TRUE]),
  date_min = min(valid_dates, na.rm = TRUE),
  date_max = max(valid_dates, na.rm = TRUE),
  data_source = basename(audit_path),
  audit_source = basename(audit_path),
  strict_analysis_source = basename(modelling_path),
  generic_disease_labels = generic_disease_labels,
  diseases = disease_counts,
  strict_diseases = strict_disease_counts,
  disease_groups = disease_group_counts,
  strict_disease_groups = strict_disease_group_counts,
  scope_counts = scope_counts,
  scope_groups = scope_group_counts,
  countries = country_counts,
  strict_countries = strict_country_counts
)

dir.create(dirname(web_data_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(web_rows, web_data_path, na = "null", auto_unbox = TRUE)
jsonlite::write_json(meta, web_meta_path, na = "null", auto_unbox = TRUE, pretty = TRUE)

if (identical(Sys.getenv("WHO_DON_V2_COPY_WEB_TO_APP"), "1")) {
  app_data_dir <- Sys.getenv("WHO_DON_APP_DATA_DIR", unset = "")
  if (!nzchar(app_data_dir)) {
    stop("WHO_DON_APP_DATA_DIR must be set when WHO_DON_V2_COPY_WEB_TO_APP=1.", call. = FALSE)
  }
  if (!dir.exists(app_data_dir)) {
    stop("WHO_DON_APP_DATA_DIR does not exist: ", app_data_dir, call. = FALSE)
  }
  file.copy(web_data_path, file.path(app_data_dir, "who_don_web.json"), overwrite = TRUE)
  file.copy(web_meta_path, file.path(app_data_dir, "who_don_meta.json"), overwrite = TRUE)
}

message(
  "Wrote v2 web exports: ",
  nrow(web_rows),
  " web rows, ",
  strict_n,
  " strict rows"
)
