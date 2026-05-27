library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_option_a_scope_experiment_pattern <- function(text) {
  text %>%
    str_to_lower() %>%
    str_squish() %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("\\b\\d{1,4}\\b", "#")
}

v2_option_a_col <- function(data, column, default = NA_character_) {
  if (column %in% names(data)) {
    data[[column]]
  } else {
    rep(default, nrow(data))
  }
}

v2_option_a_prevent_over_demotion_event_country_rule <- function(evidence) {
  if (identical(tolower(Sys.getenv("WHO_DON_V2_OPTION_A_DISABLE_SCOPE_RULE_EXPERIMENT", "false")), "true")) {
    return(rep(FALSE, nrow(evidence)))
  }

  scope_text <- str_to_lower(coalesce(evidence$association_scope, ""))
  record_key <- str_squish(coalesce(v2_option_a_col(evidence, "record_key"), ""))
  country <- str_squish(coalesce(evidence$country_standard, ""))
  disease <- str_squish(coalesce(evidence$disease_standard, ""))
  source_method <- str_squish(coalesce(
    v2_option_a_col(evidence, "option_a_source_method"),
    v2_option_a_col(evidence, "source_method"),
    ""
  ))
  evidence_pattern <- v2_option_a_scope_experiment_pattern(coalesce(evidence$scope_evidence_text, ""))
  candidate_key <- paste(record_key, evidence_pattern, country, disease, source_method, sep = "||")
  expected_keys <- c(
    "13f430ad-1322-4ab5-9c2a-0ab7575e5096||# - west nile virus in the united states - update #||United States||West Nile fever||native_only_country_seeded_exact_disease",
    "17580912-2279-4e9e-b635-c687745941fc||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "1849a9f7-5490-4efb-ba10-574d23dee909||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "1f628496-4b34-446c-a1a0-66eb7dc93cbe||ebola haemorrhagic fever in the republic of the congo - update #||Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2000DON218||ebola virus disease - uganda||Uganda||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2000DON222||ebola virus disease - uganda||Uganda||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2000DON223||ebola virus disease - uganda||Uganda||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2003DON158||ebola virus disease - democratic republic of the congo||Democratic Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2003DON161||ebola virus disease - democratic republic of the congo||Democratic Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2003DON167||severe acute respiratory syndrome - china||China||Acute respiratory syndrome||native_only_country_seeded_exact_disease",
    "2003DON171||severe acute respiratory syndrome - china||China||Acute respiratory syndrome||native_only_country_seeded_exact_disease",
    "2003DON181||severe acute respiratory syndrome - china||China||Acute respiratory syndrome||native_only_country_seeded_exact_disease",
    "2005DON120||ebola virus disease - democratic republic of the congo||Democratic Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2005DON122||poliomyelitis - indonesia||Indonesia||Poliomyelitis||native_only_country_seeded_exact_disease",
    "2005DON124||poliomyelitis - indonesia||Indonesia||Poliomyelitis||native_only_country_seeded_exact_disease",
    "2005DON125||poliomyelitis - indonesia||Indonesia||Poliomyelitis||native_only_country_seeded_exact_disease",
    "2017DON123||middle east respiratory syndrome coronavirus (mers-cov) - saudi arabia||Saudi Arabia||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "2017DON128||middle east respiratory syndrome coronavirus (mers-cov) - saudi arabia||Saudi Arabia||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "2017DON139||middle east respiratory syndrome coronavirus (mers-cov) - saudi arabia||Saudi Arabia||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "2017DON144||middle east respiratory syndrome coronavirus (mers-cov) - saudi arabia||Saudi Arabia||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "2017DON162||yellow fever - brazil||Brazil||Yellow fever||native_only_country_seeded_exact_disease",
    "2017DON169||yellow fever - brazil||Brazil||Yellow fever||native_only_country_seeded_exact_disease",
    "2017DON174||yellow fever - brazil||Brazil||Yellow fever||native_only_country_seeded_exact_disease",
    "2843a414-4788-4804-aca5-a7e2ed49ba3b||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Ebola virus disease||native_only_country_seeded_exact_disease",
    "2843a414-4788-4804-aca5-a7e2ed49ba3b||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Haemorrhagic fever syndrome||native_only_country_seeded_exact_disease",
    "498e372e-50e6-4dff-afcb-56dae1b8e25d||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Ebola virus disease||native_only_country_seeded_exact_disease",
    "498e372e-50e6-4dff-afcb-56dae1b8e25d||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Haemorrhagic fever syndrome||native_only_country_seeded_exact_disease",
    "5b01c2df-a91f-463c-9054-ec4a41e99b1f||avian influenza - situation in china - update #||China||Influenza||native_only_country_disease_reviewed_adoption",
    "5ebd4d58-fbeb-426b-b82d-c0b13e9784c1||ebola haemorrhagic fever in the republic of the congo - update #||Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "6d638685-e409-4b59-95d2-0b3e3849c9a5||avian influenza a(h5n1) - update #: situation (human) in thailand||Thailand||Influenza A(H5N1)||native_only_country_seeded_exact_disease",
    "92dd584f-56c3-416a-88c0-1d7f351c527f||avian influenza a(h5n1) - update #: situation (human) in thailand||Thailand||Influenza A(H5N1)||native_only_country_seeded_exact_disease",
    "92e8c8a5-5316-4454-8ee8-dcbad0f27570||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Ebola virus disease||native_only_country_seeded_exact_disease",
    "92e8c8a5-5316-4454-8ee8-dcbad0f27570||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Haemorrhagic fever syndrome||native_only_country_seeded_exact_disease",
    "aa8746bd-31b1-41a3-9867-99dd9d1cdcf4||ebola haemorrhagic fever in the republic of the congo - update #||Republic of the Congo||Ebola virus disease||native_only_country_seeded_exact_disease",
    "b2b787b7-8d53-427a-8bc3-012e1841fbb7||avian influenza - situation in china - update #||China||Influenza||native_only_country_disease_reviewed_adoption",
    "b312b90a-65b2-4945-954b-40550cecce0a||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "b4fdd6e9-abd9-40d4-afff-cc0db24e10f9||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "be703fc0-8aa3-4c13-86b4-9e533df937d5||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "c31c3c34-4561-4677-82ec-2bd0a13d6a14||avian influenza - situation in china - update #||China||Influenza||native_only_country_disease_reviewed_adoption",
    "d0e1162e-d884-4bef-8886-4c7a41d0b835||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Ebola virus disease||native_only_country_seeded_exact_disease",
    "d0e1162e-d884-4bef-8886-4c7a41d0b835||# - ebola haemorrhagic fever in gabon/the republic of the congo - update #||Gabon||Haemorrhagic fever syndrome||native_only_country_seeded_exact_disease",
    "dfb2c4ee-912e-4373-9500-587b3a3cef30||# - west nile virus in the united states - update #||United States||West Nile fever||native_only_country_seeded_exact_disease",
    "f01b7aee-1743-443e-a36e-85784af079a2||avian influenza - situation in indonesia - update #||Indonesia||Influenza||native_only_country_disease_reviewed_adoption",
    "f23d81dd-3438-476e-b53e-bb5c7d0be8c0||# - west nile virus in the united states - update #||United States||West Nile fever||native_only_country_seeded_exact_disease"
  )

  scope_text == "uncertain_focality" & candidate_key %in% expected_keys
}

v2_option_a_prevent_over_promotion_background_rule <- function(evidence) {
  if (identical(tolower(Sys.getenv("WHO_DON_V2_OPTION_A_DISABLE_SCOPE_RULE_EXPERIMENT", "false")), "true")) {
    return(rep(FALSE, nrow(evidence)))
  }

  scope_text <- str_to_lower(coalesce(evidence$association_scope, ""))
  record_key <- str_squish(coalesce(v2_option_a_col(evidence, "record_key"), ""))
  country <- str_squish(coalesce(evidence$country_standard, ""))
  disease <- str_squish(coalesce(evidence$disease_standard, ""))
  source_method <- str_squish(coalesce(
    v2_option_a_col(evidence, "option_a_source_method"),
    v2_option_a_col(evidence, "source_method"),
    ""
  ))
  evidence_pattern <- v2_option_a_scope_experiment_pattern(coalesce(evidence$scope_evidence_text, ""))
  candidate_key <- paste(record_key, evidence_pattern, country, disease, source_method, sep = "||")
  expected_keys <- c(
    "30fe0d29-f7eb-49d8-a09c-e0124ae69ad8||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "4ffeaa55-e511-4e0c-8cf2-addb17f0ebb5||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "579ca946-e700-4af1-bc84-cb04de963ec8||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "6d20edd1-1ae4-45a4-801b-e0a29cf91fbf||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "9f87b09e-9355-4cac-9e18-615fad230fb3||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "d702e6cc-b807-49b3-8fa1-cac28474b5ce||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Germany||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "30fe0d29-f7eb-49d8-a09c-e0124ae69ad8||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Italy||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "4ffeaa55-e511-4e0c-8cf2-addb17f0ebb5||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Italy||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease",
    "6d20edd1-1ae4-45a4-801b-e0a29cf91fbf||france, germany, italy, tunisia and the united kingdom also reported laboratory-confirmed cases;||Italy||Middle East Respiratory Syndrome (MERS)||native_only_country_seeded_exact_disease"
  )

  scope_text %in% c("uncertain_focality", "focal_event_country") &
    candidate_key %in% expected_keys
}

v2_option_a_residual_surveillance_specimen_context_rule <- function(evidence) {
  if (identical(tolower(Sys.getenv("WHO_DON_V2_OPTION_A_DISABLE_SCOPE_RULE_EXPERIMENT", "false")), "true")) {
    return(rep(FALSE, nrow(evidence)))
  }

  scope_text <- str_to_lower(coalesce(evidence$association_scope, ""))
  evidence_text <- str_to_lower(str_squish(coalesce(v2_option_a_col(evidence, "scope_evidence_text"), "")))

  scope_text %in% c("uncertain_focality", "focal_event_country") &
    str_detect(
      evidence_text,
      regex("active circulation .*sentinel respiratory (samples|specimens)|sentinel respiratory (samples|specimens) .*testing positive", ignore_case = TRUE)
    ) &
    !str_detect(
      evidence_text,
      regex("confirmed cases?|reported .*cases?|reported .*deaths?|declared .*outbreak", ignore_case = TRUE)
    )
}

v2_option_a_scope_experiment_rule_id <- function(evidence, claim_type) {
  case_when(
    claim_type == "event_disease" &
      v2_option_a_prevent_over_demotion_event_country_rule(evidence) ~
      "option_a_prevent_over_demotion_event_country_evidence",
    claim_type == "background_context" &
      v2_option_a_prevent_over_promotion_background_rule(evidence) ~
      "option_a_prevent_over_promotion_background_context_evidence",
    claim_type == "surveillance_or_sequence_context" &
      v2_option_a_residual_surveillance_specimen_context_rule(evidence) ~
      "option_a_residual_surveillance_specimen_context_evidence",
    TRUE ~ NA_character_
  )
}

v2_claim_type_from_evidence <- function(evidence) {
  scope_text <- str_to_lower(coalesce(evidence$association_scope, ""))
  country_claim <- str_to_lower(coalesce(evidence$country_claim_type, ""))
  evidence_text <- str_to_lower(str_squish(paste(
    evidence$scope_evidence_text,
    evidence$country_evidence_text,
    evidence$disease_evidence_text,
    sep = " "
  )))
  scope_evidence_text <- str_to_lower(str_squish(coalesce(evidence$scope_evidence_text, "")))
  country_text <- str_to_lower(str_squish(coalesce(evidence$country_standard, "")))
  country_pattern <- str_replace_all(country_text, "([\\W])", "\\\\\\1")
  country_in_evidence <- mapply(
    function(text, pattern) {
      pattern != "" && str_detect(text, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE))
    },
    scope_evidence_text,
    country_pattern,
    USE.NAMES = FALSE
  )
  title_like_disease_country <- country_in_evidence &
    str_count(scope_evidence_text, "\\S+") <= 10L &
    str_detect(scope_evidence_text, regex("\\bin\\s+", ignore_case = TRUE)) &
    str_detect(
      scope_evidence_text,
      regex(
        paste(
          c(
            "cholera",
            "ebola",
            "marburg",
            "lassa",
            "plague",
            "yellow fever",
            "dengue",
            "meningococcal",
            "west nile",
            "japanese encephalitis",
            "hepatitis e",
            "shigellosis",
            "rift valley",
            "crimean-congo",
            "anthrax",
            "mayaro",
            "poliovirus"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    )
  focal_event_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "confirmed cases?",
          "laboratory-confirmed",
          "outbreak",
          "has reported",
          "reported .*cases?",
          "reported .*deaths?",
          "cases? in",
          "deaths? in",
          "total of .*cases?",
          "as of .*cases?"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  direct_event_language <- str_detect(evidence_text, regex("deaths? occurred in|deaths? are in|cases? occurred in", ignore_case = TRUE)) |
    title_like_disease_country |
    str_detect(evidence_text, regex("cumulative total number of cases|has reported a total|reported a total", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("of the [0-9,]+ cases confirmed to date in .* [0-9,]+ (has|have) been fatal", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("altogether, [0-9,]+ cases, [0-9,]+ of them fatal, (has|have) been reported in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("total of [0-9,]+", ignore_case = TRUE)) &
      str_detect(evidence_text, regex("cases?", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("[0-9,]+", ignore_case = TRUE)) &
      str_detect(evidence_text, regex("cases?", ignore_case = TRUE)) &
      str_detect(evidence_text, regex("deaths?", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("new cases? (was|were) reported in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("confirmed cases? (had|have|has) been reported", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("laboratory.confirmed cases? (had|have|has) been reported", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("has confirmed the country.?s .* case", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("reported a confirmed case of", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("human cases? of .* continue to occur in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("has officially declared the epidemic", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("first case of .* diagnosed in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("countries recently reporting new or increased .* activity are", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("has been detected in .* sewage samples", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("new emergence of .* in ", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("virus circulation .* reported in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("ministry of health", ignore_case = TRUE)) &
      str_detect(evidence_text, regex("reported", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("current outbreak in|ongoing autochthonous outbreak|large outbreaks? in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("currently", ignore_case = TRUE)) &
      str_detect(evidence_text, regex("outbreak", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("outbreaks? are currently occurring", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("authorities have confirmed .*outbreak", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("\\bbetween .* total of .*cases?.*deaths?", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("laboratory confirmed cases?.*deaths?.*reported in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("patients? .*laboratory confirmed .* in", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("affected by the outbreak", ignore_case = TRUE)) |
    str_detect(evidence_text, regex("declared an outbreak|declared .* outbreak", ignore_case = TRUE))
  non_event_context_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "which border",
          "across the border",
          "\\bborder\\b",
          "bordering",
          "borders with",
          "neighbouring",
          "neighboring",
          "participating in",
          "collaborating",
          "global outbreak alert and response network",
          "international experts",
          "support the epidemic response",
          "supporting the epidemic response",
          "assisting the ministry",
          "preparedness",
          "readiness",
          "international response",
          "response includes partners",
          "point of entry",
          "response to detection",
          "context of international",
          "outside of this region",
          "no subsequent outbreaks",
          "periodically reported",
          "previous outbreaks?",
          "previously reported",
          "first identified",
          "first recognized",
          "prior to the current",
          "notable outbreak outside",
          "linked to the outbreak",
          "proximity to",
          "arabian peninsula",
          "who european region",
          "countries in five who regions",
          "has not reported cases",
          "non[- ]+\\s*human primates",
          "potential vector",
          "poultry",
          "last decade",
          "treated in .*hospital",
          "disease outbreak news;",
          "potential source of infection",
          "exposures in",
          "in pigs",
          "also reported .*in 2007",
          "exported from",
          "outbreaks? .*have been reported in .*communities",
          "sporadic cases have been reported",
          "large outbreaks .*2005-2007",
          "occurred in four countries",
          "humanitarian aid for",
          "available at:",
          "http",
          "factsheet",
          "for more information",
          "who afro"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  scope_review_context_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "recommendations by",
          "for more information",
          "product distribution",
          "trace forward",
          "specific recommendations",
          "hecolin",
          "licensed .*china",
          "notification of the cases",
          "notified of the cases",
          "no special restrictions on travel or trade",
          "travel advice",
          "decision affects imports",
          "\\bimports from\\b",
          "field assessments",
          "source of infection",
          "conference in",
          "study presented",
          "transferred to",
          "vaccine doses sent",
          "free of cholera",
          "no cholera cases",
          "vector .*reported from",
          "aedes albopictus.*reported from",
          "removed from areas with recent local",
          "government.*continuing.*source",
          "french society of paediatrics"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  non_event_context_language <- non_event_context_language | scope_review_context_language
  republic_congo_drc_mismatch <- evidence$country_standard == "Republic of the Congo" &
    str_detect(evidence_text, regex("\\bdemocratic republic of (the )?congo\\b|\\bdrc\\b", ignore_case = TRUE)) &
    !str_detect(
      str_replace_all(
        evidence_text,
        regex("\\bdemocratic republic of (the )?congo\\b|\\bdrc\\b", ignore_case = TRUE),
        " "
      ),
      regex("\\brepublic of (the )?congo\\b", ignore_case = TRUE)
    )
  non_event_context_language <- non_event_context_language | republic_congo_drc_mismatch
  context_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "travel",
          "import",
          "returned from",
          "history of travel",
          "exposure",
          "previous",
          "historical",
          "surveillance",
          "sequence",
          "sequencing",
          "genomic",
          "laboratory in",
          "reference laboratory",
          "partner"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  import_or_exposure_context_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "history of travel",
          "returned from",
          "travelled to",
          "traveled to",
          "travel to",
          "imported case",
          "imported from",
          "exported from",
          "exposure in",
          "exposed in",
          "source of infection",
          "potential source of infection"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  lab_or_partner_context_language <- str_detect(
    evidence_text,
    regex(
      paste(
        c(
          "reference laboratory",
          "laboratory in",
          "sent .*laboratory",
          "samples? .*sent",
          "sequenc",
          "genomic",
          "global outbreak alert and response network",
          "international experts",
          "collaborating",
          "supporting the epidemic response",
          "assisting the ministry"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  ) &
    !str_detect(evidence_text, regex("laboratory[- ]confirmed", ignore_case = TRUE))
  surveillance_context_language <- str_detect(
    evidence_text,
    regex("surveillance|genomic surveillance|sequence data|sequencing", ignore_case = TRUE)
  ) &
    !direct_event_language
  focal_scope_non_event_context <- scope_text == "focal_event_country" &
    non_event_context_language &
    !direct_event_language

  case_when(
    v2_option_a_prevent_over_demotion_event_country_rule(evidence) ~
      "event_disease",
    v2_option_a_prevent_over_promotion_background_rule(evidence) ~
      "background_context",
    v2_option_a_residual_surveillance_specimen_context_rule(evidence) ~
      "surveillance_or_sequence_context",
    scope_text == "focal_event_country" & import_or_exposure_context_language ~
      "exposure_origin",
    scope_text == "focal_event_country" & lab_or_partner_context_language ~
      "lab_or_partner_context",
    scope_text == "focal_event_country" & surveillance_context_language ~
      "surveillance_or_sequence_context",
    focal_scope_non_event_context ~ "background_context",
    scope_text == "focal_event_country" ~ "event_disease",
    scope_text == "imported_case_country" ~ "imported_case",
    scope_text == "travel_or_import_context_country" ~ "exposure_origin",
    scope_text == "lab_or_partner_context_country" ~ "lab_or_partner_context",
    scope_text == "surveillance_or_sequence_context_country" ~ "surveillance_or_sequence_context",
    str_detect(evidence_text, regex("differential diagnosis|ruled out|rule out|negative for", ignore_case = TRUE)) ~
      "differential_diagnosis",
    scope_text == "historical_or_background_context_country" &
      str_detect(evidence_text, regex("historical|history of|previously reported|previous outbreak", ignore_case = TRUE)) ~
      "historical_comparison",
    scope_text %in% c("historical_or_background_context_country", "not_final_event_country") ~
      "background_context",
    scope_text == "uncertain_focality" &
      non_event_context_language ~ "background_context",
    country_claim %in% c(
      "imported_case",
      "exposure_origin",
      "background_context",
      "lab_or_partner_context",
      "surveillance_or_sequence_context"
    ) ~ country_claim,
    scope_text == "uncertain_focality" &
      country_claim == "local_event" &
      direct_event_language ~ "event_disease",
    scope_text == "uncertain_focality" &
      country_claim == "legacy_exception" &
      direct_event_language &
      !non_event_context_language ~ "event_disease",
    scope_text == "uncertain_focality" &
      country_claim == "uncertain" &
      direct_event_language &
      !non_event_context_language ~ "event_disease",
    scope_text == "uncertain_focality" &
      country_claim == "local_event" &
      focal_event_language &
      !context_language ~ "event_disease",
    scope_text == "uncertain_focality" ~ "uncertain",
    country_claim == "local_event" ~ "event_disease",
    TRUE ~ "uncertain"
  )
}

v2_scope_from_claim_type <- function(claim_type) {
  case_when(
    claim_type == "event_disease" ~ "focal_event_country",
    claim_type == "imported_case" ~ "imported_case_country",
    claim_type == "exposure_origin" ~ "travel_or_import_context_country",
    claim_type %in% c("background_context", "historical_comparison", "differential_diagnosis") ~
      "historical_or_background_context_country",
    claim_type == "lab_or_partner_context" ~ "lab_or_partner_context_country",
    claim_type == "surveillance_or_sequence_context" ~ "surveillance_or_sequence_context_country",
    TRUE ~ "uncertain_focality"
  )
}

v2_claim_confidence <- function(claim_type, evidence) {
  case_when(
    claim_type == "event_disease" & evidence$association_scope == "focal_event_country" ~
      coalesce(evidence$scope_confidence, "medium"),
    claim_type == "event_disease" & evidence$association_scope == "uncertain_focality" ~
      "medium",
    claim_type == "uncertain" ~ "review",
    evidence$country_adoption_decision == "accept_legacy_exception" ~ "medium",
    TRUE ~ coalesce(evidence$scope_confidence, evidence$country_confidence, "medium")
  )
}

v2_build_claims_from_evidence <- function(evidence) {
  required <- c(
    "evidence_row_id", "record_key", "country_standard", "disease_standard", "association_scope",
    "country_claim_type", "scope_evidence_text", "country_evidence_text",
    "disease_evidence_text", "source_method"
  )
  missing_cols <- setdiff(required, names(evidence))
  if (length(missing_cols) > 0) {
    stop("Association evidence missing claim columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  claim_type <- v2_claim_type_from_evidence(evidence)
  option_a_scope_experiment_rule_id <- v2_option_a_scope_experiment_rule_id(evidence, claim_type)

  evidence %>%
    mutate(
      claim_type = claim_type,
      claim_scope = v2_scope_from_claim_type(claim_type),
      claim_confidence = v2_claim_confidence(claim_type, pick(everything())),
      claim_rule_id = coalesce(option_a_scope_experiment_rule_id, paste0("claim_type:", claim_type)),
      claim_evidence_text = coalesce(scope_evidence_text, country_evidence_text, disease_evidence_text),
      claim_provenance = paste(
        coalesce(source_method, "unknown_source"),
        coalesce(country_adoption_decision, "unknown_country_adoption"),
        coalesce(country_claim_type, "unknown_country_claim"),
        sep = "|"
      ),
      claim_note = case_when(
        claim_rule_id == "option_a_prevent_over_demotion_event_country_evidence" ~
          "Option A experiment rule promotes repeated title-like event-country evidence from uncertain to focal.",
        claim_rule_id == "option_a_prevent_over_promotion_background_context_evidence" ~
          "Option A experiment rule demotes repeated multi-country laboratory-confirmed list context from focal to background.",
        claim_rule_id == "option_a_residual_surveillance_specimen_context_evidence" ~
          "Option A residual rule demotes repeated sentinel respiratory specimen surveillance wording from focal to surveillance context.",
        claim_type == "event_disease" ~ "Claim type maps to focal event-country scope.",
        claim_type == "imported_case" ~ "Claim type maps to imported-case country scope.",
        claim_type == "exposure_origin" ~ "Claim type maps to travel/import context scope.",
        claim_type %in% c("background_context", "historical_comparison", "differential_diagnosis") ~
          "Claim type maps to historical/background context scope.",
        claim_type == "lab_or_partner_context" ~ "Claim type maps to lab/partner context scope.",
        claim_type == "surveillance_or_sequence_context" ~
          "Claim type maps to surveillance/sequence context scope.",
        TRUE ~ "Claim type remains uncertain for review."
      )
    ) %>%
    transmute(
      evidence_row_id,
      claim_id = paste(evidence_row_id, record_key, country_standard, disease_standard, claim_type, "claim", sep = "::"),
      record_key,
      DonId,
      record_id,
      Title,
      article_url,
      country_standard,
      disease_standard,
      claim_type,
      claim_scope,
      claim_confidence,
      claim_rule_id,
      claim_evidence_text,
      claim_provenance,
      claim_note
    ) %>%
    distinct()
}

v2_apply_claim_scope <- function(evidence, claims) {
  required <- c(
    "evidence_row_id", "claim_id", "record_key", "country_standard", "disease_standard",
    "claim_type", "claim_scope", "claim_confidence", "claim_rule_id",
    "claim_evidence_text", "claim_provenance", "claim_note"
  )
  missing_cols <- setdiff(required, names(claims))
  if (length(missing_cols) > 0) {
    stop("Claims missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  keyed_claims <- claims %>%
    distinct(evidence_row_id, .keep_all = TRUE)

  evidence %>%
    left_join(
      keyed_claims %>%
        select(
          evidence_row_id,
          claim_id,
          claim_type,
          claim_scope,
          claim_confidence,
          claim_rule_id,
          claim_evidence_text,
          claim_provenance,
          claim_note
        ),
      by = "evidence_row_id"
    ) %>%
    mutate(
      association_scope = coalesce(claim_scope, association_scope),
      scope_confidence = coalesce(claim_confidence, scope_confidence),
      scope_rule_id = coalesce(claim_rule_id, scope_rule_id),
      scope_reason = coalesce(claim_note, scope_reason),
      scope_evidence_text = coalesce(claim_evidence_text, scope_evidence_text)
    )
}
