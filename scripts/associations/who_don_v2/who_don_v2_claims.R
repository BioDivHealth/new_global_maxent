library(dplyr)
library(stringr)

source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_io.R"))

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

  case_when(
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

  evidence %>%
    mutate(
      claim_type = claim_type,
      claim_scope = v2_scope_from_claim_type(claim_type),
      claim_confidence = v2_claim_confidence(claim_type, pick(everything())),
      claim_rule_id = paste0("claim_type:", claim_type),
      claim_evidence_text = coalesce(scope_evidence_text, country_evidence_text, disease_evidence_text),
      claim_provenance = paste(
        coalesce(source_method, "unknown_source"),
        coalesce(country_adoption_decision, "unknown_country_adoption"),
        coalesce(country_claim_type, "unknown_country_claim"),
        sep = "|"
      ),
      claim_note = case_when(
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
