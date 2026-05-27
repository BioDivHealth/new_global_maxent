library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_extraction.R"))

v2_country_text_sections <- function(records) {
  section_cols <- c(
    title = "Title",
    summary = "summary_text",
    overview = "overview_text",
    epidemiology = "epidemiology_text",
    assessment = "assessment_text",
    response = "response_text",
    advice = "advice_text",
    further_information = "further_information_text"
  )
  section_cols <- section_cols[section_cols %in% names(records)]

  records %>%
    transmute(
      record_key,
      DonId,
      record_id,
      record_title = Title,
      publication_datetime_utc,
      article_url,
      !!!rlang::syms(unname(section_cols))
    ) %>%
    pivot_longer(
      cols = all_of(unname(section_cols)),
      names_to = "source_column",
      values_to = "section_text"
    ) %>%
    mutate(
      Title = record_title,
      country_evidence_location = names(section_cols)[match(source_column, unname(section_cols))],
      section_text = str_squish(as.character(section_text))
    ) %>%
    select(-record_title) %>%
    filter(!is.na(section_text), section_text != "")
}

v2_country_claim_type <- function(location, title, local_text, alias_is_ambiguous) {
  text <- str_to_lower(str_squish(paste(title, local_text, sep = " ")))
  title_lower <- str_to_lower(str_squish(title))

  imported_pattern <- paste(
    c("imported", "travelled from", "traveled from", "travel from", "returned from", "came from"),
    collapse = "|"
  )
  exposure_pattern <- paste(
    c("exposure in", "exposed in", "history of travel", "travel history", "visited", "originated in"),
    collapse = "|"
  )
  background_pattern <- paste(
    c("neighbouring countries", "neighboring countries", "other countries", "endemic in", "previously reported",
      "historical", "history of", "globally", "worldwide", "including countries", "countries in the region"),
    collapse = "|"
  )
  lab_pattern <- paste(c("reference laboratory", "laboratory in", "partner", "collaborating centre", "collaborating center"), collapse = "|")
  surveillance_pattern <- paste(c("surveillance", "sequence", "sequencing", "genomic", "sample sent"), collapse = "|")

  if (alias_is_ambiguous && location != "title") {
    return(tibble::tibble(
      country_claim_type = "uncertain",
      country_claim_rank = 9L,
      country_claim_reason = "Ambiguous alias outside title requires review."
    ))
  }
  if (location == "title" && str_detect(title_lower, "\\s[-–:]\\s|[-–]\\s*")) {
    return(tibble::tibble(
      country_claim_type = "local_event",
      country_claim_rank = 1L,
      country_claim_reason = "Country alias appears in title."
    ))
  }
  if (str_detect(text, regex(imported_pattern, ignore_case = TRUE))) {
    return(tibble::tibble(
      country_claim_type = "imported_case",
      country_claim_rank = 4L,
      country_claim_reason = "Country appears with imported or travel wording."
    ))
  }
  if (str_detect(text, regex(exposure_pattern, ignore_case = TRUE))) {
    return(tibble::tibble(
      country_claim_type = "exposure_origin",
      country_claim_rank = 5L,
      country_claim_reason = "Country appears with exposure-origin wording."
    ))
  }
  if (str_detect(text, regex(lab_pattern, ignore_case = TRUE))) {
    return(tibble::tibble(
      country_claim_type = "lab_or_partner_context",
      country_claim_rank = 7L,
      country_claim_reason = "Country appears in lab or partner context."
    ))
  }
  if (str_detect(text, regex(surveillance_pattern, ignore_case = TRUE))) {
    return(tibble::tibble(
      country_claim_type = "surveillance_or_sequence_context",
      country_claim_rank = 7L,
      country_claim_reason = "Country appears in surveillance or sequence context."
    ))
  }
  if (str_detect(text, regex(background_pattern, ignore_case = TRUE))) {
    return(tibble::tibble(
      country_claim_type = "background_context",
      country_claim_rank = 8L,
      country_claim_reason = "Country appears in background, regional, historical, or global context."
    ))
  }
  if (location %in% c("summary", "overview", "epidemiology")) {
    return(tibble::tibble(
      country_claim_type = "local_event",
      country_claim_rank = 2L,
      country_claim_reason = "Country appears in event summary, overview, or epidemiology text."
    ))
  }
  tibble::tibble(
    country_claim_type = "uncertain",
    country_claim_rank = 9L,
    country_claim_reason = "Country mention is outside strongest event sections."
  )
}

v2_country_false_positive_reason <- function(country_standard, country_raw, local_text) {
  text <- str_squish(coalesce(local_text, ""))
  text_without_compound_country <- function(compound_pattern) {
    str_squish(str_replace_all(text, regex(compound_pattern, ignore_case = TRUE), " "))
  }

  case_when(
    country_standard == "Guinea" &
      str_detect(text, regex("\\bPapua New Guinea\\b", ignore_case = TRUE)) &
      !str_detect(text_without_compound_country("\\bPapua New Guinea\\b"), regex("\\bGuinea\\b", ignore_case = TRUE)) ~
      "Guinea matched only inside Papua New Guinea.",
    country_standard == "Ireland" &
      str_detect(text, regex("\\bNorthern Ireland\\b", ignore_case = TRUE)) &
      !str_detect(text_without_compound_country("\\bNorthern Ireland\\b"), regex("\\bIreland\\b", ignore_case = TRUE)) ~
      "Ireland matched only inside Northern Ireland.",
    country_standard == "Mongolia" &
      str_detect(text, regex("\\bInner Mongolia\\b", ignore_case = TRUE)) &
      !str_detect(text_without_compound_country("\\bInner Mongolia\\b"), regex("\\bMongolia\\b", ignore_case = TRUE)) ~
      "Mongolia matched only inside Inner Mongolia.",
    country_standard == "Republic of the Congo" &
      str_detect(country_raw, regex("^Congo$", ignore_case = TRUE)) &
      str_detect(text, regex("Crimean[- ]Congo|CCHF", ignore_case = TRUE)) &
      !str_detect(
        text,
        regex("in\\s+(the\\s+)?Congo|Republic of the Congo", ignore_case = TRUE)
      ) ~
      "Congo matched only inside Crimean-Congo/CCHF disease wording.",
    country_standard == "Madagascar" &
      str_detect(text, regex("\\bEx[- ]+Madagascar\\b", ignore_case = TRUE)) ~
      "Madagascar appears as suspected exposure origin rather than event country.",
    TRUE ~ NA_character_
  )
}

v2_country_reported_cases_policy_hit <- function(local_text) {
  text <- str_to_lower(str_squish(coalesce(local_text, "")))
  positive_hit <- str_detect(
    text,
    regex(
      paste(
        c(
          "number of laboratory confirmed cases and deaths",
          "countries officially reporting cases",
          "following countries have reported laboratory confirmed cases",
          "reported laboratory confirmed cases with no deaths",
          "reported laboratory-confirmed cases originating in the following countries",
          "also reported laboratory-confirmed cases",
          "cases have now also been notified from",
          "has now been reported in .* in addition to",
          "laboratory[- ]confirmed cases of new influenza .* officially reported to who",
          "update on countries and cases .* cumulative total",
          "new cases were reported (in|from)",
          "new probable sars cases were also reported",
          "cases are now being reported in",
          "countries .* newly reported",
          "countries reported their first suspected cases",
          "reported their first (suspected |probable |confirmed )?cases",
          "reported their first probable cases"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  negative_hit <- str_detect(
    text,
    regex(
      paste(
        c(
          "reported no influenza activity",
          "no influenza activity was reported",
          "no similar pattern",
          "no cases",
          "without reporting cases",
          "list of definitions of qualitative indicators"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
  positive_hit & !negative_hit
}

v2_extract_native_country_candidates <- function(
  records = v2_read_records(),
  country_aliases = v2_prepare_country_aliases()
) {
  validation <- v2_validate_country_aliases(country_aliases)
  v2_write_csv(validation, who_don_v2_output_dir("qa", "v2_country_rule_validation.csv"))
  if (nrow(validation) > 0 && any(validation$severity == "blocking")) {
    stop("Blocking country rule validation issues found.", call. = FALSE)
  }

  sections <- v2_country_text_sections(records)
  aliases <- country_aliases %>%
    mutate(
      country_rule_id = paste0("country_alias:", row_number()),
      alias_pattern = vapply(alias, v2_alias_regex, character(1))
    ) %>%
    arrange(priority, desc(nchar(alias)))

  candidates <- purrr::map_dfr(seq_len(nrow(aliases)), function(i) {
    rule <- aliases[i, ]
    hit <- str_detect(sections$section_text, regex(rule$alias_pattern, ignore_case = TRUE))
    if (!any(hit)) {
      return(tibble::tibble())
    }

    sections[hit, ] %>%
      mutate(
        local_evidence_text = vapply(
          section_text,
          v2_local_evidence_text,
          character(1),
          pattern = rule$alias_pattern
        )
      ) %>%
      bind_cols(purrr::pmap_dfr(
        list(.$country_evidence_location, .$Title, .$local_evidence_text),
        ~ v2_country_claim_type(..1, ..2, ..3, rule$is_ambiguous)
      )) %>%
      transmute(
        record_key,
        DonId,
        record_id,
        Title,
        publication_datetime_utc,
        article_url,
        country_raw = v2_first_match(section_text, rule$alias_pattern),
        country_standard = rule$country_standard,
        country_evidence_text = section_text,
        local_evidence_text,
        country_evidence_location,
        country_source_method = case_when(
          country_evidence_location == "title" ~ "native_title_alias",
          country_claim_type == "local_event" ~ "native_event_text_alias",
          TRUE ~ "native_context_text_alias"
        ),
        country_rule_id = rule$country_rule_id,
        country_confidence = case_when(
          country_claim_type == "local_event" & country_evidence_location == "title" ~ "high",
          country_claim_type == "local_event" ~ "medium",
          country_claim_type == "uncertain" ~ "review",
          TRUE ~ "low"
        ),
        country_needs_review = rule$is_ambiguous | country_claim_type != "local_event",
        country_claim_type,
        country_claim_rank,
        country_claim_reason,
        alias_type = rule$alias_type,
        alias_priority = rule$priority,
        native_country_note = rule$notes
      )
  })

  candidates <- candidates %>%
    mutate(
      is_congo_republic_substring_false_positive =
        country_standard == "Republic of the Congo" &
          str_detect(
            local_evidence_text,
            regex("Democratic Republic of (the )?Congo|Congo-Kinshasa|\\bDRC\\b", ignore_case = TRUE)
          ),
      country_false_positive_reason = v2_country_false_positive_reason(
        country_standard,
        country_raw,
        local_evidence_text
      ),
      section_rank = case_when(
        country_evidence_location == "title" ~ 1L,
        country_evidence_location %in% c("summary", "overview", "epidemiology") ~ 2L,
        country_evidence_location == "assessment" ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    filter(!is_congo_republic_substring_false_positive, is.na(country_false_positive_reason)) %>%
    arrange(record_key, country_standard, country_claim_rank, section_rank, alias_priority, desc(nchar(country_raw))) %>%
    group_by(record_key, country_standard, country_claim_type) %>%
    slice(1) %>%
    ungroup() %>%
    select(-section_rank, -is_congo_republic_substring_false_positive, -country_false_positive_reason) %>%
    distinct()

  review_queue <- candidates %>%
    filter(country_needs_review) %>%
    transmute(
      review_id = paste(record_key, country_standard, country_claim_type, "country", sep = "::"),
      record_key,
      DonId,
      Title,
      article_url,
      country_standard,
      country_claim_type,
      country_claim_reason,
      local_evidence_text,
      review_task = "country_candidate_review",
      evidence_text = country_evidence_text,
      current_decision = country_claim_type,
      allowed_decisions = paste(
        c("local_event", "imported_case", "exposure_origin", "background_context",
          "lab_or_partner_context", "surveillance_or_sequence_context", "reject"),
        collapse = "|"
      ),
      reason_for_review = native_country_note
    ) %>%
    distinct()

  summary <- tibble::tibble(
    metric = c(
      "native_country_candidate_rows",
      "native_country_candidate_records",
      "native_country_standards",
      "native_country_local_event_rows",
      "native_country_review_rows"
    ),
    value = c(
      nrow(candidates),
      n_distinct(candidates$record_key),
      n_distinct(candidates$country_standard),
      sum(candidates$country_claim_type == "local_event"),
      nrow(review_queue)
    )
  )

  list(candidates = candidates, review_queue = review_queue, summary = summary)
}

v2_compare_country_candidates <- function(
  native = v2_read_csv(who_don_v2_output_dir("candidates", "who_don_country_candidates_native.csv")),
  association_contract = v2_read_association_contract()
) {
  accepted <- association_contract %>%
    transmute(
      record_key,
      country_standard = coalesce(country_standard_final, country_standard),
      accepted_country_scope = v2_scope_from_clean(pick(everything())),
      accepted_country_evidence_text = v2_first_present(
        pick(everything()),
        c("final_evidence_span", "best_evidence_span", "evidence_sentence")
      )
    ) %>%
    distinct()

  native_best <- native %>%
    arrange(record_key, country_standard, country_claim_rank, alias_priority) %>%
    group_by(record_key, country_standard) %>%
    slice(1) %>%
    ungroup()

  comparison <- full_join(
    accepted %>% mutate(present_accepted = TRUE),
    native_best %>% mutate(present_native = TRUE),
    by = c("record_key", "country_standard")
  ) %>%
    mutate(
      present_accepted = coalesce(present_accepted, FALSE),
      present_native = coalesce(present_native, FALSE),
      diff_category = case_when(
        present_accepted & present_native ~ "exact_record_country_match",
        present_accepted & !present_native ~ "accepted_missing_native",
        !present_accepted & present_native ~ "native_new_country_candidate",
        TRUE ~ "unexpected"
      )
    ) %>%
    arrange(diff_category, record_key, country_standard)

  list(
    diff = comparison,
    summary = comparison %>% count(diff_category, name = "rows"),
    summary_by_claim = comparison %>%
      count(diff_category, country_claim_type, country_evidence_location, name = "rows") %>%
      arrange(diff_category, country_claim_type, country_evidence_location),
    native_new = comparison %>%
      filter(diff_category == "native_new_country_candidate"),
    unmatched_accepted = comparison %>%
      filter(diff_category == "accepted_missing_native")
  )
}

v2_country_adoption_decisions <- function(
  comparison,
  policy_decisions = v2_prepare_country_policy_decisions()
) {
  required <- c(
    "record_key", "country_standard", "diff_category", "present_accepted",
    "present_native", "accepted_country_evidence_text", "country_raw",
    "country_evidence_text", "local_evidence_text", "country_evidence_location",
    "country_source_method", "country_rule_id", "country_confidence",
    "country_claim_type", "country_claim_reason"
  )
  missing_cols <- setdiff(required, names(comparison))
  if (length(missing_cols) > 0) {
    stop(
      "Country comparison missing columns for adoption decisions: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  missing_policy <- setdiff(unique(comparison$diff_category), policy_decisions$diff_category)
  if (length(missing_policy) > 0) {
    stop(
      "Country policy decisions missing diff categories: ",
      paste(missing_policy, collapse = ", "),
      call. = FALSE
    )
  }

  comparison %>%
    left_join(
      policy_decisions %>%
        select(diff_category, policy_decision, policy_review_priority, policy_note),
      by = "diff_category"
    ) %>%
    mutate(
      high_value_native_new_title_event =
        diff_category == "native_new_country_candidate" &
          country_claim_type == "local_event" &
          country_evidence_location == "title" &
          country_confidence == "high",
      medium_native_new_reported_cases_event =
        diff_category == "native_new_country_candidate" &
          country_claim_type == "local_event" &
          country_evidence_location %in% c("summary", "overview", "epidemiology") &
          country_confidence == "medium" &
          v2_country_reported_cases_policy_hit(local_evidence_text),
      policy_decision = case_when(
        high_value_native_new_title_event ~ "accept_native_reviewed",
        medium_native_new_reported_cases_event ~ "accept_native_reviewed",
        TRUE ~ policy_decision
      ),
      policy_review_priority = case_when(
        high_value_native_new_title_event ~ "reviewed_high_value_title_event",
        medium_native_new_reported_cases_event ~ "reviewed_medium_reported_cases_event",
        TRUE ~ policy_review_priority
      ),
      policy_note = case_when(
        high_value_native_new_title_event ~
          "Accepted by deterministic v2 policy: native-new country appears as a high-confidence local-event title country after false-positive filters.",
        medium_native_new_reported_cases_event ~
          "Accepted by deterministic v2 policy medium_native_new_reported_cases_policy: native-new country appears in medium-confidence local-event text with explicit reported-case wording.",
        TRUE ~ policy_note
      )
    ) %>%
    transmute(
      decision_id = paste("country_adoption", record_key, country_standard, sep = "::"),
      record_key,
      country_standard,
      present_accepted,
      present_native,
      diff_category,
      adoption_decision = coalesce(policy_decision, "needs_manual_review"),
      review_priority = coalesce(policy_review_priority, "manual"),
      native_country_raw = country_raw,
      native_country_evidence_text = local_evidence_text,
      native_country_evidence_full_text = country_evidence_text,
      native_country_evidence_location = country_evidence_location,
      native_country_source_method = country_source_method,
      native_country_rule_id = country_rule_id,
      native_country_confidence = country_confidence,
      native_country_claim_type = country_claim_type,
      native_country_claim_reason = country_claim_reason,
      accepted_country_evidence_text,
      adoption_source = case_when(
        high_value_native_new_title_event ~ "country_candidate_high_value_title_policy",
        medium_native_new_reported_cases_event ~ "country_candidate_medium_reported_cases_policy",
        TRUE ~ "country_candidate_policy_decisions"
      ),
      decision_note = coalesce(policy_note, "Country adoption needs manual review.")
    ) %>%
    distinct()
}

v2_apply_country_adoption_decisions <- function(evidence, country_decisions) {
  required <- c(
    "record_key", "country_standard", "adoption_decision",
    "native_country_raw", "native_country_evidence_text",
    "native_country_evidence_location", "native_country_source_method",
    "native_country_rule_id", "native_country_confidence",
    "native_country_claim_type", "native_country_claim_reason",
    "accepted_country_evidence_text"
  )
  missing_cols <- setdiff(required, names(country_decisions))
  if (length(missing_cols) > 0) {
    stop(
      "Country adoption decisions missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  native_accept_decisions <- c("accept_native", "accept_native_reviewed")
  adopted <- country_decisions %>%
    filter(adoption_decision %in% c(native_accept_decisions, "accept_legacy_exception")) %>%
    distinct(record_key, country_standard, .keep_all = TRUE)

  missing_accepted <- evidence %>%
    distinct(record_key, country_standard) %>%
    anti_join(adopted, by = c("record_key", "country_standard"))
  if (nrow(missing_accepted) > 0) {
    stop(
      "Country adoption decisions do not cover accepted evidence countries: ",
      paste(head(paste(missing_accepted$record_key, missing_accepted$country_standard, sep = "::"), 20), collapse = ", "),
      call. = FALSE
    )
  }

  reviewed_native_country_additions <- adopted %>%
    filter(adoption_decision %in% native_accept_decisions) %>%
    anti_join(
      evidence %>% distinct(record_key, country_standard),
      by = c("record_key", "country_standard")
    )

  added_country_rows <- evidence %>%
    semi_join(reviewed_native_country_additions, by = "record_key") %>%
    group_by(record_key, disease_standard, influenza_type, influenza_subtype) %>%
    slice(1) %>%
    ungroup() %>%
    select(-country_standard) %>%
    inner_join(reviewed_native_country_additions, by = "record_key", relationship = "many-to-many") %>%
    mutate(
      country_raw = native_country_raw,
      country_evidence_text = native_country_evidence_text,
      country_evidence_location = native_country_evidence_location,
      country_source_method = native_country_source_method,
      country_rule_id = native_country_rule_id,
      country_confidence = native_country_confidence,
      country_needs_review = FALSE,
      country_claim_type = native_country_claim_type,
      country_claim_reason = native_country_claim_reason,
      source_method = paste(source_method, "v2_native_new_country_reviewed_adoption", sep = "+"),
      country_adoption_decision = adoption_decision,
      country_adoption_decision_id = decision_id,
      country_adoption_note = decision_note
    )

  updated_existing <- evidence %>%
    left_join(
      adopted %>%
        select(
          record_key,
          country_standard,
          country_adoption_decision = adoption_decision,
          country_adoption_decision_id = decision_id,
          country_adoption_note = decision_note,
          native_country_raw,
          native_country_evidence_text,
          native_country_evidence_location,
          native_country_source_method,
          native_country_rule_id,
          native_country_confidence,
          native_country_claim_type,
          native_country_claim_reason,
          accepted_country_evidence_text
        ),
      by = c("record_key", "country_standard")
    ) %>%
    mutate(
      country_raw = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_raw,
        TRUE ~ country_raw
      ),
      country_evidence_text = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_evidence_text,
        country_adoption_decision == "accept_legacy_exception" ~ accepted_country_evidence_text,
        TRUE ~ country_evidence_text
      ),
      country_evidence_location = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_evidence_location,
        country_adoption_decision == "accept_legacy_exception" ~ "legacy_country_exception",
        TRUE ~ country_evidence_location
      ),
      country_source_method = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_source_method,
        country_adoption_decision == "accept_legacy_exception" ~ "legacy_country_exception",
        TRUE ~ country_source_method
      ),
      country_rule_id = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_rule_id,
        country_adoption_decision == "accept_legacy_exception" ~ "legacy_country_exception_missing_native",
        TRUE ~ country_rule_id
      ),
      country_confidence = case_when(
        country_adoption_decision %in% native_accept_decisions ~ native_country_confidence,
        country_adoption_decision == "accept_legacy_exception" ~ "medium",
        TRUE ~ country_confidence
      ),
      country_needs_review = country_adoption_decision == "accept_legacy_exception",
      country_claim_type = coalesce(native_country_claim_type, "legacy_exception"),
      country_claim_reason = coalesce(native_country_claim_reason, country_adoption_note),
      source_method = case_when(
        country_adoption_decision %in% native_accept_decisions & source_method == "clean_final_seed" ~
          "v2_native_country_seeded_disease",
        country_adoption_decision %in% native_accept_decisions & source_method == "v2_seeded_reviewed_adoption" ~
          "v2_native_country_seeded_reviewed_disease",
        country_adoption_decision %in% native_accept_decisions & source_method == "v2_native_reviewed_adoption" ~
          "v2_native_country_disease_reviewed_adoption",
        country_adoption_decision == "accept_legacy_exception" ~
          paste(source_method, "legacy_country_exception", sep = "+"),
        TRUE ~ source_method
      )
    )

  bind_rows(
    updated_existing,
    added_country_rows %>% select(all_of(names(updated_existing)))
  ) %>%
    distinct()
}
