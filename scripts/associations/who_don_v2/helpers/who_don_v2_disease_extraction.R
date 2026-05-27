library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_rules.R"))

v2_disease_text_sections <- function(records) {
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
      disease_evidence_location = names(section_cols)[match(source_column, unname(section_cols))],
      section_text = str_squish(as.character(section_text)),
      section_text_norm = str_to_lower(section_text)
    ) %>%
    select(-record_title) %>%
    filter(!is.na(section_text), section_text != "")
}

v2_alias_regex <- function(alias) {
  escaped <- stringr::str_replace_all(alias, "([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1")
  escaped <- stringr::str_replace_all(escaped, "\\s+", "\\\\s+")
  paste0("(?<![[:alnum:]])", escaped, "(?![-[:alnum:]])")
}

v2_first_match <- function(text, pattern, ignore_case = TRUE) {
  match <- str_match(text, regex(pattern, ignore_case = ignore_case))[, 1]
  ifelse(is.na(match), NA_character_, match)
}

v2_local_evidence_text <- function(text, pattern) {
  text <- str_squish(as.character(text))
  loc <- str_locate(text, regex(pattern, ignore_case = TRUE))[1, ]
  if (is.na(loc[1])) {
    return(text)
  }

  before <- str_sub(text, 1, loc[1] - 1)
  after <- str_sub(text, loc[2] + 1, str_length(text))
  left_breaks <- str_locate_all(before, "[.!?;]\\s+")[[1]]
  right_breaks <- str_locate_all(after, "[.!?;]\\s+")[[1]]

  start <- if (nrow(left_breaks) == 0) 1L else left_breaks[nrow(left_breaks), 2] + 1L
  end <- if (nrow(right_breaks) == 0) str_length(text) else loc[2] + right_breaks[1, 1]

  str_squish(str_sub(text, start, end))
}

v2_event_anchor_classify <- function(location, title, local_text, is_generic_label) {
  text <- str_to_lower(str_squish(paste(title, local_text, sep = " ")))
  background_pattern <- paste(
    c(
      "history of", "historical", "since 2003", "since 2004", "previously",
      "cumulative total", "does not change", "overall risk", "differential diagnosis",
      "advice", "preparedness", "further information", "references?", "source:",
      "centers for disease control", "who advice", "risk assessment",
      "can cause", "may cause", "known causes", "other known causes",
      "family includes", "is endemic", "are endemic", "also endemic",
      "could lead to", "in addition to", "detailed epidemiological",
      "recommendation", "recommended", "such a recommendation"
    ),
    collapse = "|"
  )
  context_pattern <- paste(
    c("rule out", "negative for", "tested for", "surveillance", "monitoring", "background"),
    collapse = "|"
  )

  has_background <- str_detect(text, regex(background_pattern, ignore_case = TRUE))
  has_context <- str_detect(text, regex(context_pattern, ignore_case = TRUE))

  if (location == "title") {
    return(tibble::tibble(
      event_anchor_class = "event_title",
      event_anchor_score = 1L,
      event_anchor_reason = "Disease alias appears in title."
    ))
  }
  if (location == "summary" && !has_background && !has_context) {
    return(tibble::tibble(
      event_anchor_class = "event_summary",
      event_anchor_score = 2L,
      event_anchor_reason = "Disease alias appears in summary without background language."
    ))
  }
  if (location %in% c("overview", "epidemiology") && !has_background && !has_context) {
    return(tibble::tibble(
      event_anchor_class = "event_overview",
      event_anchor_score = 3L,
      event_anchor_reason = "Disease alias appears in event overview or epidemiology text."
    ))
  }
  if (has_background) {
    return(tibble::tibble(
      event_anchor_class = "background",
      event_anchor_score = 7L,
      event_anchor_reason = "Disease alias appears with background, history, advice, or reference language."
    ))
  }
  if (has_context) {
    return(tibble::tibble(
      event_anchor_class = "context",
      event_anchor_score = 6L,
      event_anchor_reason = "Disease alias appears in contextual or surveillance wording."
    ))
  }
  if (is_generic_label) {
    return(tibble::tibble(
      event_anchor_class = "review",
      event_anchor_score = 5L,
      event_anchor_reason = "Generic disease alias requires review."
    ))
  }
  tibble::tibble(
    event_anchor_class = "review",
    event_anchor_score = 5L,
    event_anchor_reason = "Disease alias is outside the strongest event sections."
  )
}

v2_disease_standard_from_influenza_subtype <- function(
  subtype,
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  subtype_lookup <- stringr::str_squish(as.character(subtype))
  subtype_lookup <- dplyr::na_if(subtype_lookup, "")

  standardization <- influenza_standardization %>%
    transmute(
      influenza_subtype = stringr::str_squish(as.character(influenza_subtype)),
      canonical_disease_standard = stringr::str_squish(as.character(canonical_disease_standard))
    ) %>%
    filter(!is.na(influenza_subtype), influenza_subtype != "") %>%
    distinct(influenza_subtype, .keep_all = TRUE)

  dplyr::tibble(influenza_subtype = subtype_lookup) %>%
    left_join(standardization, by = "influenza_subtype") %>%
    mutate(
      disease_standard = case_when(
        !is.na(canonical_disease_standard) & canonical_disease_standard != "" ~ canonical_disease_standard,
        !is.na(influenza_subtype) & influenza_subtype != "" ~ paste0("Influenza A(", influenza_subtype, ")"),
        TRUE ~ NA_character_
      )
    ) %>%
    pull(disease_standard)
}

v2_extract_alias_hits <- function(sections, disease_aliases) {
  aliases <- disease_aliases %>%
    mutate(
      disease_rule_id = paste0("disease_alias:", row_number()),
      alias_pattern = vapply(alias, v2_alias_regex, character(1))
    ) %>%
    arrange(priority, desc(nchar(alias)))

  purrr::map_dfr(seq_len(nrow(aliases)), function(i) {
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
        list(.$disease_evidence_location, .$Title, .$local_evidence_text),
        ~ v2_event_anchor_classify(..1, ..2, ..3, rule$is_generic_label)
      )) %>%
      transmute(
        record_key,
        DonId,
        record_id,
        Title,
        publication_datetime_utc,
        article_url,
        disease_raw = v2_first_match(section_text, rule$alias_pattern),
        disease_standard = rule$disease_standard,
        disease_evidence_text = section_text,
        local_evidence_text,
        disease_evidence_location,
        disease_source_method = case_when(
          disease_evidence_location == "title" ~ "native_title_alias",
          rule$alias_type == "influenza_generic" ~ "native_generic_label",
          rule$alias_type == "influenza_subtype" ~ "native_influenza_subtype_rule",
          TRUE ~ "native_section_alias"
        ),
        disease_rule_id = rule$disease_rule_id,
        disease_confidence = case_when(
          disease_evidence_location == "title" & !rule$is_generic_label ~ "high",
          disease_evidence_location %in% c("overview", "epidemiology", "summary") & !rule$is_generic_label ~ "medium",
          rule$is_generic_label ~ "review",
          TRUE ~ "low"
        ),
        disease_needs_review = rule$is_generic_label |
          disease_confidence %in% c("low", "review") |
          event_anchor_class %in% c("context", "background", "review"),
        influenza_type = NA_character_,
        influenza_subtype = NA_character_,
        influenza_subtype_candidates = NA_character_,
        influenza_subtype_evidence_span = NA_character_,
        influenza_subtype_evidence_scope = NA_character_,
        event_anchor_class,
        event_anchor_score,
        event_anchor_reason,
        native_extraction_note = rule$notes,
        alias_priority = rule$priority,
        alias_type = rule$alias_type
      )
  })
}

v2_apply_influenza_subtype_rules <- function(
  candidates,
  influenza_rules,
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  subtype_rules <- influenza_rules %>%
    filter(specificity != "generic") %>%
    arrange(priority)

  find_subtype <- function(title, local_span, event_anchor_class) {
    title_hits <- subtype_rules %>%
      mutate(matched = str_detect(title, regex(pattern, ignore_case = TRUE))) %>%
      filter(matched) %>%
      arrange(priority)
    local_hits <- subtype_rules %>%
      mutate(matched = str_detect(local_span, regex(pattern, ignore_case = TRUE))) %>%
      filter(matched) %>%
      arrange(priority)
    hits <- bind_rows(title_hits, local_hits) %>%
      distinct(rule_id, .keep_all = TRUE) %>%
      arrange(priority)

    if (nrow(hits) == 0) {
      return(tibble::tibble(
        influenza_type_rule = NA_character_,
        influenza_subtype_rule = NA_character_,
        influenza_subtype_candidates_rule = NA_character_,
        influenza_title_subtype_candidates_rule = NA_character_,
        influenza_subtype_evidence_span_rule = NA_character_,
        influenza_subtype_evidence_scope_rule = "none",
        subtype_rule_id = NA_character_
      ))
    }

    best <- if (nrow(title_hits) > 0) title_hits[1, ] else hits[1, ]
    evidence_scope <- case_when(
      best$rule_id %in% title_hits$rule_id ~ "title",
      event_anchor_class %in% c("event_summary", "event_overview", "event_title") ~ "local_event",
      TRUE ~ "background_context"
    )
    tibble::tibble(
      influenza_type_rule = best$influenza_type,
      influenza_subtype_rule = best$influenza_subtype,
      influenza_subtype_candidates_rule = paste(hits$influenza_subtype, collapse = "|"),
      influenza_title_subtype_candidates_rule = if_else(
        nrow(title_hits) > 0,
        paste(title_hits$influenza_subtype, collapse = "|"),
        NA_character_
      ),
      influenza_subtype_evidence_span_rule = if_else(evidence_scope == "title", title, local_span),
      influenza_subtype_evidence_scope_rule = evidence_scope,
      subtype_rule_id = best$rule_id
    )
  }

  subtype_hits <- purrr::pmap_dfr(
    list(candidates$Title, candidates$local_evidence_text, candidates$event_anchor_class),
    find_subtype
  )

  bind_cols(candidates, subtype_hits) %>%
    mutate(
      is_influenza_candidate = str_detect(disease_standard, regex("influenza", ignore_case = TRUE)) |
        str_detect(disease_raw, regex("influenza|\\bH[0-9]", ignore_case = TRUE)),
      subtype_can_promote = is_influenza_candidate &
        !is.na(influenza_subtype_rule) &
        influenza_subtype_evidence_scope_rule %in% c("title", "local_event"),
      promoted_disease_standard = v2_disease_standard_from_influenza_subtype(
        influenza_subtype_rule,
        influenza_standardization
      ),
      disease_standard = case_when(
        subtype_can_promote ~ promoted_disease_standard,
        TRUE ~ disease_standard
      ),
      influenza_type = if_else(is_influenza_candidate, coalesce(influenza_type_rule, "influenza"), influenza_type),
      influenza_subtype = if_else(subtype_can_promote, influenza_subtype_rule, influenza_subtype),
      influenza_subtype_candidates = if_else(
        is_influenza_candidate,
        influenza_subtype_candidates_rule,
        influenza_subtype_candidates
      ),
      influenza_title_subtype_candidates = if_else(
        is_influenza_candidate,
        influenza_title_subtype_candidates_rule,
        NA_character_
      ),
      influenza_subtype_evidence_span = if_else(
        subtype_can_promote,
        influenza_subtype_evidence_span_rule,
        influenza_subtype_evidence_span
      ),
      influenza_subtype_evidence_scope = if_else(
        is_influenza_candidate,
        influenza_subtype_evidence_scope_rule,
        influenza_subtype_evidence_scope
      ),
      disease_source_method = if_else(
        subtype_can_promote,
        "native_influenza_subtype_rule",
        disease_source_method
      ),
      disease_rule_id = if_else(
        subtype_can_promote,
        paste(disease_rule_id, subtype_rule_id, sep = "|"),
        disease_rule_id
      ),
      disease_confidence = case_when(
        subtype_can_promote ~ "high",
        TRUE ~ disease_confidence
      ),
      disease_needs_review = case_when(
        is_influenza_candidate & is.na(influenza_subtype_rule) ~ TRUE,
        is_influenza_candidate & !subtype_can_promote ~ TRUE,
        subtype_can_promote ~ FALSE,
        TRUE ~ disease_needs_review
      ),
      native_extraction_note = case_when(
        subtype_can_promote ~
          paste(native_extraction_note, "Influenza subtype resolved from local/title subtype evidence."),
        is_influenza_candidate & !is.na(influenza_subtype_rule) ~
          paste(native_extraction_note, "Influenza subtype found only in background/context evidence; not promoted."),
        is_influenza_candidate & is.na(influenza_subtype_rule) ~
          paste(native_extraction_note, "Generic influenza candidate requires subtype review."),
        TRUE ~ native_extraction_note
      )
    ) %>%
    select(
      -influenza_type_rule,
      -influenza_subtype_rule,
      -influenza_subtype_candidates_rule,
      -influenza_title_subtype_candidates_rule,
      -influenza_subtype_evidence_span_rule,
      -influenza_subtype_evidence_scope_rule,
      -subtype_rule_id,
      -is_influenza_candidate,
      -subtype_can_promote,
      -promoted_disease_standard
    )
}

v2_expand_multi_title_influenza_candidates <- function(
  candidates,
  influenza_rules,
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  subtype_lookup <- influenza_rules %>%
    filter(specificity != "generic") %>%
    transmute(
      influenza_subtype_multi = influenza_subtype,
      influenza_type_multi = influenza_type,
      subtype_rule_id_multi = rule_id
    ) %>%
    distinct(influenza_subtype_multi, .keep_all = TRUE)

  expandable <- candidates %>%
    filter(
      str_detect(disease_standard, regex("influenza", ignore_case = TRUE)),
      influenza_subtype_evidence_scope == "title",
      str_detect(coalesce(influenza_title_subtype_candidates, ""), fixed("|"))
    ) %>%
    mutate(.expand_row_id = row_number())

  if (nrow(expandable) == 0) {
    return(candidates)
  }

  non_expandable <- candidates %>%
    anti_join(
      expandable %>% select(record_key, disease_standard, disease_evidence_location, disease_rule_id),
      by = c("record_key", "disease_standard", "disease_evidence_location", "disease_rule_id")
    )

  expanded <- expandable %>%
    separate_rows(influenza_title_subtype_candidates, sep = "\\|") %>%
    rename(influenza_subtype_multi = influenza_title_subtype_candidates) %>%
    left_join(subtype_lookup, by = "influenza_subtype_multi") %>%
    mutate(
      disease_standard = v2_disease_standard_from_influenza_subtype(
        influenza_subtype_multi,
        influenza_standardization
      ),
      influenza_type = coalesce(influenza_type_multi, influenza_type),
      influenza_subtype = influenza_subtype_multi,
      influenza_subtype_candidates = influenza_subtype_multi,
      influenza_title_subtype_candidates = influenza_subtype_multi,
      influenza_subtype_evidence_span = Title,
      influenza_subtype_evidence_scope = "title",
      disease_source_method = "native_influenza_subtype_rule",
      disease_rule_id = paste(disease_rule_id, "title_multi", subtype_rule_id_multi, sep = "|"),
      disease_confidence = "high",
      disease_needs_review = FALSE,
      native_extraction_note = paste(
        native_extraction_note,
        "Expanded from explicit multi-subtype influenza title evidence."
      )
    ) %>%
    select(
      -influenza_type_multi,
      -influenza_subtype_multi,
      -subtype_rule_id_multi,
      -.expand_row_id
    )

  bind_rows(non_expandable, expanded)
}

v2_collapse_generic_influenza_candidates <- function(candidates) {
  specific_influenza <- candidates %>%
    filter(str_detect(disease_standard, regex("influenza", ignore_case = TRUE)), !is.na(influenza_subtype)) %>%
    distinct(record_key)

  candidates %>%
    left_join(specific_influenza %>% mutate(has_specific_influenza = TRUE), by = "record_key") %>%
    mutate(
      has_specific_influenza = coalesce(has_specific_influenza, FALSE),
      generic_influenza = disease_standard == "Influenza" & is.na(influenza_subtype)
    ) %>%
    filter(!(generic_influenza & has_specific_influenza)) %>%
    select(-has_specific_influenza, -generic_influenza)
}

v2_annotate_disease_rule_model <- function(
  candidates,
  disease_rule_model = v2_read_csv(who_don_v2_rules_dir("disease_rule_model.csv"))
) {
  rule_policy <- disease_rule_model %>%
    group_by(standard_label) %>%
    summarise(
      disease_group_rule = first(na.omit(disease_group), default = NA_character_),
      specificity_rank = min(specificity_rank, na.rm = TRUE),
      requires_title_or_event_anchor = any(requires_title_or_event_anchor, na.rm = TRUE),
      background_exclusion_pattern = first(na.omit(background_exclusion_pattern), default = NA_character_),
      .groups = "drop"
    ) %>%
    rename(disease_standard = standard_label)

  candidates %>%
    left_join(rule_policy, by = "disease_standard") %>%
    mutate(
      specificity_rank = if_else(is.infinite(specificity_rank), NA_integer_, as.integer(specificity_rank)),
      requires_title_or_event_anchor = coalesce(requires_title_or_event_anchor, FALSE)
    )
}

v2_apply_native_candidate_promotion_gate <- function(candidates) {
  noisy_event_labels <- c(
    "Acute respiratory infection",
    "COVID-19",
    "Malaria",
    "Respiratory illness"
  )

  candidates %>%
    mutate(
      is_noisy_event_label = disease_standard %in% noisy_event_labels,
      is_influenza_specific = str_detect(disease_standard, regex("influenza", ignore_case = TRUE)) &
        !is.na(influenza_subtype),
      has_strong_event_anchor = event_anchor_class %in% c("event_title", "event_summary"),
      has_event_overview_anchor = event_anchor_class == "event_overview",
      candidate_promotion_status = case_when(
        event_anchor_class %in% c("context", "review") ~ "review_only",
        is_noisy_event_label & !has_strong_event_anchor ~ "review_only",
        is_influenza_specific & influenza_subtype_evidence_scope %in% c("title", "local_event") ~ "candidate",
        has_strong_event_anchor | has_event_overview_anchor ~ "candidate",
        TRUE ~ "review_only"
      ),
      candidate_promotion_reason = case_when(
        candidate_promotion_status == "candidate" & is_influenza_specific ~
          "Specific influenza subtype has title or local event evidence.",
        candidate_promotion_status == "candidate" & has_strong_event_anchor ~
          "Disease has title or summary event anchor.",
        candidate_promotion_status == "candidate" & has_event_overview_anchor ~
          "Disease has event overview anchor.",
        event_anchor_class %in% c("context", "review") ~
          "Disease evidence is contextual or review-only.",
        is_noisy_event_label & !has_strong_event_anchor ~
          "Background-prone disease label requires title or summary support.",
        TRUE ~ "Candidate did not meet native promotion gate."
      )
    ) %>%
    filter(candidate_promotion_status == "candidate") %>%
    select(
      -is_influenza_specific,
      -is_noisy_event_label,
      -has_strong_event_anchor,
      -has_event_overview_anchor
    )
}

v2_rank_native_disease_candidates <- function(candidates) {
  candidates %>%
    mutate(
      section_rank = case_when(
        disease_evidence_location == "title" ~ 1L,
        disease_evidence_location %in% c("overview", "epidemiology") ~ 2L,
        disease_evidence_location == "summary" ~ 3L,
        TRUE ~ 4L
      ),
      confidence_rank = case_when(
        disease_confidence == "high" ~ 1L,
        disease_confidence == "medium" ~ 2L,
        disease_confidence == "low" ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    filter(event_anchor_class != "background") %>%
    arrange(record_key, disease_standard, event_anchor_score, section_rank, confidence_rank, alias_priority, desc(nchar(disease_raw))) %>%
    group_by(record_key, disease_standard, disease_evidence_location, disease_rule_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(
      record_key,
      DonId,
      record_id,
      Title,
      publication_datetime_utc,
      article_url,
      disease_raw,
      disease_standard,
      disease_evidence_text,
      local_evidence_text,
      disease_evidence_location,
      disease_source_method,
      disease_rule_id,
      disease_confidence,
      disease_needs_review,
      influenza_type,
      influenza_subtype,
      influenza_subtype_candidates,
      influenza_title_subtype_candidates,
      influenza_subtype_evidence_span,
      influenza_subtype_evidence_scope,
      event_anchor_class,
      event_anchor_score,
      event_anchor_reason,
      candidate_promotion_status,
      candidate_promotion_reason,
      native_extraction_note
    )
}

v2_native_disease_review_queue <- function(candidates) {
  candidates %>%
    filter(disease_needs_review) %>%
    transmute(
      review_id = paste(record_key, disease_standard, disease_evidence_location, "disease", sep = "::"),
      record_key,
      DonId,
      Title,
      article_url,
      disease_standard,
      influenza_type,
      influenza_subtype,
      event_anchor_class,
      event_anchor_reason,
      candidate_promotion_status,
      candidate_promotion_reason,
      local_evidence_text,
      review_task = if_else(
        str_detect(disease_standard, regex("influenza", ignore_case = TRUE)),
        "influenza_subtype_resolution",
        "disease_candidate_review"
      ),
      evidence_text = disease_evidence_text,
      current_decision = disease_standard,
      allowed_decisions = "accept|reject|standardize|resolve_influenza_subtype",
      reason_for_review = native_extraction_note
    ) %>%
    distinct()
}

v2_extract_native_disease_candidates <- function(
  records = v2_read_records(),
  disease_aliases = v2_read_csv(who_don_v2_rules_dir("disease_aliases.csv")),
  influenza_rules = v2_read_csv(who_don_v2_rules_dir("influenza_subtype_rules.csv")),
  influenza_standardization = v2_read_csv(who_don_v2_rules_dir("influenza_label_standardization.csv"))
) {
  validation <- v2_validate_disease_aliases(disease_aliases)
  if (any(validation$severity == "blocking")) {
    v2_write_csv(validation, who_don_v2_output_dir("qa", "v2_disease_rule_validation.csv"))
    stop("Blocking disease rule validation issues found.", call. = FALSE)
  }

  sections <- v2_disease_text_sections(records)
  candidates <- sections %>%
    v2_extract_alias_hits(disease_aliases) %>%
    v2_apply_influenza_subtype_rules(influenza_rules, influenza_standardization) %>%
    v2_expand_multi_title_influenza_candidates(influenza_rules, influenza_standardization) %>%
    v2_collapse_generic_influenza_candidates() %>%
    v2_annotate_disease_rule_model() %>%
    v2_apply_native_candidate_promotion_gate() %>%
    v2_rank_native_disease_candidates()

  review_queue <- v2_native_disease_review_queue(candidates)

  summary <- tibble::tibble(
    metric = c(
      "native_disease_candidate_rows",
      "native_disease_candidate_records",
      "native_disease_standards",
      "native_disease_review_rows",
      "native_influenza_candidate_rows",
      "native_influenza_subtype_resolved_rows",
      "native_event_anchor_rows",
      "native_context_review_rows"
    ),
    value = c(
      nrow(candidates),
      n_distinct(candidates$record_key),
      n_distinct(candidates$disease_standard),
      nrow(review_queue),
      sum(str_detect(candidates$disease_standard, regex("influenza", ignore_case = TRUE))),
      sum(str_detect(candidates$disease_standard, regex("influenza", ignore_case = TRUE)) & !is.na(candidates$influenza_subtype)),
      sum(candidates$event_anchor_class %in% c("event_title", "event_summary", "event_overview")),
      sum(candidates$event_anchor_class %in% c("context", "review"))
    )
  )

  list(candidates = candidates, review_queue = review_queue, summary = summary)
}
