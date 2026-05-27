suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_compare.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_final_shaping.R"))

v2_option_a_exception_blank <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::if_else(is.na(x), "", x)
}

v2_option_a_add_compare_keys <- function(data) {
  if (!"influenza_type" %in% names(data)) {
    data$influenza_type <- NA_character_
  }
  if (!"influenza_subtype" %in% names(data)) {
    data$influenza_subtype <- NA_character_
  }

  data %>%
    v2_add_influenza_compare_keys() %>%
    mutate(
      country_key = v2_option_a_exception_blank(country_standard),
      disease_key = v2_option_a_exception_blank(disease_compare_key),
      influenza_type_key = v2_option_a_exception_blank(influenza_type_compare_key),
      influenza_subtype_key = v2_option_a_exception_blank(influenza_subtype_compare_key),
      scope_key = v2_option_a_exception_blank(association_scope),
      semantic_key = paste(
        v2_option_a_exception_blank(record_key),
        country_key,
        disease_key,
        influenza_type_key,
        influenza_subtype_key,
        sep = "::"
      ),
      row_key = paste(
        v2_option_a_exception_blank(record_key),
        country_key,
        disease_key,
        influenza_type_key,
        influenza_subtype_key,
        scope_key,
        sep = "::"
      )
    )
}

v2_option_a_keep_current_exception_path <- function() {
  who_don_v2_output_dir("review", "option_a_full_article_keep_current_exception_rows.csv")
}

v2_materialize_option_a_keep_current_exceptions <- function(
  current_audit,
  option_a_audit,
  full_article_decisions,
  full_article_workpack
) {
  decisions <- full_article_decisions %>%
    mutate(across(everything(), v2_option_a_exception_blank)) %>%
    filter(review_decision == "keep_current_exception")

  if (nrow(decisions) == 0L) {
    return(tibble::tibble())
  }

  workpack <- full_article_workpack %>%
    mutate(across(any_of(c("residual_id", "record_key", "country_standard", "disease_standard")), v2_option_a_exception_blank))

  decision_rows <- decisions %>%
    left_join(
      workpack %>%
        select(
          residual_id,
          workpack_record_key = record_key,
          workpack_influenza_type = influenza_type,
          workpack_influenza_subtype = influenza_subtype,
          review_weight,
          cluster_weighted_rows,
          cluster_pattern
        ),
      by = "residual_id"
    ) %>%
    mutate(
      record_key = workpack_record_key,
      influenza_type = workpack_influenza_type,
      influenza_subtype = workpack_influenza_subtype
    )

  missing_workpack <- decision_rows %>%
    filter(record_key == "" | is.na(record_key))
  if (nrow(missing_workpack) > 0L) {
    stop(
      "Option A keep-current exceptions failed to join workpack record keys. First residual IDs: ",
      paste(head(missing_workpack$residual_id, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  decision_option_keys <- decision_rows %>%
    mutate(association_scope = option_a_scope) %>%
    v2_option_a_add_compare_keys() %>%
    transmute(residual_id, option_row_key = row_key)

  decision_current_keys <- decision_rows %>%
    mutate(association_scope = current_scope) %>%
    v2_option_a_add_compare_keys() %>%
    transmute(residual_id, current_row_key = row_key)

  decision_rows <- decision_rows %>%
    left_join(decision_option_keys, by = "residual_id") %>%
    left_join(decision_current_keys, by = "residual_id")

  current_audit_keys <- current_audit %>% v2_option_a_add_compare_keys()
  option_a_audit_keys <- option_a_audit %>% v2_option_a_add_compare_keys()

  missing_current <- decision_rows %>%
    filter(!current_row_key %in% current_audit_keys$row_key)
  missing_option <- decision_rows %>%
    filter(!option_row_key %in% option_a_audit_keys$row_key)
  if (nrow(missing_current) > 0L || nrow(missing_option) > 0L) {
    stop(
      "Option A keep-current exceptions have unmatched current or native row keys.",
      call. = FALSE
    )
  }

  current_audit_keys %>%
    inner_join(
      decision_rows %>%
        select(
          current_row_key,
          option_row_key,
          residual_id,
          cluster_id,
          direction,
          review_decision,
          review_confidence,
          review_rationale,
          article_evidence_quote_or_excerpt,
          article_section_used,
          reviewer_id,
          reviewed_at
        ),
      by = c("row_key" = "current_row_key")
    ) %>%
    mutate(
      source_method = paste0(source_method, "+option_a_full_article_keep_current_exception"),
      scope_rule_id = "option_a_full_article_keep_current_exception",
      scope_reason = paste0("Full-article review preserves current production scope: ", review_rationale),
      scope_evidence_text = if_else(
        article_evidence_quote_or_excerpt == "",
        scope_evidence_text,
        article_evidence_quote_or_excerpt
      ),
      review_status = "option_a_full_article_keep_current_exception",
      review_decision_id = paste0("option_a_full_article:", residual_id),
      final_review_source = "option_a_full_article_scope_review",
      final_review_note = review_rationale,
      option_a_exception_current_row_key = row_key
    ) %>%
    select(
      option_a_exception_current_row_key,
      option_a_exception_option_row_key = option_row_key,
      residual_id,
      cluster_id,
      direction,
      review_decision,
      review_confidence,
      review_rationale,
      article_evidence_quote_or_excerpt,
      article_section_used,
      reviewer_id,
      reviewed_at,
      any_of(names(current_audit))
    )
}

v2_apply_option_a_keep_current_exceptions <- function(audit, exception_rows) {
  if (nrow(exception_rows) == 0L) {
    return(audit)
  }

  required_cols <- c(
    "option_a_exception_current_row_key",
    "option_a_exception_option_row_key"
  )
  missing_cols <- setdiff(required_cols, names(exception_rows))
  if (length(missing_cols) > 0L) {
    stop(
      "Option A keep-current exception rows missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  audit_keys <- audit %>% v2_option_a_add_compare_keys()
  remove_keys <- unique(c(
    exception_rows$option_a_exception_option_row_key,
    exception_rows$option_a_exception_current_row_key
  ))
  output_cols <- names(audit)

  for (col in setdiff(output_cols, names(exception_rows))) {
    exception_rows[[col]] <- NA
  }

  audit_keys %>%
    filter(!row_key %in% remove_keys) %>%
    select(all_of(output_cols)) %>%
    bind_rows(exception_rows %>% select(all_of(output_cols))) %>%
    arrange(record_key, country_standard, disease_standard, association_scope)
}

v2_apply_option_a_keep_current_exceptions_from_file <- function(audit) {
  exception_path <- v2_option_a_keep_current_exception_path()
  if (!file.exists(exception_path)) {
    return(audit)
  }

  exceptions <- v2_read_csv(exception_path)
  v2_apply_option_a_keep_current_exceptions(audit, exceptions)
}
