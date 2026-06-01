# ------------------------------------------------------------------------------
# 5_6d_QA_Vector_Competence_Unmatched.R
# ------------------------------------------------------------------------------
# Purpose: Add a lightweight QA layer to unmatched disease-vector competence
#          rows so we can separate fixable join mismatches from rows that are
#          genuinely outside the curated disease-vector inclusion gate.
#
# Inputs : vector_screening_qa_path("vector_competence_join_unmatched.csv")
# Outputs: vector_screening_qa_path(
#            "vector_competence_join_unmatched_review.csv"
#          )
#          vector_screening_qa_path(
#            "vector_competence_join_unmatched_summary.csv"
#          )
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

input_path <- vector_screening_qa_path("vector_competence_join_unmatched.csv")
review_path <- file.path(vector_screening_qa_dir, "vector_competence_join_unmatched_review.csv")
summary_path <- file.path(vector_screening_qa_dir, "vector_competence_join_unmatched_summary.csv")
dir.create(vector_screening_qa_dir, recursive = TRUE, showWarnings = FALSE)

unmatched <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    note_text = dplyr::coalesce(competence_note_examples, ""),
    unmatched_reason = case_when(
      str_detect(vector_join_key, "\\bspp\\b") |
        str_detect(vector_competence_species, "\\bspp\\b") ~ "genus_or_family_level",
      str_detect(vector_competence_species, " x ") |
        str_detect(vector_competence_species, "hybrid") ~ "hybrid_or_cross",
      str_detect(vector_competence_species, "s\\.l\\.|complex|biotype|form") ~
        "complex_or_infraspecific_label",
      str_detect(vector_competence_species, "\\(") |
        str_detect(vector_competence_species, "/") ~ "format_variant_or_authorship",
      str_detect(note_text, "spelled|spelling|OCR") ~
        "spelling_or_ocr_issue",
      vector_competence_status %in% c("competent", "mixed") ~
        "positive_competence_without_curated_match",
      TRUE ~ "not_in_curated_disease_vector_table"
    ),
    recommended_action = case_when(
      unmatched_reason %in% c("format_variant_or_authorship", "spelling_or_ocr_issue") ~
        "review_name_cleanup_or_manual_map",
      unmatched_reason == "complex_or_infraspecific_label" ~
        "review_taxonomic_grain_and_possible_manual_map",
      unmatched_reason == "genus_or_family_level" ~
        "keep_unmatched_unless_genus_level_links_are_intended",
      unmatched_reason == "hybrid_or_cross" ~
        "keep_unmatched_or_review_hybrid_handling_explicitly",
      vector_competence_status == "not_competent" ~
        "no_join_fix_needed_negative_evidence_only",
      TRUE ~ "manual_review_needed"
    ),
    qa_priority = case_when(
      vector_competence_status %in% c("competent", "mixed") ~ "high",
      vector_competence_status == "unclear" ~ "medium",
      TRUE ~ "low"
    )
  ) %>%
  select(-note_text) %>%
  arrange(desc(qa_priority), competence_disease_name, vector_competence_species)

summary_table <- unmatched %>%
  count(
    qa_priority,
    unmatched_reason,
    vector_competence_status,
    name = "row_count"
  ) %>%
  arrange(desc(qa_priority), desc(row_count), unmatched_reason)

write_csv(unmatched, review_path, na = "")
write_csv(summary_table, summary_path, na = "")

message("Wrote review file: ", review_path)
message("Wrote summary file: ", summary_path)
