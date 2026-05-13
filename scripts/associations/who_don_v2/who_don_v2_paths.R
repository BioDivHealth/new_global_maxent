library(here)

who_don_v2_scripts_dir <- function(...) {
  here::here("scripts", "associations", "who_don_v2", ...)
}

who_don_v2_rules_dir <- function(...) {
  who_don_v2_scripts_dir("rules", ...)
}

who_don_v2_output_dir <- function(...) {
  here::here("pathogen_association_data", "WHO", "disease_outbreak_news_v2", ...)
}

who_don_v2_reference_dir <- function(...) {
  who_don_v2_output_dir("reference", ...)
}

who_don_v2_qa_archive_dir <- function(...) {
  who_don_v2_output_dir("qa", "archive", ...)
}

who_don_pre_v2_archive_dir <- function(...) {
  here::here("archive", "who_don_pre_v2", ...)
}

who_don_clean_output_dir <- function(...) {
  active_dir <- here::here("pathogen_association_data", "WHO", "disease_outbreak_news_clean")
  archive_dir <- who_don_pre_v2_archive_dir("data", "disease_outbreak_news_clean")
  clean_dir <- if (dir.exists(active_dir)) active_dir else archive_dir
  file.path(clean_dir, ...)
}

who_don_v2_ensure_dirs <- function() {
  dirs <- c(
    who_don_v2_output_dir(),
    who_don_v2_output_dir("records"),
    who_don_v2_reference_dir(),
    who_don_v2_output_dir("candidates"),
    who_don_v2_output_dir("evidence"),
    who_don_v2_output_dir("review"),
    who_don_v2_output_dir("final"),
    who_don_v2_output_dir("web"),
    who_don_v2_output_dir("qa"),
    who_don_v2_qa_archive_dir("orphaned"),
    who_don_v2_qa_archive_dir("clean_audit"),
    who_don_v2_qa_archive_dir("completed_review_surfaces"),
    who_don_v2_qa_archive_dir("stage_diagnostics")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}
