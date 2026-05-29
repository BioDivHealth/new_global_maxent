library(here)

source(here::here("scripts", "associations", "working_inputs.R"))

who_don_v2_scripts_dir <- function(...) {
  here::here("scripts", "associations", "who_don_v2", ...)
}

who_don_v2_rules_dir <- function(...) {
  who_don_v2_scripts_dir("rules", ...)
}

who_don_v2_primary_layer_dir <- function(layer = NULL) {
  if (is.null(layer) || length(layer) == 0 || layer == "") {
    return(who_don_v2_evidence_dir)
  }

  switch(
    layer,
    records = who_don_v2_records_dir,
    reference = who_don_v2_reference_staged_dir,
    candidates = who_don_v2_candidates_dir,
    review = who_don_v2_review_manual_dir,
    evidence = who_don_v2_evidence_tables_dir,
    final = who_don_v2_final_dir,
    web = who_don_v2_web_dir,
    qa = who_don_v2_qa_dir,
    archive = who_don_v2_archive_dir,
    file.path(who_don_v2_evidence_dir, layer)
  )
}

who_don_v2_primary_path <- function(...) {
  parts <- c(...)

  if (length(parts) == 0) {
    return(who_don_v2_evidence_dir)
  }

  if (length(parts) == 1) {
    return(who_don_v2_primary_layer_dir(parts[[1]]))
  }

  file.path(who_don_v2_primary_layer_dir(parts[[1]]), parts[-1])
}

who_don_v2_legacy_path <- function(...) {
  file.path(who_don_v2_legacy_dir, ...)
}

who_don_v2_output_dir <- function(...) {
  prefer_existing_path(
    who_don_v2_primary_path(...),
    who_don_v2_legacy_path(...)
  )
}

who_don_v2_reference_dir <- function(...) {
  who_don_v2_output_dir("reference", ...)
}

who_don_v2_qa_archive_dir <- function(...) {
  prefer_existing_path(
    file.path(who_don_v2_qa_archive_root_dir, ...),
    who_don_v2_legacy_path("qa", "archive", ...)
  )
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
    who_don_v2_staged_dir,
    who_don_v2_manual_dir,
    who_don_v2_evidence_dir,
    who_don_v2_archive_dir,
    who_don_v2_records_dir,
    who_don_v2_reference_staged_dir,
    who_don_v2_candidates_dir,
    who_don_v2_review_manual_dir,
    who_don_v2_evidence_tables_dir,
    who_don_v2_final_dir,
    who_don_v2_web_dir,
    who_don_v2_qa_dir,
    file.path(who_don_v2_qa_archive_root_dir, "orphaned"),
    file.path(who_don_v2_qa_archive_root_dir, "clean_audit"),
    file.path(who_don_v2_qa_archive_root_dir, "completed_review_surfaces"),
    file.path(who_don_v2_qa_archive_root_dir, "stage_diagnostics")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}
