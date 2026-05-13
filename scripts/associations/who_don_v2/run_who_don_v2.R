source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_production.R"))

args <- commandArgs(trailingOnly = TRUE)
valid_args <- c("--audit-clean", "--skip-production", "--verbose-qa", "--help", "-h")
unknown_args <- setdiff(args, valid_args)
if (length(unknown_args) > 0) {
  stop("Unknown arguments: ", paste(unknown_args, collapse = ", "), call. = FALSE)
}

if (any(args %in% c("--help", "-h"))) {
  cat(
    paste(
      "WHO DON v2 runner",
      "",
      "Usage:",
      "  Rscript scripts/associations/who_don_v2/run_who_don_v2.R",
      "  Rscript scripts/associations/who_don_v2/run_who_don_v2.R --audit-clean",
      "  Rscript scripts/associations/who_don_v2/run_who_don_v2.R --skip-production --audit-clean",
      "  Rscript scripts/associations/who_don_v2/run_who_don_v2.R --verbose-qa",
      "",
      "Options:",
      "  --audit-clean       Run optional clean-vs-v2 audit after production.",
      "  --skip-production   Skip production stages; requires --audit-clean.",
      "  --verbose-qa        Write detailed stage diagnostics to qa/archive/stage_diagnostics/.",
      sep = "\n"
    ),
    "\n"
  )
  quit(status = 0)
}

audit_clean <- "--audit-clean" %in% args
skip_production <- "--skip-production" %in% args
verbose_qa <- "--verbose-qa" %in% args
if (skip_production && !audit_clean) {
  stop("--skip-production requires --audit-clean.", call. = FALSE)
}
if (verbose_qa) {
  Sys.setenv(WHO_DON_V2_VERBOSE_QA = "1")
}

production_stages <- c(
  "01_records.R",
  "02e_prepare_country_rules.R",
  "02f_extract_country_candidates_native.R",
  "02g_compare_country_candidates.R",
  "02b_prepare_disease_rules.R",
  "02c_extract_disease_candidates_native.R",
  "02d_compare_disease_candidates.R",
  "03_build_association_evidence.R",
  "04_classify_scope.R",
  "04b_write_policy_review_manifest.R",
  "05_export_final.R",
  "06_export_web.R"
)

run_stage <- function(stage_file) {
  stage_path <- normalizePath(who_don_v2_scripts_dir(stage_file), mustWork = TRUE)
  rscript <- normalizePath(file.path(R.home("bin"), "Rscript"), mustWork = TRUE)
  message("Running WHO DON v2 stage: ", stage_file)
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop("WHO DON v2 stage failed: ", stage_file, call. = FALSE)
  }
}

if (!skip_production) {
  lapply(production_stages, run_stage)
  v2_validate_production_outputs()
  manifest <- v2_write_output_manifest(v2_production_output_specs())
  message(
    "WHO DON v2 production complete: ",
    nrow(manifest),
    " manifest rows written to ",
    who_don_v2_output_dir("final", "who_don_v2_output_manifest.csv")
  )
}

if (audit_clean) {
  run_stage("06_compare_to_clean.R")
  message("WHO DON v2 clean-vs-v2 audit complete.")
}
