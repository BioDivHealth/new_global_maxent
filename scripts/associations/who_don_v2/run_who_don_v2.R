source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_production.R"))

args <- commandArgs(trailingOnly = TRUE)
flag_args <- args[grepl("^--", args)]
valid_flags <- c("--audit-clean", "--skip-production", "--verbose-qa", "--association-mode", "--help", "-h")
unknown_args <- setdiff(flag_args, valid_flags)
if (length(unknown_args) > 0) {
  stop("Unknown arguments: ", paste(unknown_args, collapse = ", "), call. = FALSE)
}

association_mode <- "native"
association_mode_index <- which(args == "--association-mode")
if (length(association_mode_index) > 1L) {
  stop("--association-mode can be supplied at most once.", call. = FALSE)
}
if (length(association_mode_index) == 1L) {
  value_index <- association_mode_index + 1L
  if (value_index > length(args) || grepl("^--", args[[value_index]])) {
    stop("--association-mode requires a value: native or contract.", call. = FALSE)
  }
  association_mode <- args[[value_index]]
}
if (!association_mode %in% c("native", "contract")) {
  stop("--association-mode must be native or contract.", call. = FALSE)
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
      "  Rscript scripts/associations/who_don_v2/run_who_don_v2.R --association-mode contract",
      "",
      "Options:",
      "  --audit-clean       Run optional clean-vs-v2 audit after production.",
      "  --skip-production   Skip production stages; requires --audit-clean.",
      "  --verbose-qa        Write detailed stage diagnostics to qa/archive/stage_diagnostics/.",
      "  --association-mode  Association evidence mode: native (default Option A) or contract (rollback/reference).",
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
Sys.setenv(WHO_DON_V2_ASSOCIATION_MODE = association_mode)

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
  message("WHO DON v2 association mode: ", association_mode)
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
