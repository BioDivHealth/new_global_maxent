# -----------------------------------------------------------------------------|
# network_wrapper_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared orchestration helpers for network-building wrapper entrypoints.
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# 1. Package and path setup ----
# -----------------------------------------------------------------------------|

require_network_wrapper_packages <- function(packages = c("here", "readr")) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Package `", package, "` is required.", call. = FALSE)
    }
  }
}

source_network_working_inputs <- function() {
  source(here::here("scripts", "associations", "working_inputs.R"))
}

network_building_script <- function(filename) {
  normalizePath(
    here::here("scripts", "associations", "network_building", filename),
    mustWork = TRUE
  )
}

# -----------------------------------------------------------------------------|
# 2. Command-line arguments ----
# -----------------------------------------------------------------------------|

parse_network_wrapper_args <- function(
  valid_flags,
  help_text,
  args = commandArgs(trailingOnly = TRUE)
) {
  unknown_args <- setdiff(args, valid_flags)
  if (length(unknown_args) > 0) {
    stop("Unknown arguments: ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }

  if (any(args %in% c("--help", "-h"))) {
    cat(help_text, "\n")
    quit(status = 0)
  }

  args
}

# -----------------------------------------------------------------------------|
# 3. Stage execution ----
# -----------------------------------------------------------------------------|

run_stage <- function(stage_file, running_label, failure_label = running_label) {
  stage_path <- network_building_script(stage_file)
  rscript <- normalizePath(file.path(R.home("bin"), "Rscript"), mustWork = TRUE)

  cat("Running ", running_label, " stage: ", stage_file, "\n", sep = "")
  status <- system2(rscript, stage_path)
  if (!identical(status, 0L)) {
    stop(failure_label, " stage failed: ", stage_file, call. = FALSE)
  }
  cat("Completed ", running_label, " stage: ", stage_file, "\n", sep = "")
}

# -----------------------------------------------------------------------------|
# 4. Output summaries ----
# -----------------------------------------------------------------------------|

summarize_csv_output <- function(name, path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Expected output is missing: ", path, call. = FALSE)
    }
    return(NULL)
  }

  data <- readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
  data.frame(
    output = name,
    rows = nrow(data),
    columns = ncol(data),
    path = path,
    check.names = FALSE
  )
}

summarize_text_output <- function(name, path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Expected output is missing: ", path, call. = FALSE)
    }
    return(NULL)
  }

  data.frame(
    output = name,
    rows = length(readLines(path, warn = FALSE)),
    columns = NA_integer_,
    path = path,
    check.names = FALSE
  )
}

summarize_csv_outputs <- function(outputs, required = TRUE) {
  summaries <- Map(
    summarize_csv_output,
    names(outputs),
    unname(outputs),
    MoreArgs = list(required = required)
  )
  summaries <- Filter(Negate(is.null), summaries)
  if (length(summaries) == 0) {
    return(NULL)
  }
  do.call(rbind, summaries)
}

summarize_text_outputs <- function(outputs, required = FALSE) {
  summaries <- Map(
    summarize_text_output,
    names(outputs),
    unname(outputs),
    MoreArgs = list(required = required)
  )
  summaries <- Filter(Negate(is.null), summaries)
  if (length(summaries) == 0) {
    return(NULL)
  }
  do.call(rbind, summaries)
}
