#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 01_prepare_gbif_occurrences.R ----
# -----------------------------------------------------------------------------|
# Compatibility wrapper. The one-species GBIF occurrence worker now lives at:
#   scripts/sdms/present/occurrences/01_prepare_one_gbif_species.R
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(
  here::here(),
  "scripts",
  "sdms",
  "present",
  "occurrences",
  "01_prepare_one_gbif_species.R"
))
