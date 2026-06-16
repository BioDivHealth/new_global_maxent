#!/usr/bin/env Rscript
# Compatibility wrapper. Prefer the scoped script path in new docs and automation.
source(here::here(
  "scripts",
  "associations",
  "role_annotation",
  "source_check",
  "02_import_source_checked_role_rows.R"
))
