#!/usr/bin/env Rscript
# Compatibility wrapper. Prefer the scoped script path in new docs and automation.
source(here::here(
  "scripts",
  "associations",
  "role_annotation",
  "source_check",
  "01_build_source_check_decision_ledger.R"
))
