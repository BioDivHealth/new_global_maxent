#!/usr/bin/env Rscript
# Compatibility wrapper. Prefer the scoped script path in new docs and automation.
source(here::here(
  "scripts",
  "associations",
  "role_annotation",
  "roster",
  "01_build_host_role_candidates.R"
))
