# -----------------------------------------------------------------------------|
# 4_CombineNetworks.R ----
# -----------------------------------------------------------------------------|
# Purpose: Combine CLOVER and VIRION WHO source-component host networks into the
#          active WHO host-pathogen backbone.
# Inputs : clover_who_network.csv and virion_who_network.csv
# Outputs: combined_who_network.csv
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# 1. Load required libraries and path helpers ----
# -----------------------------------------------------------------------------|
library(pacman)
p_load(here, tidyverse)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "legacy_who_compatibility_helpers.R"
))

source_component_columns <- c(
  "Pathogen", "Host_clean", "PathogenClass", "PathogenOrder",
  "PathogenFamily", "PathogenGenus", "HostPhylum", "HostClass",
  "HostFamily", "HostOrder", "DetectionMethod", "MainSource"
)

# -----------------------------------------------------------------------------|
# 3. Load source-component networks ----
# -----------------------------------------------------------------------------|
clover_network <- read_csv(who_network_source_component_path("clover_who_network.csv"))
virion_network <- read_csv(who_network_source_component_path("virion_who_network.csv"))

legacy_who_require_columns(
  clover_network,
  source_component_columns,
  "CLOVER source component"
)
legacy_who_require_columns(
  virion_network,
  source_component_columns,
  "VIRION source component"
)

# -----------------------------------------------------------------------------|
# 4. Harmonize schemas and bind sources ----
# -----------------------------------------------------------------------------|
clover_network$PathogenType <- "bacteria"
virion_network$PathogenType <- "virus"

clover_network <- clover_network %>%
  mutate(
    in_gibb_etal = if ("in_gibb_etal" %in% names(.)) in_gibb_etal else NA,
    in_empres_i = if ("in_empres_i" %in% names(.)) in_empres_i else NA,
    high_quality_detection = if ("high_quality_detection" %in% names(.)) high_quality_detection else DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing"),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) downstream_default_include else high_quality_detection,
    downstream_review_reason = if ("downstream_review_reason" %in% names(.)) downstream_review_reason else NA_character_
  )

virion_network <- virion_network %>%
  mutate(
    high_quality_detection = if ("high_quality_detection" %in% names(.)) high_quality_detection else DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing"),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) downstream_default_include else high_quality_detection,
    downstream_review_reason = if ("downstream_review_reason" %in% names(.)) downstream_review_reason else NA_character_
  )

clover_network <- clover_network %>%
  select(-any_of("ID"))

combined_network <- bind_rows(clover_network, virion_network)

# -----------------------------------------------------------------------------|
# 5. Standardize output fields and write backbone ----
# -----------------------------------------------------------------------------|
combined_network <- combined_network %>%
  rename(Host = Host_clean) %>%
  mutate(
    PathogenClass = tolower(PathogenClass),
    PathogenOrder = tolower(PathogenOrder),
    PathogenFamily = tolower(PathogenFamily),
    PathogenGenus = tolower(PathogenGenus),
    HostPhylum = tolower(HostPhylum),
    HostClass = tolower(HostClass),
    HostFamily = tolower(HostFamily),
    HostOrder = tolower(HostOrder)
  )

output_path <- who_raw_network_path()
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(combined_network, output_path)
