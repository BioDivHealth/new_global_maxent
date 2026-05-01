# ------------------------------------------------------------------------------
# working_inputs.R
# ------------------------------------------------------------------------------
# Purpose: Centralize the default WHO working-input paths used by the
#          downstream pathogen-association pipeline.
#
# Notes  : Keep the raw network-building stages pointed at the original WHO
#          source tables. Downstream scripts should default to the canonical
#          zoonotic working layer unless they explicitly need a broader input.
# ------------------------------------------------------------------------------

who_raw_network_path <- function() {
  here::here(
    "pathogen_association_data",
    "WHO",
    "networks",
    "combined_who_network.csv"
  )
}

who_canonical_network_path <- function() {
  here::here(
    "pathogen_association_data",
    "WHO",
    "networks",
    "combined_who_network_canonical.csv"
  )
}

who_canonical_zoonotic_network_path <- function() {
  here::here(
    "pathogen_association_data",
    "WHO",
    "networks",
    "combined_who_network_canonical_zoonotic.csv"
  )
}

who_working_network_path <- function(scope = c("zoonotic", "canonical", "raw")) {
  scope <- match.arg(scope)

  switch(
    scope,
    zoonotic = who_canonical_zoonotic_network_path(),
    canonical = who_canonical_network_path(),
    raw = who_raw_network_path()
  )
}

who_raw_pathogens_path <- function() {
  here::here(
    "pathogen_association_data",
    "WHO",
    "who_diseases",
    "who_pathogens_diseases.csv"
  )
}

who_zoonotic_pathogens_path <- function() {
  here::here(
    "pathogen_association_data",
    "WHO",
    "who_diseases",
    "who_pathogen_analysis_units_keep.csv"
  )
}

who_working_pathogens_path <- function(scope = c("zoonotic", "raw")) {
  scope <- match.arg(scope)

  switch(
    scope,
    zoonotic = who_zoonotic_pathogens_path(),
    raw = who_raw_pathogens_path()
  )
}
