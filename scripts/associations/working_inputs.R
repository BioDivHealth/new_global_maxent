# ------------------------------------------------------------------------------
# working_inputs.R
# ------------------------------------------------------------------------------
# Purpose: Centralize the default data roots and WHO working-input paths used
#          by the downstream pathogen-association pipeline.
#
# Notes  : Keep the raw network-building stages pointed at the original WHO
#          source tables. Downstream scripts should default to the canonical
#          zoonotic working layer unless they explicitly need a broader input.
# ------------------------------------------------------------------------------

pathogen_association_data_dir <- here::here("pathogen_association_data")
who_data_dir <- file.path(pathogen_association_data_dir, "WHO")
source_data_dir <- file.path(pathogen_association_data_dir, "source_data")
manual_data_dir <- file.path(pathogen_association_data_dir, "manual")
staged_data_dir <- file.path(pathogen_association_data_dir, "staged")
evidence_data_dir <- file.path(pathogen_association_data_dir, "evidence")

vectormap_raw_dir <- file.path(source_data_dir, "vectormap", "raw")
vectormap_dir <- vectormap_raw_dir
mapveu_raw_dir <- file.path(source_data_dir, "mapveu", "raw")
mapveu_dir <- mapveu_raw_dir
vector_host_dir <- file.path(evidence_data_dir, "host_vector")
readiness_dir <- file.path(pathogen_association_data_dir, "readiness")

# Current role-annotation layout. Core evidence/QA files still live under the
# existing WHO root; manual reviews/source checks live under manual/, and source
# PDFs plus extracted text live under source_data/.
role_annotation_dir <- file.path(who_data_dir, "role_annotation")
role_manual_dir <- file.path(manual_data_dir, "role_annotation")
role_reviews_dir <- file.path(role_manual_dir, "reviews")
role_deep_research_dir <- file.path(role_annotation_dir, "deep_research_inputs")
role_deep_research_consolidated_dir <- file.path(
  role_deep_research_dir,
  "consolidated"
)
role_source_check_dir <- file.path(role_manual_dir, "source_check")
role_source_check_import_dir <- file.path(role_source_check_dir, "import")
role_source_pdf_dir <- file.path(source_data_dir, "role_annotation", "papers")
role_source_pdf_text_dir <- file.path(
  source_data_dir,
  "role_annotation",
  "pdf_text"
)
role_candidates_dir <- role_annotation_dir
role_evidence_dir <- role_annotation_dir
role_roster_dir <- role_annotation_dir
role_qa_dir <- file.path(role_annotation_dir, "qa")

# Raw CLOVER checkout/vendor export. WHO-specific generated CLOVER outputs live
# under `who_clover_dir`.
clover_source_dir <- file.path(
  pathogen_association_data_dir,
  "viralemergence-clover-2604d22"
)
who_clover_dir <- file.path(who_data_dir, "clover")

vectormap_outputs_dir <- file.path(staged_data_dir, "vectormap", "outputs")
vectormap_manual_dir <- file.path(manual_data_dir, "vectormap")

mapveu_outputs_dir <- file.path(staged_data_dir, "mapveu", "outputs")
mapveu_manual_dir <- file.path(manual_data_dir, "mapveu")

vector_host_outputs_dir <- vector_host_dir

who_raw_network_path <- function() {
  file.path(
    who_data_dir,
    "networks",
    "combined_who_network.csv"
  )
}

who_canonical_network_path <- function() {
  file.path(
    who_data_dir,
    "networks",
    "combined_who_network_canonical.csv"
  )
}

who_canonical_zoonotic_network_path <- function() {
  file.path(
    who_data_dir,
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
  file.path(
    who_data_dir,
    "who_diseases",
    "who_pathogens_diseases.csv"
  )
}

who_zoonotic_pathogens_path <- function() {
  file.path(
    who_data_dir,
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
