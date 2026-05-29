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
archive_data_dir <- file.path(pathogen_association_data_dir, "archive")

vectormap_raw_dir <- file.path(source_data_dir, "vectormap", "raw")
vectormap_dir <- vectormap_raw_dir
mapveu_raw_dir <- file.path(source_data_dir, "mapveu", "raw")
mapveu_dir <- mapveu_raw_dir
vector_host_dir <- file.path(evidence_data_dir, "host_vector")
readiness_dir <- file.path(pathogen_association_data_dir, "readiness")

# Current role-annotation layout. Core evidence/QA files live under evidence/,
# manual reviews/source checks live under manual/, generated Deep Research
# prompts/reports live under staged/, and source PDFs plus extracted text live
# under source_data/.
role_annotation_dir <- file.path(evidence_data_dir, "role_annotation")
role_manual_dir <- file.path(manual_data_dir, "role_annotation")
role_reviews_dir <- file.path(role_manual_dir, "reviews")
role_deep_research_dir <- file.path(
  staged_data_dir,
  "role_annotation",
  "deep_research_inputs"
)
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

# Raw CLOVER/VIRION source exports. WHO-specific generated source outputs live
# under staged source-specific output directories.
clover_source_dir <- file.path(
  source_data_dir,
  "clover",
  "viralemergence-clover-2604d22"
)
virion_source_dir <- file.path(source_data_dir, "virion", "raw", "virion_download")
virion_source_version_dir <- file.path(virion_source_dir, "19502921")

who_clover_dir <- file.path(staged_data_dir, "clover", "outputs")
who_virion_dir <- file.path(staged_data_dir, "virion", "outputs")

# WHO diseases layout. During migration, helpers prefer the proposed lifecycle
# roots but fall back to the legacy WHO/who_diseases/ directory until files move.
who_diseases_legacy_dir <- file.path(who_data_dir, "who_diseases")
who_diseases_source_dir <- file.path(source_data_dir, "who_diseases")
who_diseases_manual_dir <- file.path(manual_data_dir, "who_diseases")
who_diseases_staged_dir <- file.path(staged_data_dir, "who_diseases")
who_diseases_evidence_dir <- file.path(evidence_data_dir, "who_diseases")
who_diseases_archive_dir <- file.path(archive_data_dir, "who_diseases")

who_diseases_regional_tables_dir <- file.path(
  who_diseases_source_dir,
  "regional_tables"
)
who_diseases_comparison_lookups_dir <- file.path(
  who_diseases_source_dir,
  "comparison_lookups"
)
who_diseases_name_resolution_dir <- file.path(
  who_diseases_manual_dir,
  "name_resolution"
)
who_diseases_transmission_rules_dir <- file.path(
  who_diseases_manual_dir,
  "transmission_rules"
)
who_diseases_pathogen_matching_manual_dir <- file.path(
  who_diseases_manual_dir,
  "pathogen_matching"
)
who_diseases_broad_taxa_manual_dir <- file.path(
  who_diseases_manual_dir,
  "broad_taxa"
)
who_diseases_staged_backbone_dir <- file.path(
  who_diseases_staged_dir,
  "backbone"
)
who_diseases_staged_master_expansion_dir <- file.path(
  who_diseases_staged_dir,
  "master_expansion"
)
who_diseases_staged_pathogen_matching_dir <- file.path(
  who_diseases_staged_dir,
  "pathogen_matching"
)
who_diseases_host_queries_dir <- file.path(
  who_diseases_staged_dir,
  "host_queries"
)
who_diseases_broad_taxa_staged_dir <- file.path(
  who_diseases_staged_dir,
  "broad_taxa"
)
who_diseases_backbone_dir <- file.path(who_diseases_evidence_dir, "backbone")
who_diseases_master_expansion_dir <- file.path(
  who_diseases_evidence_dir,
  "master_expansion"
)
who_diseases_host_species_dir <- file.path(
  who_diseases_evidence_dir,
  "host_species"
)
who_diseases_qa_dir <- file.path(who_diseases_evidence_dir, "qa")

# GenBank-simple layout. Manual query overrides live under manual/, generated
# manifests/intermediates/maps/local retrieval checkpoints live under staged/,
# and the active disease-country evidence plus QA live under evidence/.
genbank_simple_legacy_dir <- file.path(who_data_dir, "genbank_simple")
genbank_simple_evidence_dir <- file.path(evidence_data_dir, "genbank_simple")
genbank_simple_dir <- genbank_simple_evidence_dir
genbank_simple_manual_dir <- file.path(manual_data_dir, "genbank_simple")
genbank_simple_manifest_dir <- file.path(
  staged_data_dir,
  "genbank_simple",
  "manifests"
)
genbank_simple_intermediate_dir <- file.path(
  staged_data_dir,
  "genbank_simple",
  "intermediate"
)
genbank_simple_maps_dir <- file.path(staged_data_dir, "genbank_simple", "maps")
genbank_simple_standard_maps_dir <- file.path(genbank_simple_maps_dir, "standard")
genbank_simple_readiness_maps_dir <- file.path(genbank_simple_maps_dir, "readiness")
genbank_simple_local_runs_dir <- file.path(
  staged_data_dir,
  "genbank_simple",
  "local_runs"
)
genbank_simple_standard_run_dir <- file.path(
  genbank_simple_local_runs_dir,
  "pathogen_runs"
)
genbank_simple_readiness_run_dir <- file.path(
  genbank_simple_local_runs_dir,
  "pathogen_runs_readiness"
)
genbank_simple_qa_dir <- file.path(genbank_simple_evidence_dir, "qa")

# WHO DON v2 layout. The current production tree still lives under WHO/ during
# migration; proposed roots split generated fixtures/intermediates under
# staged/, durable review decisions under manual/, active evidence/final/web/QA
# under evidence/, and historical diagnostics under archive/.
who_don_v2_legacy_dir <- file.path(who_data_dir, "disease_outbreak_news_v2")
who_don_v2_staged_dir <- file.path(staged_data_dir, "who_don_v2")
who_don_v2_manual_dir <- file.path(manual_data_dir, "who_don_v2")
who_don_v2_evidence_dir <- file.path(evidence_data_dir, "who_don_v2")
who_don_v2_archive_dir <- file.path(archive_data_dir, "who_don_v2")

who_don_v2_records_dir <- file.path(who_don_v2_staged_dir, "records")
who_don_v2_reference_staged_dir <- file.path(who_don_v2_staged_dir, "reference")
who_don_v2_candidates_dir <- file.path(who_don_v2_staged_dir, "candidates")
who_don_v2_review_manual_dir <- file.path(who_don_v2_manual_dir, "review")
who_don_v2_evidence_tables_dir <- file.path(who_don_v2_evidence_dir, "evidence")
who_don_v2_final_dir <- file.path(who_don_v2_evidence_dir, "final")
who_don_v2_web_dir <- file.path(who_don_v2_evidence_dir, "web")
who_don_v2_qa_dir <- file.path(who_don_v2_evidence_dir, "qa")
who_don_v2_qa_archive_root_dir <- file.path(who_don_v2_archive_dir, "qa")

vectormap_outputs_dir <- file.path(staged_data_dir, "vectormap", "outputs")
vectormap_manual_dir <- file.path(manual_data_dir, "vectormap")

mapveu_outputs_dir <- file.path(staged_data_dir, "mapveu", "outputs")
mapveu_manual_dir <- file.path(manual_data_dir, "mapveu")

vector_host_outputs_dir <- vector_host_dir

# Vector-screening layout. Source EFSA workbooks live under source_data/,
# manual screening/crosswalk decisions under manual/, generated intermediates
# under staged/, and active vector evidence plus QA under evidence/.
vector_screening_legacy_dir <- file.path(who_data_dir, "vector_screening")
vector_screening_legacy_outputs_dir <- file.path(
  vector_screening_legacy_dir,
  "outputs"
)
vector_screening_source_dir <- file.path(source_data_dir, "vector_screening")
vector_screening_manual_dir <- file.path(manual_data_dir, "vector_screening")
vector_screening_staged_dir <- file.path(staged_data_dir, "vector_screening")
vector_screening_staged_outputs_dir <- file.path(
  vector_screening_staged_dir,
  "outputs"
)
vector_screening_evidence_dir <- file.path(evidence_data_dir, "vector_screening")
vector_screening_qa_dir <- file.path(vector_screening_evidence_dir, "qa")

vector_screening_efsa_source_dir <- file.path(
  vector_screening_source_dir,
  "efsa",
  "raw"
)
vector_screening_efsa_manual_dir <- file.path(
  vector_screening_manual_dir,
  "efsa"
)
vector_screening_efsa_outputs_dir <- file.path(
  vector_screening_staged_dir,
  "efsa",
  "outputs"
)
vector_screening_taxonomy_manual_dir <- file.path(
  vector_screening_manual_dir,
  "taxonomy"
)
vector_screening_taxonomy_review_dir <- file.path(
  vector_screening_staged_dir,
  "taxonomy_review"
)
vector_screening_vectraits_dir <- file.path(
  vector_screening_staged_dir,
  "vectraits"
)

vector_screening_legacy_input_dir <- file.path(
  vector_screening_legacy_dir,
  "inputs"
)
vector_screening_legacy_efsa_input_dir <- file.path(
  vector_screening_legacy_dir,
  "efsa",
  "inputs"
)
vector_screening_legacy_efsa_manual_dir <- file.path(
  vector_screening_legacy_dir,
  "efsa",
  "manual"
)
vector_screening_legacy_efsa_outputs_dir <- file.path(
  vector_screening_legacy_dir,
  "efsa",
  "outputs"
)
vector_screening_legacy_taxonomy_review_dir <- file.path(
  vector_screening_legacy_dir,
  "taxonomy_review"
)

prefer_existing_path <- function(primary, fallback) {
  if (file.exists(primary) || dir.exists(primary) ||
      (!file.exists(fallback) && !dir.exists(fallback))) {
    return(primary)
  }

  fallback
}

vector_screening_manual_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_manual_dir, filename),
    file.path(vector_screening_legacy_input_dir, filename)
  )
}

vector_screening_efsa_source_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_efsa_source_dir, filename),
    file.path(vector_screening_legacy_efsa_input_dir, filename)
  )
}

vector_screening_efsa_manual_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_efsa_manual_dir, filename),
    file.path(vector_screening_legacy_efsa_manual_dir, filename)
  )
}

vector_screening_efsa_staged_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_efsa_outputs_dir, filename),
    file.path(vector_screening_legacy_efsa_outputs_dir, filename)
  )
}

vector_screening_staged_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_staged_outputs_dir, filename),
    file.path(vector_screening_legacy_outputs_dir, filename)
  )
}

vector_screening_evidence_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_evidence_dir, filename),
    file.path(vector_screening_legacy_outputs_dir, filename)
  )
}

vector_screening_qa_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_qa_dir, filename),
    file.path(vector_screening_legacy_outputs_dir, filename)
  )
}

vector_screening_taxonomy_manual_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_taxonomy_manual_dir, filename),
    file.path(vector_screening_legacy_taxonomy_review_dir, filename)
  )
}

vector_screening_taxonomy_review_path <- function(filename) {
  prefer_existing_path(
    file.path(vector_screening_taxonomy_review_dir, filename),
    file.path(vector_screening_legacy_taxonomy_review_dir, filename)
  )
}

who_diseases_path <- function(filename) {
  file.path(who_diseases_legacy_dir, filename)
}

who_diseases_source_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_source_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_regional_table_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_regional_tables_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_regional_table_paths <- function(regions = c(
  "africa",
  "americas",
  "europe",
  "mediterranean",
  "se_asia",
  "western_pacific"
)) {
  vapply(
    paste0(regions, "_table.csv"),
    who_diseases_regional_table_path,
    character(1),
    USE.NAMES = FALSE
  )
}

who_diseases_comparison_lookup_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_comparison_lookups_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_name_resolution_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_name_resolution_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_transmission_rules_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_transmission_rules_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_pathogen_matching_manual_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_pathogen_matching_manual_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_broad_taxa_manual_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_broad_taxa_manual_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_staged_backbone_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_staged_backbone_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_staged_master_expansion_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_staged_master_expansion_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_staged_pathogen_matching_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_staged_pathogen_matching_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_host_query_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_host_queries_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_broad_taxa_staged_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_broad_taxa_staged_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_backbone_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_backbone_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_master_expansion_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_master_expansion_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_host_species_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_host_species_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_qa_path <- function(filename) {
  prefer_existing_path(
    file.path(who_diseases_qa_dir, filename),
    who_diseases_path(filename)
  )
}

who_diseases_translation_path <- function() {
  who_diseases_name_resolution_path("translation.csv")
}

who_disease_names_path <- function() {
  who_diseases_name_resolution_path("disease_names.csv")
}

who_diseases_gibb_lookup_path <- function() {
  who_diseases_comparison_lookup_path("diseases_in_gibb_etal.csv")
}

who_final_pathogen_data_path <- function() {
  who_diseases_staged_backbone_path("final_pathogen_data.csv")
}

who_pathogens_diseases_zoonotic_path <- function() {
  who_diseases_backbone_path("who_pathogens_diseases_zoonotic.csv")
}

who_pathogen_analysis_units_path <- function() {
  who_diseases_backbone_path("who_pathogen_analysis_units.csv")
}

who_pathogen_analysis_units_keep_path <- function() {
  who_diseases_backbone_path("who_pathogen_analysis_units_keep.csv")
}

who_master_disease_analysis_units_path <- function() {
  who_diseases_master_expansion_path("master_disease_analysis_units.csv")
}

who_master_plus_analysis_units_path <- function() {
  who_diseases_master_expansion_path("master_plus_who_analysis_units.csv")
}

who_master_pathogen_host_species_path <- function() {
  who_diseases_host_species_path("master_pathogen_host_species.csv")
}

who_master_pathogen_host_species_clean_path <- function() {
  who_diseases_host_species_path("master_pathogen_host_species_clean.csv")
}

who_master_pathogen_host_species_summary_path <- function() {
  who_diseases_qa_path("master_pathogen_host_species_summary.csv")
}

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
  who_diseases_backbone_path("who_pathogens_diseases.csv")
}

who_zoonotic_pathogens_path <- function() {
  who_pathogen_analysis_units_keep_path()
}

who_working_pathogens_path <- function(scope = c("zoonotic", "raw")) {
  scope <- match.arg(scope)

  switch(
    scope,
    zoonotic = who_zoonotic_pathogens_path(),
    raw = who_raw_pathogens_path()
  )
}
