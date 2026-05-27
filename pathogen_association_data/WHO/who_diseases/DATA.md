# WHO Disease Data

This folder holds the source and derived disease/pathogen tables used to build the WHO pathogen backbone for downstream host, vector, GenBank, and WHO DON workflows.

## WHO Source Tables

- `africa_table.csv` - WHO Africa regional priority/prototype pathogen table used by `scripts/associations/network_building/1_WHO_Diseases.R`.
- `americas_table.csv` - WHO Americas regional priority/prototype pathogen table used by `1_WHO_Diseases.R`.
- `europe_table.csv` - WHO Europe regional priority/prototype pathogen table used by `1_WHO_Diseases.R`.
- `mediterranean_table.csv` - WHO Eastern Mediterranean regional priority/prototype pathogen table used by `1_WHO_Diseases.R`.
- `se_asia_table.csv` - WHO South-East Asia regional priority/prototype pathogen table used by `1_WHO_Diseases.R`.
- `western_pacific_table.csv` - WHO Western Pacific regional priority/prototype pathogen table used by `1_WHO_Diseases.R`.
- `translation.csv` - Lookup table mapping older/common pathogen names to MSL39 viral species names; used during WHO pathogen name standardization.
- `disease_names.csv` - Manual pathogen-to-disease lookup used to attach disease names to standardized WHO pathogen rows and to seed WHO DON disease aliases.
- `diseases_in_gibb_etal.csv` - Project comparison/provenance lookup marking disease analysis units represented in Gibb et al. and/or EMPRES-i.

## Canonical WHO Backbone

- `final_pathogen_data.csv` - Intermediate standardized pathogen table written by `1_WHO_Diseases.R` before disease-name joins and comparison flags.
- `who_pathogens_diseases.csv` - Broad consolidated WHO pathogen-disease table; preserves priority/prototype status and per-region WHO provenance.
- `who_pathogens_diseases_zoonotic.csv` - Conservative zoonotic-focused subset derived from `who_pathogens_diseases.csv` by `1_1_WHO_Diseases_Zoonotic_Filter.R`; default downstream working layer.
- `who_pathogens_diseases_zoonotic.xlsx` - Excel copy of the zoonotic WHO table for manual review or sharing.
- `who_pathogen_analysis_units.csv` - Curated analysis-unit scaffold derived from the zoonotic WHO table, including broad taxa that may need splitting or review.
- `who_pathogen_analysis_units_keep.csv` - Retained modelling/query scope from the analysis-unit scaffold; used by VIRION, GenBank, and related downstream scripts.

## Disease Master List Expansion

- `master_disease_analysis_units.csv` - Additive merge of `dr/disease_master_list_v2.xlsx` with the WHO analysis-unit table; keeps master-list metadata without replacing the WHO backbone.
- `master_disease_name_resolution_review.csv` - Generated review queue for master-list diseases not already resolved to an existing WHO analysis unit.
- `master_disease_name_resolution_manual.csv` - Editable/manual resolution surface for unresolved master-list rows, including pathogen, disease, inclusion, split, and note fields.
- `master_who_pathogen_bridge.csv` - One-row-per-master-list bridge table combining master-list fields, WHO analysis-unit context, manual resolution, and VIRION/CLOVER match outputs.
- `master_plus_who_analysis_units.csv` - Compact combined analysis-unit table that lets resolved master-list rows sit alongside existing WHO analysis units.
- `master_plus_who_transmission_rules_manual.csv` - Manual transmission-rule curation scaffold used as optional input when building the combined master-plus-WHO units.
- `master_plus_who_transmission_rules_manual_completed.csv` - Later filled version of the transmission-rule scaffold retained as a manual curation snapshot.
- `master_plus_who_transmission_rules_manual_reviewed_v2.csv` - Reviewed v2 transmission-rule snapshot; use only with awareness of its manual-review status.

## Pathogen Matching And Host Evidence

- `master_pathogen_aliases.csv` - Manual aliases for matching resolved master-list pathogen names to local VIRION and CLOVER taxonomy labels.
- `master_pathogen_virion_clover_candidates.csv` - Candidate VIRION/CLOVER taxonomy matches for active resolved master-list analysis units.
- `master_pathogen_virion_clover_matches.csv` - Best source-prioritized VIRION/CLOVER matches used to build host query units.
- `master_pathogen_external_taxonomy_review.csv` - Review table for pathogen names requiring external taxonomy support when local VIRION/CLOVER evidence is absent or incomplete.
- `master_pathogen_host_query_units.csv` - Query-ready units specifying which source, pathogen names, and taxids should be used for host retrieval.
- `master_pathogen_host_species.csv` - Raw host-species associations retrieved from VIRION/CLOVER for the master-list query units.
- `master_pathogen_host_species_review.csv` - Review-bucket host-species associations (for example shared-species proxy and match-review rows) kept separate from default-clean host evidence.
- `master_pathogen_host_species_clean.csv` - QA/harmonized version of the host-species table with host-name standardization and downstream-readiness flags.
- `master_pathogen_host_species_summary.csv` - Summary counts and QA totals for the master-list host-species extraction.

## Broad-Taxa Candidate Strains And NCBI Metadata

- `who_broad_taxa_candidate_strains_seed.csv` - Manual seed inventory for candidate strains/exemplar viruses under broad WHO taxa.
- `who_broad_taxa_candidate_strains.csv` - Generated curation inventory of ICTV-supported candidate strains and examples for broad taxa under review.
- `who_broad_taxa_candidate_host_overrides.csv` - Manual host-name override notes for candidate-strain NCBI records where Datasets metadata is too broad or outdated.
- `who_broad_taxa_candidate_strains_ncbi_resolution.csv` - Accession-base to resolved accession lookup from the NCBI Datasets metadata pull.
- `who_broad_taxa_candidate_strains_ncbi_metadata.csv` - Parsed NCBI Datasets metadata for resolved candidate-strain accessions.
- `who_broad_taxa_candidate_strains_ncbi_metadata.tsv` - Raw or near-raw tabular NCBI Datasets metadata export kept for provenance.
- `who_broad_taxa_candidate_strains_ncbi_raw.jsonl` - Raw JSONL responses from the NCBI Datasets CLI.
- `who_broad_taxa_candidate_strains_ncbi_enriched.csv` - Candidate-strain table joined to NCBI metadata and host/taxonomy fields for review.
- `who_broad_taxa_candidate_strains_ncbi_enriched_slim.csv` - Slim review/export version of the enriched NCBI candidate-strain table.
- `who_broad_taxa_candidate_strains_ncbi_enriched_slim.xlsx` - Excel copy of the slim enriched candidate-strain table for manual review or sharing.

## Local Metadata

- `.DS_Store` - macOS Finder metadata; not an analysis input and safe to ignore.
