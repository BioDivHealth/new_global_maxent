## Pathogen Associations: Project Goals and Results

This directory now acts as the top-level home for the WHO pathogen workflow,
with scripts grouped by workflow block rather than kept in a single flat list.
The pipeline integrates, standardizes, and analyzes global pathogen-host
association data, with a focus on WHO priority pathogens.

### Folder Layout

- `network_building/`
  Builds the WHO disease list, matches pathogens against CLOVER and VIRION,
  cleans host taxonomy, and assembles the combined WHO host-pathogen backbone.
- `vector_screening/`
  Builds the disease/pathogen-vector workflow, including EFSA crosswalks,
  canonical disease-vector tables, taxonomy cleanup, and pathogen-vector
  backfilling. Source workbooks, manual review inputs, staged outputs, and
  active evidence are split under `pathogen_association_data/source_data/`,
  `manual/`, `staged/`, and `evidence/`.
- `host_vector_integration/`
  Joins the WHO disease/pathogen network to observational host-vector evidence
  and writes disease-level, pathogen-level, and expanded host-vector outputs.
- `role_annotation/`
  Builds conservative host/vector role candidate scaffolds and keeps final
  biological role review separate from the core network evidence tables.
- `genbank_simple/`
  Builds and runs the current GenBank-simple country-evidence workflow,
  including the expanded readiness manifest and readiness-combined summaries.

### Network Building Entrypoints

Use these wrapper scripts as the normal runnable surface for
`network_building/`. The older `1_*`, `2_*`, `3_*`, and `4_*` files remain
stage scripts called by these entrypoints.

```sh
Rscript scripts/associations/network_building/01_build_disease_scope_and_analysis_units.R
Rscript scripts/associations/network_building/02_build_master_plus_registry.R
Rscript scripts/associations/network_building/03_build_master_plus_host_network.R
Rscript scripts/associations/network_building/04_build_legacy_who_compatibility_outputs.R
Rscript scripts/associations/network_building/05_build_broad_taxa_support.R
```

- `01_build_disease_scope_and_analysis_units.R`: builds the WHO disease scope,
  zoonotic subset, analysis-unit tables, and disease-master scaffold.
- `02_build_master_plus_registry.R`: builds master-plus source matches,
  registry rows, and host-query units.
- `03_build_master_plus_host_network.R`: builds master-plus host evidence,
  host QA, and the downstream-ready master-plus host network.
- `04_build_legacy_who_compatibility_outputs.R`: rebuilds the WHO-only
  CLOVER/VIRION compatibility outputs retained for legacy provenance, contract
  checks, and a small number of explicit legacy QA scripts. By default it reuses
  existing standardized host-taxonomy CSVs; run with
  `--refresh-host-taxonomy` only when deliberately refreshing the external or
  cache-sensitive CLOVER/VIRION host-taxonomy stages.
- `05_build_broad_taxa_support.R`: rebuilds broad-taxa candidate-strain
  support outputs. By default it summarizes existing NCBI metadata files
  without refreshing them; run with `--refresh-ncbi-metadata` only when
  deliberately refreshing the external NCBI Datasets metadata stage. The older
  shell metadata workflow is retained under
  `network_building/stages/broad_taxa_support/deprecated/` for provenance only.

### Project Goals

- **Integrate Data:** Combine WHO priority pathogen lists with host-pathogen association data from the CLOVER and VIRION databases.
- **Standardize Taxonomy:** Harmonize pathogen and host names using fuzzy matching, manual curation, and taxonomic databases (GBIF, ITIS, IUCN).
- **Map Associations:** Identify which WHO pathogens are present in global databases and extract their host associations.
- **Network Analysis:** Visualize and analyze the structure of pathogen-host networks, including risk stratification, centrality, modularity, and bridge species.
- **Support Research:** Provide processed datasets and visualizations to guide research on zoonotic risk, host diversity, and sampling bias.

### Network Building Stage Notes

The stage scripts below document the internal pieces that the wrapper
entrypoints call. Run the wrappers above for routine rebuilds.

1.  **WHO Pathogen Data Processing (`network_building/stages/disease_scope/1_WHO_Diseases.R`):**
    -   Loads and standardizes WHO priority pathogen lists from various regional documents.
    -   Maps pathogens to standardized names using a translation table and fuzzy matching.
    -   Assigns risk categories and family information.
    -   Preserves whether each final pathogen is priority, prototype, or both, plus per-region WHO source status (`priority`, `prototype`, `both`, or `none`) across the WHO regions.
    -   Outputs a cleaned and consolidated pathogen list (`who_pathogens_diseases.csv`).

2.  **CLOVER Integration (Bacteria) (`network_building/stages/legacy_who_compatibility/2_1_CLOVER.R`, `network_building/stages/legacy_who_compatibility/2_2_CLOVER_Host_Clean.R`, `network_building/stages/legacy_who_compatibility/2_3_CLOVER_Network.R`):**
    -   `2_1_CLOVER.R`: Matches WHO-listed bacteria against the CLOVER database using exact, manual, and fuzzy matching. Extracts associated host species.
    -   `2_2_CLOVER_Host_Clean.R`: Takes the unique host species from CLOVER and standardizes their taxonomy using external databases (GBIF, ITIS, IUCN). Generates taxonomic summary visualizations.
    -   `2_3_CLOVER_Network.R`: Prepares the bacteria-host association data for network analysis by merging it with the cleaned host taxonomy. Outputs `clover_who_network.csv`.

3.  **VIRION Integration (Viruses) (`network_building/helpers/virion_data.R`, `network_building/stages/legacy_who_compatibility/3_1_Match_WHO_Virion.R`, `network_building/stages/legacy_who_compatibility/3_2_WHO_Virion_Hosts.R`, `network_building/stages/legacy_who_compatibility/3_3_Host_Species_Clean.R`):**
    -   `helpers/virion_data.R`: Utility script to load the comprehensive VIRION dataset.
    -   `3_1_Match_WHO_Virion.R`: Matches WHO-listed viruses against the VIRION taxonomy to find corresponding `VirusTaxID`s.
    -   `3_2_WHO_Virion_Hosts.R`: Extracts all known host associations for the matched VIRION viruses, filtering for high-quality detection methods.
    -   `3_3_Host_Species_Clean.R`: Takes the unique host species from VIRION and standardizes their taxonomy, similar to the CLOVER workflow. Generates taxonomic summary visualizations.

4.  **Network Combination and Analysis (`network_building/stages/legacy_who_compatibility/3_4_VIRION_Networks.R`, `network_building/stages/legacy_who_compatibility/4_CombineNetworks.R`):**
    -   `network_building/stages/legacy_who_compatibility/3_4_VIRION_Networks.R`: Primary VIRION network assembly and analysis script. It constructs pathogen-host networks, performs advanced analyses (centrality, modularity, bridge species, sampling bias), and writes `virion_who_network.csv` under `pathogen_association_data/staged/who_networks/source_components/`.
    -   `network_building/stages/legacy_who_compatibility/4_CombineNetworks.R`: Merges the processed network data from CLOVER (bacteria) and VIRION (viruses) into a single, comprehensive dataset for combined analysis.
    -   Deprecated exploratory VIRION visualizations live under `network_building/stages/legacy_who_compatibility/deprecated/`. They are not called by the wrapper entrypoints and should not be treated as part of the current rebuild contract unless repaired against current schemas.

### Working Input Layers

- Raw source artifacts remain in place for provenance and matching:
  - `pathogen_association_data/evidence/who_diseases/backbone/who_pathogens_diseases.csv`
  - `pathogen_association_data/source_data/who_networks/domesticated/domesticated_lab_farmed.csv`
  - `pathogen_association_data/evidence/who_networks/host_pathogen/combined_who_network.csv`
- Derived review artifact with canonical pathogen labels:
  - `pathogen_association_data/evidence/who_networks/host_pathogen/combined_who_network_canonical.csv`
- Current downstream host-pathogen working layer:
  - `pathogen_association_data/evidence/who_networks/host_pathogen/master_plus_who_host_network.csv`
- Legacy WHO-only compatibility/provenance layer:
  - `pathogen_association_data/evidence/who_diseases/backbone/who_pathogens_diseases_zoonotic.csv`
  - `pathogen_association_data/evidence/who_networks/host_pathogen/combined_who_network_canonical_zoonotic.csv`
- Shared path helpers for these layers live in `scripts/associations/working_inputs.R`.
- A separate curation layer for splitting broad pathogen taxa into narrower
  host/vector/amplifier analysis units can be generated with:
  - `scripts/associations/network_building/1_2_WHO_Pathogen_Analysis_Units.R`
  - output: `pathogen_association_data/evidence/who_diseases/backbone/who_pathogen_analysis_units.csv`
- A candidate strain inventory for ICTV-backed Sarbecovirus, Merbecovirus,
  and Vesiculovirus rows can be generated with:
  - `scripts/associations/network_building/1_3_WHO_Broad_Taxa_Candidate_Strains.R`
  - output: `pathogen_association_data/staged/who_diseases/broad_taxa/who_broad_taxa_candidate_strains.csv`

The intended workflow is:

- keep `network_building/` scripts pointed at the raw WHO files
- use `master_plus_who_host_network.csv` for active downstream host-pathogen
  consumers
- use the legacy canonical zoonotic layer only through explicit compatibility
  helpers or for frozen provenance/QA scripts

5.  **Vector Screening (`vector_screening/5_1_*` to `vector_screening/5_6_*`):**
    -   `vector_screening/5_1_Pathogen_Vector_Links_Scaffold.R` to `vector_screening/5_6_Backfill_Pathogen_Vector_Links.R`: Build, standardize, and backfill the WHO disease-pathogen-vector tables.
    -   `vector_screening/5_5b_Vector_Name_Cleanup.R` and `vector_screening/5_5c_Vector_Taxonomy_Package_Review.R`: Add conservative vector-name normalization and taxonomy-review outputs.
    -   Use `scripts/associations/working_inputs.R` helpers for all Vector Screening paths. Raw EFSA inputs live under `pathogen_association_data/source_data/vector_screening/`, manual control files under `pathogen_association_data/manual/vector_screening/`, intermediate outputs under `pathogen_association_data/staged/vector_screening/`, and active vector evidence plus QA under `pathogen_association_data/evidence/vector_screening/`.
    -   Optional VecTraits API probes live under `vector_screening/exploratory/vectraits/`; their outputs are exploratory/local and remain ignored under `pathogen_association_data/staged/vector_screening/vectraits/`.

6.  **Host-Vector Integration (`host_vector_integration/5_8_*` to `host_vector_integration/5_11_*`):**
    -   These scripts connect the master-plus WHO host network, filtered by its legacy compatibility flag where needed, to the staged VectorMap and MapVEu host-vector evidence, then write conservative disease-level, pathogen-level, expanded, and QA outputs under `pathogen_association_data/evidence/who_networks/`. Current host-vector outputs are WHO-only and live under `pathogen_association_data/evidence/who_networks/host_vector/who_only/`.

7.  **Geographic Enrichment (`genbank_simple/`):**
    -   `genbank_simple/01b_build_readiness_manifest.R`: Builds the expanded readiness manifest from the disease modelling readiness surface.
    -   `genbank_simple/02_run_genbank_full_retrieval.R`: Retrieves GenBank nuccore records for approved manifest targets. In readiness mode it writes ignored per-target checkpoints under `pathogen_association_data/staged/genbank_simple/local_runs/pathogen_runs_readiness/`.
    -   `genbank_simple/03_summarize_country_metadata.R` to `genbank_simple/06_map_disease_countries.R`: Summarize, QA, standardize, and map country evidence. With `GENBANK_SIMPLE_SUMMARY_KIND=readiness_combined`, these scripts bind the original 19-target run with the expanded readiness run. Generated manifests, intermediates, and map-control files live under `pathogen_association_data/staged/genbank_simple/`; manual query overrides live under `pathogen_association_data/manual/genbank_simple/`; active evidence and QA live under `pathogen_association_data/evidence/genbank_simple/`.
    -   Current modelling-readiness handoffs should use `pathogen_association_data/evidence/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv` when present.

8.  **Role Annotation (`role_annotation/`):**
    -   `role_annotation/roster/01_build_host_role_candidates.R`: Seeds conservative host-role candidate rows from the master-plus compatibility view for the current role-review scope. It writes generated candidate and summary tables under `pathogen_association_data/evidence/role_annotation/`.
    -   `role_annotation/roster/02_build_species_host_vector_roster.R`: Builds a collaborator-facing disease-species roster that covers both vectored and non-vectored diseases by combining host rows from the master-plus compatibility view with vector rows from the curated disease-vector table, plus host-vector observation and competence flags where available.
    -   `role_annotation/curation_inputs/01_prepare_deep_research_batches.R` to `role_annotation/curation_inputs/03_consolidate_deep_research_reports.R`: On-demand curation tooling for preparing, normalizing, and consolidating external research review batches. These are not part of the routine regeneration path unless new review batches are being prepared or imported.
    -   `role_annotation/source_check/01_build_source_check_decision_ledger.R` and `role_annotation/source_check/02_import_source_checked_role_rows.R`: Build the source-check ledger from Deep Research, manual, and generated role-gap candidates, then import accepted rows into official role evidence/assignment tables.
    -   `role_annotation/features/01_build_role_modelling_features.R`: Builds host role proxy/bucket fields plus host/vector biological evidence-tier fields consumed by readiness handoffs.
    -   `readiness/01_build_disease_modelling_readiness.R` and `readiness/02_build_modelling_evidence_tiers_handoff.R`: Assemble downstream modelling-readiness and evidence-tier handoff surfaces.
    -   The old top-level `role_annotation/6_*.R` scripts are compatibility wrappers around the scoped paths above. Prefer the scoped paths in new automation and documentation.
    -   Role annotation files are an interpretation layer. Do not treat candidate rows as final reservoir, amplifier, incidental, dead-end, or vector-role assignments without source-backed evidence review.

### Outputs

- Processed CSV files of pathogen-host associations, standardized taxonomy, and combined network data.
- Summary statistics and quality control reports printed to the console during script execution.
- Deprecated exploratory static and interactive network visualizations may exist
  under `figures/network_plots/`, but they are not regenerated by the current
  wrapper entrypoints.
