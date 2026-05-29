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

### Project Goals

- **Integrate Data:** Combine WHO priority pathogen lists with host-pathogen association data from the CLOVER and VIRION databases.
- **Standardize Taxonomy:** Harmonize pathogen and host names using fuzzy matching, manual curation, and taxonomic databases (GBIF, ITIS, IUCN).
- **Map Associations:** Identify which WHO pathogens are present in global databases and extract their host associations.
- **Network Analysis:** Visualize and analyze the structure of pathogen-host networks, including risk stratification, centrality, modularity, and bridge species.
- **Support Research:** Provide processed datasets and visualizations to guide research on zoonotic risk, host diversity, and sampling bias.

### Pipeline Overview

1.  **WHO Pathogen Data Processing (`network_building/1_WHO_Diseases.R`):**
    -   Loads and standardizes WHO priority pathogen lists from various regional documents.
    -   Maps pathogens to standardized names using a translation table and fuzzy matching.
    -   Assigns risk categories and family information.
    -   Preserves whether each final pathogen is priority, prototype, or both, plus per-region WHO source status (`priority`, `prototype`, `both`, or `none`) across the WHO regions.
    -   Outputs a cleaned and consolidated pathogen list (`who_pathogens_diseases.csv`).

2.  **CLOVER Integration (Bacteria) (`network_building/2_1_CLOVER.R`, `network_building/2_2_CLOVER_Host_Clean.R`, `network_building/2_3_CLOVER_Network.R`):**
    -   `2_1_CLOVER.R`: Matches WHO-listed bacteria against the CLOVER database using exact, manual, and fuzzy matching. Extracts associated host species.
    -   `2_2_CLOVER_Host_Clean.R`: Takes the unique host species from CLOVER and standardizes their taxonomy using external databases (GBIF, ITIS, IUCN). Generates taxonomic summary visualizations.
    -   `2_3_CLOVER_Network.R`: Prepares the bacteria-host association data for network analysis by merging it with the cleaned host taxonomy. Outputs `clover_who_network.csv`.

3.  **VIRION Integration (Viruses) (`network_building/virion_data.R`, `network_building/3_1_Match_WHO_Virion.R`, `network_building/3_2_WHO_Virion_Hosts.R`, `network_building/3_3_Host_Species_Clean.R`):**
    -   `virion_data.R`: Utility script to load the comprehensive VIRION dataset.
    -   `3_1_Match_WHO_Virion.R`: Matches WHO-listed viruses against the VIRION taxonomy to find corresponding `VirusTaxID`s.
    -   `3_2_WHO_Virion_Hosts.R`: Extracts all known host associations for the matched VIRION viruses, filtering for high-quality detection methods.
    -   `3_3_Host_Species_Clean.R`: Takes the unique host species from VIRION and standardizes their taxonomy, similar to the CLOVER workflow. Generates taxonomic summary visualizations.

4.  **Network Combination, Analysis, and Visualization (`network_building/3_4_VIRION_Networks.R`, `network_building/3_5_VIRION_Visualise_Networks.R`, `network_building/4_CombineNetworks.R`):**
    -   `network_building/3_4_VIRION_Networks.R`: Primary VIRION network assembly and analysis script. It constructs pathogen-host networks, performs advanced analyses (centrality, modularity, bridge species, sampling bias), and writes `virion_who_network.csv` under `pathogen_association_data/WHO/networks/`.
    -   `network_building/3_5_VIRION_Visualise_Networks.R`: Reads `virion_who_network.csv` and generates static (`.png`) and interactive (`.html`) visualizations, with outputs saved under `figures/network_plots/` (and related subfolders).
    -   `network_building/4_CombineNetworks.R`: Merges the processed network data from CLOVER (bacteria) and VIRION (viruses) into a single, comprehensive dataset for combined analysis.

### Working Input Layers

- Raw source artifacts remain in place for provenance and matching:
  - `pathogen_association_data/WHO/who_diseases/who_pathogens_diseases.csv`
  - `pathogen_association_data/WHO/networks/combined_who_network.csv`
- Derived review artifact with canonical pathogen labels:
  - `pathogen_association_data/WHO/networks/combined_who_network_canonical.csv`
- Default downstream working layer for the rest of `scripts/associations/`:
  - `pathogen_association_data/WHO/who_diseases/who_pathogens_diseases_zoonotic.csv`
  - `pathogen_association_data/WHO/networks/combined_who_network_canonical_zoonotic.csv`
- Shared path helpers for these layers live in `scripts/associations/working_inputs.R`.
- A separate curation layer for splitting broad pathogen taxa into narrower
  host/vector/amplifier analysis units can be generated with:
  - `scripts/associations/network_building/1_2_WHO_Pathogen_Analysis_Units.R`
  - output: `pathogen_association_data/WHO/who_diseases/who_pathogen_analysis_units.csv`
- A candidate strain inventory for ICTV-backed Sarbecovirus, Merbecovirus,
  and Vesiculovirus rows can be generated with:
  - `scripts/associations/network_building/1_3_WHO_Broad_Taxa_Candidate_Strains.R`
  - output: `pathogen_association_data/WHO/who_diseases/who_broad_taxa_candidate_strains.csv`

The intended workflow is:

- keep `network_building/` scripts pointed at the raw WHO files
- use the canonical zoonotic working layer for downstream vector-screening,
  host-vector integration, host-vector source filtering, and GenBank scripts

5.  **Vector Screening (`vector_screening/5_1_*` to `vector_screening/5_6_*`):**
    -   `vector_screening/5_1_Pathogen_Vector_Links_Scaffold.R` to `vector_screening/5_6_Backfill_Pathogen_Vector_Links.R`: Build, standardize, and backfill the WHO disease-pathogen-vector tables.
    -   `vector_screening/5_5b_Vector_Name_Cleanup.R` and `vector_screening/5_5c_Vector_Taxonomy_Package_Review.R`: Add conservative vector-name normalization and taxonomy-review outputs.
    -   Use `scripts/associations/working_inputs.R` helpers for all Vector Screening paths. Raw EFSA inputs live under `pathogen_association_data/source_data/vector_screening/`, manual control files under `pathogen_association_data/manual/vector_screening/`, intermediate outputs under `pathogen_association_data/staged/vector_screening/`, and active vector evidence plus QA under `pathogen_association_data/evidence/vector_screening/`.
    -   VecTraits API outputs are exploratory/local and remain ignored under `pathogen_association_data/staged/vector_screening/vectraits/`.

6.  **Host-Vector Integration (`host_vector_integration/5_8_*` to `host_vector_integration/5_11_*`):**
    -   These scripts connect the WHO disease/pathogen network to the staged VectorMap and MapVEu host-vector evidence, then write conservative disease-level, pathogen-level, expanded, and QA outputs under `pathogen_association_data/WHO/networks/`.

7.  **Geographic Enrichment (`genbank_simple/`):**
    -   `genbank_simple/01b_build_readiness_manifest.R`: Builds the expanded readiness manifest from the disease modelling readiness surface.
    -   `genbank_simple/02_run_genbank_full_retrieval.R`: Retrieves GenBank nuccore records for approved manifest targets. In readiness mode it writes ignored per-target checkpoints under `pathogen_association_data/staged/genbank_simple/local_runs/pathogen_runs_readiness/`.
    -   `genbank_simple/03_summarize_country_metadata.R` to `genbank_simple/06_map_disease_countries.R`: Summarize, QA, standardize, and map country evidence. With `GENBANK_SIMPLE_SUMMARY_KIND=readiness_combined`, these scripts bind the original 19-target run with the expanded readiness run. Generated manifests, intermediates, and map-control files live under `pathogen_association_data/staged/genbank_simple/`; manual query overrides live under `pathogen_association_data/manual/genbank_simple/`; active evidence and QA live under `pathogen_association_data/evidence/genbank_simple/`.
    -   Current modelling-readiness handoffs should use `pathogen_association_data/evidence/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv` when present.

8.  **Role Annotation (`role_annotation/6_1_*`):**
    -   `role_annotation/6_1_Derive_Host_Role_Candidates.R`: Seeds conservative host-role candidate rows from the canonical WHO disease-pathogen-host backbone for the current role-review scope. It writes generated candidate and summary tables under `pathogen_association_data/evidence/role_annotation/`.
    -   `role_annotation/6_2_Derive_Species_Host_Vector_Roster.R`: Builds a collaborator-facing disease-species roster that covers both vectored and non-vectored diseases by combining host rows from the canonical WHO backbone with vector rows from the curated disease-vector table, plus host-vector observation and competence flags where available.
    -   Role annotation files are an interpretation layer. Do not treat candidate rows as final reservoir, amplifier, incidental, dead-end, or vector-role assignments without source-backed evidence review.

### Outputs

- Processed CSV files of pathogen-host associations, standardized taxonomy, and combined network data.
- Summary statistics and quality control reports printed to the console during script execution.
- Static and interactive network visualizations (see `figures/network_plots/` and `figures/network_plots/advanced_analysis/`).
