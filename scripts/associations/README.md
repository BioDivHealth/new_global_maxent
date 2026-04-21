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
  backfilling.
- `host_vector_integration/`
  Joins the WHO disease/pathogen network to observational host-vector evidence
  and writes disease-level, pathogen-level, and expanded host-vector outputs.
- `genbank/`
  Adds pathogen-country enrichment from GenBank/NCBI and contains the shared
  helper layer plus source-routing rules.

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

6.  **Host-Vector Integration (`host_vector_integration/5_8_*` to `host_vector_integration/5_11_*`):**
    -   These scripts connect the WHO disease/pathogen network to the staged VectorMap and MapVEu host-vector evidence, then write conservative disease-level, pathogen-level, expanded, and QA outputs under `pathogen_association_data/WHO/networks/`.

7.  **Geographic Enrichment (`genbank/5_7_*`):**
    -   `genbank/5_7_GenBank_Pathogen_Country_Metadata.R`: Builds GenBank-ready pathogen queries from the WHO master tables, writes a per-pathogen metadata-source recommendation (`nuccore` vs `biosample`), fetches accession-level nuccore metadata, and writes pathogen-country summaries under `pathogen_association_data/WHO/genbank/`.
    -   `genbank/5_7a_GenBank_Search_Diagnostics.R`: Runs lightweight NCBI ESearch diagnostics for the same WHO pathogen manifest and writes per-query success/failure tables before any accession fetching.
    -   `genbank/genbank_metadata_helpers.R`: Shared helper layer for query building, NCBI search/fetch plumbing, taxonomy-link fallback, and metadata-source routing.
    -   `genbank/GENBANK_METADATA_SOURCE_RULES.md`: Exact routing rules for deciding when a pathogen should start in `nuccore` versus `biosample`, plus when a weak `nuccore` result should be escalated to `biosample` review.

### Outputs

- Processed CSV files of pathogen-host associations, standardized taxonomy, and combined network data.
- Summary statistics and quality control reports printed to the console during script execution.
- Static and interactive network visualizations (see `figures/network_plots/` and `figures/network_plots/advanced_analysis/`).
