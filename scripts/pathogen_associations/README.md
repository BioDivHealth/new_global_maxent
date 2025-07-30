## Pathogen Associations: Project Goals and Results

This directory contains a pipeline for integrating, standardizing, and analyzing global pathogen-host association data, with a focus on WHO priority pathogens. The workflow leverages large-scale databases (WHO, CLOVER, VIRION) to map pathogens to their hosts, standardize taxonomy, and visualize ecological networks.

### Project Goals

- **Integrate Data:** Combine WHO priority pathogen lists with host-pathogen association data from the CLOVER and VIRION databases.
- **Standardize Taxonomy:** Harmonize pathogen and host names using fuzzy matching, manual curation, and taxonomic databases (GBIF, ITIS, IUCN).
- **Map Associations:** Identify which WHO pathogens are present in global databases and extract their host associations.
- **Network Analysis:** Visualize and analyze the structure of pathogen-host networks, including risk stratification, centrality, modularity, and bridge species.
- **Support Research:** Provide processed datasets and visualizations to guide research on zoonotic risk, host diversity, and sampling bias.

### Pipeline Overview

1.  **WHO Pathogen Data Processing (`1_WHO_Diseases.R`):**
    -   Loads and standardizes WHO priority pathogen lists from various regional documents.
    -   Maps pathogens to standardized names using a translation table and fuzzy matching.
    -   Assigns risk categories and family information.
    -   Outputs a cleaned and consolidated pathogen list (`who_pathogens_diseases.csv`).

2.  **CLOVER Integration (Bacteria) (`2_1_CLOVER.R`, `2_2_CLOVER_Host_Clean.R`, `2_3_CLOVER_Network.R`):**
    -   `2_1_CLOVER.R`: Matches WHO-listed bacteria against the CLOVER database using exact, manual, and fuzzy matching. Extracts associated host species.
    -   `2_2_CLOVER_Host_Clean.R`: Takes the unique host species from CLOVER and standardizes their taxonomy using external databases (GBIF, ITIS, IUCN). Generates taxonomic summary visualizations.
    -   `2_3_CLOVER_Network.R`: Prepares the bacteria-host association data for network analysis by merging it with the cleaned host taxonomy. Outputs `clover_who_network.csv`.

3.  **VIRION Integration (Viruses) (`virion_data.R`, `3_1_Match_WHO_Virion.R`, `3_2_WHO_Virion_Hosts.R`, `3_3_Host_Species_Clean.R`):**
    -   `virion_data.R`: Utility script to load the comprehensive VIRION dataset.
    -   `3_1_Match_WHO_Virion.R`: Matches WHO-listed viruses against the VIRION taxonomy to find corresponding `VirusTaxID`s.
    -   `3_2_WHO_Virion_Hosts.R`: Extracts all known host associations for the matched VIRION viruses, filtering for high-quality detection methods.
    -   `3_3_Host_Species_Clean.R`: Takes the unique host species from VIRION and standardizes their taxonomy, similar to the CLOVER workflow. Generates taxonomic summary visualizations.

4.  **Network Combination, Analysis, and Visualization (`3_4_VIRION_Networks.R`, `4_CombineNetworks.R`):**
    -   `3_4_VIRION_Networks.R`: The primary network analysis script for VIRION data. It constructs pathogen-host networks, performs advanced analyses (centrality, modularity, bridge species, sampling bias), and generates a wide range of static (`.png`) and interactive (`.html`) visualizations. Outputs are saved in `figures/network_plots/`. This script also prepares the `virion_who_network.csv`.
    -   `4_CombineNetworks.R`: Merges the processed network data from CLOVER (bacteria) and VIRION (viruses) into a single, comprehensive dataset for combined analysis.

### Key Results

- **Comprehensive Mapping:** Most WHO priority pathogens were successfully mapped to global databases, with host associations extracted for both bacteria and viruses.
- **Standardized Datasets:** Output includes harmonized tables of pathogen-host associations and standardized host taxonomy for both CLOVER and VIRION data.
- **Network Insights:** Advanced network analysis reveals key hosts (bridge species), community structure (modularity), potential sampling biases, and high-risk pathogen ecosystems.
- **Visualizations:** The pipeline generates publication-quality figures and interactive network visualizations for further exploration.

### Outputs

- Processed CSV files of pathogen-host associations, standardized taxonomy, and combined network data.
- Summary statistics and quality control reports printed to the console during script execution.
- Static and interactive network visualizations (see `figures/network_plots/` and `figures/network_plots/advanced_analysis/`).
