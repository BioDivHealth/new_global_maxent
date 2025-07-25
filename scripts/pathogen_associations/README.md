## Pathogen Associations: Project Goals and Results

This directory contains a pipeline for integrating, standardizing, and analyzing global pathogen-host association data, with a focus on WHO priority pathogens. The workflow leverages large-scale databases (WHO, CLOVER, VIRION) to map pathogens to their hosts, standardize taxonomy, and visualize ecological networks.

### Project Goals

- **Integrate Data:** Combine WHO priority pathogen lists with host-pathogen association data from the CLOVER and VIRION databases.
- **Standardize Taxonomy:** Harmonize pathogen and host names using fuzzy matching, manual curation, and taxonomic databases (GBIF, ITIS, IUCN).
- **Map Associations:** Identify which WHO pathogens are present in global databases and extract their host associations.
- **Network Analysis:** Visualize and analyze the structure of pathogen-host networks, including risk stratification, centrality, modularity, and bridge species.
- **Support Research:** Provide processed datasets and visualizations to guide research on zoonotic risk, host diversity, and sampling bias.

### Pipeline Overview

1. **WHO Pathogen Data Processing (`1_WHO_Diseases.R`):**
   - Loads and standardizes WHO priority pathogen lists.
   - Maps pathogens to standardized names and risk categories.

2. **CLOVER Integration (`2_1_CLOVER.R`, `2_2_CLOVER_Host_Clean.R`):**
   - Matches WHO bacteria to CLOVER database entries.
   - Extracts host associations and standardizes host taxonomy.
   - Summarizes host diversity and taxonomic breakdown.

3. **VIRION Integration (`virion_data.R`, `2_Match_WHO_Virion.R`, `3_WHO_Virion_Hosts.R`):**
   - Loads VIRION host-virus interaction data.
   - Matches WHO viral pathogens to VIRION taxonomy.
   - Extracts and summarizes host associations for WHO viruses.

4. **Host Taxonomy Standardization (`4_Host_Species_Clean.R`):**
   - Standardizes host species names for downstream analysis and visualization.

5. **Network Visualization and Analysis (`5_Network_Visualization.R`, `6_Network_Examples.R`):**
   - Constructs and visualizes pathogen-host networks.
   - Analyzes network properties (centrality, modularity, bridge species, sampling bias).
   - Produces static and interactive network plots.

### Key Results

- **Comprehensive Mapping:** Most WHO priority pathogens were successfully mapped to global databases, with host associations extracted for both bacteria and viruses.
- **Standardized Datasets:** Output includes harmonized tables of pathogen-host associations and standardized host taxonomy.
- **Network Insights:** Network analysis reveals key hosts (bridge species), community structure, and potential sampling biases.
- **Visualizations:** The pipeline generates publication-quality figures and interactive network visualizations for further exploration.

### Outputs

- Processed CSV files of pathogen-host associations and standardized taxonomy.
- Summary statistics and quality control reports.
- Static and interactive network visualizations (see `figures/network_plots/`).
