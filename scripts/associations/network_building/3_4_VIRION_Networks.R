# ------------------------------------------------------------------------------|
# 5_Network_Visualization.R
# ------------------------------------------------------------------------------|
# Purpose: Create network visualizations of WHO pathogen-host associations
#          using standardized data from previous processing steps
#
# Input:   who_pathogens_virion_hosts_summary.csv (from 3_WHO_Virion_Hosts.R)
#          who_host_species_standardized.csv (from 4_Host_Species_Clean.R)
#
# Output:  Interactive and static network plots
# ------------------------------------------------------------------------------|

# ------------------------------| Load libraries |------------------------------
library(pacman)
p_load(here, tidyverse, igraph, ggraph, networkD3, visNetwork, 
       plotly, RColorBrewer, viridis, cowplot, scales)

# Additional network packages
if (!require(tidygraph)) install.packages("tidygraph")
library(tidygraph)
library(magrittr)

source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------| Load data |--------------------------------
cat("Loading pathogen-host association data...\n")
# Load the main association data
host_associations <- read_csv(file.path(who_virion_dir, "who_pathogens_virion_hosts_summary.csv"))
host_detection_methods_keep <- c("Isolation/Observation", "PCR/Sequencing")

# Harmonize Virus names to standardized taxonomy
synonyms <- c(
  "influenza a virus"                       = "alphainfluenzavirus influenzae",
  "chikungunya virus"                       = "alphavirus chikungunya",
  "dengue virus"                            = "orthoflavivirus denguei",
  "zika virus"                              = "orthoflavivirus zikaense",
  "west nile virus"                         = "orthoflavivirus nilense",
  "yellow fever virus"                      = "orthoflavivirus flavi",
  "monkeypox virus"                         = "orthopoxvirus monkeypox",
  "henipavirus nipahense"                   = "henipavirus nipahense",
  "middle east respiratory syndrome-related coronavirus" =
    "betacoronavirus cameli"
  # Leave SARS rows as they are unless you deliberately merge them
)
 host_associations$Virus_std <- unname(synonyms[ tolower(host_associations$Virus) ])
 host_associations$Virus_std[ is.na(host_associations$Virus_std) ] <- host_associations$Virus[ is.na(host_associations$Virus_std) ]

# Rename columns
host_associations$Virus_og = host_associations$Virus  # Keep original names for reference
host_associations$Virus = host_associations$Virus_std  # Use standardized names for analysis

# Load standardized host taxonomy  
host_taxonomy <- read_csv(file.path(who_virion_dir, "who_host_species_standardized.csv"))
host_taxonomy$Host_lower = str_to_lower(host_taxonomy$Host)

# Clean and prepare data for network analysis
network_data <- host_associations %>%
  # Create lowercase host names for matching
  mutate(Host_lower = str_to_lower(Host)) %>%
  # Use standardized host names if available (match on lowercase)
  left_join(host_taxonomy %>% select(Host, Host_lower, correct_name,Phylum, Class, Family, Order), 
            by = "Host_lower") %>%
  mutate(
    # Use standardized name if available, otherwise original
    Host_clean = coalesce(correct_name, Host.x),  # Host.x is from host_associations
    Pathogen_clean = Pathogens,
    high_quality_detection = if ("high_quality_detection" %in% names(.)) {
      coalesce(high_quality_detection, FALSE)
    } else {
      DetectionMethod %in% host_detection_methods_keep
    },
    host_taxonomy_ready = !is.na(Host_clean) & !is.na(HostTaxID),
    host_flag_review = coalesce(HostFlagID, FALSE),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) {
      coalesce(downstream_default_include, FALSE) & host_taxonomy_ready
    } else {
      high_quality_detection & host_taxonomy_ready & !host_flag_review
    },
    downstream_review_reason = case_when(
      downstream_default_include ~ NA_character_,
      !high_quality_detection ~ paste0("detection_method=", DetectionMethod),
      host_flag_review ~ "virion_host_flag_review",
      !host_taxonomy_ready ~ "host_taxonomy_not_ready",
      TRUE ~ "manual_review"
    ),
    # Create risk categories
    Risk_category = case_when(
      str_detect(`PHEIC risk`, "High") ~ "High Risk",
      str_detect(`PHEIC risk`, "Medium") ~ "Medium Risk", 
      str_detect(`PHEIC risk`, "Low") ~ "Low Risk",
      TRUE ~ "Unknown Risk"
    )
  ) %>%
  # Filter for high-quality detections
  #filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  # Remove uncertain host identifications if desired
  select(Pathogen = Pathogen_clean, Host_clean, Disease_name, HostTaxID, PathogenTaxID = VirusTaxID, 
         PathogenGenus = VirusGenus, PathogenFamily = VirusFamily, 
         PathogenOrder = VirusOrder, PathogenClass = VirusClass, HostPhylum = Phylum,
         HostClass = Class, HostFamily = Family, HostOrder = Order, DetectionMethod,
         high_quality_detection, downstream_default_include, downstream_review_reason,
         in_gibb_etal, in_empres_i,
         `PHEIC risk`) %>%
  distinct() %>%
  mutate(MainSource = "VIRION") %>%
  filter(!is.na(Pathogen), !is.na(Host_clean))

cat("Prepared", nrow(network_data), "pathogen-host associations for visualization\n")
output_path <- who_network_source_component_path("virion_who_network.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(network_data, output_path)
