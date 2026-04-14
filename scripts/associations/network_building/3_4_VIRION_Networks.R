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

# ------------------------------| Load data |--------------------------------
cat("Loading pathogen-host association data...\n")
# Load the main association data
host_associations <- read_csv(here("pathogen_association_data", "WHO", "virion", "who_pathogens_virion_hosts_summary.csv"))
host_associations_long <- read_csv(here("pathogen_association_data", "WHO", "virion", "who_pathogens_virion_hosts_long.csv"))
disease_names = host_associations_long %>% select(Disease_name, Virion_VirusName) %>% distinct() %>% filter(!is.na(Virion_VirusName) & !is.na(Disease_name))
disease_names$Virus = disease_names$Virion_VirusName
disease_names$Disease_name[disease_names$Virus=="alphainfluenzavirus influenzae"] = "Influenza"
disease_names = unique(disease_names)

host_associations = host_associations %>%
  left_join(disease_names %>% select(Virus, Disease_name), by = "Virus")


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
host_taxonomy <- read_csv(here("pathogen_association_data", "WHO", "virion", "who_host_species_standardized.csv"))
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
    # Simplify virus names for better visualization
    Virus_clean = str_remove(Virus, " sp\\.$|strain.*$"),
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
  filter(!HostFlagID | is.na(HostFlagID)) %>%
  select(Pathogen = Virus_clean, Host_clean, Disease_name, HostTaxID, PathogenTaxID = VirusTaxID, 
         PathogenGenus = VirusGenus, PathogenFamily = VirusFamily, 
         PathogenOrder = VirusOrder, PathogenClass = VirusClass, HostPhylum = Phylum,
         HostClass = Class, HostFamily = Family, HostOrder = Order, DetectionMethod,
         `PHEIC risk`) %>%
  distinct() %>%
  mutate(MainSource = "VIRION") %>%
  filter(!is.na(Pathogen), !is.na(Host_clean))

cat("Prepared", nrow(network_data), "pathogen-host associations for visualization\n")
dir.create(here("pathogen_association_data", "WHO", "networks"), showWarnings = FALSE)
write_csv(network_data, here("pathogen_association_data", "WHO", "networks", "virion_who_network.csv"))