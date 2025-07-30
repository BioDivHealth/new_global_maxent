# ------------------------------------------------------------------------------|
# 3_WHO_Virion_Hosts.R
# ------------------------------------------------------------------------------|
# Purpose: Extract host associations for WHO priority pathogens from VIRION
#          database using the VirusTaxID matches from 2_Match_WHO_Virion.R
#
# Input:   who_pathogens_virion_taxid.csv (from 2_Match_WHO_Virion.R)
#          virion_data object (from virion_data.R)
#
# Output:  who_pathogens_hosts_long.csv - detailed host-pathogen interactions
#          who_pathogens_hosts_summary.csv - summarised host counts per pathogen
# ------------------------------------------------------------------------------|

# ------------------------------| Load libraries |------------------------------
library(tidyverse)
library(here)
library(pacman)
library(magrittr)
library(frictionless)
# install.packages("devtools")
#devtools::install_github("frictionlessdata/frictionless-r")

# ------------------------------| Helper paths  |------------------------------
input_csv_path <- here("pathogen_association_data", "WHO", "virion", "who_pathogens_virion_taxid.csv")
output_long_path <- here("pathogen_association_data", "WHO", "virion", "who_pathogens_virion_hosts_long.csv")
output_summary_path <- here("pathogen_association_data", "WHO", "virion", "who_pathogens_virion_hosts_summary.csv")

# ----------------------------- Load datasets ------------------------------
## 1. WHO pathogen-virion matches -----------------------------------------------
cat("Loading WHO-VIRION matched data...\n")
who_virion <- read_csv(input_csv_path, show_col_types = FALSE)

## 2. VIRION database -----------------------------------------------------------
cat("Loading VIRION database...\n")
if (!exists("virion_data")) {
  source(here("scripts", "pathogen_associations", "virion_data.R"))
}

dictionaries = virionData::get_data_dictionary(datapackage_json = here("pathogen_association_data","virion_download",
                                                                       "15896981","datapackage.json"))

print(dictionaries)
dictionaries$virion_csv

virion = virion_data$virion

# ----------------------------- Data processing ---------------------------
cat("Processing pathogen-host associations...\n")

## 1: Normalise VirusTaxID column (handle multiple IDs per row) -------
# Handle cases where columns have different numbers of semicolon-separated values
who_virion_long <- who_virion %>%
  # Only process rows that have VirusTaxID matches
  filter(!is.na(VirusTaxID)) %>%
  # Split semicolon-separated values into lists
  mutate(
    VirusTaxID_list = str_split(VirusTaxID, ";\\s*"),
    Virion_VirusName_list = str_split(Virion_VirusName, ";\\s*"),
    Virion_VirusFamily_list = str_split(Virion_VirusFamily, ";\\s*"),
    Virion_Database_list = str_split(Virion_Database, ";\\s*")
  ) %>%
  # Remove original columns
  select(-VirusTaxID, -Virion_VirusName, -Virion_VirusFamily, -Virion_Database) %>%
  # Unnest VirusTaxID first
  unnest(VirusTaxID_list) %>%
  # Create position index for matching other columns
  group_by(ID) %>%
  mutate(virus_position = row_number()) %>%
  ungroup() %>%
  # Extract corresponding values from other list columns by position
  mutate(
    VirusTaxID = as.integer(str_trim(VirusTaxID_list)),
    Virion_VirusName = map2_chr(Virion_VirusName_list, virus_position, 
                                ~if(length(.x) >= .y) str_trim(.x[.y]) else NA_character_),
    Virion_VirusFamily = map2_chr(Virion_VirusFamily_list, virus_position,
                                  ~if(length(.x) >= .y) str_trim(.x[.y]) else NA_character_),
    Virion_Database = map2_chr(Virion_Database_list, virus_position,
                               ~if(length(.x) >= .y) str_trim(.x[.y]) else NA_character_)) %>%
  # Clean up temporary columns
  select(-VirusTaxID_list, -Virion_VirusName_list, -Virion_VirusFamily_list, -Virion_Database_list, -virus_position) %>%
  # Remove any rows where VirusTaxID conversion failed
  filter(!is.na(VirusTaxID))

cat("  - Expanded to", nrow(who_virion_long), "WHO pathogen-VirusTaxID pairs\n")

# Filter selected viruses from Virion  ----------------------------
virion_who_taxa = virion %>% 
  filter(VirusTaxID %in% unique(who_virion_long$VirusTaxID)) %>%
  filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  select(Host, Virus, 
         HostTaxID,HostGenus, HostFamily, HostOrder, HostClass, 
         VirusTaxID,VirusGenus,VirusFamily,VirusOrder,VirusClass,
         DetectionMethod) %>%
  distinct() %>% 
  filter(!is.na(Host))

# This includes all types fo risk + all a lot of metadata
virion_who_taxa_detailed = virion %>% 
  filter(VirusTaxID %in% unique(who_virion_long$VirusTaxID)) %>%
  filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  select(VirusTaxID, HostTaxID, Host, Virus, 
         HostGenus, HostFamily, HostOrder, HostClass, 
         VirusGenus, VirusFamily, VirusOrder, VirusClass,
         DetectionMethod, DetectionOriginal, HostFlagID,
         # Provenance columns
         Database, DatabaseVersion, PublicationYear, ReferenceText, PMID,
         ReleaseYear, CollectionYear, AssocID, NCBIAccession) %>%
  distinct() %>% 
  filter(!is.na(Host)) %>%
  # Group by unique virus-host pairs and collapse metadata columns
  group_by(VirusTaxID, HostTaxID, Host, Virus, 
           HostGenus, HostFamily, HostOrder, HostClass, 
           VirusGenus, VirusFamily, VirusOrder, VirusClass, 
           DetectionMethod, DetectionOriginal, HostFlagID) %>%
  summarise(
    Databases = paste(sort(unique(na.omit(Database))), collapse = "; "),
    DatabaseVersions = paste(sort(unique(na.omit(DatabaseVersion))), collapse = "; "),
    PublicationYears = paste(sort(unique(na.omit(PublicationYear))), collapse = "; "),
    ReferenceTexts = paste(unique(na.omit(ReferenceText)), collapse = " | "),
    PMIDs = paste(sort(unique(na.omit(PMID))), collapse = "; "),
    ReleaseYears = paste(sort(unique(na.omit(ReleaseYear))), collapse = "; "),
    CollectionYears = paste(sort(unique(na.omit(CollectionYear))), collapse = "; "),
    AssocIDs = paste(sort(unique(na.omit(AssocID))), collapse = "; "),
    NCBIAccessions = paste(unique(na.omit(NCBIAccession)), collapse = "; "),
    n_records = n(),  # Number of original records collapsed
    .groups = "drop"
  )

# Join Virion with full WHO -------------------------------
who_virion_hosts_complete <- who_virion_long %>%
  left_join(virion_who_taxa_detailed, by = "VirusTaxID", relationship = "many-to-many", suffix = c("_WHO", "_VIRION")) %>%
  # Filter to only include rows where host data was found
  filter(!is.na(Host)) %>%
  # Reorder columns for better readability
  select(
    # WHO pathogen information
    ID, Family, `PHEIC risk`, Pathogens, Disease_name, previous_name, msl39_viral_name,
    # VIRION virus information
    VirusTaxID, Virion_VirusName, Virion_VirusFamily, Virion_Database,
    # Host information  
    HostTaxID, Host, HostGenus, HostFamily, HostOrder, HostClass,
    # Virus taxonomy from VIRION (more detailed)
    Virus, VirusGenus, VirusFamily, VirusOrder, VirusClass,
    # Detection evidence and quality
    DetectionMethod, DetectionOriginal, HostFlagID,
    # Provenance and metadata
    Databases, DatabaseVersions, PublicationYears, ReferenceTexts, PMIDs,
    ReleaseYears, CollectionYears, AssocIDs, NCBIAccessions, n_records,
    # Match metadata
    matched_name_type, matched_virus_name, match_source, fuzzy_scores, num_virion_matches
  )

cat("Combined WHO-VIRION dataset contains", nrow(who_virion_hosts_complete), "pathogen-host associations\n")
cat("Covering", n_distinct(who_virion_hosts_complete$Pathogens), "WHO pathogens\n")
cat("Covering", n_distinct(who_virion_hosts_complete$Virus), "Virion pathogens\n")
cat("With", n_distinct(who_virion_hosts_complete$Host), "unique host species\n")


write_csv(who_virion_hosts_complete, output_long_path)

# ----------------------------- Summarise host associations -------------------
who_virion_hosts_short = who_virion_hosts_complete %>% select(Host, 
                                                              Virus, 
                                                              `PHEIC risk`, 
                                                              #Pathogens, 
                                                              #Disease_name,
                                                              HostTaxID,
                                                              HostGenus, 
                                                              HostFamily, 
                                                              HostOrder, 
                                                              HostClass, 
                                                              VirusTaxID,
                                                              VirusGenus,
                                                              VirusFamily,
                                                              VirusOrder,
                                                              VirusClass,
                                                              DetectionMethod,
                                                              HostFlagID) %>% filter(!is.na(Host)) %>% distinct()

write_csv(who_virion_hosts_short, output_summary_path)

host_species = who_virion_hosts_short %>%
  select(Host, HostGenus, HostFamily, HostFlagID) %>%
  mutate(
    Host = str_to_sentence(Host),
    HostGenus = str_to_sentence(HostGenus),
    HostFamily = str_to_sentence(HostFamily)
  ) %>%
  distinct()



dim(host_species)

# ----------------------------- Quality control -----------------------------
cat("\n=== QUALITY CONTROL SUMMARY ===\n")

# Check for WHO pathogens without host data
unmatched_pathogens <- who_virion %>%
  filter(is.na(VirusTaxID)) %>%
  distinct(Pathogens)

cat("WHO pathogens without VIRION matches:", nrow(unmatched_pathogens), "\n")
if (nrow(unmatched_pathogens) > 0) {
  cat("Unmatched pathogens:\n")
  print(unmatched_pathogens$Pathogens)
}

# Check for uncertain host identifications
uncertain_hosts <- virus_host_final %>%
  filter(HostFlagID == TRUE) %>%
  distinct(Pathogens, Host) %>%
  group_by(Pathogens) %>%
  summarise(uncertain_hosts = paste(Host, collapse = "; "), .groups = "drop")

cat("Pathogens with uncertain host IDs:", nrow(uncertain_hosts), "\n")

# Summary statistics
cat("\n=== FINAL SUMMARY ===\n")
cat("Total WHO pathogens processed:", n_distinct(virus_host_final$Pathogens), "\n")
cat("Total unique hosts found:", n_distinct(virus_host_final$Host), "\n")
cat("Total pathogen-host interactions:", nrow(virus_host_final), "\n")
cat("Host classes represented:", n_distinct(virus_host_final$HostClass, na.rm = TRUE), "\n")
cat("Detection methods available:", n_distinct(virus_host_final$DetectionMethod, na.rm = TRUE), "\n")

# Show top 10 pathogens by host diversity
cat("\nTop 10 pathogens by host diversity:\n")
print(hosts_per_pathogen %>% 
      select(Pathogens, n_unique_hosts, n_host_families) %>% 
      head(10))

cat("\nHost class distribution:\n")
print(host_class_summary)

cat("\nProcessing complete!\n") 