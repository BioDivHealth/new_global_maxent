library(pacman)
p_load(here, tidyverse, igraph, ggraph, networkD3, visNetwork, 
       plotly, RColorBrewer, viridis, cowplot, scales, magrittr,dplyr)

source(here("scripts", "associations", "working_inputs.R"))

host_taxonomy = read_csv(file.path(who_clover_dir, "clover_host_species_standardized.csv"))
host_taxonomy$Host_lower = str_to_lower(host_taxonomy$Host)
# names(host_taxonomy)
#  [1] "Host"                  "HostTaxID"             "correct_name"          "type"                  "taxon_level"           "Genus"                
#  [7] "Family"                "Order"                 "Class"                 "Phylum"                "host_species"          "Spp_syn"              
# [13] "IUCN_spp"              "Or_name"               "IUCN_Present"          "IUCN_id"               "IUCN_name"             "IUCN_latest"          
# [19] "IUCN_date"             "IUCN_Category"         "IUCN_N_syn"            "IUCN_syn"              "IUCN_status"           "IUCN_Phylum"          
# [25] "IUCN_Class"            "IUCN_Order"            "IUCN_Family"           "ITIS_Present"          "ITIS_is_valid"         "ITIS_id"              
# [31] "ITIS_name"             "ITIS_Phylum"           "ITIS_N_syn"            "ITIS_syn"              "ITIS_Class"            "ITIS_Order"           
# [37] "ITIS_Family"           "ITIS_species_in_genus" "GBIF_Present"          "GBIF_id"               "GBIF_name"             "GBIF_N_syn"           
# [43] "GBIF_syn"              "GBIF_Phylum"           "GBIF_Status"           "GBIF_Class"            "GBIF_Order"            "GBIF_Family"  

unique(host_taxonomy$Host)
unique(host_taxonomy$correct_name)
unique(host_taxonomy$host_species)

disease_names = read_csv(file.path(who_clover_dir, "who_bacteria_clover_taxid.csv"))
host_associations = read_csv(file.path(who_clover_dir, "who_bacteria_clover_hosts.csv"))
host_detection_methods_keep <- c("Isolation/Observation", "PCR/Sequencing")
names(host_associations)
#  [1] "ID"                      "bacteria_name"           "name_type"               "match_source"            "dist"                    "PathogenTaxID"          
#  [7] "Pathogen"                "PathogenType"            "PathogenClass"           "PathogenOrder"           "PathogenFamily"          "PathogenGenus"          
# [13] "PathogenNCBIResolved"    "Host"                    "HostTaxID"               "HostGenus"               "HostFamily"              "HostOrder"              
# [19] "HostClass"               "HostNCBIResolved"        "DetectionMethod"         "DetectionMethodOriginal" "ICTVRatified"            "Database"               
# [25] "DatabaseVersion"         "DatabaseDOI"             "PublicationYear"         "ReferenceText"           "PMID"                    "ReleaseYear"            
# [31] "AssocID"                 "NCBIAccession"         

# Match disease names to host_associations
host_associations = host_associations %>%
  left_join(disease_names %>% select(ID, Disease_name), by = "ID")

# $bacteria_name is the original name
# $Pathogen is the standardized name
# $Host is the CLOVER host name
 unique(host_associations$bacteria_name)
 unique(host_associations$Pathogen)
 table(unique(host_associations$Host) %in% host_taxonomy$Host_lower)
 
# Clean and prepare data for network analysis
network_data <- host_associations %>%
  # Create lowercase host names for matching
  mutate(Host_lower = str_to_lower(Host)) %>%
  # Use standardized host names if available (match on lowercase)
  left_join(host_taxonomy %>% select(Host, HostTaxID, Host_lower, correct_name,Spp_syn,Phylum, Class, Family, Order), 
            by = "Host_lower") %>%
  mutate(
    # Use standardized name if available, otherwise original
    Host_clean = coalesce(correct_name, Host.x),  # Host.x is from host_associations
    high_quality_detection = if ("high_quality_detection" %in% names(.)) {
      coalesce(high_quality_detection, FALSE)
    } else {
      DetectionMethod %in% host_detection_methods_keep
    },
    host_taxonomy_ready = !is.na(Host_clean) & !is.na(HostTaxID.x),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) {
      coalesce(downstream_default_include, FALSE)
    } else {
      high_quality_detection & host_taxonomy_ready
    },
    downstream_review_reason = case_when(
      downstream_default_include ~ NA_character_,
      !high_quality_detection ~ paste0("detection_method=", DetectionMethod),
      !host_taxonomy_ready ~ "host_taxonomy_not_ready",
      TRUE ~ "manual_review"
    )
  ) %>%
  # Filter for high-quality detections
  #filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  # Remove uncertain host identifications if desired
  select(ID, Pathogen, PathogenTaxID, `PHEIC risk`, Disease_name, HostTaxID = HostTaxID.x,-HostTaxID.y,Host_clean, PathogenClass,PathogenOrder,PathogenFamily, PathogenGenus, HostPhylum = Phylum, HostClass = Class, HostFamily = Family, HostOrder = Order, DetectionMethod, high_quality_detection, downstream_default_include, downstream_review_reason) %>%
  distinct() %>%
  filter(!is.na(Host_clean)) %>% 
  mutate(MainSource = "CLOVER")

cat("Prepared", nrow(network_data), "pathogen-host associations for visualization\n")

output_path <- who_network_source_component_path("clover_who_network.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(network_data, output_path)
