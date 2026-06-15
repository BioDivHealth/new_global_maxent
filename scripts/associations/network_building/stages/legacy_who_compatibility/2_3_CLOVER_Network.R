# -----------------------------------------------------------------------------|
# 2_3_CLOVER_Network.R ----
# -----------------------------------------------------------------------------|
# Purpose: Join WHO-CLOVER bacteria host associations to standardized host
#          taxonomy and write the CLOVER source-component network table.
# Inputs : clover_host_species_standardized.csv, who_bacteria_clover_taxid.csv,
#          and who_bacteria_clover_hosts.csv
# Outputs: clover_who_network.csv
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# 1. Load required libraries and path helpers ----
# -----------------------------------------------------------------------------|
library(pacman)
p_load(here, tidyverse)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "legacy_who_compatibility_helpers.R"
))

# -----------------------------------------------------------------------------|
# 2. Load and validate inputs ----
# -----------------------------------------------------------------------------|
host_taxonomy = read_csv(file.path(who_clover_dir, "clover_host_species_standardized.csv"))
disease_names = read_csv(file.path(who_clover_dir, "who_bacteria_clover_taxid.csv"))
host_associations = read_csv(file.path(who_clover_dir, "who_bacteria_clover_hosts.csv"))
host_detection_methods_keep <- c("Isolation/Observation", "PCR/Sequencing")

legacy_who_require_columns(
  host_taxonomy,
  c("Host", "HostTaxID", "correct_name", "Spp_syn", "Phylum", "Class", "Family", "Order"),
  "CLOVER host taxonomy"
)

legacy_who_require_columns(
  disease_names,
  c("ID", "Disease_name"),
  "CLOVER disease-name lookup"
)

legacy_who_require_columns(
  host_associations,
  c(
    "ID", "Pathogen", "PathogenTaxID", "PHEIC risk", "Host", "HostTaxID",
    "PathogenClass", "PathogenOrder", "PathogenFamily", "PathogenGenus",
    "DetectionMethod"
  ),
  "CLOVER host associations"
)

host_taxonomy$Host_lower = str_to_lower(host_taxonomy$Host)

# -----------------------------------------------------------------------------|
# 3. Build CLOVER source-component network table ----
# -----------------------------------------------------------------------------|
# Match disease names to host_associations
host_associations = host_associations %>%
  left_join(disease_names %>% select(ID, Disease_name), by = "ID")

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

# -----------------------------------------------------------------------------|
# 4. Write source-component output ----
# -----------------------------------------------------------------------------|
output_path <- who_network_source_component_path("clover_who_network.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(network_data, output_path)
