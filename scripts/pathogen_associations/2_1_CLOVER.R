# ------------------------------------------------------------------------------
# 2_1_CLOVER.R
# ------------------------------------------------------------------------------
# Purpose: Identify which WHO priority bacteria are present in the CLOVER
#          database and retrieve their host associations directly.
#
# Output : CSV files containing WHO bacteria matches and their host associations
# ------------------------------------------------------------------------------|

# ------------------------------| Load libraries |------------------------------
library(tidyverse)
library(here)
library(stringdist)
library(fuzzyjoin)
library(magrittr)
library(dplyr)
# ------------------------------| Helper paths  |------------------------------
who_csv_path   <- file.path("data_artur", "WHO", "who_diseases", "who_pathogens_diseases.csv")
output_csv_path <- file.path("data_artur", "WHO", "clover", "who_bacteria_clover_taxid.csv")
output_hosts_path <- file.path("data_artur", "WHO", "clover", "who_bacteria_clover_hosts.csv")
output_unique_hosts_path <- file.path("data_artur", "WHO", "clover", "who_bacteria_clover_unique_hosts.csv")

# ------------------------------| Load datasets |------------------------------
# 1. WHO pathogen list ---------------------------------------------------------
who_df <- read_csv(who_csv_path, show_col_types = FALSE)
who_df %<>% filter(Family=="Bacteria") # Include only bacteria
who_df$ID = paste0("B",1:nrow(who_df))

unique(who_df$Pathogens)
# [1] "Klebsiella pneumoniae"                     
# [2] "Salmonella enterica non typhoidal serovars"
# [3] "Shigella dysenteriae serotype 1"           
# [4] "Vibrio cholerae serogroup 0139"            
# [5] "Yersinia pestis"

# 2. CLOVER bacteria database --------------------------------------------------
# Read column descriptions
clover_col_desc <- read_csv(here("data_artur","viralemergence-clover-2604d22",
                                 "clover","clover_1.0_allpathogens",
                                 "CLOVER_ColumnDescriptions.csv"))

# Read bacteria dataset
clover_bacteria <- read_csv(here("data_artur","viralemergence-clover-2604d22",
                                 "clover","clover_1.0_allpathogens",
                                 "CLOVER_1.0_Bacteria_AssociationsFlatFile.csv"))

# ------------------------------| Pre-processing |-----------------------------
# Create a lowercase, trimmed helper column for safer joins --------------------
who_long <- who_df %>%
  # For bacteria, we only have the main Pathogens column (no previous_name, msl39_viral_name)
  dplyr::select(ID, Pathogens, Family) %>%
  mutate(
    bacteria_name = Pathogens,
    name_type = "Pathogens",
    bacteria_key = str_to_lower(str_trim(Pathogens))
  )

# Process CLOVER bacteria data - keep all information including hosts
clover_bacteria_proc <- clover_bacteria %>%
  # Filter for bacteria and apply quality filters
  filter(PathogenType == "bacteria/rickettsia") %>%
  filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  filter(!is.na(Host)) %>%
  mutate(bacteria_key = str_to_lower(str_trim(Pathogen))) %>%
  # Keep all relevant columns since each row is already a pathogen-host association
  dplyr::select(bacteria_key, Host, HostTaxID, HostGenus, HostFamily, HostOrder, HostClass, HostNCBIResolved,
         Pathogen, PathogenTaxID, PathogenType, PathogenClass, PathogenOrder, 
         PathogenFamily, PathogenGenus, PathogenNCBIResolved,
         DetectionMethod, DetectionMethodOriginal, ICTVRatified,
         Database, DatabaseVersion, DatabaseDOI, PublicationYear, ReferenceText, PMID,
         ReleaseYear, AssocID, NCBIAccession)

# ------------------------------| Manual matches |-----------------------------
# Define manual mappings for specific cases where fuzzy matching needs correction
manual_mappings <- tribble(
  ~who_name_lower, ~clover_name_lower,
  # Klebsiella pneumoniae variants
  "klebsiella pneumoniae", "klebsiella pneumoniae",
  # Salmonella variants
  "salmonella enterica non typhoidal serovars", "salmonella enterica",
  # Shigella variants  
  "shigella dysenteriae serotype 1", "shigella dysenteriae",
  # Vibrio variants
  "vibrio cholerae serogroup 0139", "vibrio cholerae",
  # Yersinia variants
  "yersinia pestis", "yersinia pestis"
)

# Apply manual mappings - directly get host associations
manual_matches <- who_long %>%
  inner_join(manual_mappings, by = c("bacteria_key" = "who_name_lower")) %>%
  inner_join(clover_bacteria_proc, by = c("clover_name_lower" = "bacteria_key")) %>%
  mutate(dist = 0.5, match_source = "manual") %>%
  dplyr::select(ID, name_type, bacteria_name, all_of(names(clover_bacteria_proc)), dist, match_source)

# ------------------------------| Exact matches |------------------------------
exact_matches <- who_long %>%
  inner_join(clover_bacteria_proc, by = "bacteria_key") %>%
  mutate(dist = 0, match_source = "exact") %>%
  dplyr::select(ID, name_type, bacteria_name, all_of(names(clover_bacteria_proc)), dist, match_source)

# ------------------------------| Fuzzy matches |------------------------------
# Attempt fuzzy matching only for those still unmatched ------------------------
unmatched <- who_long %>%
  filter(!ID %in% c(exact_matches$ID, manual_matches$ID))

if (nrow(unmatched) > 0) {
      fuzzy_candidates <- stringdist_left_join(
      unmatched %>% dplyr::select(bacteria_key, ID, name_type, bacteria_name),
      clover_bacteria_proc,
      by = "bacteria_key",
      method = "jw",   # Jaro-Winkler distance
      max_dist = 0.2,   # Slightly higher threshold for bacteria names
      distance_col = "dist"
    ) %>%
      group_by(bacteria_key.x, ID, name_type) %>%
      slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      dplyr::select(ID, name_type, bacteria_name, Host, HostTaxID, HostGenus, HostFamily, HostOrder, HostClass, HostNCBIResolved,
             Pathogen, PathogenTaxID, PathogenType, PathogenClass, PathogenOrder, 
             PathogenFamily, PathogenGenus, PathogenNCBIResolved,
             DetectionMethod, DetectionMethodOriginal, ICTVRatified,
             Database, DatabaseVersion, DatabaseDOI, PublicationYear, ReferenceText, PMID,
             ReleaseYear, AssocID, NCBIAccession, dist) %>%
      mutate(match_source = "fuzzy")
  
  fuzzy_matches <- fuzzy_candidates
} else {
  fuzzy_matches <- tibble()
}

# ------------------------------| Combine results |----------------------------
# Combine all matches - each row is now a WHO bacteria - CLOVER host association
all_host_associations <- bind_rows(exact_matches, manual_matches, fuzzy_matches)
all_host_associations = unique(all_host_associations)
names(all_host_associations)
unique(all_host_associations$match_source)

# Process matches: prioritize exact over fuzzy
final_host_associations_raw <- all_host_associations %>%
  arrange(ID, match_source, dist) %>%
  group_by(ID) %>%
  filter(
    # If we have exact matches for this ID, only keep exact matches
    if(any(match_source == "exact")) {
      match_source == "exact"
    } else {
      # Otherwise, keep manual and fuzzy matches
      TRUE
    }
  ) %>%
  ungroup() %>%
  # Reorder columns for better readability
  dplyr::select(
    # WHO pathogen information
    ID, bacteria_name, name_type, match_source, dist,
    # CLOVER pathogen information  
    PathogenTaxID, Pathogen, PathogenType, PathogenClass, PathogenOrder, 
    PathogenFamily, PathogenGenus, PathogenNCBIResolved,
    # Host information
    Host, HostTaxID, HostGenus, HostFamily, HostOrder, HostClass, HostNCBIResolved,
    # Detection and provenance
    DetectionMethod, DetectionMethodOriginal, ICTVRatified,
    Database, DatabaseVersion, DatabaseDOI, PublicationYear, ReferenceText, PMID,
    ReleaseYear, AssocID, NCBIAccession
  )

# Collapse to one row per pathogen/host association
final_host_associations <- final_host_associations_raw %>%
  group_by(ID, bacteria_name, PathogenTaxID, Pathogen, Host, HostTaxID) %>%
  summarise(
    # Keep first values for pathogen info
    name_type = first(name_type),
    match_source = paste(unique(match_source), collapse = "; "),
    dist = paste(unique(dist), collapse = "; "),
    PathogenType = paste(unique(PathogenType), collapse = "; "),
    PathogenClass = paste(unique(PathogenClass), collapse = "; "),
    PathogenOrder = paste(unique(PathogenOrder), collapse = "; "),
    PathogenFamily = paste(unique(PathogenFamily), collapse = "; "),
    PathogenGenus = paste(unique(PathogenGenus), collapse = "; "),
    PathogenNCBIResolved = paste(unique(PathogenNCBIResolved), collapse = "; "),
    # Keep first values for host taxonomy
    HostGenus = first(HostGenus),
    HostFamily = first(HostFamily),
    HostOrder = first(HostOrder),
    HostClass = first(HostClass),
    HostNCBIResolved = first(HostNCBIResolved),
    # Collapse metadata columns with semicolons
    DetectionMethod = paste(unique(DetectionMethod), collapse = "; "),
    DetectionMethodOriginal = paste(unique(DetectionMethodOriginal), collapse = "; "),
    ICTVRatified = paste(unique(ICTVRatified), collapse = "; "),
    Database = paste(unique(Database), collapse = "; "),
    DatabaseVersion = paste(unique(DatabaseVersion), collapse = "; "),
    DatabaseDOI = paste(unique(DatabaseDOI), collapse = "; "),
    PublicationYear = paste(unique(PublicationYear), collapse = "; "),
    ReferenceText = paste(unique(ReferenceText), collapse = "; "),
    PMID = paste(unique(PMID), collapse = "; "),
    ReleaseYear = paste(unique(ReleaseYear), collapse = "; "),
    AssocID = paste(unique(AssocID), collapse = "; "),
    NCBIAccession = paste(unique(NCBIAccession), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Reorder columns for better readability
  dplyr::select(
    # WHO pathogen information
    ID, bacteria_name, name_type, match_source, dist,
    # CLOVER pathogen information  
    PathogenTaxID, Pathogen, PathogenType, PathogenClass, PathogenOrder, 
    PathogenFamily, PathogenGenus, PathogenNCBIResolved,
    # Host information
    Host, HostTaxID, HostGenus, HostFamily, HostOrder, HostClass, HostNCBIResolved,
    # Detection and provenance
    DetectionMethod, DetectionMethodOriginal, ICTVRatified,
    Database, DatabaseVersion, DatabaseDOI, PublicationYear, ReferenceText, PMID,
    ReleaseYear, AssocID, NCBIAccession
  )

# Join back with WHO metadata
final_output <- who_df %>%
  left_join(
    final_host_associations %>%
      group_by(ID) %>%
      summarise(
        PathogenTaxID = paste(unique(PathogenTaxID), collapse = "; "),
        Clover_PathogenName = paste(unique(Pathogen), collapse = "; "),
        Clover_PathogenClass = paste(unique(PathogenClass), collapse = "; "),
        match_source = paste(unique(match_source), collapse = "; "),
        num_clover_matches = n_distinct(PathogenTaxID),
        num_host_associations = n(),
        .groups = "drop"
      ),
    by = "ID"
  )


# Save CSV files -------------------------------------------------------------
# Create output directory if it doesn't exist
dir.create(dirname(output_csv_path), recursive = TRUE, showWarnings = FALSE)

# Save WHO-CLOVER matches summary
write_csv(final_output, output_csv_path)

# Save detailed host associations
write_csv(final_host_associations, output_hosts_path)

# Create host species summary
host_species <- final_host_associations %>%
  dplyr::select(Host, HostGenus, HostFamily, HostNCBIResolved) %>%
  mutate(
    Host = str_to_sentence(Host),
    HostGenus = str_to_sentence(HostGenus),
    HostFamily = str_to_sentence(HostFamily)
  ) %>%
  distinct()

# Save detailed host associations
write_csv(host_species, output_unique_hosts_path)

# ------------------------------| Console summary |---------------------------
cat("Match summary:\n")
cat(" - Total WHO bacteria records: ", nrow(who_df), "\n", sep = "")
cat(" - Successfully matched: ", sum(!is.na(final_output$PathogenTaxID)), "\n", sep = "")
cat(" - Unmatched: ", sum(is.na(final_output$PathogenTaxID)), "\n", sep = "")
cat(" - Total pathogen-host associations found: ", nrow(final_host_associations), "\n", sep = "")
cat(" - Unique host species: ", n_distinct(final_host_associations$Host), "\n", sep = "")

if (any(is.na(final_output$PathogenTaxID))) {
  cat("\nUnmatched WHO bacteria names:\n")
  print(final_output$Pathogens[is.na(final_output$PathogenTaxID)])
}

# Show successful matches
if (any(!is.na(final_output$PathogenTaxID))) {
  cat("\nSuccessful matches:\n")
  successful_matches <- final_output %>%
    filter(!is.na(PathogenTaxID)) %>%
    select(Pathogens, Clover_PathogenName, match_source, num_host_associations)
  print(successful_matches)
}

# Summary statistics for host associations
if (nrow(final_host_associations) > 0) {
  hosts_per_pathogen <- final_host_associations %>%
    group_by(bacteria_name) %>%
    summarise(
      n_unique_hosts = n_distinct(Host),
      n_host_families = n_distinct(HostFamily),
      n_host_orders = n_distinct(HostOrder),
      .groups = "drop"
    ) %>%
    arrange(desc(n_unique_hosts))
  
  # Host class distribution
  host_class_summary <- final_host_associations %>%
    count(HostClass, sort = TRUE, name = "n_associations")
  
  cat("\nTop bacteria by host diversity:\n")
  print(head(hosts_per_pathogen, 10))
  
  cat("\nHost class distribution:\n")
  print(host_class_summary)
}

cat("\nProcessing complete!\n")
cat("Files saved:\n")
cat("  - WHO-CLOVER summary:", output_csv_path, "\n")
cat("  - Detailed host associations:", output_hosts_path, "\n")

# Display column descriptions for reference
cat("\nCLOVER Column Descriptions (first 10):\n")
print(head(clover_col_desc, 10))
