# ------------------------------------------------------------------------------
# 2_Match_WHO_Virion.R
# ------------------------------------------------------------------------------
# Purpose: Identify which WHO priority pathogens are present in the VIRION
#          database and retrieve their corresponding VirusTaxID identifiers.
#
# Output : A CSV file `who_pathogens_virion_taxid.csv` saved into the WHO
#          documents folder containing all exact (and selected fuzzy) matches.
# ------------------------------------------------------------------------------|

# ------------------------------| Load libraries |------------------------------
library(tidyverse)
library(stringdist)
library(fuzzyjoin)
library(magrittr)

# ------------------------------| Helper paths  |------------------------------
who_csv_path   <- file.path("data_artur", "WHO", "who_diseases", "who_pathogens_diseases.csv")
output_csv_path <- file.path("data_artur", "WHO", "virion", "who_pathogens_virion_taxid.csv")

# ------------------------------| Load datasets |------------------------------
# 1. WHO pathogen list ---------------------------------------------------------
who_df <- read_csv(who_csv_path, show_col_types = FALSE)
who_df %<>% filter(Family!="Bacteria") # Exclude bacteria as per original script
who_df$ID = 1:nrow(who_df)

# 2. VIRION taxonomy -----------------------------------------------------------
if (!exists("virion_data")) {
  source(file.path("scripts", "pathogen_associations", "virion_data.R"))
  # `virion_data` object is created inside the sourced script
}

taxonomy_virus <- virion_data$taxonomy_virus

# ------------------------------| Pre-processing |-----------------------------
# Create a lowercase, trimmed helper column for safer joins --------------------
who_long <- who_df %>%
  # Gather all possible virus name columns into one
  select(ID, Pathogens, previous_name, msl39_viral_name, Family) %>%
  pivot_longer(cols = c(Pathogens, previous_name, msl39_viral_name), names_to = "name_type", values_to = "virus_name") %>%
  filter(!is.na(virus_name) & virus_name != "") %>%
  mutate(virus_key = str_to_lower(str_trim(virus_name)))
  # Remove distinct() to preserve all IDs with the same virus name

taxonomy_virus_proc <- taxonomy_virus %>%
  mutate(virus_key = str_to_lower(str_trim(Virus))) %>%
  select(VirusTaxID, Virus, virus_key, VirusFamily, Database)

# ------------------------------| Manual matches |-----------------------------
# Define manual mappings for specific cases where fuzzy matching needs correction
manual_mappings <- tribble(
  ~who_name_lower, ~virion_name_lower,
  # SARS-CoV related
  "subgenus sarbecovirus", "severe acute respiratory syndrome-related coronavirus",
  "subgenus sarbecovirus", "betacoronavirus pandemicum", 
  "severe acute respiratory syndrome coronavirus", "severe acute respiratory syndrome-related coronavirus",
  "severe acute respiratory syndrome coronavirus", "betacoronavirus pandemicum",
  # MERS-CoV related  
  "subgenus merbecovirus", "middle east respiratory syndrome-related coronavirus",
  "subgenus merbecovirus", "betacoronavirus cameli",
  "middle east respiratory syndrome coronavirus", "middle east respiratory syndrome-related coronavirus",
  "middle east respiratory syndrome coronavirus", "betacoronavirus cameli",
  # Add more manual mappings as needed
  # Human Mastadenovirus B
  "human mastadenovirus b", "mastadenovirus blackbeardi",
  # human polioviruses
  "human polioviruses", "enterovirus c",
  "human polioviruses", "human poliovirus sp.",
  # carnivore parvoviruses
  "protoparvovirus carnivoran", "protoparvovirus carnivoran1",
  "protoparvovirus carnivoran", "protoparvovirus carnivoran2",
  "protoparvovirus carnivoran", "protoparvovirus carnivoran3",
  "protoparvovirus carnivoran", "protoparvovirus carnivoran4",
  "protoparvovirus carnivoran", "protoparvovirus carnivoran5"
  # MISSING (not found through fuzzy matching)
  
)

# Apply manual mappings
manual_matches <- who_long %>%
  inner_join(manual_mappings, by = c("virus_key" = "who_name_lower")) %>%
  inner_join(taxonomy_virus_proc, by = c("virion_name_lower" = "virus_key")) %>%
  mutate(dist = 0.5, match_source = "manual") %>%  # Flag as manual matches with intermediate distance
  select(ID, name_type, virus_name, VirusTaxID, Virus, VirusFamily, Database, dist, match_source)

# ------------------------------| Exact matches |------------------------------
exact_matches <- who_long %>%
  inner_join(taxonomy_virus_proc, by = "virus_key") %>%
  mutate(dist = 0, match_source = "exact") %>%  # Add distance column and source for exact matches
  select(ID, name_type, virus_name, VirusTaxID, Virus, VirusFamily, Database, dist, match_source)

# ------------------------------| Fuzzy matches |------------------------------
# Attempt fuzzy matching only for those still unmatched ------------------------
unmatched <- who_long %>%
  filter(!ID %in% c(exact_matches$ID, manual_matches$ID))

if (nrow(unmatched) > 0) {
  fuzzy_candidates <- stringdist_left_join(
    unmatched %>% select(virus_key, ID, name_type),
    taxonomy_virus_proc %>% select(virus_key, VirusTaxID, Virus, VirusFamily, Database),
    by = "virus_key",
    method = "jw",   # Jaro-Winkler distance
    max_dist = 0.1,   # adjust threshold as needed
    distance_col = "dist"
  ) %>%
    group_by(virus_key.x, ID, name_type) %>%
    slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(virus_key = virus_key.x, ID, name_type, VirusTaxID, Virus, VirusFamily, Database, dist)

  fuzzy_matches <- unmatched %>%
    left_join(fuzzy_candidates, by = c("virus_key", "ID", "name_type")) %>%
    filter(!is.na(VirusTaxID)) %>%
    mutate(match_source = "fuzzy") %>%
    select(ID, name_type, virus_name, VirusTaxID, Virus, VirusFamily, Database, dist, match_source)
} else {
  fuzzy_matches <- tibble()
}

# ------------------------------| Combine & save |-----------------------------
all_matches <- bind_rows(exact_matches, manual_matches, fuzzy_matches)

# Process matches: prioritize exact over fuzzy, then by name_type, but keep all VirusTaxIDs per ID --
processed_matches <- all_matches %>%
  mutate(
    #match_type = ifelse(ID %in% exact_matches$ID, "exact", "fuzzy"),
    name_priority = case_when(
      name_type == "Pathogens" ~ 1,
      name_type == "msl39_viral_name" ~ 2,
      name_type == "previous_name" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(ID, match_source, name_priority)

# For each ID, if we have exact matches, only keep exact matches
# If we only have fuzzy matches, keep the best fuzzy matches (by name_priority)
final_matches <- processed_matches %>%
  group_by(ID) %>%
  filter(
    # If we have exact matches for this ID, only keep exact matches
    if(any(match_source == "exact")) {
      match_source == "exact"
    } else {
      # Otherwise, keep only the best name_priority fuzzy matches
      name_priority == min(name_priority)
    }
  ) %>%
  ungroup()

# Collapse multiple matches per ID into comma-separated strings ---------------
collapsed_matches <- final_matches %>%
  group_by(ID) %>%
  summarise(
    VirusTaxID = paste(unique(VirusTaxID), collapse = "; "),
    Virion_VirusName = paste(unique(Virus), collapse = "; "),
    Virion_VirusFamily = paste(unique(VirusFamily), collapse = "; "),
    Virion_Database = paste(unique(Database), collapse = "; "),
    matched_name_type = paste(unique(name_type), collapse = "; "),
    matched_virus_name = paste(unique(virus_name), collapse = "; "),
    #match_type = paste(unique(match_type), collapse = "; "),
    match_source = paste(unique(match_source), collapse = "; "),  # Add this line
    fuzzy_scores = if_else(
      any(match_source == "fuzzy"),
      paste(unique(round(dist[match_source == "fuzzy"], 3)), collapse = "; "),
      NA_character_
    ),
    num_virion_matches = n(),
    .groups = "drop"
  )


# Attach metadata to the original WHO table -----------------------------------
final_output <- who_df %>%
  left_join(collapsed_matches, by = "ID")


View(final_output %>% select(Pathogens, previous_name, msl39_viral_name, ID, VirusTaxID, Virion_VirusName, match_type))


# Save CSV --------------------------------------------------------------------
write_csv(final_output, output_csv_path)


# ------------------------------| Console summary |---------------------------
cat("Match summary:\n")
cat(" - Total WHO pathogen records: ", nrow(who_df), "\n", sep = "")
cat(" - Successfully matched (exact + fuzzy): ", sum(!is.na(final_output$VirusTaxID)), "\n", sep = "")
cat(" - Unmatched: ", sum(is.na(final_output$VirusTaxID)), "\n", sep = "")

if (any(is.na(final_output$VirusTaxID))) {
  cat("\nUnmatched WHO pathogen names (first 20 shown):\n")
  print(head(final_output$Pathogens[is.na(final_output$VirusTaxID)], 20))
} 