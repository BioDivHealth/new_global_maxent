# -----------------------------------------------------------------------------|
# 1_WHO_Diseases.R ----
# -----------------------------------------------------------------------------|
# Purpose: Build the consolidated WHO pathogen backbone from regional tables,
#          translation lookups, disease names, and source-presence flags.
# Inputs : WHO regional pathogen tables, translation table, disease-name lookup,
#          and source-presence lookup.
# Outputs: who_pathogens_diseases.csv
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# 1. Load required libraries and path helpers ----
# -----------------------------------------------------------------------------|
library(tidyverse)
library(here)
library(pacman)
p_load(fuzzyjoin, stringdist)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "disease_scope_helpers.R"
))

region_levels <- c(
  "africa",
  "americas",
  "europe",
  "mediterranean",
  "se_asia",
  "western_pacific"
)

# -----------------------------------------------------------------------------|
# 2. Load and combine WHO disease data ----
# -----------------------------------------------------------------------------|
csv_files <- who_diseases_regional_table_paths(region_levels)

# Read and combine all tables
document_tables <- csv_files %>%
  set_names(~basename(.) %>% tools::file_path_sans_ext()) %>%
  map(read_csv) %>%
  # Standardize column names: harmonize 'Family Risk' to 'PHEIC risk'
  map(~{
    df <- .x
    if ("Family Risk" %in% names(df)) {
      names(df)[names(df) == "Family Risk"] <- "PHEIC risk"
    }
    df
  })

# Combine all tables and add source-region column
who_diseases_all <- imap_dfr(document_tables, ~mutate(.x, source_region = str_remove(.y, "_table$"))) %>%
  # Remove rows where both pathogen columns are empty
  filter(
    !(is.na(`Priority Pathogens`) | `Priority Pathogens` == "") |
    !(is.na(`Prototype Pathogens`) | `Prototype Pathogens` == "")
  )

# -----------------------------------------------------------------------------|
# 3. Standardize pathogen names and regional provenance ----
# -----------------------------------------------------------------------------|
# Function to standardize pathogen names
standardize_pathogen_name <- function(x) {
  x %>%
    str_trim() %>% # Remove leading/trailing spaces
    str_replace_all(" +", " ") %>% # Remove double spaces
    str_replace_all("\\( ", "(") %>% # Remove space after (
    str_replace_all(" \\)", ")") %>% # Remove space before )
    str_replace_all("[lI]nfluenzae H[lI]", "Influenzae H1") %>% # Fix H1 typos
    str_replace_all("[lI]nfluenzae Hl", "Influenzae H1") %>% # Another H1 typo
    str_replace_all("[lI]nfluenzae HIN1", "Influenzae H1N1") %>% # Fix H1N1 typo
    str_replace_all("[lI]nfluenzae H10Nx", "Influenzae H10Nx") %>% # Standardize H10Nx
    str_replace_all("Orthoeb olavirus", "Orthoebolavirus") %>% # Remove space typo
    str_replace_all("Orthopicobimavirus", "Orthopicobirnavirus") %>% # Fix typo
    str_replace_all("Orthonairovirus haemorhagiae", "Orthonairovirus haemorrhagiae") %>% # Fix typo
    str_replace_all("Paslahepevirus balayani, genotype 3", "Paslahepevirus balayani genotype 3") %>% # Remove comma
    str_replace_all("Vibrio cholera \\(0139\\)", "Vibrio cholerae serogroup 0139") %>% # Standardize Vibrio cholerae
    str_replace_all("Subgenus Sarbecoviruses", "Subgenus Sarbecovirus") %>% # Singular
    str_replace_all("Lentivirus humimdef[ 1lI]", "Lentivirus humimdef1") %>% # Standardize Lentivirus
    str_replace_all("Mammarenavirus lassa ense", "Mammarenavirus lassaense") %>% # Remove space
    str_replace_all("Influenzae h", "Influenzae H") %>% # Capitalize H after Influenzae
    str_replace_all("\u2013", "-") %>% # Replace en-dash with hyphen if present
    # Standardize influenza subtypes
    str_replace_all("Alphainfluenzavirus Influenzae H([0-9]+)", "Alphainfluenzavirus influenzae (H\\1N1)") %>%
    str_replace_all("Alphainfluenzavirus influenzae \\(H([0-9]+)N([0-9xX]+)\\)", "Alphainfluenzavirus influenzae (H\\1N\\2)") %>%
    # Capitalize abbreviations
    str_replace_all("hiv-1", "HIV-1") %>%
    str_replace_all("ev-a71", "EV-A71") %>%
    str_replace_all("ev-d68", "EV-D68") %>%
    str_replace_all("gil.b-human", "GII.B-human") %>%
    str_replace_all("Yersinia Pestis", "Yersinia pestis") %>%
    # Fix standalone 'encephalitidis'
    {ifelse(tolower(.) == "encephalitidis", NA, .)}
}

region_status_for <- function(source_region, source_pathogen_type, region_name) {
  has_priority <- any(
    source_region == region_name & source_pathogen_type == "priority",
    na.rm = TRUE
  )
  has_prototype <- any(
    source_region == region_name & source_pathogen_type == "prototype",
    na.rm = TRUE
  )

  case_when(
    has_priority & has_prototype ~ "both",
    has_priority ~ "priority",
    has_prototype ~ "prototype",
    TRUE ~ "none"
  )
}

# Apply standardization to pathogen columns
who_diseases_all <- who_diseases_all %>%
  mutate(
    `Priority Pathogens` = ifelse(!is.na(`Priority Pathogens`), standardize_pathogen_name(`Priority Pathogens`), NA),
    `Prototype Pathogens` = ifelse(!is.na(`Prototype Pathogens`), standardize_pathogen_name(`Prototype Pathogens`), NA)
  )

who_diseases_long <- who_diseases_all %>%
  pivot_longer(
    cols = c(`Priority Pathogens`, `Prototype Pathogens`),
    names_to = "source_pathogen_type",
    values_to = "Pathogens"
  ) %>%
  filter(!is.na(Pathogens) & Pathogens != "") %>%
  mutate(
    source_pathogen_type = recode(
      source_pathogen_type,
      `Priority Pathogens` = "priority",
      `Prototype Pathogens` = "prototype"
    ),
    source_region = factor(source_region, levels = region_levels)
  )

# Create list of all unique pathogens
pathogens_all <- who_diseases_long %>%
  pull(Pathogens) %>%
  unique()

# Create pathogen-family-risk mapping plus provenance
pathogens_with_family_risk <- who_diseases_long %>%
  group_by(Pathogens) %>%
  summarise(
    Family = disease_scope_first_non_missing(Family),
    `PHEIC risk` = disease_scope_first_non_missing(`PHEIC risk`),
    is_priority_pathogen = any(source_pathogen_type == "priority"),
    is_prototype_pathogen = any(source_pathogen_type == "prototype"),
    region_africa = region_status_for(source_region, source_pathogen_type, "africa"),
    region_americas = region_status_for(source_region, source_pathogen_type, "americas"),
    region_europe = region_status_for(source_region, source_pathogen_type, "europe"),
    region_mediterranean = region_status_for(source_region, source_pathogen_type, "mediterranean"),
    region_se_asia = region_status_for(source_region, source_pathogen_type, "se_asia"),
    region_western_pacific = region_status_for(source_region, source_pathogen_type, "western_pacific"),
    .groups = "drop"
  ) %>%
  mutate(
    priority_prototype_status = case_when(
      is_priority_pathogen & is_prototype_pathogen ~ "both",
      is_priority_pathogen ~ "priority",
      is_prototype_pathogen ~ "prototype",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Family,
    `PHEIC risk`,
    Pathogens,
    is_priority_pathogen,
    is_prototype_pathogen,
    priority_prototype_status,
    region_africa,
    region_americas,
    region_europe,
    region_mediterranean,
    region_se_asia,
    region_western_pacific
  )

# -----------------------------------------------------------------------------|
# 4. Map WHO names to translation-table names ----
# -----------------------------------------------------------------------------|
translation <- read_csv(who_diseases_translation_path())
names(translation) <- c("Family", "Previous_Name", "MSL39_Viral_Species_Name")

# Define manual fuzzy matches
manual_fuzzy_matches <- tribble(
  ~pathogen, ~matched_column,
  "Carnivore protoparvoviruses (CPV)", "Previous_Name",
  "Mamastrovirus 9 (GIl.B-human)", "Previous_Name",
  "Mastadenovirus blackbeardi serotype 14", "MSL39_Viral_Species_Name",
  "Paslahepevirus balayani genotype 3", "MSL39_Viral_Species_Name"
)

# Helper function to check if pathogen is Alphainfluenzavirus variant
is_alpha_influenza <- function(pathogen) {
  str_detect(pathogen, "^Alphainfluenzavirus influenzae \\(H[0-9]+N[0-9xX]+\\)$")
}

# Create comprehensive pathogen mapping
pathogen_mapping <- tibble(pathogen = pathogens_all) %>%
  rowwise() %>%
  mutate(
    # Check exact matches
    exact_prev = pathogen %in% translation$Previous_Name | 
      (is_alpha_influenza(pathogen) & "Alphainfluenzavirus influenzae" %in% translation$Previous_Name),
    exact_msl39 = pathogen %in% translation$MSL39_Viral_Species_Name | 
      (is_alpha_influenza(pathogen) & "Alphainfluenzavirus influenzae" %in% translation$MSL39_Viral_Species_Name),
    
    # Get fuzzy matches for non-exact matches
    fuzzy_prev = if (!exact_prev && !is.na(pathogen)) {
      dists <- stringdist(pathogen, translation$Previous_Name, method = "jw")
      translation$Previous_Name[which.min(dists)]
    } else NA_character_,
    
    fuzzy_msl39 = if (!exact_msl39 && !is.na(pathogen)) {
      dists <- stringdist(pathogen, translation$MSL39_Viral_Species_Name, method = "jw")
      translation$MSL39_Viral_Species_Name[which.min(dists)]
    } else NA_character_,
    
    # Apply manual matches
    manual_match = manual_fuzzy_matches$matched_column[match(pathogen, manual_fuzzy_matches$pathogen)],
    
    # Final mapping logic
    previous_name = case_when(
      exact_prev ~ if (is_alpha_influenza(pathogen)) "Alphainfluenzavirus influenzae" else pathogen,
      !is.na(manual_match) && manual_match == "Previous_Name" ~ fuzzy_prev,
      !is.na(manual_match) && manual_match == "MSL39_Viral_Species_Name" ~ {
        idx <- match(fuzzy_msl39, translation$MSL39_Viral_Species_Name)
        if (!is.na(idx)) translation$Previous_Name[idx] else NA_character_
      },
      exact_msl39 ~ {
        base_name <- if (is_alpha_influenza(pathogen)) "Alphainfluenzavirus influenzae" else pathogen
        idx <- match(base_name, translation$MSL39_Viral_Species_Name)
        if (!is.na(idx)) translation$Previous_Name[idx] else NA_character_
      },
      TRUE ~ NA_character_
    ),
    
    msl39_viral_name = case_when(
      exact_msl39 ~ if (is_alpha_influenza(pathogen)) "Alphainfluenzavirus influenzae" else pathogen,
      !is.na(manual_match) && manual_match == "MSL39_Viral_Species_Name" ~ fuzzy_msl39,
      !is.na(manual_match) && manual_match == "Previous_Name" ~ {
        idx <- match(fuzzy_prev, translation$Previous_Name)
        if (!is.na(idx)) translation$MSL39_Viral_Species_Name[idx] else NA_character_
      },
      exact_prev ~ {
        base_name <- if (is_alpha_influenza(pathogen)) "Alphainfluenzavirus influenzae" else pathogen
        idx <- match(base_name, translation$Previous_Name)
        if (!is.na(idx)) translation$MSL39_Viral_Species_Name[idx] else NA_character_
      },
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(pathogen, previous_name, msl39_viral_name)

# -----------------------------------------------------------------------------|
# 5. Write interim pathogen mapping and summary ----
# -----------------------------------------------------------------------------|
# Merge pathogen mapping with family and risk data
final_pathogen_data <- pathogens_with_family_risk %>%
  left_join(pathogen_mapping, by = c("Pathogens" = "pathogen"))

# Save final_pathogen_data to csv
write_csv(final_pathogen_data, who_final_pathogen_data_path())

# Summary statistics
mapped_count <- sum(!is.na(pathogen_mapping$previous_name) | !is.na(pathogen_mapping$msl39_viral_name))
total_count <- nrow(pathogen_mapping)

cat("Pathogen mapping summary:\n")
cat("- Total pathogens:", total_count, "\n")
cat("- Successfully mapped:", mapped_count, "\n")
cat("- Unmapped:", total_count - mapped_count, "\n")

# Show unmapped pathogens for review
unmapped_pathogens <- pathogen_mapping %>% 
  filter(is.na(previous_name) & is.na(msl39_viral_name))

if (nrow(unmapped_pathogens) > 0) {
  cat("\nUnmapped pathogens requiring attention:\n")
  print(unmapped_pathogens$pathogen)
}

# -----------------------------------------------------------------------------|
# 6. Add disease names and source-presence flags ----
# -----------------------------------------------------------------------------|

# Read final pathogen data
final_pathogen_data = read_csv(who_final_pathogen_data_path())

# Read disease names
diseases = read_csv(who_disease_names_path())
diseases = diseases %>% distinct()

point_data_lookup <- read_csv(
  who_diseases_gibb_lookup_path(),
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  transmute(
    Disease_name = source_disease_name,
    in_gibb_etal = coalesce(in_gibb_etal, FALSE),
    in_empres_i = coalesce(in_empres_i, FALSE)
  ) %>%
  distinct(Disease_name, .keep_all = TRUE)

# Add disease names to final_pathogen_data
final_pathogen_data = final_pathogen_data %>%
  left_join(diseases, by = c("Pathogens" = "Pathogens")) %>%
  left_join(point_data_lookup, by = "Disease_name") %>%
  mutate(
    in_gibb_etal = coalesce(in_gibb_etal, FALSE),
    in_empres_i = coalesce(in_empres_i, FALSE)
  )

final_pathogen_data <- final_pathogen_data %>%
  select(
    Family,
    `PHEIC risk`,
    Pathogens,
    previous_name,
    msl39_viral_name,
    Disease_name,
    in_gibb_etal,
    in_empres_i,
    is_priority_pathogen,
    is_prototype_pathogen,
    priority_prototype_status,
    region_africa,
    region_americas,
    region_europe,
    region_mediterranean,
    region_se_asia,
    region_western_pacific
  )

# Save final_pathogen_data to csv
write_csv(final_pathogen_data, who_raw_pathogens_path())
