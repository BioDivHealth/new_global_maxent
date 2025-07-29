# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(tidyverse)
library(here)
library(pacman)
p_load(fuzzyjoin, stringdist)

# ------------------------------------------------------------------------------|
#      Load and combine WHO disease data --------------------------------------
# ------------------------------------------------------------------------------|
# Read all WHO CSV files (excluding translation.csv)
csv_files <- list.files(path = here("data_artur","WHO","who_diseases"), pattern = "\\.csv$", full.names = TRUE) %>%
  .[!grepl("translation.csv|final_pathogen_data.csv|disease_names.csv", .)]

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

# Combine all tables and add region column
who_diseases_all <- imap_dfr(document_tables, ~mutate(.x, Region = str_remove(.y, "_table$"))) %>%
  # Remove rows where both pathogen columns are empty
  filter(
    !(is.na(`Priority Pathogens`) | `Priority Pathogens` == "") |
    !(is.na(`Prototype Pathogens`) | `Prototype Pathogens` == "")
  )

# ------------------------------------------------------------------------------|
#      Pathogen name standardization ------------------------------------------
# ------------------------------------------------------------------------------|
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

# Apply standardization to pathogen columns
who_diseases_all <- who_diseases_all %>%
  mutate(
    `Priority Pathogens` = ifelse(!is.na(`Priority Pathogens`), standardize_pathogen_name(`Priority Pathogens`), NA),
    `Prototype Pathogens` = ifelse(!is.na(`Prototype Pathogens`), standardize_pathogen_name(`Prototype Pathogens`), NA)
  )

# Create list of all unique pathogens
pathogens_all <- who_diseases_all %>%
  select(`Priority Pathogens`, `Prototype Pathogens`) %>%
  pivot_longer(everything(), values_to = "Pathogens") %>%
  filter(!is.na(Pathogens) & Pathogens != "") %>%
  pull(Pathogens) %>%
  unique()

# Create pathogen-family-risk mapping
pathogens_with_family_risk <- who_diseases_all %>%
  pivot_longer(
    cols = c(`Priority Pathogens`, `Prototype Pathogens`),
    names_to = "Pathogen_Type",
    values_to = "Pathogens"
  ) %>%
  filter(!is.na(Pathogens) & Pathogens != "") %>%
  select(Family, `PHEIC risk`, Pathogens) %>%
  distinct()

# ------------------------------------------------------------------------------|
#      Load translation data and create mapping -------------------------------
# ------------------------------------------------------------------------------|
translation <- read_csv(here("data_artur","WHO","who_diseases", "translation.csv"))
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

# ------------------------------------------------------------------------------|
#      Final results and summary -----------------------------------------------
# ------------------------------------------------------------------------------|
# Merge pathogen mapping with family and risk data
final_pathogen_data <- pathogens_with_family_risk %>%
  left_join(pathogen_mapping, by = c("Pathogens" = "pathogen"))

# Save final_pathogen_data to csv
write_csv(final_pathogen_data, here("data_artur","WHO","who_diseases","final_pathogen_data.csv"))


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

# ------------------------------------------------------------------------------|
#      Add disease names to final_pathogen_data -------------------------------
# ------------------------------------------------------------------------------|

# Read final pathogen data
final_pathogen_data = read_csv(here("data_artur","WHO","who_diseases","final_pathogen_data.csv"))

# Read disease names
diseases = read_csv(here("data_artur","WHO","who_diseases","disease_names.csv"))
diseases = diseases %>% distinct()
# Check which missing pathogens are in the disease names
missing_pathogens = diseases$Pathogens[!diseases$Pathogens %in% final_pathogen_data$Pathogens]

# Add disease names to final_pathogen_data
final_pathogen_data = final_pathogen_data %>%
  left_join(diseases, by = c("Pathogens" = "Pathogens"))

dim(final_pathogen_data)

# Save final_pathogen_data to csv
write_csv(final_pathogen_data, here("data_artur","WHO","who_diseases","who_pathogens_diseases.csv"))