#!/usr/bin/env Rscript
################################################################################
# 1_1_WHO_Diseases_Zoonotic_Filter.R
################################################################################
# Purpose: Derive a conservative zoonotic-focused WHO pathogen table from the
#          broader `who_pathogens_diseases.csv` source.
# Input  : `who_raw_pathogens_path()`
#          `who_disease_names_path()`
# Output : `who_zoonotic_pathogens_path()`
#
# Notes  : The filter is intentionally conservative. It removes:
#          - clear duplicate / case-only duplicate rows
#          - obvious human-only, healthcare-associated, waterborne, or
#            vaccine/laboratory-associated entries
#          - broad influenza subtype umbrella rows when a more specific row is
#            already present in the source table
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "disease_scope_helpers.R"
))

character_cols <- c(
  "Family",
  "PHEIC risk",
  "Pathogens",
  "previous_name",
  "msl39_viral_name",
  "Disease_name",
  "priority_prototype_status"
)

provenance_cols <- disease_scope_provenance_cols()

# ------------------------------------------------------------------------------|
#      Paths and lookup tables -------------------------------------------------|
# ------------------------------------------------------------------------------|
source_path <- who_raw_pathogens_path()
disease_lookup_path <- who_disease_names_path()
output_path <- who_pathogens_diseases_zoonotic_path()

excluded_pathogens <- c(
  "Klebsiella pneumoniae",
  "Shigella dysenteriae serotype 1",
  "Vibrio cholerae serogroup 0139",
  "Enterovirus coxsackiepol",
  "Human polioviruses",
  "Enterovirus alphacoxsackie 71",
  "Enterovirus A71 (EV-A71)",
  "Enterovirus deconjucti 68",
  "Enterovirus D68 (EV-D68)",
  "Lentivirus humimdef1",
  "Human immunodeficiency virus 1 (HIV-1)",
  "Recombinant mastadenovirus",
  "Mastadenovirus blackbeardi serotype 14",
  "Human mastadenovirus B",
  "Metapneumovirus hominis",
  "Mamastrovirus virginiaense",
  "Mamastrovirus 9 (GII.B-human)",
  "Orthopicobirnavirus hominis",
  "Human picobirnavirus",
  "Genus Rotavirus",
  "Orthopoxvirus vaccinia",
  "Orthohepadnavirus hominoidei",
  "Protoparvovirus carnivoran",
  "Carivore protoparvoviruses (CPV)"
)

overlapping_influenza_rows <- c(
  "Alphainfluenzavirus influenzae (H2Nx)",
  "Alphainfluenzavirus influenzae (H5Nx)",
  "Alphainfluenzavirus influenzae (H6Nx)",
  "Alphainfluenzavirus influenzae (H7Nx)",
  "Alphainfluenzavirus influenzae (H10Nx)"
)

disease_lookup <- readr::read_csv(disease_lookup_path, show_col_types = FALSE) %>%
  mutate(
    Pathogens = str_squish(Pathogens),
    pathogen_key = str_to_lower(Pathogens)
  ) %>%
  distinct(pathogen_key, .keep_all = TRUE) %>%
  transmute(
    pathogen_key,
    disease_name_lookup = na_if(Disease_name, "NA")
  )

# ------------------------------------------------------------------------------|
#      Clean and filter the WHO source table -----------------------------------|
# ------------------------------------------------------------------------------|
who_pathogens_zoonotic <- readr::read_csv(source_path, show_col_types = FALSE) %>%
  mutate(
    source_row = row_number(),
    across(any_of(character_cols), ~na_if(.x, "NA")),
    across(any_of(character_cols), ~ifelse(is.na(.x), NA_character_, str_squish(.x))),
    Pathogens = case_when(
      str_to_lower(Pathogens) == "subgenus sarbecovirus" ~ "Subgenus Sarbecovirus",
      str_to_lower(Pathogens) == "subgenus merbecovirus" ~ "Subgenus Merbecovirus",
      TRUE ~ Pathogens
    ),
    pathogen_key = str_to_lower(Pathogens)
  ) %>%
  left_join(disease_lookup, by = "pathogen_key") %>%
  mutate(Disease_name = coalesce(Disease_name, disease_name_lookup)) %>%
  select(-disease_name_lookup) %>%
  filter(
    !Pathogens %in% excluded_pathogens,
    !Pathogens %in% overlapping_influenza_rows,
    !is.na(Disease_name)
  ) %>%
  arrange(source_row) %>%
  group_by(pathogen_key) %>%
  arrange(
    desc(!is.na(Disease_name)),
    desc(!is.na(msl39_viral_name)),
    desc(!is.na(previous_name)),
    source_row,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(source_row) %>%
  select(
    Family,
    `PHEIC risk`,
    Pathogens,
    previous_name,
    msl39_viral_name,
    Disease_name,
    any_of(provenance_cols)
  )

# ------------------------------------------------------------------------------|
#      QA checks ----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
duplicate_pathogens <- who_pathogens_zoonotic %>%
  count(Pathogens, sort = TRUE) %>%
  filter(n > 1)

duplicate_diseases <- who_pathogens_zoonotic %>%
  count(Disease_name, sort = TRUE) %>%
  filter(n > 1)

if (nrow(duplicate_pathogens) > 0) {
  stop(
    "Duplicate pathogen rows remain in `who_pathogens_diseases_zoonotic.csv`: ",
    paste(duplicate_pathogens$Pathogens, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(duplicate_diseases) > 0) {
  stop(
    "Duplicate disease rows remain in `who_pathogens_diseases_zoonotic.csv`: ",
    paste(duplicate_diseases$Disease_name, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------|
#      Write output -------------------------------------------------------------|
# ------------------------------------------------------------------------------|
readr::write_csv(who_pathogens_zoonotic, output_path, na = "")

cat("Wrote zoonotic WHO pathogen table to:\n")
cat(output_path, "\n\n")
cat("Rows kept:", nrow(who_pathogens_zoonotic), "\n")
