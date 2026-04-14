#!/usr/bin/env Rscript
################################################################################
# 4_1_Canonicalize_Combined_WHO_Network.R
################################################################################
# Purpose: Canonicalize pathogen names in the merged WHO host-pathogen network,
#          add zoonotic status, and preserve the raw labels for provenance.
#
# Inputs : pathogen_association_data/WHO/networks/combined_who_network.csv
#          pathogen_association_data/WHO/who_diseases/who_pathogens_diseases.csv
#          pathogen_association_data/WHO/who_diseases/disease_names.csv
#          pathogen_association_data/WHO/who_diseases/who_pathogens_diseases_zoonotic.csv
#
# Outputs: pathogen_association_data/WHO/networks/combined_who_pathogen_canonical_lookup.csv
#          pathogen_association_data/WHO/networks/combined_who_network_canonical.csv
#          pathogen_association_data/WHO/networks/combined_who_network_canonical_zoonotic.csv
#
# Notes  : This script does not overwrite the existing raw combined network.
#          It creates a derived canonical artifact so downstream scripts can
#          switch deliberately rather than changing behavior silently.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr, tidyr)

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- ifelse(is.na(x), NA_character_, x)
  x <- na_if(x, "")
  x <- na_if(x, "NA")
  ifelse(is.na(x), NA_character_, str_squish(x))
}

first_non_missing <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = "; ")
}

safe_lower <- function(x) {
  ifelse(is.na(x), NA_character_, str_to_lower(clean_text(x)))
}

normalize_who_pathogen <- function(x) {
  case_when(
    safe_lower(x) == "subgenus sarbecovirus" ~ "Subgenus Sarbecovirus",
    safe_lower(x) == "subgenus merbecovirus" ~ "Subgenus Merbecovirus",
    TRUE ~ clean_text(x)
  )
}

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
who_dir <- here::here("pathogen_association_data", "WHO")

network_path <- file.path(who_dir, "networks", "combined_who_network.csv")
who_path <- file.path(who_dir, "who_diseases", "who_pathogens_diseases.csv")
disease_names_path <- file.path(who_dir, "who_diseases", "disease_names.csv")
zoonotic_path <- file.path(who_dir, "who_diseases", "who_pathogens_diseases_zoonotic.csv")

lookup_output_path <- file.path(who_dir, "networks", "combined_who_pathogen_canonical_lookup.csv")
canonical_output_path <- file.path(who_dir, "networks", "combined_who_network_canonical.csv")
zoonotic_network_output_path <- file.path(
  who_dir, "networks", "combined_who_network_canonical_zoonotic.csv"
)

# ------------------------------------------------------------------------------|
#      Manual overrides for known synonym / specificity cases -------------------|
# ------------------------------------------------------------------------------|
manual_pathogen_map <- tibble::tribble(
  ~Pathogen_raw, ~Pathogen_canonical, ~canonicalization_status,
  "Salmonella enterica", "Salmonella enterica non typhoidal serovars", "manual_specificity_map",
  "Shigella dysenteriae", "Shigella dysenteriae serotype 1", "manual_specificity_map",
  "Vibrio cholerae", "Vibrio cholerae serogroup 0139", "manual_specificity_map",
  "Betacoronavirus pandemicum", "Subgenus Sarbecovirus", "manual_synonym_map",
  "Severe acute respiratory syndrome-related coronavirus", "Subgenus Sarbecovirus", "manual_synonym_map",
  "Betacoronavirus cameli", "Subgenus Merbecovirus", "manual_synonym_map",
  "Zaire ebolavirus", "Orthoebolavirus zairense", "manual_synonym_map",
  "Enterovirus c", "Enterovirus coxsackiepol", "manual_synonym_map",
  "Human poliovirus", "Enterovirus coxsackiepol", "manual_synonym_map",
  "Enterovirus alphacoxsackie", "Enterovirus alphacoxsackie 71", "manual_synonym_map",
  "Enterovirus a", "Enterovirus alphacoxsackie 71", "manual_synonym_map",
  "Enterovirus deconjuncti", "Enterovirus deconjucti 68", "manual_synonym_map",
  "Alphainfluenzavirus influenzae", "Alphainfluenzavirus influenzae", "manual_group_retained",
  "Protoparvovirus carnivoran1", "Protoparvovirus carnivoran", "manual_group_map",
  "Protoparvovirus carnivoran3", "Protoparvovirus carnivoran", "manual_group_map",
  "Protoparvovirus carnivoran4", "Protoparvovirus carnivoran", "manual_group_map",
  "Protoparvovirus carnivoran5", "Protoparvovirus carnivoran", "manual_group_map"
) %>%
  mutate(Pathogen_raw_key = safe_lower(Pathogen_raw))

zoonotic_override <- tibble::tribble(
  ~Pathogen_canonical, ~is_zoonotic_override, ~zoonotic_status_override,
  "Alphainfluenzavirus influenzae", TRUE, "zoonotic_group_retained"
) %>%
  mutate(Pathogen_canonical_key = safe_lower(Pathogen_canonical))

# ------------------------------------------------------------------------------|
#      Load and prepare WHO lookup tables --------------------------------------|
# ------------------------------------------------------------------------------|
disease_names <- read_csv(disease_names_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    Pathogens = normalize_who_pathogen(Pathogens),
    Disease_name = clean_text(Disease_name),
    pathogen_key = safe_lower(Pathogens)
  ) %>%
  distinct(pathogen_key, .keep_all = TRUE) %>%
  select(pathogen_key, disease_name_lookup = Disease_name)

who_pathogens <- read_csv(who_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    Pathogens = normalize_who_pathogen(Pathogens),
    pathogen_key = safe_lower(Pathogens)
  ) %>%
  left_join(disease_names, by = "pathogen_key") %>%
  mutate(Disease_name = coalesce(Disease_name, disease_name_lookup)) %>%
  select(-disease_name_lookup)

who_canonical <- who_pathogens %>%
  transmute(
    Pathogen_canonical = Pathogens,
    pathogen_canonical_key = safe_lower(Pathogens),
    Disease_name_canonical = Disease_name,
    Family_canonical = Family,
    PHEIC_risk_canonical = `PHEIC risk`,
    previous_name_canonical = previous_name,
    msl39_viral_name_canonical = msl39_viral_name
  ) %>%
  group_by(pathogen_canonical_key, Pathogen_canonical) %>%
  summarise(
    Disease_name_canonical = if (n_distinct(Disease_name_canonical, na.rm = TRUE) == 1) {
      first_non_missing(Disease_name_canonical)
    } else {
      NA_character_
    },
    Family_canonical = first_non_missing(Family_canonical),
    PHEIC_risk_canonical = first_non_missing(PHEIC_risk_canonical),
    previous_name_canonical = collapse_unique(previous_name_canonical),
    msl39_viral_name_canonical = collapse_unique(msl39_viral_name_canonical),
    .groups = "drop"
  )

who_alias_lookup <- who_pathogens %>%
  transmute(
    Pathogen_canonical = Pathogens,
    Disease_name_canonical = Disease_name,
    alias_1 = Pathogens,
    alias_2 = previous_name,
    alias_3 = msl39_viral_name
  ) %>%
  pivot_longer(
    cols = starts_with("alias_"),
    names_to = "alias_type",
    values_to = "alias"
  ) %>%
  mutate(alias = clean_text(alias)) %>%
  filter(!is.na(alias)) %>%
  mutate(alias_key = safe_lower(alias)) %>%
  distinct(alias_key, Pathogen_canonical, Disease_name_canonical)

who_alias_resolved <- who_alias_lookup %>%
  group_by(alias_key) %>%
  summarise(
    Pathogen_canonical = if (n_distinct(Pathogen_canonical) == 1) first(Pathogen_canonical) else NA_character_,
    Disease_name_canonical = if (n_distinct(Disease_name_canonical) == 1) first(Disease_name_canonical) else NA_character_,
    alias_candidate_count = n_distinct(Pathogen_canonical),
    .groups = "drop"
  )

zoonotic_lookup <- read_csv(zoonotic_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    Pathogens = normalize_who_pathogen(Pathogens),
    pathogen_canonical_key = safe_lower(Pathogens)
  ) %>%
  distinct(pathogen_canonical_key, .keep_all = TRUE) %>%
  transmute(
    pathogen_canonical_key,
    is_zoonotic_lookup = TRUE,
    zoonotic_status_lookup = "zoonotic"
  )

# ------------------------------------------------------------------------------|
#      Build a canonical lookup for raw network pathogens -----------------------|
# ------------------------------------------------------------------------------|
network_targets <- read_csv(network_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text)) %>%
  distinct(Pathogen, PathogenTaxID, Disease_name) %>%
  mutate(
    Pathogen_raw = Pathogen,
    Disease_name_raw = Disease_name,
    Pathogen_raw_key = safe_lower(Pathogen_raw)
  )

canonical_lookup <- network_targets %>%
  left_join(
    manual_pathogen_map %>%
      select(Pathogen_raw_key, Pathogen_canonical_manual = Pathogen_canonical, canonicalization_status_manual = canonicalization_status),
    by = "Pathogen_raw_key"
  ) %>%
  left_join(
    who_alias_resolved %>%
      rename(
        Pathogen_canonical_alias = Pathogen_canonical,
        Disease_name_canonical_alias = Disease_name_canonical
      ),
    by = c("Pathogen_raw_key" = "alias_key")
  ) %>%
  mutate(
    Pathogen_canonical = coalesce(Pathogen_canonical_manual, Pathogen_canonical_alias, Pathogen_raw),
    canonicalization_status = case_when(
      !is.na(Pathogen_canonical_manual) ~ canonicalization_status_manual,
      !is.na(Pathogen_canonical_alias) & Pathogen_canonical_alias == Pathogen_raw ~ "identity",
      !is.na(Pathogen_canonical_alias) ~ "who_alias_match",
      TRUE ~ "raw_retained_no_match"
    ),
    pathogen_canonical_key = safe_lower(Pathogen_canonical)
  ) %>%
  left_join(who_canonical, by = c("pathogen_canonical_key", "Pathogen_canonical")) %>%
  mutate(
    canonicalization_status = case_when(
      canonicalization_status == "raw_retained_no_match" &
        !is.na(Family_canonical) &
        Pathogen_canonical == Pathogen_raw ~ "identity",
      TRUE ~ canonicalization_status
    ),
    Disease_name_canonical = coalesce(Disease_name_raw, Disease_name_canonical_alias, Disease_name_canonical)
  ) %>%
  left_join(zoonotic_lookup, by = "pathogen_canonical_key") %>%
  left_join(
    zoonotic_override %>%
      select(Pathogen_canonical_key, is_zoonotic_override, zoonotic_status_override),
    by = c("pathogen_canonical_key" = "Pathogen_canonical_key")
  ) %>%
  mutate(
    is_zoonotic = case_when(
      !is.na(is_zoonotic_override) ~ is_zoonotic_override,
      !is.na(is_zoonotic_lookup) ~ is_zoonotic_lookup,
      !is.na(Family_canonical) ~ FALSE,
      TRUE ~ NA
    ),
    zoonotic_status = case_when(
      !is.na(zoonotic_status_override) ~ zoonotic_status_override,
      !is.na(zoonotic_status_lookup) ~ zoonotic_status_lookup,
      !is.na(Family_canonical) ~ "not_zoonotic_or_out_of_scope",
      TRUE ~ "unknown"
    )
  ) %>%
  transmute(
    Pathogen_raw,
    PathogenTaxID,
    Disease_name_raw,
    Pathogen_canonical,
    Disease_name_canonical,
    canonicalization_status,
    Family_canonical,
    PHEIC_risk_canonical,
    previous_name_canonical,
    msl39_viral_name_canonical,
    is_zoonotic,
    zoonotic_status
  ) %>%
  distinct()

# ------------------------------------------------------------------------------|
#      Apply the canonical lookup to the full network ---------------------------|
# ------------------------------------------------------------------------------|
combined_network <- read_csv(network_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    Pathogen_raw = Pathogen,
    Disease_name_raw = Disease_name
  )

combined_network_canonical <- combined_network %>%
  left_join(
    canonical_lookup,
    by = c(
      "Pathogen_raw" = "Pathogen_raw",
      "PathogenTaxID" = "PathogenTaxID",
      "Disease_name_raw" = "Disease_name_raw"
    )
  ) %>%
  mutate(
    Pathogen = coalesce(Pathogen_canonical, Pathogen_raw),
    Disease_name = coalesce(Disease_name_raw, Disease_name_canonical),
    `PHEIC risk` = coalesce(`PHEIC risk`, PHEIC_risk_canonical)
  ) %>%
  group_by(Pathogen, PathogenTaxID, Disease_name, Host, HostTaxID) %>%
  summarise(
    Pathogen_raw_examples = collapse_unique(Pathogen_raw),
    Disease_name_raw_examples = collapse_unique(Disease_name_raw),
    canonicalization_status = collapse_unique(canonicalization_status),
    is_zoonotic = dplyr::first(is_zoonotic),
    zoonotic_status = first_non_missing(zoonotic_status),
    `PHEIC risk` = first_non_missing(`PHEIC risk`),
    PathogenClass = first_non_missing(PathogenClass),
    PathogenOrder = first_non_missing(PathogenOrder),
    PathogenFamily = first_non_missing(PathogenFamily),
    PathogenGenus = first_non_missing(PathogenGenus),
    HostPhylum = first_non_missing(HostPhylum),
    HostClass = first_non_missing(HostClass),
    HostFamily = first_non_missing(HostFamily),
    HostOrder = first_non_missing(HostOrder),
    DetectionMethod = collapse_unique(DetectionMethod),
    MainSource = collapse_unique(MainSource),
    PathogenType = first_non_missing(PathogenType),
    .groups = "drop"
  ) %>%
  select(
    Pathogen,
    Pathogen_raw_examples,
    PathogenTaxID,
    `PHEIC risk`,
    Disease_name,
    Disease_name_raw_examples,
    HostTaxID,
    Host,
    PathogenClass,
    PathogenOrder,
    PathogenFamily,
    PathogenGenus,
    HostPhylum,
    HostClass,
    HostFamily,
    HostOrder,
    DetectionMethod,
    MainSource,
    PathogenType,
    is_zoonotic,
    zoonotic_status,
    canonicalization_status
  ) %>%
  arrange(Disease_name, Pathogen, Host)

combined_network_canonical_zoonotic <- combined_network_canonical %>%
  filter(!is.na(is_zoonotic), is_zoonotic)

# ------------------------------------------------------------------------------|
#      Write outputs and report ------------------------------------------------|
# ------------------------------------------------------------------------------|
write_csv(canonical_lookup, lookup_output_path, na = "")
write_csv(combined_network_canonical, canonical_output_path, na = "")
write_csv(combined_network_canonical_zoonotic, zoonotic_network_output_path, na = "")

cat("Wrote pathogen lookup to:\n")
cat(lookup_output_path, "\n\n")
cat("Wrote canonical network to:\n")
cat(canonical_output_path, "\n\n")
cat("Wrote zoonotic-only canonical network to:\n")
cat(zoonotic_network_output_path, "\n\n")
cat("Distinct raw pathogen labels:", nrow(network_targets), "\n")
cat("Distinct canonical pathogen labels:", n_distinct(combined_network_canonical$Pathogen), "\n")
cat(
  "Distinct canonical zoonotic pathogen labels:",
  combined_network_canonical_zoonotic %>% pull(Pathogen) %>% n_distinct(),
  "\n"
)

unmatched_kept <- canonical_lookup %>%
  filter(canonicalization_status == "raw_retained_no_match")

if (nrow(unmatched_kept) > 0) {
  cat("\nPathogens retained without a WHO canonical match:\n")
  print(unmatched_kept$Pathogen_raw)
}
