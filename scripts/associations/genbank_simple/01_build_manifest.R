# ------------------------------------------------------------------------------|
#      01_build_manifest.R -----------------------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Build the strict GenBank-simple manifest from current WHO/network
#          zoonotic targets with Gibb et al. or EMPRES-i point-data support.
# Inputs : who_pathogens_diseases_zoonotic.csv
#          combined_who_network_canonical_zoonotic.csv
# Outputs: genbank_simple_manifest.csv
#          excluded_targets.csv
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))
source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Define input and output paths ------------------------------------------
# ------------------------------------------------------------------------------|
who_path <- who_zoonotic_pathogens_path()
network_path <- who_canonical_zoonotic_network_path()
output_dir <- genbank_simple_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Load WHO point-supported targets ---------------------------------------
# ------------------------------------------------------------------------------|
# Limit the old 19-target manifest to zoonotic WHO rows with Gibb et al. or
# EMPRES-i point-data support.
who_pathogens <- read_csv(who_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    Pathogens = clean_text(Pathogens),
    Disease_name = clean_text(Disease_name),
    in_gibb_etal = as_logical_flag(in_gibb_etal),
    in_empres_i = as_logical_flag(in_empres_i)
  ) %>%
  filter(in_gibb_etal | in_empres_i) %>%
  distinct(
    Family,
    `PHEIC risk`,
    Pathogens,
    previous_name,
    msl39_viral_name,
    Disease_name,
    in_gibb_etal,
    in_empres_i,
    .keep_all = TRUE
  )

# ------------------------------------------------------------------------------|
#      Load network taxonomy --------------------------------------------------
# ------------------------------------------------------------------------------|
network_targets <- read_csv(network_path, show_col_types = FALSE, na = c("", "NA")) %>%
  transmute(
    Pathogens = clean_text(Pathogen),
    Disease_name = clean_text(Disease_name),
    PathogenTaxID = clean_text(PathogenTaxID),
    network_pathogen_type = clean_text(PathogenType),
    network_zoonotic_status = clean_text(zoonotic_status),
    network_canonicalization_status = clean_text(canonicalization_status)
  ) %>%
  distinct()

# ------------------------------------------------------------------------------|
#      Join targets and apply scope guardrails --------------------------------
# ------------------------------------------------------------------------------|
joined_targets <- who_pathogens %>%
  inner_join(
    network_targets,
    by = c("Pathogens", "Disease_name"),
    relationship = "many-to-many"
  )

target_summary <- joined_targets %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    Family = collapse_unique(Family),
    `PHEIC risk` = collapse_unique(`PHEIC risk`),
    previous_name = collapse_unique(previous_name),
    msl39_viral_name = collapse_unique(msl39_viral_name),
    in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
    in_empres_i = any(in_empres_i, na.rm = TRUE),
    PathogenTaxID = collapse_unique(PathogenTaxID),
    network_pathogen_type = collapse_unique(network_pathogen_type),
    network_zoonotic_status = collapse_unique(network_zoonotic_status),
    network_canonicalization_status = collapse_unique(network_canonicalization_status),
    .groups = "drop"
  ) %>%
  mutate(
    target_id = make_target_id(Pathogens, Disease_name),
    excluded_coronavirus = is_coronavirus_excluded(Pathogens, Disease_name),
    excluded_influenza = is_broad_or_unwanted_influenza(Pathogens),
    exclusion_reason = case_when(
      excluded_coronavirus ~ "coronavirus_scope_deferred",
      excluded_influenza ~ "broad_or_unwanted_influenza",
      TRUE ~ NA_character_
    )
  )

excluded_targets <- target_summary %>%
  filter(!is.na(exclusion_reason)) %>%
  select(
    target_id,
    Pathogens,
    Disease_name,
    PathogenTaxID,
    in_gibb_etal,
    in_empres_i,
    exclusion_reason
  ) %>%
  arrange(exclusion_reason, Pathogens, Disease_name)

# ------------------------------------------------------------------------------|
#      Build strict GenBank-simple manifest -----------------------------------
# ------------------------------------------------------------------------------|
manifest <- target_summary %>%
  filter(is.na(exclusion_reason)) %>%
  rowwise() %>%
  mutate(
    simple_query = build_simple_query(
      pathogen = Pathogens,
      tax_ids = unlist(strsplit(dplyr::coalesce(PathogenTaxID, ""), ";\\s*"))
    ),
    query_used = simple_query,
    source_db = "nuccore",
    query_strategy = case_when(
      is_allowed_influenza_target(Pathogens) ~ "influenza_subtype_constrained_full_retrieval",
      !is.na(PathogenTaxID) ~ "taxid_full_retrieval",
      TRUE ~ "organism_name_full_retrieval"
    ),
    query_source = case_when(
      is_allowed_influenza_target(Pathogens) ~ "simple_subtype_guardrail",
      TRUE ~ "simple_generated_query"
    ),
    old_query_profile = NA_character_,
    old_taxid_query = NA_character_,
    old_organism_query = NA_character_,
    retrieval_policy = "full_deterministic_pagination",
    allow_sampling = FALSE
  ) %>%
  ungroup() %>%
  select(
    target_id,
    Pathogens,
    Disease_name,
    Family,
    `PHEIC risk`,
    previous_name,
    msl39_viral_name,
    in_gibb_etal,
    in_empres_i,
    PathogenTaxID,
    network_pathogen_type,
    network_zoonotic_status,
    network_canonicalization_status,
    source_db,
    query_strategy,
    query_source,
    old_query_profile,
    old_taxid_query,
    old_organism_query,
    simple_query,
    query_used,
    retrieval_policy,
    allow_sampling
  ) %>%
  arrange(Pathogens, Disease_name)

# ------------------------------------------------------------------------------|
#      Write outputs -----------------------------------------------------------
# ------------------------------------------------------------------------------|
write_csv(
  manifest,
  genbank_simple_file_path(output_dir, "genbank_simple_manifest.csv", create_parent = TRUE)
)
write_csv(
  excluded_targets,
  genbank_simple_file_path(output_dir, "excluded_targets.csv", create_parent = TRUE)
)

message("Wrote manifest rows: ", nrow(manifest))
message("Wrote excluded rows: ", nrow(excluded_targets))
