# ------------------------------------------------------------------------------
# 1_2e_Master_Host_Species.R
# ------------------------------------------------------------------------------
# Purpose: Retrieve host species associations for disease master list query units
#          using source-prioritized VIRION/CLOVER mappings from
#          master_pathogen_host_query_units.csv.
#
# Inputs : who_diseases_host_query_path(
#            "master_pathogen_host_query_units.csv"
#          )
#          local VIRION and CLOVER source tables
#
# Outputs: who_master_pathogen_host_species_path()
#          who_master_pathogen_host_species_summary_path()
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_host_network_helpers.R"
))

host_query_path <- who_diseases_host_query_path(
  "master_pathogen_host_query_units.csv"
)

host_output_path <- who_master_pathogen_host_species_path()
summary_output_path <- who_master_pathogen_host_species_summary_path()

clover_dir <- file.path(
  clover_source_dir,
  "clover", "clover_1.0_allpathogens"
)

clover_paths <- file.path(
  clover_dir,
  c(
    "CLOVER_1.0_Bacteria_AssociationsFlatFile.csv",
    "CLOVER_1.0_Viruses_AssociationsFlatFile.csv",
    "CLOVER_1.0_HelminthProtozoaFungi_AssociationsFlatFile.csv"
  )
)

required_inputs <- c(host_query_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input files: ", paste(missing_inputs, collapse = "; "))
}

missing_clover <- clover_paths[!file.exists(clover_paths)]
if (length(missing_clover) > 0) {
  stop("Missing CLOVER input files: ", paste(missing_clover, collapse = "; "))
}

host_queries <- read_csv(host_query_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), host_network_clean_text))

required_query_cols <- c(
  "analysis_unit_id",
  "master_row",
  "disease_master_name",
  "resolved_disease_name",
  "resolved_pathogen_name",
  "resolved_pathogen_rank",
  "host_query_include_default",
  "host_query_bucket",
  "host_query_source",
  "host_query_pathogen_names",
  "host_query_taxids",
  "match_review_flag",
  "shared_species_proxy_flag",
  "match_review_notes"
)

missing_query_cols <- setdiff(required_query_cols, names(host_queries))
if (length(missing_query_cols) > 0) {
  stop(
    "master_pathogen_host_query_units.csv missing required columns: ",
    paste(missing_query_cols, collapse = ", ")
  )
}

active_queries <- host_queries %>%
  filter(!is.na(host_query_source), host_query_source %in% c("virion", "clover")) %>%
  mutate(
    host_query_include_default = coalesce(host_query_include_default, FALSE),
    match_review_flag = coalesce(match_review_flag, FALSE),
    shared_species_proxy_flag = coalesce(shared_species_proxy_flag, FALSE),
    query_pathogen_keys = map(
      host_query_pathogen_names,
      ~ host_network_clean_key(host_network_split_semicolon_values(.x))
    ),
    query_taxids = map(host_query_taxids, host_network_split_semicolon_values),
    query_has_taxids = map_int(query_taxids, length) > 0,
    query_has_names = map_int(query_pathogen_keys, length) > 0
  )

if (!exists("virion_data")) {
  source(file.path("scripts", "associations", "network_building", "helpers", "virion_loaders.R"))
  virion_data <- load_virion_data()
}

virion_links <- virion_data$virion %>%
  transmute(
    source = "virion",
    source_pathogen_name = host_network_clean_text(Virus),
    source_pathogen_taxid = host_network_clean_text(VirusTaxID),
    source_pathogen_type = "virus",
    source_pathogen_family = host_network_clean_text(VirusFamily),
    source_pathogen_order = host_network_clean_text(VirusOrder),
    source_pathogen_class = host_network_clean_text(VirusClass),
    source_host_name = host_network_clean_text(Host),
    source_host_taxid = host_network_clean_text(HostTaxID),
    source_host_genus = host_network_clean_text(HostGenus),
    source_host_family = host_network_clean_text(HostFamily),
    source_host_order = host_network_clean_text(HostOrder),
    source_host_class = host_network_clean_text(HostClass),
    source_database = host_network_clean_text(Database),
    source_assoc_id = host_network_clean_text(AssocID),
    source_detection_method = host_network_clean_text(DetectionMethod),
    source_host_flag_id = HostFlagID
  ) %>%
  mutate(
    source_pathogen_key = host_network_clean_key(source_pathogen_name),
    source_pathogen_taxid = str_remove(source_pathogen_taxid, "\\.0+$")
  )

clover_links <- map_dfr(clover_paths, ~ read_csv(.x, show_col_types = FALSE, na = c("", "NA"))) %>%
  transmute(
    source = "clover",
    source_pathogen_name = host_network_clean_text(Pathogen),
    source_pathogen_taxid = host_network_clean_text(PathogenTaxID),
    source_pathogen_type = host_network_clean_text(PathogenType),
    source_pathogen_family = host_network_clean_text(PathogenFamily),
    source_pathogen_order = host_network_clean_text(PathogenOrder),
    source_pathogen_class = host_network_clean_text(PathogenClass),
    source_host_name = host_network_clean_text(Host),
    source_host_taxid = host_network_clean_text(HostTaxID),
    source_host_genus = host_network_clean_text(HostGenus),
    source_host_family = host_network_clean_text(HostFamily),
    source_host_order = host_network_clean_text(HostOrder),
    source_host_class = host_network_clean_text(HostClass),
    source_database = host_network_clean_text(Database),
    source_assoc_id = host_network_clean_text(AssocID),
    source_detection_method = host_network_clean_text(DetectionMethod),
    source_host_flag_id = NA
  ) %>%
  mutate(
    source_pathogen_key = host_network_clean_key(source_pathogen_name),
    source_pathogen_taxid = str_remove(source_pathogen_taxid, "\\.0+$")
  )

all_source_links <- bind_rows(virion_links, clover_links) %>%
  filter(
    !is.na(source_pathogen_name),
    !is.na(source_host_name)
  )

matched_rows <- map_dfr(
  seq_len(nrow(active_queries)),
  ~ host_network_match_one_query(active_queries[.x, , drop = FALSE], all_source_links)
)

matched_rows <- matched_rows %>%
  mutate(
    host_record_key = paste(
      analysis_unit_id,
      source,
      source_pathogen_taxid,
      source_pathogen_name,
      source_host_taxid,
      source_host_name,
      source_detection_method,
      sep = "|"
    )
  ) %>%
  distinct(host_record_key, .keep_all = TRUE) %>%
  select(
    analysis_unit_id,
    master_row,
    disease_master_name,
    resolved_disease_name,
    resolved_pathogen_name,
    resolved_pathogen_rank,
    preferred_match_source,
    host_query_include_default,
    host_query_bucket,
    host_query_source,
    host_query_pathogen_names,
    host_query_taxids,
    match_method,
    match_review_flag,
    shared_species_proxy_flag,
    match_review_notes,
    source_database,
    source_assoc_id,
    source_detection_method,
    source_host_flag_id,
    source_pathogen_name,
    source_pathogen_taxid,
    source_pathogen_type,
    source_pathogen_family,
    source_pathogen_order,
    source_pathogen_class,
    source_host_name,
    source_host_taxid,
    source_host_genus,
    source_host_family,
    source_host_order,
    source_host_class
  ) %>%
  arrange(master_row, source_pathogen_name, source_host_name)

default_output <- matched_rows %>%
  filter(host_query_bucket == "default_clean", host_query_include_default)

review_output <- matched_rows %>%
  filter(host_query_bucket != "default_clean" | !host_query_include_default)

host_species_output <- matched_rows %>%
  transmute(
    Pathogen = source_pathogen_name,
    PathogenTaxID = source_pathogen_taxid,
    Disease_name = resolved_disease_name,
    HostTaxID = source_host_taxid,
    Host = source_host_name,
    PathogenClass = source_pathogen_class,
    PathogenOrder = source_pathogen_order,
    PathogenFamily = source_pathogen_family,
    HostClass = source_host_class,
    HostFamily = source_host_family,
    HostOrder = source_host_order,
    DetectionMethod = source_detection_method,
    MainSource = str_to_upper(host_query_source),
    PathogenType = source_pathogen_type,
    analysis_unit_id,
    master_row,
    disease_master_name,
    resolved_pathogen_name,
    host_query_bucket,
    host_query_include_default,
    match_review_flag,
    shared_species_proxy_flag,
    match_review_notes,
    match_method,
    source_database,
    source_assoc_id,
    source_host_flag_id,
    host_query_pathogen_names,
    host_query_taxids
  ) %>%
  arrange(master_row, Pathogen, Host)

default_query_keys <- active_queries %>%
  filter(host_query_bucket == "default_clean", host_query_include_default) %>%
  distinct(analysis_unit_id, disease_master_name, resolved_disease_name, host_query_source)

default_matches_by_unit <- default_output %>%
  distinct(analysis_unit_id) %>%
  mutate(has_host_matches = TRUE)

default_zero_match_diseases <- default_query_keys %>%
  left_join(default_matches_by_unit, by = "analysis_unit_id") %>%
  mutate(has_host_matches = coalesce(has_host_matches, FALSE)) %>%
  filter(!has_host_matches) %>%
  select(analysis_unit_id, disease_master_name, resolved_disease_name, host_query_source)

summary_counts <- bind_rows(
  default_output %>%
    summarise(
      table_name = "default_host_species",
      rows = n(),
      distinct_analysis_units = n_distinct(analysis_unit_id),
      distinct_diseases = n_distinct(disease_master_name),
      distinct_pathogens = n_distinct(source_pathogen_name),
      distinct_hosts = n_distinct(source_host_name),
      host_taxid_missing_rows = sum(is.na(source_host_taxid)),
      pathogen_taxid_missing_rows = sum(is.na(source_pathogen_taxid))
    ),
  review_output %>%
    summarise(
      table_name = "review_host_species",
      rows = n(),
      distinct_analysis_units = n_distinct(analysis_unit_id),
      distinct_diseases = n_distinct(disease_master_name),
      distinct_pathogens = n_distinct(source_pathogen_name),
      distinct_hosts = n_distinct(source_host_name),
      host_taxid_missing_rows = sum(is.na(source_host_taxid)),
      pathogen_taxid_missing_rows = sum(is.na(source_pathogen_taxid))
    )
)

summary_by_bucket_source <- matched_rows %>%
  count(host_query_bucket, source = host_query_source, match_method, name = "rows") %>%
  mutate(table_name = "bucket_source_match_method")

summary_by_host_rank <- matched_rows %>%
  count(host_query_bucket, source = host_query_source, source_host_class, name = "rows") %>%
  mutate(table_name = "bucket_source_host_class")

summary_zero_match <- default_zero_match_diseases %>%
  mutate(
    table_name = "default_clean_zero_matches",
    rows = 0L,
    distinct_analysis_units = 1L,
    distinct_diseases = 1L,
    distinct_pathogens = NA_integer_,
    distinct_hosts = NA_integer_,
    host_taxid_missing_rows = NA_integer_,
    pathogen_taxid_missing_rows = NA_integer_,
    host_query_bucket = "default_clean",
    source = host_query_source,
    match_method = NA_character_,
    source_host_class = NA_character_
  ) %>%
  select(
    table_name,
    rows,
    distinct_analysis_units,
    distinct_diseases,
    distinct_pathogens,
    distinct_hosts,
    host_taxid_missing_rows,
    pathogen_taxid_missing_rows,
    host_query_bucket,
    source,
    match_method,
    source_host_class,
    analysis_unit_id,
    disease_master_name,
    resolved_disease_name
  )

summary_output <- summary_counts %>%
  mutate(
    host_query_bucket = NA_character_,
    source = NA_character_,
    match_method = NA_character_,
    source_host_class = NA_character_,
    analysis_unit_id = NA_character_,
    disease_master_name = NA_character_,
    resolved_disease_name = NA_character_
  ) %>%
  bind_rows(
    summary_by_bucket_source %>%
      mutate(
        distinct_analysis_units = NA_integer_,
        distinct_diseases = NA_integer_,
        distinct_pathogens = NA_integer_,
        distinct_hosts = NA_integer_,
        host_taxid_missing_rows = NA_integer_,
        pathogen_taxid_missing_rows = NA_integer_,
        source_host_class = NA_character_,
        analysis_unit_id = NA_character_,
        disease_master_name = NA_character_,
        resolved_disease_name = NA_character_
      ) %>%
      select(names(summary_zero_match)),
    summary_by_host_rank %>%
      mutate(
        distinct_analysis_units = NA_integer_,
        distinct_diseases = NA_integer_,
        distinct_pathogens = NA_integer_,
        distinct_hosts = NA_integer_,
        host_taxid_missing_rows = NA_integer_,
        pathogen_taxid_missing_rows = NA_integer_,
        match_method = NA_character_,
        analysis_unit_id = NA_character_,
        disease_master_name = NA_character_,
        resolved_disease_name = NA_character_
      ) %>%
      select(names(summary_zero_match)),
    summary_zero_match
  ) %>%
  select(names(summary_zero_match))

write_csv(host_species_output, host_output_path, na = "")
write_csv(summary_output, summary_output_path, na = "")

cat("Active host-query rows:", nrow(active_queries), "\n")
cat("Default-clean query rows:", nrow(default_query_keys), "\n")
cat("Host-species rows written:", nrow(host_species_output), "\n")
cat("Default host-species rows:", nrow(default_output), "\n")
cat("Review host-species rows retained:", nrow(review_output), "\n")
cat("Distinct default-clean diseases with host matches:", n_distinct(default_output$disease_master_name), "\n")
cat("Default-clean diseases with zero matches:", nrow(default_zero_match_diseases), "\n")
cat("Wrote:", host_output_path, "\n")
cat("Wrote:", summary_output_path, "\n")
