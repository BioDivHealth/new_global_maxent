# ------------------------------------------------------------------------------
# 1_2f_Master_Host_Species_QA_Clean.R
# ------------------------------------------------------------------------------
# Purpose: Add conservative QA, host-name harmonization, and downstream-readiness
#          flags to the disease master host-species table without dropping
#          evidence rows.
#
# Input : who_master_pathogen_host_species_path()
# Output: who_master_pathogen_host_species_clean_path()
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

large_host_threshold <- 75L
narrow_host_threshold <- 3L
host_detection_methods_keep <- c("Isolation/Observation", "PCR/Sequencing")

host_input_path <- who_master_pathogen_host_species_path()
host_clean_output_path <- who_master_pathogen_host_species_clean_path()
analysis_units_path <- who_master_plus_analysis_units_path()

virion_host_standardized_path <- file.path(
  who_virion_dir,
  "who_host_species_standardized.csv"
)
clover_host_standardized_path <- file.path(
  who_clover_dir,
  "clover_host_species_standardized.csv"
)

required_input_cols <- c(
  "Pathogen",
  "PathogenTaxID",
  "Disease_name",
  "HostTaxID",
  "Host",
  "DetectionMethod",
  "MainSource",
  "analysis_unit_id",
  "disease_master_name",
  "resolved_pathogen_name",
  "host_query_bucket",
  "host_query_include_default",
  "match_review_flag",
  "shared_species_proxy_flag"
)

required_paths <- c(
  host_input_path,
  virion_host_standardized_path,
  clover_host_standardized_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = "; "))
}

host_rows_raw <- read_csv(host_input_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(
    across(where(is.character), host_network_clean_text),
    HostTaxID = host_network_clean_text(HostTaxID),
    PathogenTaxID = host_network_clean_text(PathogenTaxID),
    host_query_include_default = coalesce(host_query_include_default, FALSE),
    match_review_flag = coalesce(match_review_flag, FALSE),
    shared_species_proxy_flag = coalesce(shared_species_proxy_flag, FALSE)
  )

missing_cols <- setdiff(required_input_cols, names(host_rows_raw))
if (length(missing_cols) > 0) {
  stop(
    "master_pathogen_host_species.csv missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

analysis_unit_metadata <- tibble(
  analysis_unit_id = character(),
  `PHEIC risk` = character(),
  in_gibb_etal = logical(),
  in_empres_i = logical()
)

if (file.exists(analysis_units_path)) {
  analysis_unit_metadata <- read_csv(analysis_units_path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), host_network_clean_text)) %>%
    transmute(
      analysis_unit_id,
      `PHEIC risk` = host_network_clean_text(pheic_risk),
      in_gibb_etal = if ("in_gibb_etal" %in% names(.)) in_gibb_etal else NA,
      in_empres_i = if ("in_empres_i" %in% names(.)) in_empres_i else NA
    ) %>%
    distinct(analysis_unit_id, .keep_all = TRUE)
}

standard_hosts <- bind_rows(
  host_network_read_host_standardization(virion_host_standardized_path, "VIRION"),
  host_network_read_host_standardization(clover_host_standardized_path, "CLOVER")
)

host_lookup_specs <- tribble(
  ~lookup_name, ~group_cols, ~method_name, ~suffix,
  "source_taxid_host", list(c("host_source", "HostTaxID", "raw_host_key")), "source_taxid_host", "source_taxid_host",
  "source_taxid", list(c("host_source", "HostTaxID")), "source_taxid_unique", "source_taxid_unique",
  "source_name", list(c("host_source", "raw_host_key")), "source_name", "source_name",
  "cross_taxid", list(c("HostTaxID")), "cross_source_taxid_unique", "cross_source_taxid_unique"
)

host_lookups <- host_lookup_specs %>%
  mutate(
    lookup_data = pmap(
      list(group_cols, method_name, suffix),
      ~ host_network_summarise_host_lookup(standard_hosts, unlist(..1), ..2, ..3)
    )
  ) %>%
  select(lookup_name, lookup_data) %>%
  deframe()

host_rows <- host_rows_raw %>%
  mutate(
    Host_raw = Host,
    host_source = str_to_upper(MainSource),
    raw_host_key = host_network_clean_key(Host_raw)
  ) %>%
  left_join(
    host_lookups$source_taxid_host,
    by = c("host_source", "HostTaxID", "raw_host_key")
  ) %>%
  left_join(
    host_lookups$source_taxid,
    by = c("host_source", "HostTaxID")
  ) %>%
  left_join(
    host_lookups$source_name,
    by = c("host_source", "raw_host_key")
  ) %>%
  left_join(
    host_lookups$cross_taxid,
    by = "HostTaxID"
  ) %>%
  left_join(analysis_unit_metadata, by = "analysis_unit_id") %>%
  mutate(
    Host = coalesce(
      clean_host_source_taxid_host,
      clean_host_source_taxid_unique,
      clean_host_source_name,
      clean_host_cross_source_taxid_unique,
      Host_raw
    ),
    host_name_cleaning_method = case_when(
      !is.na(clean_host_source_taxid_host) ~ method_source_taxid_host,
      !is.na(clean_host_source_taxid_unique) ~ method_source_taxid_unique,
      !is.na(clean_host_source_name) ~ method_source_name,
      !is.na(clean_host_cross_source_taxid_unique) ~ method_cross_source_taxid_unique,
      TRUE ~ "source_name_fallback"
    ),
    HostPhylum = coalesce(
      HostPhylum_source_taxid_host,
      HostPhylum_source_taxid_unique,
      HostPhylum_source_name,
      HostPhylum_cross_source_taxid_unique
    ),
    HostClass = coalesce(
      HostClass_source_taxid_host,
      HostClass_source_taxid_unique,
      HostClass_source_name,
      HostClass_cross_source_taxid_unique,
      HostClass
    ),
    HostFamily = coalesce(
      HostFamily_source_taxid_host,
      HostFamily_source_taxid_unique,
      HostFamily_source_name,
      HostFamily_cross_source_taxid_unique,
      HostFamily
    ),
    HostOrder = coalesce(
      HostOrder_source_taxid_host,
      HostOrder_source_taxid_unique,
      HostOrder_source_name,
      HostOrder_cross_source_taxid_unique,
      HostOrder
    ),
    PathogenGenus = NA_character_
  ) %>%
  select(
    -host_source,
    -raw_host_key,
    -starts_with("clean_host_"),
    -starts_with("HostPhylum_"),
    -starts_with("HostClass_"),
    -starts_with("HostFamily_"),
    -starts_with("HostOrder_"),
    -starts_with("method_")
  )

default_host_counts <- host_rows %>%
  filter(
    host_query_include_default,
    host_query_bucket == "default_clean",
    DetectionMethod %in% host_detection_methods_keep
  ) %>%
  distinct(analysis_unit_id, HostTaxID, Host) %>%
  count(analysis_unit_id, name = "host_count_for_analysis_unit")

all_method_host_counts <- host_rows %>%
  filter(host_query_include_default, host_query_bucket == "default_clean") %>%
  distinct(analysis_unit_id, HostTaxID, Host) %>%
  count(analysis_unit_id, name = "all_method_host_count_for_analysis_unit")

host_rows_with_counts <- host_rows %>%
  left_join(default_host_counts, by = "analysis_unit_id") %>%
  left_join(all_method_host_counts, by = "analysis_unit_id") %>%
  mutate(
    host_count_for_analysis_unit = coalesce(host_count_for_analysis_unit, 0L),
    all_method_host_count_for_analysis_unit = coalesce(all_method_host_count_for_analysis_unit, 0L),
    high_quality_detection = DetectionMethod %in% host_detection_methods_keep,
    host_count_flag = case_when(
      host_query_bucket != "default_clean" | !host_query_include_default ~ "not_default_clean",
      host_count_for_analysis_unit >= large_host_threshold ~ "very_large",
      host_count_for_analysis_unit <= narrow_host_threshold ~ "very_narrow",
      TRUE ~ "expected_range"
    )
  )

model_or_lab_pattern <- paste(
  c(
    "^homo sapiens$",
    "^mus musculus$",
    "^rattus norvegicus$",
    "^rattus rattus$",
    "^cavia porcellus$",
    "^mesocricetus auratus$",
    "^oryctolagus cuniculus$",
    "^macaca\\b",
    "^chlorocebus\\b",
    "^callithrix\\b",
    "^gallus gallus$"
  ),
  collapse = "|"
)

domestic_or_livestock_pattern <- paste(
  c(
    "^bos taurus$",
    "^bos indicus$",
    "^bubalus bubalis$",
    "^ovis aries$",
    "^capra hircus$",
    "^sus scrofa$",
    "^equus caballus$",
    "^equus asinus$",
    "^camelus\\b",
    "^lama glama$",
    "^alpaca$",
    "^vicugna pacos$",
    "^gallus gallus$",
    "^meleagris gallopavo$",
    "^anas platyrhynchos$",
    "^anas platyrhynchos domesticus$",
    "^canis lupus familiaris$",
    "^felis catus$"
  ),
  collapse = "|"
)

host_clean <- host_rows_with_counts %>%
  mutate(
    host_name_key = str_to_lower(Host),
    host_taxonomy_flag = case_when(
      is.na(Host) ~ "missing_name",
      is.na(HostTaxID) ~ "missing_taxid",
      str_detect(host_name_key, "\\b(sp|spp|species|unidentified|unknown|uncultured)\\b\\.?") ~ "unresolved_sp",
      str_detect(host_name_key, "^[a-z][a-z-]+\\s+[a-z][a-z.-]+(\\s+[a-z][a-z.-]+)?$") ~ "species_like",
      TRUE ~ "unresolved_sp"
    ),
    host_taxonomy_ready = host_taxonomy_flag == "species_like" & !is.na(HostTaxID),
    is_human_host = host_name_key == "homo sapiens" | HostTaxID == "9606",
    is_model_or_lab_host = str_detect(host_name_key, model_or_lab_pattern),
    is_domestic_or_livestock_hint = str_detect(host_name_key, domestic_or_livestock_pattern),
    downstream_default_include = host_query_include_default &
      host_query_bucket == "default_clean" &
      !shared_species_proxy_flag &
      !match_review_flag &
      high_quality_detection &
      host_taxonomy_ready
  ) %>%
  rowwise() %>%
  mutate(
    downstream_review_reason = host_network_collapse_reasons(
      if_else(!host_query_include_default, "not_default_host_query", NA_character_),
      if_else(host_query_bucket != "default_clean", paste0("host_query_bucket=", host_query_bucket), NA_character_),
      if_else(shared_species_proxy_flag, "shared_species_proxy", NA_character_),
      if_else(match_review_flag, "match_review", NA_character_),
      if_else(!high_quality_detection, paste0("detection_method=", DetectionMethod), NA_character_),
      if_else(!host_taxonomy_ready, paste0("host_taxonomy_flag=", host_taxonomy_flag), NA_character_),
      if_else(host_count_flag == "very_large", "large_host_list_review", NA_character_),
      if_else(host_count_flag == "very_narrow", "narrow_host_list_review", NA_character_),
      if_else(is_human_host, "human_host_flag", NA_character_),
      if_else(is_model_or_lab_host, "model_or_lab_host_flag", NA_character_),
      if_else(is_domestic_or_livestock_hint, "domestic_or_livestock_hint", NA_character_)
    )
  ) %>%
  ungroup() %>%
  select(
    Pathogen,
    PathogenTaxID,
    `PHEIC risk`,
    Disease_name,
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
    high_quality_detection,
    downstream_default_include,
    downstream_review_reason,
    MainSource,
    PathogenType,
    in_gibb_etal,
    in_empres_i,
    Host_raw,
    host_name_cleaning_method,
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
    host_query_taxids,
    host_count_for_analysis_unit,
    all_method_host_count_for_analysis_unit,
    host_count_flag,
    host_taxonomy_ready,
    host_taxonomy_flag,
    is_human_host,
    is_model_or_lab_host,
    is_domestic_or_livestock_hint
  )

stopifnot(nrow(host_clean) == nrow(host_rows_raw))
stopifnot(!any(is.na(host_clean$downstream_default_include)))
stopifnot(!any(is.na(host_clean$Host_raw)))
stopifnot(!any(is.na(host_clean$Host[host_clean$downstream_default_include])))

write_csv(host_clean, host_clean_output_path, na = "")

default_clean_hosts <- host_clean %>%
  filter(host_query_include_default, host_query_bucket == "default_clean")

default_clean_host_counts <- default_clean_hosts %>%
  distinct(
    analysis_unit_id,
    disease_master_name,
    resolved_pathogen_name,
    host_count_for_analysis_unit,
    all_method_host_count_for_analysis_unit,
    host_count_flag
  )

default_query_units <- host_rows_raw %>%
  filter(host_query_include_default, host_query_bucket == "default_clean") %>%
  distinct(analysis_unit_id, disease_master_name, Disease_name)

default_units_with_hosts <- default_clean_hosts %>%
  distinct(analysis_unit_id)

default_zero_match <- default_query_units %>%
  anti_join(default_units_with_hosts, by = "analysis_unit_id")

cat("Host rows read:", nrow(host_rows_raw), "\n")
cat("Host rows written:", nrow(host_clean), "\n")
cat("Downstream default include rows:", sum(host_clean$downstream_default_include), "\n")
cat("All-method default-clean rows:", sum(host_clean$host_query_include_default & host_clean$host_query_bucket == "default_clean"), "\n")
cat("High-quality default-clean rows:", sum(host_clean$host_query_include_default & host_clean$host_query_bucket == "default_clean" & host_clean$high_quality_detection), "\n")
cat("Default-clean diseases with zero matches:", nrow(default_zero_match), "\n")

qa_tables <- list(
  "Detection methods:" = count(host_clean, high_quality_detection, DetectionMethod),
  "Host name cleaning methods:" = count(host_clean, host_name_cleaning_method),
  "Host count flags:" = count(host_clean, host_count_flag),
  "Host taxonomy flags:" = count(host_clean, host_taxonomy_flag),
  "Human/model/livestock flags:" = host_clean %>%
    summarise(
      human_rows = sum(is_human_host),
      model_or_lab_rows = sum(is_model_or_lab_host),
      domestic_or_livestock_hint_rows = sum(is_domestic_or_livestock_hint)
    ),
  "Large default-clean host lists:" = default_clean_host_counts %>%
    select(-all_method_host_count_for_analysis_unit) %>%
    filter(host_count_flag == "very_large") %>%
    arrange(desc(host_count_for_analysis_unit), disease_master_name),
  "Narrow default-clean host lists:" = default_clean_host_counts %>%
    filter(host_count_flag == "very_narrow") %>%
    arrange(host_count_for_analysis_unit, disease_master_name),
  "Narrow clean lists with broader all-method evidence:" = default_clean_host_counts %>%
    filter(
      host_count_flag == "very_narrow",
      all_method_host_count_for_analysis_unit > host_count_for_analysis_unit
    ) %>%
    arrange(
      desc(all_method_host_count_for_analysis_unit - host_count_for_analysis_unit),
      disease_master_name
    ),
  "Review/proxy buckets:" = host_clean %>%
    filter(host_query_bucket != "default_clean" | !host_query_include_default) %>%
    distinct(disease_master_name, resolved_pathogen_name, host_query_bucket, match_review_flag, shared_species_proxy_flag) %>%
    arrange(host_query_bucket, disease_master_name)
)

iwalk(qa_tables, ~ {
  cat(.y, "\n", sep = "")
  print(.x, n = Inf)
})
cat("Wrote:", host_clean_output_path, "\n")
