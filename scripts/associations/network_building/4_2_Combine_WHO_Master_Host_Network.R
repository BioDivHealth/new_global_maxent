# ------------------------------------------------------------------------------
# 4_2_Combine_WHO_Master_Host_Network.R
# ------------------------------------------------------------------------------
# Purpose: Append the WHO host network and disease-master host network into one
#          downstream-ready host table without dropping review evidence rows.
#
# Inputs : pathogen_association_data/WHO/networks/combined_who_network.csv
#          pathogen_association_data/WHO/who_diseases/
#            master_pathogen_host_species_clean.csv
#            master_plus_who_analysis_units.csv
#
# Output: pathogen_association_data/WHO/networks/
#           master_plus_who_host_network.csv
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

source(here("scripts", "associations", "working_inputs.R"))

network_dir <- here("pathogen_association_data", "WHO", "networks")
who_network_path <- who_raw_network_path()
master_host_path <- who_master_pathogen_host_species_clean_path()
analysis_units_path <- who_master_plus_analysis_units_path()
who_keep_path <- who_pathogen_analysis_units_keep_path()
combined_output_path <- file.path(network_dir, "master_plus_who_host_network.csv")

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

common_columns <- c(
  "Pathogen",
  "PathogenTaxID",
  "PHEIC risk",
  "Disease_name",
  "HostTaxID",
  "Host",
  "PathogenClass",
  "PathogenOrder",
  "PathogenFamily",
  "PathogenGenus",
  "HostPhylum",
  "HostClass",
  "HostFamily",
  "HostOrder",
  "DetectionMethod",
  "high_quality_detection",
  "downstream_default_include",
  "downstream_review_reason",
  "MainSource",
  "PathogenType",
  "in_gibb_etal",
  "in_empres_i",
  "modelling_scope_status",
  "modelling_scope_reason",
  "Host_raw",
  "host_name_cleaning_method",
  "source_database",
  "source_assoc_id",
  "source_host_flag_id",
  "host_taxonomy_ready",
  "host_taxonomy_flag",
  "is_human_host",
  "is_model_or_lab_host",
  "is_domestic_or_livestock_hint"
)

provenance_columns <- c(
  "host_network_source",
  "source_table",
  "possible_cross_source_duplicate_flag"
)

required_common <- c(
  "Pathogen",
  "PathogenTaxID",
  "Disease_name",
  "HostTaxID",
  "Host",
  "DetectionMethod",
  "MainSource",
  "high_quality_detection",
  "downstream_default_include"
)

master_audit_only_columns <- c(
  "host_query_include_default",
  "match_method",
  "host_query_pathogen_names",
  "host_query_taxids",
  "all_method_host_count_for_analysis_unit"
)

broad_source_pathogens <- c(
  "Genus Vesiculovirus",
  "Subgenus Merbecovirus",
  "Subgenus Sarbecovirus"
)

required_paths <- c(who_network_path, master_host_path, analysis_units_path, who_keep_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = "; "))
}

add_missing_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  for (col in missing) {
    data[[col]] <- NA
  }
  data
}

scope_key <- function(x) {
  x %>%
    clean_text() %>%
    str_to_lower() %>%
    str_replace_all("&", " and ") %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
}

read_network <- function(path, host_network_source, source_table) {
  data <- read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(
      across(where(is.character), clean_text),
      PathogenTaxID = clean_text(PathogenTaxID),
      HostTaxID = clean_text(HostTaxID),
      high_quality_detection = coalesce(high_quality_detection, FALSE),
      downstream_default_include = coalesce(downstream_default_include, FALSE),
      host_network_source = host_network_source,
      source_table = source_table
    )

  missing_required <- setdiff(required_common, names(data))
  if (length(missing_required) > 0) {
    stop(
      source_table,
      " missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  data <- add_missing_columns(
    data,
    c(
      "Host_raw",
      "host_name_cleaning_method",
      "source_database",
      "source_assoc_id",
      "source_host_flag_id",
      "host_taxonomy_ready",
      "host_taxonomy_flag",
      "is_human_host",
      "is_model_or_lab_host",
      "is_domestic_or_livestock_hint"
    )
  )

  data %>%
    mutate(
      Host_raw = coalesce(Host_raw, Host),
      host_name_cleaning_method = coalesce(host_name_cleaning_method, "existing_who_network_host"),
      source_database = coalesce(source_database, MainSource),
      host_taxonomy_flag = case_when(
        !is.na(host_taxonomy_flag) ~ host_taxonomy_flag,
        is.na(Host) ~ "missing_name",
        is.na(HostTaxID) ~ "missing_taxid",
        str_detect(str_to_lower(Host), "\\b(sp|spp|species|unidentified|unknown|uncultured)\\b\\.?") ~ "unresolved_sp",
        str_detect(str_to_lower(Host), "^[a-z][a-z-]+\\s+[a-z][a-z.-]+(\\s+[a-z][a-z.-]+)?$") ~ "species_like",
        TRUE ~ "unresolved_sp"
      ),
      host_taxonomy_ready = coalesce(
        host_taxonomy_ready,
        host_taxonomy_flag == "species_like" & !is.na(HostTaxID)
      ),
      is_human_host = coalesce(
        is_human_host,
        str_to_lower(Host) == "homo sapiens" | HostTaxID == "9606"
      ),
      is_model_or_lab_host = coalesce(
        is_model_or_lab_host,
        str_detect(
          str_to_lower(Host),
          paste(
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
        )
      ),
      is_domestic_or_livestock_hint = coalesce(
        is_domestic_or_livestock_hint,
        str_detect(
          str_to_lower(Host),
          paste(
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
        )
      )
    )
}

analysis_units <- read_csv(analysis_units_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

who_keep_units <- read_csv(who_keep_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    modelling_scope_status = case_when(
      source_pathogen %in% broad_source_pathogens ~ "defer_broad_or_aggregate_unit",
      TRUE ~ "include"
    ),
    modelling_scope_reason = case_when(
      source_pathogen %in% broad_source_pathogens ~ paste0(
        "Source pathogen ", source_pathogen,
        " is a broad genus/subgenus scope that has been deferred in later role/modelling work; retain WHO rows for audit only."
      ),
      TRUE ~ "Included because this WHO-only host-network row is present in who_pathogen_analysis_units_keep.csv and has no stricter manual master-plus scope."
    )
  )

scope_by_analysis_unit_id <- analysis_units %>%
  filter(!is.na(analysis_unit_id)) %>%
  transmute(
    analysis_unit_id,
    modelling_scope_status,
    modelling_scope_reason
  ) %>%
  distinct(analysis_unit_id, .keep_all = TRUE)

make_scope_aliases <- function(data, source_priority) {
  alias_cols <- intersect(
    c("analysis_unit", "analysis_unit_label", "source_pathogen", "source_previous_name", "source_msl39_viral_name"),
    names(data)
  )

  data %>%
    select(
      any_of(c("source_disease_name", "modelling_scope_status", "modelling_scope_reason")),
      all_of(alias_cols)
    ) %>%
    pivot_longer(
      cols = all_of(alias_cols),
      names_to = "pathogen_alias_source",
      values_to = "pathogen_alias"
    ) %>%
    transmute(
      scope_priority = source_priority,
      disease_key = scope_key(source_disease_name),
      pathogen_key = scope_key(pathogen_alias),
      modelling_scope_status,
      modelling_scope_reason
    )
}

scope_by_who_key <- bind_rows(
  make_scope_aliases(analysis_units, 1L),
  make_scope_aliases(who_keep_units, 2L)
) %>%
  filter(!is.na(disease_key), !is.na(pathogen_key)) %>%
  arrange(scope_priority) %>%
  distinct(disease_key, pathogen_key, .keep_all = TRUE) %>%
  select(-scope_priority)

scope_by_who_disease <- bind_rows(
  analysis_units %>%
    transmute(
      scope_priority = 1L,
      disease_key = scope_key(source_disease_name),
      modelling_scope_status,
      modelling_scope_reason
    ),
  who_keep_units %>%
    transmute(
      scope_priority = 2L,
      disease_key = scope_key(source_disease_name),
      modelling_scope_status,
      modelling_scope_reason
    )
) %>%
  filter(!is.na(disease_key), !is.na(modelling_scope_status)) %>%
  group_by(disease_key) %>%
  arrange(scope_priority, .by_group = TRUE) %>%
  summarise(
    modelling_scope_status = if_else(
      n_distinct(modelling_scope_status) == 1,
      first(modelling_scope_status),
      "review_before_modelling"
    ),
    modelling_scope_reason = if_else(
      n_distinct(modelling_scope_status) == 1,
      first(modelling_scope_reason),
      "Disease has multiple analysis-unit scope statuses; review before using as a default modelling row."
    ),
    .groups = "drop"
  )

who_network <- read_network(
  who_network_path,
  host_network_source = "who",
  source_table = "combined_who_network.csv"
) %>%
  mutate(
    disease_key = scope_key(Disease_name),
    pathogen_key = scope_key(Pathogen)
  ) %>%
  left_join(scope_by_who_key, by = c("disease_key", "pathogen_key")) %>%
  left_join(
    scope_by_who_disease,
    by = "disease_key",
    suffix = c("", "_disease")
  ) %>%
  mutate(
    modelling_scope_status = coalesce(
      modelling_scope_status,
      modelling_scope_status_disease,
      "include"
    ),
    modelling_scope_reason = coalesce(
      modelling_scope_reason,
      modelling_scope_reason_disease,
      "Included as a legacy WHO host-network row with no explicit master-plus disease scope match."
    )
  ) %>%
  select(
    -disease_key,
    -pathogen_key,
    -ends_with("_disease")
  )

master_network <- read_network(
  master_host_path,
  host_network_source = "master",
  source_table = "master_pathogen_host_species_clean.csv"
) %>%
  select(-any_of(c("modelling_scope_status", "modelling_scope_reason"))) %>%
  left_join(scope_by_analysis_unit_id, by = "analysis_unit_id")

all_columns <- unique(c(
  common_columns,
  provenance_columns,
  setdiff(names(master_network), c(common_columns, provenance_columns, master_audit_only_columns)),
  setdiff(names(who_network), c(common_columns, provenance_columns, master_audit_only_columns))
))

combined_network <- bind_rows(
  add_missing_columns(who_network, all_columns),
  add_missing_columns(master_network, all_columns)
) %>%
  group_by(Disease_name, PathogenTaxID, HostTaxID, DetectionMethod, MainSource) %>%
  mutate(
    possible_cross_source_duplicate_flag = n_distinct(host_network_source) > 1
  ) %>%
  ungroup() %>%
  select(all_of(all_columns))

stopifnot(nrow(combined_network) == nrow(who_network) + nrow(master_network))
stopifnot(!any(is.na(combined_network$downstream_default_include)))
stopifnot(!any(is.na(combined_network$Host[combined_network$downstream_default_include])))
stopifnot(!any(is.na(combined_network$Host_raw)))
stopifnot(!any(is.na(combined_network$host_taxonomy_ready)))
stopifnot(!any(is.na(combined_network$modelling_scope_status)))
stopifnot(!any(is.na(combined_network$modelling_scope_reason)))

write_csv(combined_network, combined_output_path, na = "")

cat("WHO network rows:", nrow(who_network), "\n")
cat("Master host rows:", nrow(master_network), "\n")
cat("Combined host-network rows:", nrow(combined_network), "\n")
cat("Downstream default include rows:", sum(combined_network$downstream_default_include), "\n")
cat("Possible cross-source duplicate rows:", sum(combined_network$possible_cross_source_duplicate_flag), "\n")
cat("Detection methods:\n")
print(count(combined_network, host_network_source, DetectionMethod, downstream_default_include), n = Inf)
cat("Rows by source:\n")
print(count(combined_network, host_network_source, source_table), n = Inf)
cat("Rows by modelling scope:\n")
print(
  count(
    combined_network,
    host_network_source,
    modelling_scope_status,
    downstream_default_include
  ),
  n = Inf
)
cat("Wrote:", combined_output_path, "\n")
