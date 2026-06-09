# ------------------------------------------------------------------------------|
#      05_standardize_countries.R ---------------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Add a conservative country cleanup layer to GenBank-simple records.
# Inputs : genbank_readiness_country_records.csv
# Outputs: genbank_readiness_country_records_standardized.csv
#          genbank_readiness_pathogen_country_summary_standardized.csv
#          genbank_readiness_disease_country_summary_standardized.csv
#          genbank_readiness_country_standardization_qa.csv
#
# Notes  : Raw provenance columns (`country_raw`, `geo_loc_name_raw`, `lat_lon`)
#          are left untouched. The existing `country` column is treated as a
#          first-pass parsed value. This script adds standardized companion
#          fields rather than overwriting prior columns.
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, purrr, readr, rnaturalearth, sf, stringr, tibble, tidyr)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))
source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Resolve run mode and input records -------------------------------------
# ------------------------------------------------------------------------------|
output_dir <- genbank_simple_dir
summary_kind <- Sys.getenv("GENBANK_SIMPLE_SUMMARY_KIND", unset = "readiness_combined") %>%
  clean_text() %>%
  stringr::str_to_lower()

summary_kind <- case_when(
  summary_kind %in% c("standard", "simple", "current") ~ "standard",
  summary_kind %in% c("readiness", "readiness_combined", "expanded_readiness") ~
    "readiness_combined",
  TRUE ~ NA_character_
)

if (is.na(summary_kind)) {
  stop(
    "GENBANK_SIMPLE_SUMMARY_KIND must be `standard` or `readiness_combined`.",
    call. = FALSE
  )
}

output_prefix <- if_else(summary_kind == "readiness_combined", "genbank_readiness_", "genbank_")
records_path <- genbank_simple_existing_file_path(output_dir, paste0(output_prefix, "country_records.csv"))

records <- read_csv(
  records_path,
  col_types = cols(.default = col_character()),
  na = c("", "NA")
)

# ------------------------------------------------------------------------------|
#      Latitude/longitude country lookup helpers ------------------------------
# ------------------------------------------------------------------------------|
parse_lat_lon <- function(lat_lon) {
  lat_lon <- clean_text(lat_lon)

  if (is.na(lat_lon)) {
    return(tibble(latitude = NA_real_, longitude = NA_real_))
  }

  parts <- stringr::str_match(
    lat_lon,
    stringr::regex(
      "^\\s*([0-9.]+)\\s*([NS])\\s+([0-9.]+)\\s*([EW])\\s*$",
      ignore_case = TRUE
    )
  )

  if (all(is.na(parts))) {
    return(tibble(latitude = NA_real_, longitude = NA_real_))
  }

  latitude <- suppressWarnings(as.numeric(parts[, 2]))
  longitude <- suppressWarnings(as.numeric(parts[, 4]))

  if (toupper(parts[, 3]) == "S") {
    latitude <- -latitude
  }

  if (toupper(parts[, 5]) == "W") {
    longitude <- -longitude
  }

  tibble(latitude = latitude, longitude = longitude)
}

lookup_country_from_lat_lon <- function(records_to_lookup) {
  if (nrow(records_to_lookup) == 0) {
    return(tibble(accession_key = character(), country_from_lat_lon = character()))
  }

  parsed_points <- records_to_lookup %>%
    select(accession_key, lat_lon) %>%
    mutate(parsed = purrr::map(lat_lon, parse_lat_lon)) %>%
    tidyr::unnest(parsed) %>%
    filter(!is.na(latitude), !is.na(longitude))

  if (nrow(parsed_points) == 0) {
    return(tibble(accession_key = character(), country_from_lat_lon = character()))
  }

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    select(country_from_lat_lon = name_long, geometry)

  point_sf <- sf::st_as_sf(
    parsed_points,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

  joined <- suppressWarnings(sf::st_join(point_sf, world, left = TRUE))

  joined %>%
    sf::st_drop_geometry() %>%
    transmute(
      accession_key = accession_key,
      country_from_lat_lon = clean_text(country_from_lat_lon)
    )
}

# ------------------------------------------------------------------------------|
#      Country standardization rules ------------------------------------------
# ------------------------------------------------------------------------------|
standardize_country_value <- function(country) {
  country <- clean_text(country)

  dplyr::case_when(
    is.na(country) ~ NA_character_,
    country %in% c("missing", "unknown", "Unknown", "not provided") ~ NA_character_,
    country == "United States of America" ~ "United States",
    country == "USA" ~ "United States",
    country == "U.S.A." ~ "United States",
    country == "UK" ~ "United Kingdom",
    country == "U.K." ~ "United Kingdom",
    country == "Viet Nam" ~ "Vietnam",
    country == "Russian Federation" ~ "Russia",
    country == "Czech Republic" ~ "Czechia",
    country == "Cote D'Ivoire" ~ "Cote d'Ivoire",
    country == "Côte d'Ivoire" ~ "Cote d'Ivoire",
    country == "Lao People's Democratic Republic" ~ "Laos",
    country == "LAO PEOPLE'S DEMOCRATIC REPUBLIC" ~ "Laos",
    country == "Tanzania, United Republic of" ~ "Tanzania",
    country == "Virgin Islands, U.S" ~ "U.S. Virgin Islands",
    country == "Zaire" ~ "Democratic Republic of the Congo",
    TRUE ~ country
  )
}

territory_values <- c(
  "American Samoa",
  "Anguilla",
  "Antarctica",
  "Aruba",
  "Bermuda",
  "British Virgin Islands",
  "Cayman Islands",
  "Cook Islands",
  "Curacao",
  "French Guiana",
  "French Polynesia",
  "Greenland",
  "Guadeloupe",
  "Guam",
  "Martinique",
  "Mayotte",
  "Montserrat",
  "New Caledonia",
  "Northern Mariana Islands",
  "Puerto Rico",
  "Reunion",
  "Saint Barthelemy",
  "Saint Martin",
  "Sint Maarten",
  "South Georgia and the South Sandwich Islands",
  "Tokelau",
  "U.S. Virgin Islands",
  "Wallis and Futuna"
)

ocean_values <- c("Atlantic Ocean", "Pacific Ocean", "Indian Ocean", "Southern Ocean", "Arctic Ocean")
historical_values <- c("USSR", "Yugoslavia", "Czechoslovakia", "Netherlands Antilles")
review_values <- c("Borneo")

# ------------------------------------------------------------------------------|
#      Add standardized country fields ----------------------------------------
# ------------------------------------------------------------------------------|
records_with_keys <- records %>%
  mutate(
    accession_key = dplyr::coalesce(accession_key, accession_version, primary_accession),
    country = clean_text(country),
    lat_lon = clean_text(lat_lon)
  )

lat_lon_lookup <- records_with_keys %>%
  filter(is.na(country), !is.na(lat_lon)) %>%
  lookup_country_from_lat_lon()

standardized_records <- records_with_keys %>%
  left_join(lat_lon_lookup, by = "accession_key") %>%
  mutate(
    country_standardized_candidate = dplyr::coalesce(country, country_from_lat_lon),
    country_standardized = standardize_country_value(country_standardized_candidate),
    country_source = case_when(
      !is.na(country) ~ "parsed_geo_loc_name",
      is.na(country) & !is.na(country_from_lat_lon) ~ "lat_lon_polygon_lookup",
      TRUE ~ NA_character_
    ),
    country_status = case_when(
      is.na(country_standardized_candidate) ~ "missing",
      country_standardized_candidate %in% c("missing", "unknown", "Unknown", "not provided") ~ "missing",
      country_standardized_candidate %in% ocean_values ~ "ocean_or_non_country",
      country_standardized_candidate %in% historical_values ~ "historical",
      country_standardized %in% territory_values ~ "territory",
      country_standardized %in% review_values ~ "review",
      is.na(country_standardized) ~ "missing",
      TRUE ~ "country"
    )
  ) %>%
  select(
    everything(),
    country_standardized,
    country_status,
    country_source,
    country_from_lat_lon
  )

# ------------------------------------------------------------------------------|
#      Summarize standardized evidence ----------------------------------------
# ------------------------------------------------------------------------------|
standardized_records <- standardized_records %>%
  mutate(
    country_standardized = if_else(
      country_status %in% c("missing", "ocean_or_non_country"),
      NA_character_,
      country_standardized
    )
  )

pathogen_country_summary <- standardized_records %>%
  filter(country_status %in% c("country", "territory", "historical", "review")) %>%
  group_by(target_id, Pathogens, Disease_name, country_standardized, country_status) %>%
  summarise(
    records_with_country = n(),
    accessions = collapse_unique(accession_key),
    organisms = collapse_unique(organism),
    hosts = collapse_unique(host),
    in_gibb_etal = any(as_logical_flag(in_gibb_etal), na.rm = TRUE),
    in_empres_i = any(as_logical_flag(in_empres_i), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Pathogens, Disease_name, country_standardized)

disease_country_summary <- pathogen_country_summary %>%
  group_by(Disease_name, country_standardized, country_status) %>%
  summarise(
    records_with_country = sum(records_with_country, na.rm = TRUE),
    pathogens = collapse_unique(Pathogens),
    target_ids = collapse_unique(target_id),
    in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
    in_empres_i = any(in_empres_i, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Disease_name, country_standardized)

qa_summary <- bind_rows(
  standardized_records %>%
    count(country_status, name = "records") %>%
    mutate(metric = paste0("status_", country_status)) %>%
    select(metric, value = records),
  tibble(
    metric = c(
      "records_total",
      "records_with_standardized_country",
      "records_recovered_from_lat_lon",
      "unique_country_standardized"
    ),
    value = c(
      nrow(standardized_records),
      sum(!is.na(standardized_records$country_standardized)),
      sum(standardized_records$country_source == "lat_lon_polygon_lookup", na.rm = TRUE),
      n_distinct(standardized_records$country_standardized, na.rm = TRUE)
    )
  )
)

# ------------------------------------------------------------------------------|
#      Write outputs -----------------------------------------------------------
# ------------------------------------------------------------------------------|
write_csv(
  standardized_records,
  genbank_simple_file_path(
    output_dir,
    paste0(output_prefix, "country_records_standardized.csv"),
    create_parent = TRUE
  )
)
write_csv(
  pathogen_country_summary,
  genbank_simple_file_path(
    output_dir,
    paste0(output_prefix, "pathogen_country_summary_standardized.csv"),
    create_parent = TRUE
  )
)
write_csv(
  disease_country_summary,
  genbank_simple_file_path(
    output_dir,
    paste0(output_prefix, "disease_country_summary_standardized.csv"),
    create_parent = TRUE
  )
)
write_csv(
  qa_summary,
  genbank_simple_file_path(
    output_dir,
    paste0(output_prefix, "country_standardization_qa.csv"),
    create_parent = TRUE
  )
)

message("Summary kind: ", summary_kind)
message("Wrote standardized records: ", nrow(standardized_records))
message("Recovered country from lat_lon: ", sum(standardized_records$country_source == "lat_lon_polygon_lookup", na.rm = TRUE))
message("Wrote standardized pathogen-country rows: ", nrow(pathogen_country_summary))
message("Wrote standardized disease-country rows: ", nrow(disease_country_summary))
