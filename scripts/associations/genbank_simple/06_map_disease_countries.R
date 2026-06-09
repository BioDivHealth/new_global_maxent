# ------------------------------------------------------------------------------|
#      06_map_disease_countries.R ---------------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Map per-disease countries recovered by GenBank-simple country runs.
# Inputs : genbank_readiness_disease_country_summary_standardized.csv
# Outputs: maps/disease_country_records/*.png
#          maps/genbank_disease_country_map_countries.csv
#          maps/genbank_disease_country_map_unmatched.csv
#
# Notes  : This script maps countries from the additive standardized country
#          layer. It does not reinterpret raw GenBank locations.
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, ggplot2, here, purrr, readr, rnaturalearth, sf, stringr, tibble, tidyr)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))
source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Resolve run mode and map paths -----------------------------------------
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

map_dir <- genbank_simple_map_dir(summary_kind, output_dir)
disease_map_dir <- file.path(map_dir, "disease_country_records")
dir.create(disease_map_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- if_else(
  summary_kind == "readiness_combined",
  "genbank_readiness_disease_country_summary_standardized.csv",
  "genbank_disease_country_summary_standardized.csv"
)
summary_path <- genbank_simple_existing_file_path(output_dir, summary_file)

country_summary <- read_csv(
  summary_path,
  col_types = cols(
    Disease_name = col_character(),
    country_standardized = col_character(),
    country_status = col_character(),
    records_with_country = col_double(),
    pathogens = col_character(),
    target_ids = col_character(),
    in_gibb_etal = col_logical(),
    in_empres_i = col_logical()
  ),
  na = c("", "NA")
) %>%
  mutate(
    Disease_name = clean_text(Disease_name),
    country_standardized = clean_text(country_standardized),
    country_status = clean_text(country_status)
  ) %>%
  filter(
    !is.na(Disease_name),
    !is.na(country_standardized),
    country_status %in% c("country", "territory", "historical", "review")
  )

if (nrow(country_summary) == 0) {
  stop("No disease-country rows available to map in: ", summary_path)
}

# ------------------------------------------------------------------------------|
#      Prepare country-name map joins -----------------------------------------
# ------------------------------------------------------------------------------|
country_name_overrides <- tibble::tribble(
  ~country_standardized, ~map_country,
  "Brunei", "Brunei Darussalam",
  "Cape Verde", "Cabo Verde",
  "Curacao", "Curaçao",
  "Cote d'Ivoire", "Côte d'Ivoire",
  "Eswatini", "eSwatini",
  "Micronesia", "Federated States of Micronesia",
  "Laos", "Lao PDR",
  "South Georgia and the South Sandwich Islands", "South Georgia and the Islands",
  "United States", "United States of America",
  "Virgin Islands", "United States Virgin Islands"
)

plot_country_summary <- country_summary %>%
  group_by(Disease_name, country_standardized, country_status) %>%
  summarise(
    records_with_country = sum(records_with_country, na.rm = TRUE),
    pathogens = collapse_unique(pathogens),
    target_ids = collapse_unique(target_ids),
    in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
    in_empres_i = any(in_empres_i, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(country_name_overrides, by = "country_standardized") %>%
  mutate(map_country = dplyr::coalesce(map_country, country_standardized))

# ------------------------------------------------------------------------------|
#      Join GenBank countries to world geometry -------------------------------
# ------------------------------------------------------------------------------|
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  select(
    map_country = name_long,
    map_admin = admin,
    map_sovereignt = sovereignt,
    geometry
  ) %>%
  mutate(
    across(c(map_country, map_admin, map_sovereignt), clean_text),
    map_join_names = purrr::pmap(
      list(map_country, map_admin, map_sovereignt),
      ~ unique(stats::na.omit(c(...)))
    )
  ) %>%
  select(map_country, map_join_names, geometry) %>%
  tidyr::unnest(map_join_names) %>%
  distinct(map_join_names, .keep_all = TRUE)

mapped_countries <- plot_country_summary %>%
  left_join(world, by = c("map_country" = "map_join_names")) %>%
  mutate(map_matched = !sf::st_is_empty(geometry) & !is.na(sf::st_dimension(geometry))) %>%
  sf::st_as_sf()

unmatched_countries <- mapped_countries %>%
  filter(!map_matched) %>%
  sf::st_drop_geometry() %>%
  select(
    Disease_name,
    country_standardized,
    country_status,
    records_with_country,
    pathogens,
    target_ids
  ) %>%
  arrange(Disease_name, country_standardized)

# ------------------------------------------------------------------------------|
#      Write map tables and per-disease PNGs ----------------------------------
# ------------------------------------------------------------------------------|
map_country_records <- mapped_countries %>%
  filter(map_matched) %>%
  sf::st_drop_geometry() %>%
  select(
    Disease_name,
    country_standardized,
    country_status,
    records_with_country,
    pathogens,
    target_ids,
    in_gibb_etal,
    in_empres_i
  ) %>%
  arrange(Disease_name, country_standardized)

write_csv(map_country_records, file.path(map_dir, "genbank_disease_country_map_countries.csv"))
write_csv(unmatched_countries, file.path(map_dir, "genbank_disease_country_map_unmatched.csv"))

world_base <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

plot_one_disease <- function(disease_name) {
  disease_countries <- mapped_countries %>%
    filter(Disease_name == disease_name, map_matched)

  output_path <- file.path(
    disease_map_dir,
    paste0(sanitize_filename(disease_name), ".png")
  )

  if (nrow(disease_countries) == 0) {
    return(tibble(Disease_name = disease_name, output_path = output_path, mapped_countries = 0L))
  }

  subtitle <- paste0(
    nrow(disease_countries),
    " mapped countries/territories; fill is log10 GenBank records with standardized country metadata"
  )

  plot <- ggplot() +
    geom_sf(data = world_base, fill = "grey94", color = "white", linewidth = 0.12) +
    geom_sf(
      data = disease_countries,
      aes(fill = records_with_country),
      color = "grey25",
      linewidth = 0.08
    ) +
    scale_fill_viridis_c(
      trans = "log10",
      option = "magma",
      direction = -1,
      na.value = "grey94",
      name = "Records"
    ) +
    coord_sf(crs = "+proj=robin", datum = NA) +
    labs(
      title = disease_name,
      subtitle = subtitle,
      caption = "Source: GenBank-simple standardized country summaries"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30"),
      plot.caption = element_text(size = 7, color = "grey45"),
      legend.position = "bottom",
      legend.key.width = unit(1.6, "cm")
    )

  ggsave(
    filename = output_path,
    plot = plot,
    width = 11,
    height = 6.2,
    dpi = 300,
    bg = "white"
  )

  tibble(
    Disease_name = disease_name,
    output_path = output_path,
    mapped_countries = nrow(disease_countries)
  )
}

map_manifest <- purrr::map_dfr(
  sort(unique(plot_country_summary$Disease_name)),
  plot_one_disease
)

write_csv(map_manifest, file.path(map_dir, "genbank_disease_country_map_manifest.csv"))

message("Summary kind: ", summary_kind)
message("Wrote disease maps: ", nrow(map_manifest))
message("Wrote mapped country rows: ", nrow(map_country_records))
message("Wrote unmatched country rows: ", nrow(unmatched_countries))
