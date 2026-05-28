# ------------------------------------------------------------------------------
# Build a first-pass VectorMap vector-host table from direct host-linked layers
# ------------------------------------------------------------------------------

library(pacman)
p_load(data.table, dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))

# Normalize text while keeping this export close to the source tables.
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

clean_names_bom <- function(dt) {
  data.table::setnames(dt, names(dt), str_replace(names(dt), "^\ufeff", ""))
  dt
}

output_dir <- vectormap_outputs_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bloodmeal_path <- file.path(vectormap_dir, "BloodMealMap_Layer_-3496204453665016601.csv")
tick_path <- file.path(vectormap_dir, "TickMap_4464597498443279194.csv")
flea_path <- file.path(vectormap_dir, "FleaMap_-6875364799851429947.csv")
mite_path <- file.path(vectormap_dir, "MiteMap_-6324776740768397246.csv")

output_path <- file.path(output_dir, "vectormap_vector_host_links_raw.csv")

target_cols <- c(
  "source_dataset",
  "record_id",
  "interaction_type",
  "basis_of_record",
  "source_citation",
  "submitter_organization",
  "submitter_person",
  "collector",
  "institution_code",
  "collection_code",
  "catalog_number",
  "vector_family",
  "vector_genus",
  "vector_subgenus",
  "vector_species",
  "vector_scientific_name",
  "vector_common_name",
  "vector_identification_method",
  "host_family",
  "host_genus",
  "host_species",
  "host_scientific_name",
  "host_common_name",
  "host_identification_method",
  "verbatim_host_name",
  "host_individual_count",
  "host_count_infested",
  "associated_pathogen",
  "pathogen_identification_method",
  "pathogen_tester_organization",
  "remarks_about_pathogen",
  "associated_parasite",
  "parasite_identification_method",
  "parasite_tested_by",
  "number_parasite_positive",
  "number_specimens_parasite_tested",
  "country",
  "state_province",
  "county",
  "locality",
  "latitude",
  "longitude",
  "coordinate_uncertainty_m",
  "verbatim_collecting_date",
  "earliest_date_collected",
  "latest_date_collected",
  "time_collected",
  "collecting_method",
  "individual_count",
  "remarks"
)

add_missing_target_cols <- function(df) {
  missing_cols <- setdiff(target_cols, names(df))

  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      df[[col]] <- NA_character_
    }
  }

  df %>% select(all_of(target_cols))
}

read_vectormap_csv <- function(path, select_cols = NULL) {
  # Each VectorMap export uses a slightly different schema, so only request
  # columns that are actually present in the current file.
  header_dt <- fread(
    path,
    encoding = "UTF-8",
    nrows = 0,
    fill = TRUE,
    showProgress = FALSE
  )

  available_cols <- names(clean_names_bom(header_dt))
  select_cols <- if (is.null(select_cols)) {
    NULL
  } else {
    intersect(select_cols, available_cols)
  }

  dt <- fread(
    path,
    encoding = "UTF-8",
    na.strings = c("", "NA", "NaN", "No data", "null", "Null"),
    select = select_cols,
    fill = TRUE,
    showProgress = FALSE
  )

  clean_names_bom(dt)
}

read_vectormap_csv_readr <- function(path, select_cols = NULL) {
  header <- readr::read_csv(
    path,
    n_max = 0,
    show_col_types = FALSE,
    progress = FALSE
  )

  select_cols <- if (is.null(select_cols)) {
    NULL
  } else {
    intersect(select_cols, names(header))
  }

  readr::read_csv(
    path,
    col_select = dplyr::any_of(select_cols),
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 1000
  )
}

bloodmeal <- read_vectormap_csv(
  bloodmeal_path,
  select_cols = c(
    "FID",
    "SubmitterOrganization",
    "SubmitterPerson",
    "Collector",
    "InstitutionCode",
    "CollectionCode",
    "CatalogNumber",
    "BasisOfRecord",
    "Citation",
    "PublicationStaus",
    "HematophageFamily",
    "HematophageGenus",
    "HematophageSubGenus",
    "HematophageSpecies",
    "HematopageScientificName",
    "HematopageCommonName",
    "HematophageIdentificationMethod",
    "HostFamily",
    "HostGenus",
    "HostSpecies",
    "HostScientificName",
    "HostCommonName",
    "HostIdentificationMethod",
    "VerbatimCollectingDate",
    "EarliestDateCollected",
    "LatestDateCollected",
    "TimeCollected",
    "Country",
    "StateProvince",
    "County",
    "Locality",
    "DecimalLongitude",
    "DecimalLatitude",
    "CoordinateUncertaintyInMeters",
    "HematophageCollectingMethod",
    "HematophageIndividualCount",
    "HostIndividualCount",
    "Remarks",
    "HabitatType",
    "Elevation",
    "AssociatedPathogen",
    "IdentificationMethodForPathogen",
    "PathogenTesterOrganization",
    "RemarksAboutPathogen"
  )
) %>%
  as_tibble() %>%
  mutate(across(where(is.character), clean_text)) %>%
  transmute(
    source_dataset = "BloodMealMap",
    record_id = FID,
    interaction_type = "blood_meal",
    basis_of_record = BasisOfRecord,
    source_citation = Citation,
    submitter_organization = SubmitterOrganization,
    submitter_person = SubmitterPerson,
    collector = Collector,
    institution_code = InstitutionCode,
    collection_code = CollectionCode,
    catalog_number = CatalogNumber,
    vector_family = HematophageFamily,
    vector_genus = HematophageGenus,
    vector_subgenus = HematophageSubGenus,
    vector_species = HematophageSpecies,
    vector_scientific_name = HematopageScientificName,
    vector_common_name = HematopageCommonName,
    vector_identification_method = HematophageIdentificationMethod,
    host_family = HostFamily,
    host_genus = HostGenus,
    host_species = HostSpecies,
    host_scientific_name = HostScientificName,
    host_common_name = HostCommonName,
    host_identification_method = HostIdentificationMethod,
    verbatim_host_name = NA_character_,
    host_individual_count = as.character(HostIndividualCount),
    host_count_infested = NA_character_,
    associated_pathogen = AssociatedPathogen,
    pathogen_identification_method = IdentificationMethodForPathogen,
    pathogen_tester_organization = PathogenTesterOrganization,
    remarks_about_pathogen = RemarksAboutPathogen,
    associated_parasite = NA_character_,
    parasite_identification_method = NA_character_,
    parasite_tested_by = NA_character_,
    number_parasite_positive = NA_character_,
    number_specimens_parasite_tested = NA_character_,
    country = Country,
    state_province = StateProvince,
    county = County,
    locality = Locality,
    latitude = as.character(DecimalLatitude),
    longitude = as.character(DecimalLongitude),
    coordinate_uncertainty_m = as.character(CoordinateUncertaintyInMeters),
    verbatim_collecting_date = VerbatimCollectingDate,
    earliest_date_collected = EarliestDateCollected,
    latest_date_collected = LatestDateCollected,
    time_collected = TimeCollected,
    collecting_method = HematophageCollectingMethod,
    individual_count = as.character(HematophageIndividualCount),
    remarks = Remarks
  ) %>%
  add_missing_target_cols()

read_host_attached_layer <- function(path, source_name) {
  host_layer_cols <- c(
    "ObjectID",
    "OBJECTID",
    "OBJECTID_1",
    "SubmitterOrganization",
    "SubmitterPerson",
    "Collector",
    "InstitutionCode",
    "CollectionCode",
    "CatalogNumber",
    "BasisOfRecord",
    "Source",
    "Family",
    "Genus",
    "SubGenus",
    "Species",
    "ScientificName",
    "IdentificationMethod",
    "EarliestDateCollected",
    "LatestDateCollected",
    "VerbatimCollectingDate",
    "TimeCollected",
    "Country",
    "StateProvince",
    "County",
    "Locality",
    "DecimalLongitude",
    "DecimalLatitude",
    "CoordinateUncertaintyInMeters",
    "IndividualCount",
    "Remarks",
    "CollectingMethod",
    "VerbatimHostName",
    "HostFamily",
    "HostScientificName",
    "HostCommonName",
    "HostIndividualCount",
    "HostCountInfested",
    "AssociatedParasite",
    "NumberParasitePositive",
    "IdentificationMethodForParasite",
    "ParasiteTestedBy",
    "NumberSpecimensParasiteTested"
  )

  # TickMap has quoting/row-width irregularities that `readr` handles more
  # gracefully than `fread`, so keep a tolerant path here.
  dt <- if (identical(source_name, "TickMap")) {
    read_vectormap_csv_readr(path, select_cols = host_layer_cols)
  } else {
    read_vectormap_csv(path, select_cols = host_layer_cols)
  }

  record_id_col <- c("ObjectID", "OBJECTID", "OBJECTID_1")
  record_id_col <- record_id_col[record_id_col %in% names(dt)][1]
  vector_subgenus_col <- if ("SubGenus" %in% names(dt)) "SubGenus" else NULL

  tibble::as_tibble(dt) %>%
    mutate(across(where(is.character), clean_text)) %>%
    transmute(
      source_dataset = source_name,
      record_id = .data[[record_id_col]],
      interaction_type = "on_host_occurrence",
      basis_of_record = BasisOfRecord,
      source_citation = Source,
      submitter_organization = SubmitterOrganization,
      submitter_person = SubmitterPerson,
      collector = Collector,
      institution_code = InstitutionCode,
      collection_code = CollectionCode,
      catalog_number = CatalogNumber,
      vector_family = Family,
      vector_genus = Genus,
      vector_subgenus = if (!is.null(vector_subgenus_col)) .data[[vector_subgenus_col]] else NA_character_,
      vector_species = Species,
      vector_scientific_name = ScientificName,
      vector_common_name = NA_character_,
      vector_identification_method = IdentificationMethod,
      host_family = if ("HostFamily" %in% names(dt)) HostFamily else NA_character_,
      host_genus = NA_character_,
      host_species = NA_character_,
      host_scientific_name = if ("HostScientificName" %in% names(dt)) HostScientificName else NA_character_,
      host_common_name = if ("HostCommonName" %in% names(dt)) HostCommonName else NA_character_,
      host_identification_method = NA_character_,
      verbatim_host_name = if ("VerbatimHostName" %in% names(dt)) VerbatimHostName else NA_character_,
      host_individual_count = if ("HostIndividualCount" %in% names(dt)) as.character(HostIndividualCount) else NA_character_,
      host_count_infested = if ("HostCountInfested" %in% names(dt)) as.character(HostCountInfested) else NA_character_,
      associated_pathogen = NA_character_,
      pathogen_identification_method = NA_character_,
      pathogen_tester_organization = NA_character_,
      remarks_about_pathogen = NA_character_,
      associated_parasite = if ("AssociatedParasite" %in% names(dt)) AssociatedParasite else NA_character_,
      parasite_identification_method = if ("IdentificationMethodForParasite" %in% names(dt)) IdentificationMethodForParasite else NA_character_,
      parasite_tested_by = if ("ParasiteTestedBy" %in% names(dt)) ParasiteTestedBy else NA_character_,
      number_parasite_positive = if ("NumberParasitePositive" %in% names(dt)) as.character(NumberParasitePositive) else NA_character_,
      number_specimens_parasite_tested = if ("NumberSpecimensParasiteTested" %in% names(dt)) as.character(NumberSpecimensParasiteTested) else NA_character_,
      country = Country,
      state_province = StateProvince,
      county = County,
      locality = Locality,
      latitude = as.character(DecimalLatitude),
      longitude = as.character(DecimalLongitude),
      coordinate_uncertainty_m = as.character(CoordinateUncertaintyInMeters),
      verbatim_collecting_date = VerbatimCollectingDate,
      earliest_date_collected = EarliestDateCollected,
      latest_date_collected = LatestDateCollected,
      time_collected = TimeCollected,
      collecting_method = CollectingMethod,
      individual_count = if ("IndividualCount" %in% names(dt)) as.character(IndividualCount) else NA_character_,
      remarks = Remarks
    ) %>%
    add_missing_target_cols()
}

tickmap <- read_host_attached_layer(tick_path, "TickMap")
fleamap <- read_host_attached_layer(flea_path, "FleaMap")
mitemap <- read_host_attached_layer(mite_path, "MiteMap")

vectormap_vector_host_links_raw <- bind_rows(
  bloodmeal,
  tickmap,
  fleamap,
  mitemap
) %>%
  mutate(across(where(is.character), clean_text))

write_csv(vectormap_vector_host_links_raw, output_path, na = "")

cat("Rows written:", nrow(vectormap_vector_host_links_raw), "\n")
print(vectormap_vector_host_links_raw %>% count(source_dataset, sort = TRUE))
cat("Unique vector scientific names:", n_distinct(vectormap_vector_host_links_raw$vector_scientific_name, na.rm = TRUE), "\n")
cat("Unique host scientific names:", n_distinct(vectormap_vector_host_links_raw$host_scientific_name, na.rm = TRUE), "\n")
cat("Wrote raw links to", output_path, "\n")
