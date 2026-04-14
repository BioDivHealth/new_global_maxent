# ------------------------------------------------------------------------------
# 5_7_GenBank_Pathogen_Country_Metadata.R
# ------------------------------------------------------------------------------
# Purpose: Build GenBank-ready queries for the disease-pathogen targets in the
#          combined WHO network, fetch accession/sample-level country metadata
#          from NCBI, and collapse the results into pathogen-country and
#          disease-country summary tables for downstream geographic review.
#
# Inputs : combined_who_network_canonical_zoonotic.csv
#          who_pathogens_diseases_zoonotic.csv
#          who_pathogens_virion_taxid.csv
#          who_bacteria_clover_taxid.csv
# Outputs: genbank_pathogen_query_manifest.csv
#          genbank_nuccore_search_log.csv
#          genbank_nuccore_country_records.csv
#          genbank_nuccore_country_summary.csv
#          genbank_network_disease_country_summary.csv
#
# Notes  : Set `NCBI_API_KEY` or `ENTREZ_KEY` before running for faster,
#          higher-volume requests. Optional environment variables:
#          `GENBANK_PATHOGEN_FILTER`, `GENBANK_MAX_PATHOGENS`,
#          `GENBANK_RETMAX`, `GENBANK_BATCH_SIZE`,
#          `GENBANK_PLATEAU_BATCHES`, `GENBANK_MIN_NEW_COUNTRIES`,
#          `GENBANK_MIN_BATCHES_BEFORE_STOP`,
#          `GENBANK_MIN_RECORDS_BEFORE_STOP`,
#          `GENBANK_INTERVAL_SAMPLE_THRESHOLD`, `GENBANK_QUERY_ONLY`,
#          `GENBANK_RESUME`, `GENBANK_SECOND_PASS_ONLY`,
#          `GENBANK_SECOND_PASS_PATH`, `GENBANK_FORCE_RERUN`,
#          `GENBANK_USE_RECOMMENDED_RETMAX`.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, rentrez, stringr, tibble, tidyr, xml2)

source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------
# Shared helpers ---------------------------------------------------------------
# ------------------------------------------------------------------------------
source(here("scripts", "associations", "genbank", "genbank_metadata_helpers.R"))

sanitize_filename <- function(x) {
  x <- clean_text(x)

  if (is.na(x)) {
    return("missing")
  }

  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_{2,}", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x <- clean_text(x)

  if (is.na(x)) {
    return("missing")
  }

  x
}

read_checkpoint_csv <- function(path, template) {
  if (!file.exists(path)) {
    return(template %>% slice(0))
  }

  checkpoint <- read_csv(path, show_col_types = FALSE, na = c("", "NA"))

  if (nrow(checkpoint) == 0) {
    return(template %>% slice(0))
  }

  missing_cols <- setdiff(names(template), names(checkpoint))
  extra_cols <- setdiff(names(checkpoint), names(template))

  if (length(missing_cols) > 0) {
    for (col_name in missing_cols) {
      checkpoint[[col_name]] <- NA
    }
  }

  if (length(extra_cols) > 0) {
    checkpoint <- checkpoint %>% select(-all_of(extra_cols))
  }

  checkpoint <- checkpoint %>% select(all_of(names(template)))

  for (col_name in names(template)) {
    template_col <- template[[col_name]]

    if (is.character(template_col)) {
      checkpoint[[col_name]] <- as.character(checkpoint[[col_name]])
    } else if (is.integer(template_col)) {
      checkpoint[[col_name]] <- suppressWarnings(as.integer(checkpoint[[col_name]]))
    } else if (is.double(template_col)) {
      checkpoint[[col_name]] <- suppressWarnings(as.double(checkpoint[[col_name]]))
    } else if (is.logical(template_col)) {
      checkpoint[[col_name]] <- as.logical(checkpoint[[col_name]])
    }
  }

  checkpoint
}

checkpoint_is_reusable <- function(log_path, records_path, log_template) {
  if (!file.exists(log_path) || !file.exists(records_path)) {
    return(FALSE)
  }

  checkpoint_log <- read_checkpoint_csv(log_path, log_template)

  if (nrow(checkpoint_log) == 0) {
    return(FALSE)
  }

  statuses <- clean_text(checkpoint_log$status)

  if (all(is.na(statuses))) {
    return(FALSE)
  }

  !any(stringr::str_detect(statuses, "search_failed"), na.rm = TRUE)
}

build_archive_dir <- function(pathogen_runs_dir) {
  base_name <- paste0("archive_", Sys.Date())
  archive_dir <- file.path(pathogen_runs_dir, base_name)
  suffix <- 1L

  while (file.exists(archive_dir)) {
    suffix <- suffix + 1L
    archive_dir <- file.path(pathogen_runs_dir, paste0(base_name, "_", suffix))
  }

  archive_dir
}

archive_stale_checkpoint_artifacts <- function(
    pathogen_runs_dir,
    pathogen_log_dir,
    pathogen_records_dir,
    active_log_paths,
    active_record_paths,
    verbose = TRUE
) {
  active_log_paths <- normalizePath(active_log_paths, winslash = "/", mustWork = FALSE)
  active_record_paths <- normalizePath(active_record_paths, winslash = "/", mustWork = FALSE)

  log_entries <- list.files(pathogen_log_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  record_entries <- list.files(pathogen_records_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)

  keep_top_level_entry <- function(path) {
    stringr::str_detect(basename(path), "^archive_")
  }

  orphan_log_entries <- log_entries[
    !normalizePath(log_entries, winslash = "/", mustWork = FALSE) %in% active_log_paths &
      !vapply(log_entries, keep_top_level_entry, logical(1))
  ]

  orphan_record_entries <- record_entries[
    !normalizePath(record_entries, winslash = "/", mustWork = FALSE) %in% active_record_paths &
      !vapply(record_entries, keep_top_level_entry, logical(1))
  ]

  if (length(orphan_log_entries) == 0 && length(orphan_record_entries) == 0) {
    return(invisible(NULL))
  }

  archive_dir <- build_archive_dir(pathogen_runs_dir)
  archive_log_dir <- file.path(archive_dir, "search_logs")
  archive_record_dir <- file.path(archive_dir, "country_records")

  dir.create(archive_log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(archive_record_dir, recursive = TRUE, showWarnings = FALSE)

  move_entry <- function(entry_path, destination_dir) {
    target_path <- file.path(destination_dir, basename(entry_path))

    if (file.exists(target_path) || dir.exists(target_path)) {
      stop("Archive target already exists: ", target_path)
    }

    moved <- file.rename(entry_path, target_path)

    if (!isTRUE(moved)) {
      stop("Failed to archive stale checkpoint artifact: ", entry_path)
    }
  }

  purrr::walk(orphan_log_entries, move_entry, destination_dir = archive_log_dir)
  purrr::walk(orphan_record_entries, move_entry, destination_dir = archive_record_dir)

  log_progress(
    "Archived ",
    length(orphan_log_entries),
    " stale search logs and ",
    length(orphan_record_entries),
    " stale record artifacts to ",
    archive_dir,
    verbose = verbose
  )

  invisible(archive_dir)
}

# ------------------------------------------------------------------------------
# Parameters and paths ---------------------------------------------------------
# ------------------------------------------------------------------------------
max_pathogens <- parse_env_integer("GENBANK_MAX_PATHOGENS", default = NA_integer_)
retmax <- parse_env_integer("GENBANK_RETMAX", default = 15000L)
batch_size <- parse_env_integer("GENBANK_BATCH_SIZE", default = 250L)
plateau_batches <- parse_env_integer("GENBANK_PLATEAU_BATCHES", default = 3L)
min_new_countries <- parse_env_integer("GENBANK_MIN_NEW_COUNTRIES", default = 0L)
min_batches_before_stop <- parse_env_integer("GENBANK_MIN_BATCHES_BEFORE_STOP", default = 3L)
min_records_before_stop <- parse_env_integer("GENBANK_MIN_RECORDS_BEFORE_STOP", default = 7500L)
strict_plateau_batches <- parse_env_integer("GENBANK_STRICT_PLATEAU_BATCHES", default = 5L)
strict_min_batches_before_stop <- parse_env_integer("GENBANK_STRICT_MIN_BATCHES_BEFORE_STOP", default = 5L)
strict_min_records_before_stop <- parse_env_integer("GENBANK_STRICT_MIN_RECORDS_BEFORE_STOP", default = 15000L)
interval_sample_threshold <- parse_env_integer("GENBANK_INTERVAL_SAMPLE_THRESHOLD", default = 50000L)
query_only <- parse_env_flag("GENBANK_QUERY_ONLY", default = FALSE)
resume_mode <- parse_env_flag("GENBANK_RESUME", default = TRUE)
verbose <- parse_env_flag("GENBANK_VERBOSE", default = TRUE)
second_pass_only <- parse_env_flag("GENBANK_SECOND_PASS_ONLY", default = FALSE)
pathogen_filter <- clean_text(Sys.getenv("GENBANK_PATHOGEN_FILTER", unset = NA_character_))
#pathogen_filter <- 'Human mastadenovirus B|Mastadenovirus blackbeardi serotype 14|Mammarenavirus lassaense|Mammarenavirus lujoense|Mamastrovirus virginiaense|Mamastrovirus 9|Vibrio cholerae serogroup 0139|Klebsiella pneumoniae|Salmonella enterica non typhoidal serovars|Yersinia pestis|Shigella dysenteriae serotype 1|Orthobornavirus bornaense|Subgenus Sarbecovirus|subgenus Merbecovirus|Orthoebolavirus zairense|Orthomarburgvirus marburgense|Orthoflavivirus denguei|Orthoflavivirus encephalitidis|Orthoflavivirus nilense|Orthoflavivirus flavi'
if (retmax <= 0) {
  stop("GENBANK_RETMAX must be greater than 0.")
}

if (batch_size <= 0) {
  stop("GENBANK_BATCH_SIZE must be greater than 0.")
}

if (plateau_batches <= 0) {
  stop("GENBANK_PLATEAU_BATCHES must be greater than 0.")
}

if (min_batches_before_stop <= 0) {
  stop("GENBANK_MIN_BATCHES_BEFORE_STOP must be greater than 0.")
}

if (min_records_before_stop <= 0) {
  stop("GENBANK_MIN_RECORDS_BEFORE_STOP must be greater than 0.")
}

if (interval_sample_threshold <= 0) {
  stop("GENBANK_INTERVAL_SAMPLE_THRESHOLD must be greater than 0.")
}

if (strict_plateau_batches <= 0) {
  stop("GENBANK_STRICT_PLATEAU_BATCHES must be greater than 0.")
}

if (strict_min_batches_before_stop <= 0) {
  stop("GENBANK_STRICT_MIN_BATCHES_BEFORE_STOP must be greater than 0.")
}

if (strict_min_records_before_stop <= 0) {
  stop("GENBANK_STRICT_MIN_RECORDS_BEFORE_STOP must be greater than 0.")
}

if (min_new_countries < 0) {
  stop("GENBANK_MIN_NEW_COUNTRIES cannot be negative.")
}

api_key <- clean_text(Sys.getenv("NCBI_API_KEY", unset = Sys.getenv("ENTREZ_KEY", unset = NA_character_)))
if (is.na(api_key)) {
  dotenv_path <- here(".env")
  api_key <- dplyr::coalesce(
    read_dotenv_value(dotenv_path, "NCBI_API_KEY"),
    read_dotenv_value(dotenv_path, "ENTREZ_KEY"),
    read_dotenv_value(dotenv_path, "ncbi_api_key"),
    read_dotenv_value(dotenv_path, "entrez_key")
  )
}

if (!is.na(api_key)) {
  rentrez::set_entrez_key(api_key)
}

who_dir <- here("pathogen_association_data", "WHO")
genbank_dir <- file.path(who_dir, "genbank")
dir.create(genbank_dir, recursive = TRUE, showWarnings = FALSE)

who_path <- who_working_pathogens_path()

virion_path <- file.path(who_dir, "virion", "who_pathogens_virion_taxid.csv")
clover_path <- file.path(who_dir, "clover", "who_bacteria_clover_taxid.csv")
network_path <- who_working_network_path()

manifest_output_path <- file.path(genbank_dir, "genbank_pathogen_query_manifest.csv")
search_log_output_path <- file.path(genbank_dir, "genbank_nuccore_search_log.csv")
records_output_path <- file.path(genbank_dir, "genbank_nuccore_country_records.csv")
summary_output_path <- file.path(genbank_dir, "genbank_nuccore_country_summary.csv")
disease_summary_output_path <- file.path(genbank_dir, "genbank_network_disease_country_summary.csv")
second_pass_output_path <- file.path(genbank_dir, "genbank_nuccore_second_pass_candidates.csv")
specificity_review_output_path <- file.path(genbank_dir, "genbank_query_specificity_review.csv")
second_pass_input_path <- clean_text(Sys.getenv("GENBANK_SECOND_PASS_PATH", unset = second_pass_output_path))
force_rerun <- parse_env_flag("GENBANK_FORCE_RERUN", default = second_pass_only)
use_recommended_retmax <- parse_env_flag("GENBANK_USE_RECOMMENDED_RETMAX", default = TRUE)
pathogen_runs_dir <- file.path(genbank_dir, "pathogen_runs")
pathogen_log_dir <- file.path(pathogen_runs_dir, "search_logs")
pathogen_records_dir <- file.path(pathogen_runs_dir, "country_records")

dir.create(pathogen_runs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathogen_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathogen_records_dir, recursive = TRUE, showWarnings = FALSE)

empty_search_log <- tibble(
  Pathogens = character(),
  Disease_name = character(),
  network_pathogen = character(),
  network_disease_name = character(),
  Family = character(),
  who_match_status = character(),
  who_match_count = integer(),
  who_pathogen_candidates = character(),
  target_resolution_type = character(),
  query_profile = character(),
  influenza_search_mode = character(),
  preferred_metadata_source = character(),
  metadata_source_reason = character(),
  query_strategy = character(),
  query_used = character(),
  records_found = integer(),
  records_fetched = integer(),
  batches_fetched = integer(),
  countries_observed = integer(),
  organism_count_observed = integer(),
  last_batch_new_countries = integer(),
  query_specificity_status = character(),
  query_specificity_reason = character(),
  status = character(),
  note = character()
)

empty_country_records <- tibble(
  Pathogens = character(),
  Disease_name = character(),
  network_pathogen = character(),
  network_disease_name = character(),
  Family = character(),
  who_match_status = character(),
  who_match_count = integer(),
  who_pathogen_candidates = character(),
  target_resolution_type = character(),
  query_profile = character(),
  influenza_search_mode = character(),
  influenza_query_label = character(),
  influenza_query_class = character(),
  preferred_metadata_source = character(),
  metadata_source_reason = character(),
  query_strategy = character(),
  query_used = character(),
  record_source_db = character(),
  biosample_accession = character(),
  biosample_id = character(),
  taxid_query = character(),
  organism_query = character(),
  all_fields_query = character(),
  name_query = character(),
  accession_version = character(),
  primary_accession = character(),
  definition = character(),
  organism = character(),
  taxonomy = character(),
  sequence_length = integer(),
  country_raw = character(),
  geo_loc_name_raw = character(),
  country = character(),
  lat_lon = character(),
  collection_date = character(),
  host = character(),
  isolate = character(),
  strain = character(),
  isolate_source = character(),
  db_xref = character(),
  fetch_error = character()
)

empty_country_summary <- tibble(
  Pathogens = character(),
  Disease_name = character(),
  Family = character(),
  preferred_metadata_source = character(),
  country = character(),
  accession_count = integer(),
  countries_raw = character(),
  organisms = character(),
  collection_dates = character(),
  hosts = character()
)

empty_disease_country_summary <- tibble(
  Disease_name = character(),
  country = character(),
  pathogen_count = integer(),
  record_count = integer(),
  pathogens = character(),
  source_dbs = character(),
  organisms = character(),
  collection_dates = character(),
  hosts = character()
)

empty_specificity_review <- tibble(
  Pathogens = character(),
  Disease_name = character(),
  network_pathogen = character(),
  network_disease_name = character(),
  who_match_status = character(),
  who_match_count = integer(),
  who_pathogen_candidates = character(),
  target_resolution_type = character(),
  query_profile = character(),
  influenza_search_mode = character(),
  query_used = character(),
  preferred_metadata_source = character(),
  records_found = integer(),
  records_fetched = integer(),
  countries_observed = integer(),
  organism_count_observed = integer(),
  query_specificity_status = character(),
  query_specificity_reason = character(),
  note = character()
)

# ------------------------------------------------------------------------------
# Load and harmonize WHO pathogen inputs ---------------------------------------
# ------------------------------------------------------------------------------
who_pathogens <- read_csv(
  who_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

network_data <- read_csv(
  network_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

virion_taxids <- read_csv(
  virion_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

clover_taxids <- read_csv(
  clover_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

who_manifest <- who_pathogens %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    Family = collapse_unique(Family),
    PHEIC_risk = collapse_unique(`PHEIC risk`),
    previous_name = collapse_unique(previous_name),
    msl39_viral_name = collapse_unique(msl39_viral_name),
    .groups = "drop"
  )

virion_manifest <- virion_taxids %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    virion_tax_ids = collapse_unique(VirusTaxID),
    virion_names = collapse_unique(Virion_VirusName),
    matched_name_types = collapse_unique(matched_name_type),
    .groups = "drop"
  )

clover_manifest <- clover_taxids %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    clover_tax_ids = collapse_unique(PathogenTaxID),
    clover_names = collapse_unique(Clover_PathogenName),
    .groups = "drop"
  )

query_manifest <- build_network_query_manifest(
  network_data = network_data,
  who_manifest = who_manifest,
  virion_manifest = virion_manifest,
  clover_manifest = clover_manifest
) %>%
  arrange(Family, Disease_name, Pathogens)

if (!is.na(pathogen_filter)) {
  query_manifest <- query_manifest %>%
    filter(
      stringr::str_detect(Pathogens, stringr::regex(pathogen_filter, ignore_case = TRUE)) |
        stringr::str_detect(Disease_name, stringr::regex(pathogen_filter, ignore_case = TRUE))
    )
}

if (!is.na(max_pathogens)) {
  query_manifest <- query_manifest %>%
    slice_head(n = max_pathogens)
}

query_manifest <- query_manifest %>%
  mutate(
    pathogen_run_id = sprintf(
      "network_v2_%03d_%s__%s",
      dplyr::row_number(),
      purrr::map_chr(Pathogens, sanitize_filename),
      purrr::map_chr(Disease_name, sanitize_filename)
    ),
    pathogen_log_path = file.path(pathogen_log_dir, paste0(pathogen_run_id, ".csv")),
    pathogen_records_path = file.path(pathogen_records_dir, paste0(pathogen_run_id, ".csv"))
  )

write_csv(query_manifest, manifest_output_path, na = "")

query_manifest <- query_manifest %>%
  mutate(
    retmax_override = NA_integer_,
    run_selected = TRUE
  )

if (second_pass_only) {
  if (!file.exists(second_pass_input_path)) {
    stop("GENBANK_SECOND_PASS_ONLY is TRUE but candidate file was not found: ", second_pass_input_path)
  }

  second_pass_input <- read_csv(
    second_pass_input_path,
    show_col_types = FALSE,
    na = c("", "NA")
  )

  if (!"Pathogens" %in% names(second_pass_input)) {
    stop("Second-pass candidate file is missing required column `Pathogens`: ", second_pass_input_path)
  }

  if (!"Disease_name" %in% names(second_pass_input)) {
    second_pass_input$Disease_name <- NA_character_
  }

  if (!"recommended_retmax" %in% names(second_pass_input)) {
    second_pass_input$recommended_retmax <- NA_integer_
  }

  second_pass_input <- second_pass_input %>%
    mutate(
      across(where(is.character), clean_text),
      recommended_retmax = suppressWarnings(as.integer(recommended_retmax)),
      disease_key = dplyr::coalesce(Disease_name, "__NA__")
    ) %>%
    filter(!is.na(Pathogens)) %>%
    group_by(Pathogens, disease_key) %>%
    summarise(
      recommended_retmax = suppressWarnings(max(recommended_retmax, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      recommended_retmax = dplyr::if_else(
        is.infinite(recommended_retmax),
        NA_integer_,
        as.integer(recommended_retmax)
      )
    )

  query_manifest <- query_manifest %>%
    mutate(disease_key = dplyr::coalesce(Disease_name, "__NA__")) %>%
    left_join(
      second_pass_input %>%
        transmute(
          Pathogens,
          disease_key,
          second_pass_recommended_retmax = recommended_retmax
        ),
      by = c("Pathogens", "disease_key")
    ) %>%
    mutate(
      run_selected = !is.na(second_pass_recommended_retmax),
      retmax_override = dplyr::if_else(
        run_selected & use_recommended_retmax,
        second_pass_recommended_retmax,
        retmax_override
      )
    ) %>%
    select(-disease_key, -second_pass_recommended_retmax)

  selected_n <- sum(query_manifest$run_selected, na.rm = TRUE)

  log_progress(
    "Second-pass mode enabled | candidate_file=",
    second_pass_input_path,
    " | rows_selected=",
    selected_n,
    " | force_rerun=",
    force_rerun,
    " | use_recommended_retmax=",
    use_recommended_retmax,
    verbose = verbose
  )

  if (selected_n == 0) {
    stop(
      "GENBANK_SECOND_PASS_ONLY is TRUE but no manifest rows matched the second-pass candidate file: ",
      second_pass_input_path
    )
  }
}

if (nrow(query_manifest) == 0) {
  write_csv(empty_search_log, search_log_output_path, na = "")
  write_csv(empty_country_records, records_output_path, na = "")
  write_csv(empty_country_summary, summary_output_path, na = "")
  write_csv(empty_disease_country_summary, disease_summary_output_path, na = "")
  write_csv(empty_specificity_review, specificity_review_output_path, na = "")

  cat("Wrote query manifest to", manifest_output_path, "\n")
  cat("No pathogens matched the current filters; no NCBI requests were made.\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# Query GenBank and extract country metadata -----------------------------------
# ------------------------------------------------------------------------------
if (query_only) {
  empty_log <- query_manifest %>%
    transmute(
      Pathogens,
      Disease_name,
      network_pathogen,
      network_disease_name,
      Family,
      who_match_status,
      who_match_count,
      who_pathogen_candidates,
      target_resolution_type,
      query_profile,
      influenza_search_mode,
      preferred_metadata_source,
      metadata_source_reason,
      query_strategy,
      query_used = search_query,
      records_found = NA_integer_,
      records_fetched = 0L,
      batches_fetched = 0L,
      countries_observed = 0L,
      organism_count_observed = 0L,
      last_batch_new_countries = 0L,
      query_specificity_status,
      query_specificity_reason,
      status = "query_only",
      note = NA_character_
    )

  write_csv(empty_log, search_log_output_path, na = "")
  write_csv(empty_country_records, records_output_path, na = "")
  write_csv(empty_country_summary, summary_output_path, na = "")
  write_csv(empty_disease_country_summary, disease_summary_output_path, na = "")
  write_csv(empty_specificity_review, specificity_review_output_path, na = "")
  write_csv(tibble(), second_pass_output_path, na = "")

  cat("Wrote query manifest to", manifest_output_path, "\n")
  cat("Query-only mode enabled; no NCBI requests were made.\n")
  quit(save = "no", status = 0)
}

query_rows <- split(query_manifest, seq_len(nrow(query_manifest)))

for (idx in seq_along(query_rows)) {
  manifest_row <- tibble::as_tibble(query_rows[[idx]])
  run_selected <- isTRUE(manifest_row$run_selected[[1]])
  row_retmax_override <- suppressWarnings(as.integer(manifest_row$retmax_override[[1]]))
  row_retmax <- retmax

  if (!is.na(row_retmax_override) && row_retmax_override > 0L) {
    row_retmax <- max(retmax, row_retmax_override)
  }

  if (!run_selected) {
    next
  }

  log_path <- manifest_row$pathogen_log_path[[1]]
  records_path <- manifest_row$pathogen_records_path[[1]]

  if (
    resume_mode &&
      !force_rerun &&
      checkpoint_is_reusable(
        log_path = log_path,
        records_path = records_path,
        log_template = empty_search_log
      )
  ) {
    log_progress(
      "[",
      idx,
      "/",
      length(query_rows),
      "] Reusing checkpoint for ",
      manifest_row$Pathogens,
      verbose = verbose
    )
    next
  }

  query_result <- run_manifest_query(
    manifest_row = manifest_row,
    retmax = row_retmax,
    batch_size = batch_size,
    plateau_batches = plateau_batches,
    min_new_countries = min_new_countries,
    min_batches_before_stop = min_batches_before_stop,
    min_records_before_stop = min_records_before_stop,
    strict_plateau_batches = strict_plateau_batches,
    strict_min_batches_before_stop = strict_min_batches_before_stop,
    strict_min_records_before_stop = strict_min_records_before_stop,
    interval_sample_threshold = interval_sample_threshold,
    manifest_index = as.integer(idx),
    manifest_total = length(query_rows),
    verbose = verbose
  )

  log_to_write <- query_result$log
  records_to_write <- query_result$records

  if (nrow(log_to_write) == 0) {
    log_to_write <- empty_search_log %>% slice(0)
  }

  if (nrow(records_to_write) == 0) {
    records_to_write <- empty_country_records %>% slice(0)
  }

  write_csv(log_to_write, log_path, na = "")
  write_csv(records_to_write, records_path, na = "")
}

search_log <- purrr::map_dfr(
  query_manifest$pathogen_log_path,
  ~read_checkpoint_csv(.x, empty_search_log)
)

country_records <- purrr::map_dfr(
  query_manifest$pathogen_records_path,
  ~read_checkpoint_csv(.x, empty_country_records)
) %>%
  mutate(across(where(is.character), clean_text))

if (nrow(search_log) == 0) {
  search_log <- empty_search_log
}

if (nrow(country_records) == 0) {
  country_records <- empty_country_records
  country_summary <- empty_country_summary
  disease_country_summary <- empty_disease_country_summary
  specificity_review <- empty_specificity_review
} else {
  country_summary <- country_records %>%
    filter(!is.na(country)) %>%
    group_by(Pathogens, Disease_name, Family, preferred_metadata_source, country) %>%
    summarise(
      accession_count = dplyr::n_distinct(
        dplyr::coalesce(biosample_accession, accession_version, primary_accession),
        na.rm = TRUE
      ),
      countries_raw = collapse_unique(dplyr::coalesce(country_raw, geo_loc_name_raw)),
      organisms = collapse_unique(organism),
      collection_dates = collapse_unique(collection_date),
      hosts = collapse_unique(host),
      .groups = "drop"
    ) %>%
    arrange(Pathogens, country)

  if (nrow(country_summary) == 0) {
    country_summary <- empty_country_summary
  }

  disease_country_summary <- country_records %>%
    filter(!is.na(Disease_name), !is.na(country)) %>%
    group_by(Disease_name, country) %>%
    summarise(
      pathogen_count = dplyr::n_distinct(Pathogens, na.rm = TRUE),
      record_count = dplyr::n_distinct(
        dplyr::coalesce(biosample_accession, accession_version, primary_accession),
        na.rm = TRUE
      ),
      pathogens = collapse_unique(Pathogens),
      source_dbs = collapse_unique(record_source_db),
      organisms = collapse_unique(organism),
      collection_dates = collapse_unique(collection_date),
      hosts = collapse_unique(host),
      .groups = "drop"
    ) %>%
    arrange(Disease_name, country)

  if (nrow(disease_country_summary) == 0) {
    disease_country_summary <- empty_disease_country_summary
  }

  specificity_review <- search_log %>%
    filter(query_specificity_status %in% c("too_broad_review", "review_needed")) %>%
    arrange(desc(records_found), Pathogens)

  if (nrow(specificity_review) == 0) {
    specificity_review <- empty_specificity_review
  }
}

second_pass_candidates <- search_log %>%
  filter(
    (note == "truncated_at_retmax" & last_batch_new_countries > 0) |
      (note == "stopped_after_country_plateau" & last_batch_new_countries > 0) |
      note %in% c("fetch_error", "search_failed", "search_id_batch_error")
  ) %>%
  mutate(
    followup_metadata_source = purrr::pmap_chr(
      list(preferred_metadata_source, records_found, countries_observed, note),
      ~recommend_followup_metadata_source(..1, ..2, ..3, ..4)
    ),
    followup_metadata_reason = purrr::pmap_chr(
      list(preferred_metadata_source, records_found, countries_observed, note),
      ~followup_metadata_source_reason(..1, ..2, ..3, ..4)
    ),
    recommended_retmax = dplyr::case_when(
      note %in% c("truncated_at_retmax", "stopped_after_country_plateau") ~ pmax(records_fetched * 2L, 2000L),
      TRUE ~ pmax(records_fetched, 1000L)
    )
  ) %>%
  arrange(desc(records_found), Pathogens)

write_csv(search_log, search_log_output_path, na = "")
write_csv(country_records, records_output_path, na = "")
write_csv(country_summary, summary_output_path, na = "")
write_csv(disease_country_summary, disease_summary_output_path, na = "")
write_csv(second_pass_candidates, second_pass_output_path, na = "")
write_csv(specificity_review, specificity_review_output_path, na = "")

if (!any(stringr::str_detect(search_log$status, "search_failed"), na.rm = TRUE)) {
  archive_stale_checkpoint_artifacts(
    pathogen_runs_dir = pathogen_runs_dir,
    pathogen_log_dir = pathogen_log_dir,
    pathogen_records_dir = pathogen_records_dir,
    active_log_paths = query_manifest$pathogen_log_path,
    active_record_paths = query_manifest$pathogen_records_path,
    verbose = verbose
  )
} else {
  log_progress(
    "Skipping checkpoint archive because one or more current manifest rows still have search_failed status.",
    verbose = verbose
  )
}

# ------------------------------------------------------------------------------
# Console summary --------------------------------------------------------------
# ------------------------------------------------------------------------------
cat("Pathogens in manifest:", nrow(query_manifest), "\n")
cat("Search log rows:", nrow(search_log), "\n")
cat("Accession-level records written:", nrow(country_records), "\n")
cat("Pathogen-country rows written:", nrow(country_summary), "\n")
cat("Disease-country rows written:", nrow(disease_country_summary), "\n")
cat("Manifest path:", manifest_output_path, "\n")
cat("Search log path:", search_log_output_path, "\n")
cat("Country records path:", records_output_path, "\n")
cat("Country summary path:", summary_output_path, "\n")
cat("Disease summary path:", disease_summary_output_path, "\n")
cat("Second-pass candidate path:", second_pass_output_path, "\n")
cat("Specificity review path:", specificity_review_output_path, "\n")
