# ------------------------------------------------------------------------------|
#      02_run_genbank_full_retrieval.R ----------------------------------------
# ------------------------------------------------------------------------------|
# Purpose: Retrieve all nuccore records for the approved GenBank-simple manifest
#          through deterministic pagination and per-target checkpoints.
# Inputs : genbank_simple_readiness_manifest.csv
# Outputs: pathogen_runs_readiness/search_logs/*.csv
#          pathogen_runs_readiness/country_records/*.csv
#
# Notes  : This script intentionally does not use interval/random sampling or
#          adaptive country-plateau stopping. Set `NCBI_API_KEY` or `ENTREZ_KEY`
#          for higher request limits. Useful optional environment variables:
#          `GENBANK_SIMPLE_TARGET_FILTER`, `GENBANK_SIMPLE_MAX_TARGETS`,
#          `GENBANK_SIMPLE_SEARCH_PAGE_SIZE`, `GENBANK_SIMPLE_FETCH_BATCH_SIZE`,
#          `GENBANK_SIMPLE_RESUME`, `GENBANK_SIMPLE_FORCE_RERUN`.
#          Readiness mode is the default. Set
#          `GENBANK_SIMPLE_MANIFEST_KIND=standard` to use the frozen 19-target
#          manifest and write checkpoints under `pathogen_runs/`.
#          In readiness mode, `GENBANK_SIMPLE_READINESS_ONLY_NEW` defaults to
#          FALSE so the expanded manifest is the authoritative target surface.
#          `GENBANK_SIMPLE_DRY_RUN=TRUE` validates target selection without
#          contacting NCBI.
#          `GENBANK_SIMPLE_MAX_RECORDS_FOUND` skips targets above a count-only
#          eSearch threshold before collecting IDs/fetching XML.
#          Readiness mode always skips the broad Salmonella target for GenBank.
#          `GENBANK_SIMPLE_MAX_RECORDS_PER_TARGET` is a smoke-test/debug cap;
#          leave it unset for full retrieval.
#          Large targets stream parsed rows to disk. Tune with
#          `GENBANK_SIMPLE_STREAM_THRESHOLD` and
#          `GENBANK_SIMPLE_RECORD_FLUSH_SIZE`.
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Load required libraries -------------------------------------------------
# ------------------------------------------------------------------------------|
library(pacman)
p_load(dplyr, here, purrr, readr, rentrez, stringr, tibble, xml2)

source(here("scripts", "associations", "genbank_simple", "genbank_simple_helpers.R"))
source(here("scripts", "associations", "working_inputs.R"))

configure_entrez_key(here(".env"))

# ------------------------------------------------------------------------------|
#      Resolve manifest and run directories -----------------------------------
# ------------------------------------------------------------------------------|
output_dir <- genbank_simple_dir
standard_manifest_path <- genbank_simple_existing_file_path(
  output_dir,
  "genbank_simple_manifest.csv"
)
readiness_manifest_path <- genbank_simple_existing_file_path(
  output_dir,
  "genbank_simple_readiness_manifest.csv"
)

manifest_kind <- Sys.getenv("GENBANK_SIMPLE_MANIFEST_KIND", unset = "readiness") %>%
  clean_text() %>%
  stringr::str_to_lower()

manifest_kind <- case_when(
  manifest_kind %in% c("standard", "simple", "current", "19_target") ~ "standard",
  manifest_kind %in% c("readiness", "expanded_readiness") ~ "readiness",
  TRUE ~ NA_character_
)

if (is.na(manifest_kind)) {
  stop(
    "GENBANK_SIMPLE_MANIFEST_KIND must be `standard` or `readiness`.",
    call. = FALSE
  )
}

manifest_path_override <- clean_text(Sys.getenv("GENBANK_SIMPLE_MANIFEST_PATH", unset = NA_character_))
manifest_path <- dplyr::coalesce(
  manifest_path_override,
  if_else(manifest_kind == "readiness", readiness_manifest_path, standard_manifest_path)
)

run_dir <- if_else(
  manifest_kind == "readiness",
  genbank_simple_readiness_run_dir,
  genbank_simple_standard_run_dir
)
search_log_dir <- file.path(run_dir, "search_logs")
country_record_dir <- file.path(run_dir, "country_records")

dir.create(search_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(country_record_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------|
#      Parse run controls ------------------------------------------------------
# ------------------------------------------------------------------------------|
search_page_size <- parse_env_integer("GENBANK_SIMPLE_SEARCH_PAGE_SIZE", default = 5000L)
fetch_batch_size <- parse_env_integer("GENBANK_SIMPLE_FETCH_BATCH_SIZE", default = 200L)
max_targets <- parse_env_integer("GENBANK_SIMPLE_MAX_TARGETS", default = NA_integer_)
max_records_per_target <- parse_env_integer("GENBANK_SIMPLE_MAX_RECORDS_PER_TARGET", default = NA_integer_)
max_records_found <- parse_env_integer("GENBANK_SIMPLE_MAX_RECORDS_FOUND", default = NA_integer_)
stream_threshold <- parse_env_integer("GENBANK_SIMPLE_STREAM_THRESHOLD", default = 100000L)
record_flush_size <- parse_env_integer("GENBANK_SIMPLE_RECORD_FLUSH_SIZE", default = 10000L)
target_filter <- Sys.getenv("GENBANK_SIMPLE_TARGET_FILTER", unset = "")
resume <- parse_env_flag("GENBANK_SIMPLE_RESUME", default = TRUE)
force_rerun <- parse_env_flag("GENBANK_SIMPLE_FORCE_RERUN", default = FALSE)
dry_run <- parse_env_flag("GENBANK_SIMPLE_DRY_RUN", default = FALSE)
readiness_only_new <- parse_env_flag(
  "GENBANK_SIMPLE_READINESS_ONLY_NEW",
  default = FALSE
)

# ------------------------------------------------------------------------------|
#      Manifest readers --------------------------------------------------------
# ------------------------------------------------------------------------------|
stop_if_missing <- function(data, cols, label) {
  missing_cols <- setdiff(cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(data)
}

read_standard_manifest <- function(path) {
  manifest <- read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text))

  stop_if_missing(
    manifest,
    c("target_id", "Pathogens", "Disease_name", "query_used", "source_db"),
    basename(path)
  )

  manifest %>%
    mutate(
      target_id = clean_text(target_id),
      Pathogens = clean_text(Pathogens),
      Disease_name = clean_text(Disease_name),
      PathogenTaxID = clean_text(PathogenTaxID),
      query_used = clean_text(query_used),
      source_db = dplyr::coalesce(clean_text(source_db), "nuccore"),
      current_target_id = NA_character_
    )
}

read_readiness_manifest <- function(path) {
  manifest <- read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text))

  stop_if_missing(
    manifest,
    c(
      "target_id",
      "query_pathogen_label",
      "readiness_disease_names",
      "pathogen_taxid",
      "query_used",
      "source_db",
      "manifest_status"
    ),
    basename(path)
  )

  if (!"current_target_id" %in% names(manifest)) {
    manifest <- manifest %>%
      mutate(current_target_id = NA_character_)
  }

  non_ready <- manifest %>%
    filter(manifest_status != "ready_for_future_retrieval")

  if (nrow(non_ready) > 0) {
    stop(
      "Readiness manifest contains non-ready rows: ",
      paste(non_ready$target_id, collapse = ", "),
      call. = FALSE
    )
  }

  manifest %>%
    filter(!readiness_only_new | is.na(current_target_id)) %>%
    transmute(
      target_id = clean_text(target_id),
      Pathogens = clean_text(query_pathogen_label),
      Disease_name = clean_text(readiness_disease_names),
      PathogenTaxID = clean_text(pathogen_taxid),
      query_used = clean_text(query_used),
      source_db = dplyr::coalesce(clean_text(source_db), "nuccore"),
      current_target_id = clean_text(current_target_id)
  )
}

# ------------------------------------------------------------------------------|
#      Select targets ----------------------------------------------------------
# ------------------------------------------------------------------------------|
manifest <- if (manifest_kind == "readiness") {
  read_readiness_manifest(manifest_path)
} else {
  read_standard_manifest(manifest_path)
}

if (manifest_kind == "readiness") {
  manifest <- manifest %>%
    filter(!str_detect(
      paste(Pathogens, Disease_name, target_id, sep = " | "),
      regex("Salmonella", ignore_case = TRUE)
    ))
}

if (nzchar(target_filter)) {
  manifest <- manifest %>%
    filter(str_detect(paste(Pathogens, Disease_name, target_id, sep = " | "), regex(target_filter, ignore_case = TRUE)))
}

if (!is.na(max_targets) && max_targets > 0) {
  manifest <- manifest %>% slice_head(n = max_targets)
}

if (nrow(manifest) == 0) {
  stop("No GenBank-simple retrieval targets selected.", call. = FALSE)
}

message("Manifest kind: ", manifest_kind)
message("Manifest path: ", manifest_path)
message("Run directory: ", run_dir)
message("Targets selected: ", nrow(manifest))

# ------------------------------------------------------------------------------|
#      Empty checkpoint schemas ------------------------------------------------
# ------------------------------------------------------------------------------|
empty_records <- tibble(
  target_id = character(),
  Pathogens = character(),
  Disease_name = character(),
  PathogenTaxID = character(),
  query_used = character(),
  source_db = character(),
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
  db_xref = character()
)

empty_log <- tibble(
  target_id = character(),
  Pathogens = character(),
  Disease_name = character(),
  source_db = character(),
  query_used = character(),
  status = character(),
  records_found = integer(),
  ids_collected = integer(),
  records_parsed = integer(),
  countries_observed = integer(),
  started_at = character(),
  finished_at = character(),
  note = character()
)

# ------------------------------------------------------------------------------|
#      Checkpoint helpers ------------------------------------------------------
# ------------------------------------------------------------------------------|
write_skip_checkpoint <- function(manifest_row, status, records_found = NA_integer_, note = NA_character_) {
  log_path <- file.path(search_log_dir, paste0(manifest_row$target_id, ".csv"))
  records_path <- file.path(country_record_dir, paste0(manifest_row$target_id, ".csv"))

  log_row <- empty_log %>%
    add_row(
      target_id = manifest_row$target_id,
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      source_db = manifest_row$source_db,
      query_used = manifest_row$query_used,
      status = status,
      records_found = records_found,
      ids_collected = 0L,
      records_parsed = 0L,
      countries_observed = 0L,
      started_at = as.character(Sys.time()),
      finished_at = as.character(Sys.time()),
      note = note
    )

  write_csv(log_row, log_path)

  if (!file.exists(records_path)) {
    write_csv(empty_records, records_path)
  }

  invisible(NULL)
}

# ------------------------------------------------------------------------------|
#      Retrieve one target -----------------------------------------------------
# ------------------------------------------------------------------------------|
run_target <- function(row_index) {
  manifest_row <- manifest[row_index, ]
  target_id <- manifest_row$target_id
  log_path <- file.path(search_log_dir, paste0(target_id, ".csv"))
  records_path <- file.path(country_record_dir, paste0(target_id, ".csv"))

  if (resume && !force_rerun && file.exists(log_path) && file.exists(records_path)) {
    message("Skipping existing checkpoint: ", target_id)
    return(invisible(NULL))
  }

  if (force_rerun) {
    unlink(c(log_path, records_path), force = TRUE)
  }

  started_at <- as.character(Sys.time())
  status <- "success"
  note <- NA_character_
  records_found <- NA_integer_
  ids <- character(0)
  records <- empty_records
  records_written <- 0L
  countries_seen <- character(0)
  record_buffer <- list()
  streaming_records <- FALSE
  page_size_for_target <- if (!is.na(max_records_per_target) && max_records_per_target > 0) {
    min(search_page_size, max_records_per_target)
  } else {
    search_page_size
  }

  message("Running target ", row_index, "/", nrow(manifest), ": ", manifest_row$Pathogens)

  # Large targets are written in chunks so the script does not need to hold every
  # parsed GenBank record in memory at once.
  flush_record_buffer <- function(force = FALSE) {
    if (length(record_buffer) == 0) {
      return(invisible(NULL))
    }

    buffered_records <- bind_rows(record_buffer)

    if (!force && nrow(buffered_records) < record_flush_size) {
      return(invisible(NULL))
    }

    write_csv(
      buffered_records,
      records_path,
      append = file.exists(records_path),
      col_names = !file.exists(records_path)
    )

    records_written <<- records_written + nrow(buffered_records)
    countries_seen <<- unique(c(countries_seen, clean_text(buffered_records$country)))
    record_buffer <<- list()

    message("  wrote checkpoint rows: ", records_written)
    invisible(NULL)
  }

  try_result <- tryCatch({
    first_page <- search_ids_page(
      query = manifest_row$query_used,
      retmax = page_size_for_target,
      retstart = 0L,
      db = manifest_row$source_db
    )

    records_found <- first_page$count

    if (!is.na(max_records_found) && records_found > max_records_found) {
      status <- "skipped_records_found_exceeds_limit"
      note <- paste0(
        "records_found_exceeds_limit:",
        records_found,
        ">",
        max_records_found
      )

      write_skip_checkpoint(
        manifest_row = manifest_row,
        status = status,
        records_found = records_found,
        note = note
      )

      message("  skipped: ", note)
      return(invisible(NULL))
    }

    ids <- first_page$ids[!is.na(first_page$ids)]

    target_id_count <- if (!is.na(max_records_per_target) && max_records_per_target > 0) {
      min(records_found, max_records_per_target, na.rm = TRUE)
    } else {
      records_found
    }

    if (!is.na(target_id_count) && target_id_count > length(ids)) {
      starts <- seq.int(from = page_size_for_target, to = target_id_count - 1L, by = page_size_for_target)

      for (retstart in starts) {
        page <- search_ids_page(
          query = manifest_row$query_used,
          retmax = min(page_size_for_target, target_id_count - retstart),
          retstart = retstart,
          db = manifest_row$source_db
        )
        ids <- unique(c(ids, page$ids[!is.na(page$ids)]))
        Sys.sleep(0.12)
      }
    }

    if (!is.na(max_records_per_target) && max_records_per_target > 0 && length(ids) > max_records_per_target) {
      ids <- utils::head(ids, max_records_per_target)
      note <- paste0("debug_record_cap_applied:", max_records_per_target)
    }

    if (length(ids) == 0) {
      status <- "no_records"
    } else {
      streaming_records <- !is.na(records_found) &&
        records_found >= stream_threshold &&
        is.na(max_records_per_target)

      if (streaming_records) {
        message(
          "  streaming enabled | records_found=",
          records_found,
          " | flush_size=",
          record_flush_size
        )
      }

      id_batches <- split(ids, ceiling(seq_along(ids) / fetch_batch_size))

      record_batches <- purrr::map(seq_along(id_batches), function(batch_index) {
        batch_ids <- id_batches[[batch_index]]
        batch_xml <- fetch_nuccore_batch_with_retry(batch_ids)

        if (inherits(batch_xml, "error")) {
          stop(conditionMessage(batch_xml))
        }

        Sys.sleep(0.12)
        parsed_records <- parse_nuccore_records(batch_xml, manifest_row)

        if (streaming_records) {
          record_buffer <<- c(record_buffer, list(parsed_records))
          flush_record_buffer(force = FALSE)

          if (batch_index %% 50 == 0) {
            message(
              "  fetched batch ",
              batch_index,
              "/",
              length(id_batches),
              " | checkpoint_rows=",
              records_written
            )
          }

          return(empty_records)
        }

        parsed_records
      })

      if (streaming_records) {
        flush_record_buffer(force = TRUE)
        records <- empty_records
      } else {
        records <- bind_rows(record_batches)
        countries_seen <- unique(clean_text(records$country))
        countries_seen <- countries_seen[!is.na(countries_seen)]
      }
    }
    invisible(NULL)
  }, error = function(e) {
    status <<- "search_or_fetch_failed"
    note <<- redact_sensitive_text(conditionMessage(e))
    message("  failed: ", note)
    invisible(NULL)
  })

  if (nrow(records) == 0) {
    records <- empty_records
  }

  records_parsed <- if (streaming_records || records_written > 0L) {
    records_written
  } else {
    nrow(records)
  }

  countries_observed <- if (streaming_records || length(countries_seen) > 0) {
    length(countries_seen[!is.na(countries_seen)])
  } else {
    n_distinct(records$country, na.rm = TRUE)
  }

  log_row <- empty_log %>%
    add_row(
      target_id = target_id,
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      source_db = manifest_row$source_db,
      query_used = manifest_row$query_used,
      status = status,
      records_found = records_found,
      ids_collected = length(ids),
      records_parsed = records_parsed,
      countries_observed = countries_observed,
      started_at = started_at,
      finished_at = as.character(Sys.time()),
      note = note
    )

  write_csv(log_row, log_path)

  if (!file.exists(records_path)) {
    write_csv(records, records_path)
  } else if (!streaming_records && records_written == 0L) {
    write_csv(records, records_path)
  }

  message(
    "  status=",
    status,
    " | records_found=",
    dplyr::coalesce(as.character(records_found), "NA"),
    " | ids_collected=",
    length(ids),
    " | records_parsed=",
    records_parsed
  )

  invisible(try_result)
}

# ------------------------------------------------------------------------------|
#      Execute retrieval -------------------------------------------------------
# ------------------------------------------------------------------------------|
if (dry_run) {
  message("Dry run requested; no NCBI search/fetch calls were made.")
  print(
    manifest %>%
      select(target_id, Pathogens, Disease_name, PathogenTaxID, source_db, query_used) %>%
      slice_head(n = 20),
    n = 20
  )
} else {
  purrr::walk(seq_len(nrow(manifest)), run_target)

  message("Finished GenBank-simple retrieval targets: ", nrow(manifest))
}
