#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 09_audit_gbif_synonyms.R ----
# -----------------------------------------------------------------------------|
# Purpose: Audit whether GBIF vector downloads likely cover older names and
#          synonyms for each target species.
#
# This is a diagnostic script only. It does not submit GBIF downloads, rewrite
# occurrence files, update request ledgers, or change model inputs.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# 1. RStudio config: edit this block before sourcing the script ----
# -----------------------------------------------------------------------------|

if (!exists("batch_config", inherits = FALSE)) {
  batch_config <- list(
    roles = "vector",
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y")),
    species_filter = character(),
    max_species = Inf,
    query_gbif_api = FALSE,
    query_gbif_counts = FALSE,
    gbif_api_timeout_seconds = 20,
    dry_run = FALSE
  )
}

# -----------------------------------------------------------------------------|
# 2. Internal defaults ----
# -----------------------------------------------------------------------------|

lacie_vector_root <- "/Volumes/LaCie/new_global_maxent/sdms/runs_artur/vector_sdm_push"

default_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "vector_species_sdm_targets.csv"),
  occurrence_root = file.path(lacie_vector_root, "occurrences"),
  request_manifest_path = file.path(lacie_vector_root, "gbif_download_requests.csv"),
  audit_run_root = file.path(lacie_vector_root, "gbif_synonym_audit_runs"),
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = TRUE,
  species_filter = character(),
  max_species = Inf,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  query_gbif_api = FALSE,
  query_gbif_counts = FALSE,
  gbif_api_timeout_seconds = 20,
  dry_run = FALSE
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# 3. Config and preflight ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
}

target_manifest_path <- config_arg("target-manifest-path")
occurrence_root <- config_arg("occurrence-root")
request_manifest_path <- config_arg("request-manifest-path")
audit_run_root <- config_arg("audit-run-root")
roles <- split_arg(config_arg("roles"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
species_filter <- split_arg(config_arg("species-filter"))
max_species <- as.numeric(config_arg("max-species"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
query_gbif_api <- as_logical_arg(config_arg("query-gbif-api"))
query_gbif_counts <- as_logical_arg(config_arg("query-gbif-counts"))
gbif_api_timeout_seconds <- as.numeric(config_arg("gbif-api-timeout-seconds"))
dry_run <- as_logical_arg(config_arg("dry-run")) || has_flag(args, "dry-run")

if (!file.exists(target_manifest_path)) {
  stop("Missing vector target manifest: ", target_manifest_path, call. = FALSE)
}

if (!dir.exists(occurrence_root)) {
  stop(
    "LaCie occurrence root not found: ", occurrence_root, "\n",
    "Mount the external drive or set `batch_config$occurrence_root` explicitly.",
    call. = FALSE
  )
}

if (query_gbif_api && !requireNamespace("curl", quietly = TRUE)) {
  stop("Package `curl` is required when `query_gbif_api = TRUE`.", call. = FALSE)
}

if (query_gbif_api && !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package `jsonlite` is required when `query_gbif_api = TRUE`.", call. = FALSE)
}

# -----------------------------------------------------------------------------|
# 4. Small helpers ----
# -----------------------------------------------------------------------------|

utc_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

read_csv_if_exists <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }

  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

value_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }

  as.character(x[[1]])
}

int_or_na <- function(x) {
  out <- suppressWarnings(as.integer(value_or_na(x)))
  if (length(out) == 0) {
    return(NA_integer_)
  }

  out
}

first_present_column <- function(data, candidates) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

column_or_na <- function(data, col) {
  if (is.na(col) || !col %in% names(data)) {
    return(rep(NA_character_, nrow(data)))
  }

  as.character(data[[col]])
}

gbif_api_get <- function(path, query = list()) {
  query <- query[!vapply(query, function(x) length(x) == 0 || is.na(x), logical(1))]
  if (length(query) > 0) {
    query_values <- vapply(query, function(x) {
      if (is.logical(x)) {
        x <- ifelse(x, "true", "false")
      }
      utils::URLencode(as.character(x), reserved = TRUE)
    }, character(1))
    query_string <- paste(names(query_values), query_values, sep = "=", collapse = "&")
  } else {
    query_string <- ""
  }

  url <- paste0("https://api.gbif.org/v1/", path)
  if (nzchar(query_string)) {
    url <- paste0(url, "?", query_string)
  }

  handle <- curl::new_handle(
    timeout = gbif_api_timeout_seconds,
    connecttimeout = gbif_api_timeout_seconds
  )
  response <- curl::curl_fetch_memory(url, handle = handle)
  jsonlite::fromJSON(rawToChar(response$content), flatten = TRUE)
}

gbif_species_match <- function(name) {
  gbif_api_get("species/match", list(name = name, rank = "SPECIES"))
}

gbif_species_synonyms <- function(primary_key) {
  out <- gbif_api_get(paste0("species/", primary_key, "/synonyms"), list(limit = 1000))
  if (!"results" %in% names(out) || is.null(out$results) || nrow(out$results) == 0) {
    return(data.frame())
  }

  as.data.frame(out$results, stringsAsFactors = FALSE)
}

gbif_occurrence_count_api <- function(candidate_name, year_arg) {
  out <- gbif_api_get(
    "occurrence/search",
    list(
      scientificName = candidate_name,
      hasCoordinate = TRUE,
      hasGeospatialIssue = FALSE,
      occurrenceStatus = "PRESENT",
      year = year_arg,
      limit = 0
    )
  )
  if (!"count" %in% names(out)) {
    return(NA_integer_)
  }

  suppressWarnings(as.integer(out$count))
}

collapse_unique <- function(x, max_values = 50) {
  x <- unique(as.character(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0) {
    return(NA_character_)
  }

  shown <- head(sort(x), max_values)
  suffix <- if (length(x) > max_values) {
    paste0(" ... (+", length(x) - max_values, " more)")
  } else {
    ""
  }
  paste0(paste(shown, collapse = "; "), suffix)
}

normalize_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

split_candidate_values <- function(x) {
  x <- as.character(x)
  x <- x[nonempty(x)]
  if (length(x) == 0) {
    return(character())
  }

  pieces <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  pieces <- trimws(pieces)
  unique(pieces[nonempty(pieces)])
}

gbif_source_path <- function(species) {
  species_safe <- safe_species_name(species)
  species_dir <- file.path(occurrence_root, species_safe, "gbif-download")
  candidates <- c(
    file.path(species_dir, "deduplicated", paste0(species_safe, "_gbif_key_deduplicated.csv")),
    file.path(species_dir, "cleaned", paste0(species_safe, "_cleaned.csv")),
    file.path(species_dir, "raw", paste0(species, ".csv"))
  )

  hits <- candidates[file.exists(candidates) & !startsWith(basename(candidates), "._")]
  if (length(hits) > 0) {
    return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
  }

  raw_dir <- file.path(species_dir, "raw")
  if (!dir.exists(raw_dir)) {
    return(NA_character_)
  }

  raw_hits <- list.files(raw_dir, pattern = "[.]csv$", full.names = TRUE)
  raw_hits <- raw_hits[
    !startsWith(basename(raw_hits), "._") &
      basename(raw_hits) != "raw_download_manifest.csv"
  ]
  if (length(raw_hits) == 0) {
    return(NA_character_)
  }

  normalizePath(raw_hits[[1]], winslash = "/", mustWork = TRUE)
}

source_file_type <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    return("none")
  }
  if (grepl("/deduplicated/", path, fixed = TRUE)) {
    return("deduplicated")
  }
  if (grepl("/cleaned/", path, fixed = TRUE)) {
    return("cleaned")
  }
  if (grepl("/raw/", path, fixed = TRUE)) {
    return("raw")
  }

  "unknown"
}

empty_local_name_summary <- function(species, source_path) {
  data.frame(
    species_name = species,
    source_path = source_path,
    scientificName = NA_character_,
    verbatimScientificName = NA_character_,
    gbif_species = NA_character_,
    taxonKey = NA_integer_,
    speciesKey = NA_integer_,
    taxonRank = NA_character_,
    rows = 0L,
    present_rows = 0L,
    individual_count_zero_rows = 0L,
    unique_gbif_ids = 0L,
    unique_coordinate_rows = 0L,
    stringsAsFactors = FALSE
  )
}

summarize_local_names <- function(data, species, source_path) {
  if (is.null(data) || nrow(data) == 0) {
    return(empty_local_name_summary(species, source_path))
  }

  dt <- data.table::as.data.table(data)
  required <- c(
    "scientificName",
    "verbatimScientificName",
    "species",
    "taxonKey",
    "speciesKey",
    "taxonRank",
    "occurrenceStatus",
    "individualCount",
    "gbifID",
    "key",
    "decimalLongitude",
    "decimalLatitude"
  )
  for (col in setdiff(required, names(dt))) {
    dt[, (col) := NA]
  }

  dt[, present_record := toupper(trimws(as.character(occurrenceStatus))) == "PRESENT"]
  dt[is.na(present_record), present_record := FALSE]
  dt[, zero_count_record := suppressWarnings(as.numeric(individualCount)) == 0]
  dt[is.na(zero_count_record), zero_count_record := FALSE]
  dt[, record_id := data.table::fifelse(nonempty(gbifID), as.character(gbifID), as.character(key))]
  dt[, coordinate_key := data.table::fifelse(
    !is.na(decimalLongitude) & !is.na(decimalLatitude),
    paste(decimalLongitude, decimalLatitude, sep = "|"),
    NA_character_
  )]

  out <- dt[
    ,
    .(
      rows = .N,
      present_rows = sum(present_record),
      individual_count_zero_rows = sum(zero_count_record),
      unique_gbif_ids = data.table::uniqueN(record_id[nonempty(record_id)]),
      unique_coordinate_rows = data.table::uniqueN(coordinate_key[nonempty(coordinate_key)])
    ),
    by = .(scientificName, verbatimScientificName, gbif_species = species, taxonKey, speciesKey, taxonRank)
  ]
  data.table::setorder(out, -rows, scientificName, verbatimScientificName)

  out[, species_name := species]
  out[, source_path := source_path]
  out <- out[
    ,
    .(
      species_name,
      source_path,
      scientificName,
      verbatimScientificName,
      gbif_species,
      taxonKey,
      speciesKey,
      taxonRank,
      rows,
      present_rows,
      individual_count_zero_rows,
      unique_gbif_ids,
      unique_coordinate_rows
    )
  ]
  as.data.frame(out)
}

local_name_values <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(character())
  }

  values <- c(
    column_or_na(data, "scientificName"),
    column_or_na(data, "verbatimScientificName"),
    column_or_na(data, "species")
  )
  normalize_name(values[nonempty(values)])
}

local_taxon_keys <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(integer())
  }

  keys <- unique(c(
    suppressWarnings(as.integer(column_or_na(data, "taxonKey"))),
    suppressWarnings(as.integer(column_or_na(data, "speciesKey")))
  ))
  keys[!is.na(keys)]
}

species_request_row <- function(requests, species) {
  if (is.null(requests) || nrow(requests) == 0) {
    return(data.frame())
  }

  requests$species_key <- canonical_species_name(if ("species_name_canonical" %in% names(requests)) {
    requests$species_name_canonical
  } else {
    requests$species_name
  })
  species_key <- canonical_species_name(species)

  hits <- requests[
    requests$species_key == species_key &
      requests$occurrence_method == "gbif-download",
    ,
    drop = FALSE
  ]
  if (nrow(hits) == 0) {
    return(data.frame())
  }

  hits$start_year_int <- suppressWarnings(as.integer(hits$start_year))
  hits$end_year_int <- suppressWarnings(as.integer(hits$end_year))
  exact <- hits$start_year_int == start_year & hits$end_year_int == end_year
  if (any(exact, na.rm = TRUE)) {
    hits <- hits[exact, , drop = FALSE]
  }

  has_key <- !is.na(hits$taxon_key) & nzchar(as.character(hits$taxon_key))
  has_download <- !is.na(hits$gbif_download_key) & nzchar(as.character(hits$gbif_download_key))
  hits <- hits[order(!has_key, !has_download, -hits$end_year_int, -hits$start_year_int), , drop = FALSE]
  hits[1, , drop = FALSE]
}

gbif_backbone_row <- function(name) {
  if (!query_gbif_api) {
    return(data.frame())
  }

  backbone <- tryCatch(
    gbif_species_match(name),
    error = function(err) {
      structure(list(error = conditionMessage(err)), class = "gbif_audit_error")
    }
  )
  if (inherits(backbone, "gbif_audit_error")) {
    return(data.frame(
      queried_name = name,
      usageKey = NA_integer_,
      acceptedUsageKey = NA_integer_,
      scientificName = NA_character_,
      canonicalName = NA_character_,
      rank = NA_character_,
      status = NA_character_,
      matchType = NA_character_,
      confidence = NA_real_,
      note = NA_character_,
      error = backbone$error,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    queried_name = name,
    usageKey = int_or_na(backbone$usageKey),
    acceptedUsageKey = int_or_na(backbone$acceptedUsageKey),
    scientificName = value_or_na(backbone$scientificName),
    canonicalName = value_or_na(backbone$canonicalName),
    rank = value_or_na(backbone$rank),
    status = value_or_na(backbone$status),
    matchType = value_or_na(backbone$matchType),
    confidence = suppressWarnings(as.numeric(value_or_na(backbone$confidence))),
    note = value_or_na(backbone$note),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

accepted_key_from_backbone <- function(row) {
  if (nrow(row) == 0) {
    return(NA_integer_)
  }
  accepted <- suppressWarnings(as.integer(row$acceptedUsageKey[[1]]))
  usage <- suppressWarnings(as.integer(row$usageKey[[1]]))
  if (!is.na(accepted)) {
    return(accepted)
  }
  usage
}

gbif_synonym_names <- function(primary_key) {
  if (!query_gbif_api || is.na(primary_key)) {
    return(data.frame())
  }

  synonyms <- tryCatch(
    gbif_species_synonyms(primary_key),
    error = function(err) data.frame(error = conditionMessage(err))
  )
  if (nrow(synonyms) == 0) {
    return(data.frame())
  }

  if ("error" %in% names(synonyms)) {
    return(data.frame(
      candidate_name = NA_character_,
      candidate_source = "gbif_synonyms",
      source_usage_key = NA_integer_,
      source_taxonomic_status = NA_character_,
      source_rank = NA_character_,
      source_error = synonyms$error[[1]],
      stringsAsFactors = FALSE
    ))
  }

  name_col <- first_present_column(synonyms, c("canonicalName", "scientificName", "species"))
  if (is.na(name_col)) {
    return(data.frame())
  }

  data.frame(
    candidate_name = column_or_na(synonyms, name_col),
    candidate_source = "gbif_synonyms",
    source_usage_key = suppressWarnings(as.integer(column_or_na(synonyms, first_present_column(synonyms, c("key", "usageKey"))))),
    source_taxonomic_status = column_or_na(synonyms, first_present_column(synonyms, c("taxonomicStatus", "status"))),
    source_rank = column_or_na(synonyms, first_present_column(synonyms, c("rank", "taxonRank"))),
    source_error = NA_character_,
    stringsAsFactors = FALSE
  )
}

spatial_spp_candidate_names <- function(species) {
  taxonomy_dir <- file.path(occurrence_root, safe_species_name(species), "spatial-spp", "taxonomy")
  if (!dir.exists(taxonomy_dir)) {
    return(character())
  }

  paths <- list.files(taxonomy_dir, pattern = "[.]csv$", full.names = TRUE)
  paths <- paths[!startsWith(basename(paths), "._")]
  if (length(paths) == 0) {
    return(character())
  }

  values <- character()
  for (path in paths) {
    taxonomy <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(err) NULL)
    if (is.null(taxonomy) || nrow(taxonomy) == 0) {
      next
    }

    candidate_cols <- intersect(
      c("Spp_syn", "IUCN_syn", "ITIS_syn", "Synonym_names_used", "GBIF_backbone_names_used"),
      names(taxonomy)
    )
    for (col in candidate_cols) {
      values <- c(values, split_candidate_values(taxonomy[[col]]))
    }
  }

  unique(values[nonempty(values)])
}

gbif_occurrence_count <- function(candidate_name) {
  if (!query_gbif_api || !query_gbif_counts || is.na(candidate_name) || !nzchar(candidate_name)) {
    return(NA_integer_)
  }

  year_arg <- if (is.na(end_year)) {
    paste0(start_year, ",")
  } else {
    paste0(start_year, ",", end_year)
  }

  out <- tryCatch(
    gbif_occurrence_count_api(candidate_name, year_arg),
    error = function(err) NA_integer_
  )

  suppressWarnings(as.integer(out))
}

candidate_table <- function(species, primary_key, local_names) {
  candidates <- data.frame(
    candidate_name = species,
    candidate_source = "requested_name",
    source_usage_key = NA_integer_,
    source_taxonomic_status = NA_character_,
    source_rank = NA_character_,
    source_error = NA_character_,
    stringsAsFactors = FALSE
  )
  spatial_candidates <- spatial_spp_candidate_names(species)
  spatial_candidates <- if (length(spatial_candidates) == 0) {
    data.frame()
  } else {
    data.frame(
      candidate_name = spatial_candidates,
      candidate_source = "spatial_spp_taxonomy",
      source_usage_key = NA_integer_,
      source_taxonomic_status = NA_character_,
      source_rank = NA_character_,
      source_error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  candidates <- data.table::rbindlist(
    list(
      candidates,
      gbif_synonym_names(primary_key),
      spatial_candidates
    ),
    use.names = TRUE,
    fill = TRUE
  )
  candidates <- as.data.frame(candidates)
  candidates <- candidates[nonempty(candidates$candidate_name), , drop = FALSE]
  candidates$canonical_candidate_key <- normalize_name(candidates$candidate_name)
  candidates <- candidates[!duplicated(paste(candidates$canonical_candidate_key, candidates$candidate_source)), , drop = FALSE]

  if (nrow(candidates) == 0) {
    return(data.frame())
  }

  if (query_gbif_api) {
    backbone_rows <- lapply(candidates$candidate_name, gbif_backbone_row)
    backbone <- data.table::rbindlist(backbone_rows, use.names = TRUE, fill = TRUE)
  } else {
    backbone <- data.frame(
      queried_name = candidates$candidate_name,
      usageKey = NA_integer_,
      acceptedUsageKey = NA_integer_,
      scientificName = NA_character_,
      canonicalName = NA_character_,
      rank = NA_character_,
      status = NA_character_,
      matchType = NA_character_,
      confidence = NA_real_,
      note = NA_character_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  names(backbone) <- paste0("candidate_backbone_", names(backbone))
  candidates <- cbind(candidates, as.data.frame(backbone))

  candidate_accepted <- suppressWarnings(as.integer(candidates$candidate_backbone_acceptedUsageKey))
  candidates$candidate_usage_key <- suppressWarnings(as.integer(candidates$candidate_backbone_usageKey))
  candidates$candidate_accepted_key <- ifelse(
    is.na(candidate_accepted),
    candidates$candidate_usage_key,
    candidate_accepted
  )
  candidates$covered_by_primary_taxon_key <- !is.na(primary_key) & (
    candidates$candidate_usage_key == primary_key |
      candidates$candidate_accepted_key == primary_key
  )
  candidates$name_found_in_local_records <- vapply(
    candidates$candidate_name,
    function(name) {
      any(grepl(normalize_name(name), local_names, fixed = TRUE))
    },
    logical(1)
  )
  candidates$gbif_present_coordinate_count <- vapply(candidates$candidate_name, gbif_occurrence_count, integer(1))
  candidates$candidate_classification <- ifelse(
    is.na(candidates$candidate_usage_key),
    "unmatched",
    ifelse(
      candidates$covered_by_primary_taxon_key,
      "covered_by_primary_taxon_key",
      "possibly_different_taxon_concept"
    )
  )

  uncertain_match <- !is.na(candidates$candidate_backbone_matchType) &
    !candidates$candidate_backbone_matchType %in% c("EXACT", "NONE")
  low_confidence <- !is.na(candidates$candidate_backbone_confidence) &
    candidates$candidate_backbone_confidence < 90
  non_species_rank <- !is.na(candidates$candidate_backbone_rank) &
    toupper(candidates$candidate_backbone_rank) != "SPECIES"
  candidates$needs_review <- candidates$candidate_classification != "covered_by_primary_taxon_key" |
    uncertain_match |
    low_confidence |
    non_species_rank
  candidates$species_name <- species
  candidates$primary_taxon_key <- primary_key

  candidates[
    ,
    c(
      "species_name",
      "primary_taxon_key",
      "candidate_name",
      "candidate_source",
      "candidate_classification",
      "covered_by_primary_taxon_key",
      "name_found_in_local_records",
      "gbif_present_coordinate_count",
      "candidate_usage_key",
      "candidate_accepted_key",
      "candidate_backbone_canonicalName",
      "candidate_backbone_scientificName",
      "candidate_backbone_rank",
      "candidate_backbone_status",
      "candidate_backbone_matchType",
      "candidate_backbone_confidence",
      "source_usage_key",
      "source_taxonomic_status",
      "source_rank",
      "source_error",
      "needs_review"
    ),
    drop = FALSE
  ]
}

flag_row <- function(species, flag, detail, severity = "review") {
  data.frame(
    species_name = species,
    flag = flag,
    severity = severity,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

process_species <- function(species, requests) {
  source_path <- gbif_source_path(species)
  request_row <- species_request_row(requests, species)
  requested_taxon_key <- if (nrow(request_row) == 0) NA_integer_ else suppressWarnings(as.integer(request_row$taxon_key[[1]]))
  requested_matched_name <- if (nrow(request_row) == 0) NA_character_ else value_or_na(request_row$gbif_matched_name)

  data <- read_csv_if_exists(source_path)
  local_names_summary <- summarize_local_names(data, species, source_path)
  local_names <- local_name_values(data)
  local_keys <- local_taxon_keys(data)

  backbone <- gbif_backbone_row(species)
  backbone_key <- accepted_key_from_backbone(backbone)
  primary_key <- if (!is.na(requested_taxon_key)) requested_taxon_key else backbone_key
  candidates <- candidate_table(species, primary_key, local_names)

  flags <- list()
  if (is.na(source_path) || is.null(data) || nrow(data) == 0) {
    flags <- c(flags, list(flag_row(species, "no_local_gbif_download", "No local gbif-download source file found.", "warning")))
  }

  if (nrow(backbone) > 0) {
    backbone_match <- value_or_na(backbone$matchType)
    backbone_rank <- value_or_na(backbone$rank)
    backbone_confidence <- suppressWarnings(as.numeric(value_or_na(backbone$confidence)))
    review_backbone <- (!is.na(backbone_match) && !backbone_match %in% c("EXACT", "NONE")) ||
      (!is.na(backbone_rank) && toupper(backbone_rank) != "SPECIES") ||
      (!is.na(backbone_confidence) && backbone_confidence < 90) ||
      (!is.na(requested_taxon_key) && !is.na(backbone_key) && requested_taxon_key != backbone_key)
    if (review_backbone) {
      flags <- c(flags, list(flag_row(
        species,
        "review_backbone_match",
        paste0(
          "Backbone matchType=", backbone_match,
          "; rank=", backbone_rank,
          "; confidence=", backbone_confidence,
          "; ledger_taxon_key=", requested_taxon_key,
          "; backbone_accepted_key=", backbone_key
        ),
        "review"
      )))
    }
  }

  non_species_rows <- local_names_summary[
    !is.na(local_names_summary$taxonRank) &
      nzchar(local_names_summary$taxonRank) &
      !toupper(local_names_summary$taxonRank) %in% c("SPECIES", "SUBSPECIES"),
    ,
    drop = FALSE
  ]
  bold_rows <- local_names_summary[
    grepl("^BOLD:", local_names_summary$scientificName) |
      grepl("^BOLD:", local_names_summary$verbatimScientificName),
    ,
    drop = FALSE
  ]
  unexpected_key_rows <- local_names_summary[
    !is.na(primary_key) &
      ((!is.na(local_names_summary$taxonKey) & local_names_summary$taxonKey != primary_key) |
        (!is.na(local_names_summary$speciesKey) & local_names_summary$speciesKey != primary_key)),
    ,
    drop = FALSE
  ]
  if (nrow(non_species_rows) > 0 || nrow(bold_rows) > 0 || nrow(unexpected_key_rows) > 0) {
    flags <- c(flags, list(flag_row(
      species,
      "possible_non_species_records",
      paste0(
        "non_species_name_groups=", nrow(non_species_rows),
        "; bold_name_groups=", nrow(bold_rows),
        "; unexpected_key_groups=", nrow(unexpected_key_rows)
      ),
      "review"
    )))
  }

  candidate_count_signal <- !is.na(candidates$gbif_present_coordinate_count) &
    candidates$gbif_present_coordinate_count > 0
  missing_candidates <- candidates[
    candidates$candidate_source != "requested_name" &
      candidates$candidate_classification != "covered_by_primary_taxon_key" &
      (candidate_count_signal | is.na(candidates$gbif_present_coordinate_count)),
    ,
    drop = FALSE
  ]
  if (nrow(missing_candidates) > 0) {
    flags <- c(flags, list(flag_row(
      species,
      "review_possible_missing_synonym_records",
      paste0(
        "Candidate synonym names needing review: ",
        collapse_unique(missing_candidates$candidate_name, max_values = 20)
      ),
      "review"
    )))
  }

  covered_synonyms <- candidates[
    candidates$candidate_source != "requested_name" &
      candidates$covered_by_primary_taxon_key,
    ,
    drop = FALSE
  ]
  if (nrow(covered_synonyms) > 0 && nrow(missing_candidates) == 0) {
    flags <- c(flags, list(flag_row(
      species,
      "ok_primary_taxon_key_likely_covers_synonyms",
      paste0("Covered candidate names: ", collapse_unique(covered_synonyms$candidate_name, max_values = 20)),
      "ok"
    )))
  }

  if (length(flags) == 0) {
    flags <- list(flag_row(species, "ok_no_synonym_issues_detected", "No review flags produced by this audit.", "ok"))
  }

  status <- if (any(vapply(flags, function(x) x$severity[[1]] == "warning", logical(1)))) {
    "warning"
  } else if (any(vapply(flags, function(x) x$severity[[1]] == "review", logical(1)))) {
    "review"
  } else {
    "ok"
  }

  total_rows <- if (is.null(data)) 0L else nrow(data)
  present_rows <- if (is.null(data) || !"occurrenceStatus" %in% names(data)) {
    NA_integer_
  } else {
    sum(toupper(trimws(as.character(data$occurrenceStatus))) == "PRESENT", na.rm = TRUE)
  }
  zero_count_rows <- if (is.null(data) || !"individualCount" %in% names(data)) {
    NA_integer_
  } else {
    sum(suppressWarnings(as.numeric(data$individualCount)) == 0, na.rm = TRUE)
  }

  record_id_col <- if (is.null(data)) NA_character_ else first_present_column(data, c("gbifID", "key"))
  unique_ids <- if (is.null(data) || is.na(record_id_col)) {
    NA_integer_
  } else {
    length(unique(column_or_na(data, record_id_col)[nonempty(column_or_na(data, record_id_col))]))
  }
  unique_coords <- if (
    is.null(data) ||
      !all(c("decimalLongitude", "decimalLatitude") %in% names(data))
  ) {
    NA_integer_
  } else {
    coords <- paste(data$decimalLongitude, data$decimalLatitude, sep = "|")
    length(unique(coords[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude)]))
  }

  species_summary <- data.frame(
    species_name = species,
    audit_status = status,
    source_path = source_path,
    source_file_type = source_file_type(source_path),
    requested_start_year = start_year,
    requested_end_year = end_year,
    gbif_download_rows = total_rows,
    present_rows = present_rows,
    individual_count_zero_rows = zero_count_rows,
    unique_gbif_ids = unique_ids,
    unique_coordinate_rows = unique_coords,
    local_scientific_name_groups = nrow(local_names_summary[local_names_summary$rows > 0, , drop = FALSE]),
    local_taxon_keys = collapse_unique(as.character(local_keys)),
    ledger_taxon_key = requested_taxon_key,
    ledger_matched_name = requested_matched_name,
    backbone_usage_key = if (nrow(backbone) == 0) NA_integer_ else suppressWarnings(as.integer(backbone$usageKey[[1]])),
    backbone_accepted_key = backbone_key,
    backbone_canonical_name = if (nrow(backbone) == 0) NA_character_ else value_or_na(backbone$canonicalName),
    backbone_scientific_name = if (nrow(backbone) == 0) NA_character_ else value_or_na(backbone$scientificName),
    backbone_rank = if (nrow(backbone) == 0) NA_character_ else value_or_na(backbone$rank),
    backbone_status = if (nrow(backbone) == 0) NA_character_ else value_or_na(backbone$status),
    backbone_match_type = if (nrow(backbone) == 0) NA_character_ else value_or_na(backbone$matchType),
    backbone_confidence = if (nrow(backbone) == 0) NA_real_ else suppressWarnings(as.numeric(value_or_na(backbone$confidence))),
    candidate_name_count = nrow(candidates),
    candidate_names_needing_review = nrow(candidates[candidates$needs_review, , drop = FALSE]),
    flags = paste(vapply(flags, function(x) x$flag[[1]], character(1)), collapse = "; "),
    audited_at = utc_now(),
    stringsAsFactors = FALSE
  )

  list(
    species_summary = species_summary,
    candidate_names = candidates,
    local_record_names = local_names_summary,
    flags = data.table::rbindlist(flags, use.names = TRUE, fill = TRUE)
  )
}

# -----------------------------------------------------------------------------|
# 5. Select targets and run audit ----
# -----------------------------------------------------------------------------|

target_manifest <- read.csv(target_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
targets <- select_sdm_targets(
  target_manifest = target_manifest,
  roles = roles,
  species_filter = species_filter,
  include_not_needed = include_not_needed,
  include_already_available = include_already_available,
  max_species = max_species
)

requests <- if (file.exists(request_manifest_path)) {
  read.csv(request_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  warning("GBIF request manifest not found: ", request_manifest_path, call. = FALSE)
  data.frame()
}

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
run_dir <- file.path(audit_run_root, timestamp)

cat("Selected vector species:", nrow(targets), "\n")
cat("Occurrence root:", occurrence_root, "\n")
if (dry_run) {
  cat("Dry run: audit CSVs will not be written.\n")
}
if (!query_gbif_api) {
  cat("GBIF API lookups disabled; local files will still be summarized.\n")
} else if (!query_gbif_counts) {
  cat("GBIF occurrence-count lookups disabled; taxonomy lookups will still run.\n")
}
if (query_gbif_api) {
  cat("GBIF API timeout per call:", gbif_api_timeout_seconds, "seconds\n")
}

if (nrow(targets) == 0) {
  warning("No target species selected.", call. = FALSE)
  audit_results <- list()
} else {
  audit_results <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    species <- targets$species_name_canonical[[i]]
    audit_results[[i]] <- tryCatch(
      process_species(species, requests),
      error = function(err) {
        list(
          species_summary = data.frame(
            species_name = species,
            audit_status = "failed_error",
            source_path = NA_character_,
            source_file_type = "none",
            requested_start_year = start_year,
            requested_end_year = end_year,
            gbif_download_rows = NA_integer_,
            present_rows = NA_integer_,
            individual_count_zero_rows = NA_integer_,
            unique_gbif_ids = NA_integer_,
            unique_coordinate_rows = NA_integer_,
            local_scientific_name_groups = NA_integer_,
            local_taxon_keys = NA_character_,
            ledger_taxon_key = NA_integer_,
            ledger_matched_name = NA_character_,
            backbone_usage_key = NA_integer_,
            backbone_accepted_key = NA_integer_,
            backbone_canonical_name = NA_character_,
            backbone_scientific_name = NA_character_,
            backbone_rank = NA_character_,
            backbone_status = NA_character_,
            backbone_match_type = NA_character_,
            backbone_confidence = NA_real_,
            candidate_name_count = NA_integer_,
            candidate_names_needing_review = NA_integer_,
            flags = "failed_error",
            audited_at = utc_now(),
            notes = conditionMessage(err),
            stringsAsFactors = FALSE
          ),
          candidate_names = data.frame(),
          local_record_names = data.frame(),
          flags = flag_row(species, "failed_error", conditionMessage(err), "warning")
        )
      }
    )

    cat(
      "[",
      i,
      "/",
      nrow(targets),
      "] ",
      species,
      ": ",
      audit_results[[i]]$species_summary$audit_status[[1]],
      " (",
      audit_results[[i]]$species_summary$flags[[1]],
      ")\n",
      sep = ""
    )
  }
}

species_summary <- data.table::rbindlist(lapply(audit_results, `[[`, "species_summary"), use.names = TRUE, fill = TRUE)
candidate_names <- data.table::rbindlist(lapply(audit_results, `[[`, "candidate_names"), use.names = TRUE, fill = TRUE)
local_record_names <- data.table::rbindlist(lapply(audit_results, `[[`, "local_record_names"), use.names = TRUE, fill = TRUE)
flags <- data.table::rbindlist(lapply(audit_results, `[[`, "flags"), use.names = TRUE, fill = TRUE)

if (!dry_run) {
  ensure_dir(run_dir)
  data.table::fwrite(species_summary, file.path(run_dir, "gbif_synonym_audit_species_summary.csv"), na = "")
  data.table::fwrite(candidate_names, file.path(run_dir, "gbif_synonym_audit_candidate_names.csv"), na = "")
  data.table::fwrite(local_record_names, file.path(run_dir, "gbif_synonym_audit_existing_record_names.csv"), na = "")
  data.table::fwrite(flags, file.path(run_dir, "gbif_synonym_audit_flags.csv"), na = "")
  cat("Wrote GBIF synonym audit outputs:", run_dir, "\n")
}

cat("Audit complete.\n")
