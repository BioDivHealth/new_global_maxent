# -----------------------------------------------------------------------------|
# scripts/sdms/present/utils.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared helpers for present-day SDM fitting and calibration scripts.
# -----------------------------------------------------------------------------|

repo_root <- function() {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# -----------------------------------------------------------------------------|
# Command-line helpers ----
# -----------------------------------------------------------------------------|

parse_cli_args <- function(args) {
  parsed <- list()
  i <- 1

  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }

    key <- sub("^--", "", key)
    next_value <- if (i < length(args)) args[[i + 1]] else NULL

    if (is.null(next_value) || grepl("^--", next_value)) {
      parsed[[key]] <- TRUE
      i <- i + 1
    } else {
      parsed[[key]] <- next_value
      i <- i + 2
    }
  }

  parsed
}

get_arg <- function(parsed, key, default = NULL) {
  value <- parsed[[key]]

  if (is.null(value)) {
    return(default)
  }

  value
}

has_flag <- function(parsed, key) {
  isTRUE(parsed[[key]])
}

# -----------------------------------------------------------------------------|
# Text and scalar helpers ----
# -----------------------------------------------------------------------------|

safe_species_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(sort(unique(x)), collapse = "; ")
}

coalesce_scalar <- function(x, default = NA_character_) {
  if (length(x) == 0 || is.null(x) || all(is.na(x))) {
    return(default)
  }

  x[[1]]
}

split_arg <- function(value) {
  if (length(value) == 0 || is.null(value)) {
    return(character())
  }
  if (length(value) > 1) {
    return(value)
  }

  value <- as.character(value)
  if (!nzchar(value)) {
    return(character())
  }

  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

as_logical_arg <- function(value) {
  if (is.logical(value)) {
    return(value)
  }

  value <- as.character(value)
  if (length(value) == 0 || is.na(value[[1]]) || !nzchar(value[[1]])) {
    return(FALSE)
  }

  tolower(value[[1]]) %in% c("true", "t", "yes", "y", "1")
}

canonical_species_name <- function(x) {
  x <- trimws(as.character(x))
  vapply(strsplit(tolower(x), "\\s+"), function(parts) {
    if (length(parts) == 0 || is.na(parts[[1]]) || !nzchar(parts[[1]])) {
      return(NA_character_)
    }

    parts[[1]] <- paste0(toupper(substr(parts[[1]], 1, 1)), substr(parts[[1]], 2, nchar(parts[[1]])))
    paste(parts, collapse = " ")
  }, character(1))
}

# -----------------------------------------------------------------------------|
# Chikungunya target-manifest helpers ----
# -----------------------------------------------------------------------------|

select_sdm_targets <- function(target_manifest,
                               roles = "vector",
                               species_filter = character(),
                               include_not_needed = FALSE,
                               include_already_available = FALSE,
                               max_species = Inf) {
  required_cols <- c(
    "species_name",
    "species_role",
    "sdm_needed_for_disease",
    "run_priority",
    "run_status",
    "sdm_available"
  )
  missing_cols <- setdiff(required_cols, names(target_manifest))
  if (length(missing_cols) > 0) {
    stop("Target manifest is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  roles <- split_arg(roles)
  species_filter <- canonical_species_name(split_arg(species_filter))

  targets <- target_manifest
  if (!include_not_needed && length(species_filter) == 0) {
    targets <- targets[targets$sdm_needed_for_disease == "yes", , drop = FALSE]
  }

  if (length(roles) > 0 && !any(roles %in% c("all", "ALL"))) {
    targets <- targets[targets$species_role %in% roles, , drop = FALSE]
  }

  targets$species_name_canonical <- canonical_species_name(targets$species_name)
  if (length(species_filter) > 0) {
    targets <- targets[targets$species_name_canonical %in% species_filter, , drop = FALSE]
  }

  if (!include_already_available) {
    available <- tolower(as.character(targets$sdm_available)) %in% c("true", "yes", "1")
    already_available <- available | targets$run_status == "already_available"
    targets <- targets[!already_available, , drop = FALSE]
  }

  targets <- targets[
    order(targets$run_priority, targets$species_role, targets$species_name_canonical),
    ,
    drop = FALSE
  ]
  targets <- targets[!duplicated(targets$species_name_canonical), , drop = FALSE]

  max_species <- as.numeric(max_species)
  if (is.finite(max_species)) {
    targets <- head(targets, max_species)
  }

  targets
}

################################################################################
# GBIF download helpers
################################################################################

read_dotenv_value <- function(keys, env_path = file.path(repo_root(), ".env")) {
  if (!file.exists(env_path)) {
    return("")
  }

  lines <- readLines(env_path, warn = FALSE)
  for (key in keys) {
    pattern <- paste0("^[[:space:]]*", key, "[[:space:]]*=")
    hit <- grep(pattern, lines, value = TRUE)
    if (length(hit) == 0) {
      next
    }

    value <- sub(pattern, "", hit[[1]])
    value <- sub("[[:space:]]+#.*$", "", value)
    value <- trimws(value)
    value <- sub("^\"", "", value)
    value <- sub("\"$", "", value)
    value <- sub("^'", "", value)
    value <- sub("'$", "", value)
    if (nzchar(value)) {
      return(value)
    }
  }

  ""
}

gbif_credential <- function(keys) {
  for (key in keys) {
    value <- Sys.getenv(key, unset = "")
    if (nzchar(value)) {
      return(value)
    }
  }

  read_dotenv_value(keys)
}

gbif_credentials <- function() {
  gbif_user <- gbif_credential(c("GBIF_USER", "GBIF_USERNAME", "gbif_user", "gbif_username"))
  gbif_password <- gbif_credential(c("GBIF_PASSWORD", "GBIF_PWD", "gbif_password", "gbif_pwd"))
  gbif_email <- gbif_credential(c("GBIF_EMAIL", "gbif_email"))
  if (!nzchar(gbif_email) && grepl("@", gbif_user, fixed = TRUE)) {
    gbif_email <- gbif_user
  }

  if (!nzchar(gbif_user) || !nzchar(gbif_password) || !nzchar(gbif_email)) {
    stop(
      "GBIF downloads require GBIF_USER/GBIF_USERNAME/gbif_username, ",
      "GBIF_PASSWORD/gbif_password, and GBIF_EMAIL/gbif_email. ",
      "If gbif_username is an email address it will be reused as the email.",
      call. = FALSE
    )
  }

  list(user = gbif_user, password = gbif_password, email = gbif_email)
}

gbif_download_status_row <- function(download_key) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package `rgbif` is required.", call. = FALSE)
  }

  meta <- rgbif::occ_download_meta(download_key)
  data.frame(
    gbif_download_key = as.character(download_key),
    gbif_status = coalesce_scalar(meta$status),
    gbif_doi = coalesce_scalar(meta$doi),
    gbif_download_link = coalesce_scalar(meta$downloadLink),
    gbif_total_records = suppressWarnings(as.integer(coalesce_scalar(meta$totalRecords))),
    gbif_size = suppressWarnings(as.integer(coalesce_scalar(meta$size))),
    gbif_created = coalesce_scalar(meta$created),
    gbif_modified = coalesce_scalar(meta$modified),
    status_checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )
}

submit_gbif_occurrence_download <- function(species, start_year, end_year) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package `rgbif` is required.", call. = FALSE)
  }

  credentials <- gbif_credentials()
  backbone <- rgbif::name_backbone(name = species, rank = "species")
  taxon_key <- backbone$usageKey
  if (is.null(taxon_key) || is.na(taxon_key)) {
    stop("GBIF backbone lookup did not return a usageKey for ", species, call. = FALSE)
  }

  predicates <- list(
    rgbif::pred("taxonKey", taxon_key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    rgbif::pred_gte("year", start_year)
  )
  if (!is.na(end_year)) {
    predicates <- c(predicates, list(rgbif::pred_lte("year", end_year)))
  }

  download_key <- do.call(
    rgbif::occ_download,
    c(
      predicates,
      list(
        format = "SIMPLE_CSV",
        user = credentials$user,
        pwd = credentials$password,
        email = credentials$email
      )
    )
  )

  list(
    gbif_download_key = as.character(download_key),
    taxon_key = as.integer(taxon_key),
    matched_name = coalesce_scalar(backbone$canonicalName, default = species)
  )
}
