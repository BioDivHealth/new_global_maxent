# ------------------------------------------------------------------------------|
#      genbank_simple_helpers.R ------------------------------------------------
# ------------------------------------------------------------------------------|
# Small helper set for the simplified GenBank pathogen-country workflow.
# ------------------------------------------------------------------------------|

# ------------------------------------------------------------------------------|
#      Text, environment, and provenance helpers ------------------------------
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

parse_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  tolower(value) %in% c("1", "true", "yes", "y")
}

parse_env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))

  if (is.na(parsed)) {
    return(default)
  }

  parsed
}

read_dotenv_value <- function(path, key) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  lines <- readLines(path, warn = FALSE)
  match <- grep(paste0("^\\s*", key, "\\s*="), lines, value = TRUE)

  if (length(match) == 0) {
    return(NA_character_)
  }

  value <- sub(paste0("^\\s*", key, "\\s*=\\s*"), "", match[[1]])
  value <- sub("\\s+#.*$", "", value)
  value <- stringr::str_remove_all(value, "^['\"]|['\"]$")
  clean_text(value)
}

configure_entrez_key <- function(dotenv_path = ".env") {
  api_key <- clean_text(Sys.getenv("NCBI_API_KEY", unset = Sys.getenv("ENTREZ_KEY", unset = NA_character_)))

  if (is.na(api_key)) {
    api_key <- dplyr::coalesce(
      read_dotenv_value(dotenv_path, "NCBI_API_KEY"),
      read_dotenv_value(dotenv_path, "ENTREZ_KEY"),
      read_dotenv_value(dotenv_path, "ncbi_api_key"),
      read_dotenv_value(dotenv_path, "entrez_key")
    )
  }

  if (!is.na(api_key)) {
    rentrez::set_entrez_key(api_key)
    Sys.setenv(NCBI_API_KEY = api_key)
  }

  invisible(api_key)
}

redact_sensitive_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "([?&]api_key=)[^&'\"\\s]+", "\\1[REDACTED]")
  x <- stringr::str_replace_all(x, "(api_key=)[^&'\"\\s]+", "\\1[REDACTED]")
  x
}

collapse_unique <- function(x, sep = "; ") {
  x <- clean_text(x)
  x <- unique(x[!is.na(x)])

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = sep)
}

# ------------------------------------------------------------------------------|
#      GenBank-simple output path helpers -------------------------------------
# ------------------------------------------------------------------------------|
genbank_simple_qa_files <- c(
  "genbank_simple_readiness_manifest_qa.csv",
  "genbank_readiness_search_logs.csv",
  "genbank_readiness_qa_summary.csv",
  "genbank_readiness_target_qa.csv",
  "genbank_readiness_country_standardization_qa.csv",
  "genbank_search_logs.csv",
  "genbank_simple_qa_summary.csv",
  "genbank_simple_target_qa.csv",
  "genbank_country_standardization_qa.csv"
)

genbank_simple_intermediate_files <- c(
  "genbank_readiness_country_records.csv",
  "genbank_readiness_country_records_standardized.csv",
  "genbank_readiness_pathogen_country_summary.csv",
  "genbank_readiness_pathogen_country_summary_standardized.csv",
  "genbank_readiness_disease_country_summary.csv",
  "genbank_country_records.csv",
  "genbank_country_records_standardized.csv",
  "genbank_pathogen_country_summary.csv",
  "genbank_pathogen_country_summary_standardized.csv",
  "genbank_disease_country_summary.csv"
)

genbank_simple_manifest_files <- c(
  "genbank_simple_manifest.csv",
  "genbank_simple_readiness_manifest.csv",
  "excluded_targets.csv"
)

genbank_simple_manual_files <- c(
  "genbank_readiness_query_overrides.csv"
)

genbank_simple_evidence_files <- c(
  "genbank_readiness_disease_country_summary_standardized.csv",
  "genbank_disease_country_summary_standardized.csv"
)

genbank_simple_file_path <- function(output_dir, file_name, create_parent = FALSE) {
  directory <- dplyr::case_when(
    file_name %in% genbank_simple_qa_files &&
      exists("genbank_simple_qa_dir", inherits = TRUE) ~
      get("genbank_simple_qa_dir", inherits = TRUE),
    file_name %in% genbank_simple_intermediate_files &&
      exists("genbank_simple_intermediate_dir", inherits = TRUE) ~
      get("genbank_simple_intermediate_dir", inherits = TRUE),
    file_name %in% genbank_simple_manifest_files &&
      exists("genbank_simple_manifest_dir", inherits = TRUE) ~
      get("genbank_simple_manifest_dir", inherits = TRUE),
    file_name %in% genbank_simple_manual_files &&
      exists("genbank_simple_manual_dir", inherits = TRUE) ~
      get("genbank_simple_manual_dir", inherits = TRUE),
    file_name %in% genbank_simple_evidence_files &&
      exists("genbank_simple_evidence_dir", inherits = TRUE) ~
      get("genbank_simple_evidence_dir", inherits = TRUE),
    TRUE ~ NA_character_
  )

  if (is.na(directory)) {
    subdir <- dplyr::case_when(
      file_name %in% genbank_simple_qa_files ~ "qa",
      file_name %in% genbank_simple_intermediate_files ~ "intermediate",
      TRUE ~ NA_character_
    )

    directory <- if (is.na(subdir)) {
      output_dir
    } else {
      file.path(output_dir, subdir)
    }
  }

  path <- file.path(directory, file_name)

  if (create_parent) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }

  path
}

genbank_simple_legacy_file_candidates <- function(output_dir, file_name) {
  candidate_roots <- output_dir

  if (exists("genbank_simple_legacy_dir", inherits = TRUE)) {
    candidate_roots <- c(
      candidate_roots,
      get("genbank_simple_legacy_dir", inherits = TRUE)
    )
  }

  unique(c(
    file.path(candidate_roots, file_name),
    file.path(candidate_roots, "qa", file_name),
    file.path(candidate_roots, "intermediate", file_name)
  ))
}

genbank_simple_existing_dir <- function(preferred_dir, legacy_dir = NULL) {
  if (dir.exists(preferred_dir) || is.null(legacy_dir) || !dir.exists(legacy_dir)) {
    return(preferred_dir)
  }

  legacy_dir
}

genbank_simple_map_dir <- function(summary_kind, output_dir) {
  summary_kind <- clean_text(summary_kind)

  if (
    summary_kind == "readiness_combined" &&
      exists("genbank_simple_readiness_maps_dir", inherits = TRUE)
  ) {
    return(get("genbank_simple_readiness_maps_dir", inherits = TRUE))
  }

  if (
    summary_kind != "readiness_combined" &&
      exists("genbank_simple_standard_maps_dir", inherits = TRUE)
  ) {
    return(get("genbank_simple_standard_maps_dir", inherits = TRUE))
  }

  if (summary_kind == "readiness_combined") {
    file.path(output_dir, "maps_readiness")
  } else {
    file.path(output_dir, "maps")
  }
}

genbank_simple_existing_map_dir <- function(summary_kind, output_dir) {
  preferred_dir <- genbank_simple_map_dir(summary_kind, output_dir)
  legacy_subdir <- if (summary_kind == "readiness_combined") {
    "maps_readiness"
  } else {
    "maps"
  }

  legacy_dir <- if (exists("genbank_simple_legacy_dir", inherits = TRUE)) {
    file.path(get("genbank_simple_legacy_dir", inherits = TRUE), legacy_subdir)
  } else {
    file.path(output_dir, legacy_subdir)
  }

  genbank_simple_existing_dir(preferred_dir, legacy_dir)
}

genbank_simple_existing_file_path <- function(output_dir, file_name) {
  preferred_path <- genbank_simple_file_path(output_dir, file_name)
  legacy_paths <- genbank_simple_legacy_file_candidates(output_dir, file_name)
  existing_legacy_path <- legacy_paths[file.exists(legacy_paths)]

  if (file.exists(preferred_path) || length(existing_legacy_path) == 0) {
    return(preferred_path)
  }

  existing_legacy_path[[1]]
}

# ------------------------------------------------------------------------------|
#      Manifest target and query helpers --------------------------------------
# ------------------------------------------------------------------------------|
as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(dplyr::coalesce(x, FALSE))
  }

  x <- tolower(clean_text(x))
  dplyr::coalesce(x %in% c("true", "t", "1", "yes", "y"), FALSE)
}

sanitize_filename <- function(x) {
  x <- clean_text(x)
  x <- dplyr::coalesce(x, "missing")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_{2,}", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x <- clean_text(x)
  dplyr::coalesce(x, "missing")
}

make_target_id <- function(pathogen, disease_name) {
  paste(
    sanitize_filename(pathogen),
    sanitize_filename(disease_name),
    sep = "__"
  )
}

coronavirus_exclusion_pattern <- stringr::regex(
  paste(
    c(
      "sarbecovirus",
      "merbecovirus",
      "\\bsars\\b",
      "\\bmers\\b",
      "sars-like",
      "mers-like",
      "covid",
      "sars-cov",
      "sars-cov-2",
      "mers-cov"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

is_coronavirus_excluded <- function(pathogen, disease_name) {
  target_text <- paste(clean_text(pathogen), clean_text(disease_name), sep = " | ")
  stringr::str_detect(target_text, coronavirus_exclusion_pattern)
}

is_allowed_influenza_target <- function(pathogen) {
  pathogen <- clean_text(pathogen)
  pathogen %in% c(
    "Alphainfluenzavirus influenzae (H5N1)",
    "Alphainfluenzavirus influenzae (H7N9)"
  )
}

is_broad_or_unwanted_influenza <- function(pathogen) {
  pathogen <- clean_text(pathogen)

  stringr::str_detect(
    pathogen,
    stringr::regex("^Alphainfluenzavirus influenzae( \\(|$)|Influenza A virus|H\\d+Nx", ignore_case = TRUE)
  ) & !is_allowed_influenza_target(pathogen)
}

extract_influenza_subtype <- function(pathogen) {
  pathogen <- clean_text(pathogen)
  subtype <- stringr::str_match(pathogen, "\\((H\\d+N\\d+)\\)")[, 2]
  clean_text(subtype)
}

build_simple_query <- function(pathogen, tax_ids) {
  pathogen <- clean_text(pathogen)
  tax_ids <- clean_text(tax_ids)
  tax_ids <- unique(tax_ids[!is.na(tax_ids)])

  subtype <- extract_influenza_subtype(pathogen)

  if (length(tax_ids) > 0) {
    taxid_terms <- paste0("txid", tax_ids, "[Organism:exp]")
    taxid_query <- if (length(taxid_terms) == 1) {
      taxid_terms
    } else {
      paste0("(", paste(taxid_terms, collapse = " OR "), ")")
    }

    if (!is.na(subtype)) {
      return(paste0("(", taxid_query, ") AND \"", subtype, "\"[All Fields]"))
    }

    return(taxid_query)
  }

  if (!is.na(subtype)) {
    return(paste0("\"Influenza A virus\"[Organism] AND \"", subtype, "\"[All Fields]"))
  }

  paste0("\"", pathogen, "\"[Organism]")
}

# ------------------------------------------------------------------------------|
#      GenBank XML parsing helpers --------------------------------------------
# ------------------------------------------------------------------------------|
standardize_country_name <- function(country_raw, geo_loc_name_raw) {
  location_raw <- dplyr::coalesce(clean_text(country_raw), clean_text(geo_loc_name_raw))

  if (is.na(location_raw)) {
    return(NA_character_)
  }

  country <- stringr::str_split_fixed(location_raw, ":", n = 2)[, 1]
  country <- stringr::str_split_fixed(country, ",", n = 2)[, 1]
  country <- clean_text(country)

  dplyr::case_when(
    country %in% c("USA", "U.S.A.", "United States of America") ~ "United States",
    country %in% c("UK", "U.K.") ~ "United Kingdom",
    country == "Viet Nam" ~ "Vietnam",
    country == "Russian Federation" ~ "Russia",
    country == "Czech Republic" ~ "Czechia",
    TRUE ~ country
  )
}

extract_gbqual_value <- function(source_feature, qualifier_name) {
  if (length(source_feature) == 0 || inherits(source_feature, "xml_missing")) {
    return(NA_character_)
  }

  qualifier_paths <- c(
    paste0(
      ".//GBQualifier[GBQualifier_name='",
      qualifier_name,
      "']/GBQualifier_value"
    ),
    paste0(
      ".//INSDQualifier[INSDQualifier_name='",
      qualifier_name,
      "']/INSDQualifier_value"
    )
  )

  value <- xml2::xml_text(
    xml2::xml_find_first(source_feature, paste(qualifier_paths, collapse = " | "))
  )

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

extract_seq_value <- function(seq_node, field_name) {
  field_paths <- c(
    paste0("./GBSeq_", field_name),
    paste0("./INSDSeq_", field_name)
  )

  value <- xml2::xml_text(
    xml2::xml_find_first(seq_node, paste(field_paths, collapse = " | "))
  )

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

parse_nuccore_records <- function(xml_text, manifest_row) {
  if (is.na(xml_text) || !nzchar(xml_text)) {
    return(tibble::tibble())
  }

  doc <- xml2::read_xml(xml_text)
  seq_nodes <- xml2::xml_find_all(doc, ".//GBSeq | .//INSDSeq")

  if (length(seq_nodes) == 0) {
    return(tibble::tibble())
  }

  purrr::map_dfr(seq_nodes, function(seq_node) {
    source_feature <- xml2::xml_find_first(
      seq_node,
      ".//GBFeature[GBFeature_key='source'] | .//INSDFeature[INSDFeature_key='source']"
    )

    country_raw <- extract_gbqual_value(source_feature, "country")
    geo_loc_name_raw <- extract_gbqual_value(source_feature, "geo_loc_name")

    tibble::tibble(
      target_id = manifest_row$target_id,
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      PathogenTaxID = manifest_row$PathogenTaxID,
      query_used = manifest_row$query_used,
      source_db = manifest_row$source_db,
      accession_version = extract_seq_value(seq_node, "accession-version"),
      primary_accession = extract_seq_value(seq_node, "primary-accession"),
      definition = extract_seq_value(seq_node, "definition"),
      organism = extract_seq_value(seq_node, "organism"),
      taxonomy = extract_seq_value(seq_node, "taxonomy"),
      sequence_length = suppressWarnings(as.integer(extract_seq_value(seq_node, "length"))),
      country_raw = country_raw,
      geo_loc_name_raw = geo_loc_name_raw,
      country = standardize_country_name(country_raw, geo_loc_name_raw),
      lat_lon = extract_gbqual_value(source_feature, "lat_lon"),
      collection_date = extract_gbqual_value(source_feature, "collection_date"),
      host = extract_gbqual_value(source_feature, "host"),
      isolate = extract_gbqual_value(source_feature, "isolate"),
      strain = extract_gbqual_value(source_feature, "strain"),
      isolate_source = extract_gbqual_value(source_feature, "isolation_source"),
      db_xref = extract_gbqual_value(source_feature, "db_xref")
    )
  })
}

# ------------------------------------------------------------------------------|
#      NCBI E-utilities search and fetch helpers ------------------------------
# ------------------------------------------------------------------------------|
build_esearch_url <- function(query, retmax = 0L, retstart = 0L, db = "nuccore") {
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=", db,
    "&retmode=xml",
    "&usehistory=n",
    "&retmax=", as.integer(retmax),
    "&retstart=", as.integer(retstart),
    "&term=", utils::URLencode(query, reserved = FALSE)
  )

  api_key <- Sys.getenv("NCBI_API_KEY", unset = Sys.getenv("ENTREZ_KEY", unset = ""))

  if (nzchar(api_key)) {
    url <- paste0(url, "&api_key=", utils::URLencode(api_key, reserved = FALSE))
  }

  url
}

read_url_text <- function(url) {
  con <- base::url(url, open = "rb")
  on.exit(close(con), add = TRUE)
  paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

parse_esearch_response <- function(xml_text) {
  doc <- xml2::read_xml(xml_text)
  error_nodes <- xml2::xml_find_all(doc, ".//ERROR | .//Error")

  if (length(error_nodes) > 0) {
    stop("NCBI ESearch error: ", collapse_unique(xml2::xml_text(error_nodes)))
  }

  list(
    count = suppressWarnings(as.integer(clean_text(xml2::xml_text(xml2::xml_find_first(doc, ".//Count"))))),
    ids = clean_text(xml2::xml_text(xml2::xml_find_all(doc, ".//IdList/Id")))
  )
}

search_ids_page <- function(query, retmax, retstart, db = "nuccore") {
  entrez_result <- tryCatch(
    rentrez::entrez_search(
      db = db,
      term = query,
      retmax = retmax,
      retstart = retstart,
      use_history = FALSE
    ),
    error = function(e) e
  )

  if (!inherits(entrez_result, "error")) {
    return(list(
      count = suppressWarnings(as.integer(entrez_result$count)),
      ids = clean_text(entrez_result$ids),
      backend = paste0("rentrez:", db)
    ))
  }

  url <- build_esearch_url(query = query, retmax = retmax, retstart = retstart, db = db)
  http_result <- suppressWarnings(parse_esearch_response(read_url_text(url)))
  http_result$backend <- paste0("http:", db)
  http_result
}

fetch_nuccore_batch <- function(ids) {
  rentrez::entrez_fetch(
    db = "nuccore",
    id = ids,
    rettype = "gbc",
    retmode = "xml"
  )
}

strip_xml_declaration <- function(x) {
  x <- as.character(x)
  stringr::str_replace(x, "^\\s*<\\?xml[^>]*>\\s*", "")
}

combine_gbset_xml <- function(left, right) {
  left <- strip_xml_declaration(left)
  right <- strip_xml_declaration(right)

  left_inner <- stringr::str_replace(left, "^\\s*<GBSet>\\s*", "")
  left_inner <- stringr::str_replace(left_inner, "\\s*</GBSet>\\s*$", "")
  right_inner <- stringr::str_replace(right, "^\\s*<GBSet>\\s*", "")
  right_inner <- stringr::str_replace(right_inner, "\\s*</GBSet>\\s*$", "")

  paste0("<GBSet>", left_inner, right_inner, "</GBSet>")
}

fetch_nuccore_batch_with_retry <- function(ids, max_attempts = 3L, retry_wait_seconds = 1) {
  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      fetch_nuccore_batch(ids = ids),
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- result

    if (attempt < max_attempts) {
      Sys.sleep(retry_wait_seconds * attempt)
    }
  }

  if (length(ids) > 1) {
    split_index <- floor(length(ids) / 2)

    left <- fetch_nuccore_batch_with_retry(
      ids = ids[seq_len(split_index)],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )
    right <- fetch_nuccore_batch_with_retry(
      ids = ids[seq.int(split_index + 1L, length(ids))],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )

    if (!inherits(left, "error") && !inherits(right, "error")) {
      return(combine_gbset_xml(left, right))
    }
  }

  last_error
}
