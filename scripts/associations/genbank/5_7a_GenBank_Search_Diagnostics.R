# ------------------------------------------------------------------------------
# 5_7a_GenBank_Search_Diagnostics.R
# ------------------------------------------------------------------------------
# Purpose: Diagnose which GenBank/NCBI search queries succeed or fail for WHO
#          pathogens before attempting accession-level metadata fetching.
#
# Inputs : who_pathogen_analysis_units_keep.csv
#          who_pathogens_virion_taxid.csv
#          who_bacteria_clover_taxid.csv
# Outputs: genbank_search_diagnostic_manifest.csv
#          genbank_search_diagnostic_results.csv
#          genbank_search_diagnostic_summary.csv
#
# Notes  : Optional environment variables:
#          `GENBANK_PATHOGEN_FILTER`, `GENBANK_MAX_PATHOGENS`,
#          `GENBANK_DIAGNOSTIC_RETMAX`, `GENBANK_VERBOSE`.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, rentrez, stringr, tibble, xml2)

source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------
# Shared helpers ---------------------------------------------------------------
# ------------------------------------------------------------------------------
source(here("scripts", "associations", "genbank", "genbank_metadata_helpers.R"))

normalize_who_pathogen_input <- function(who_pathogens) {
  if ("analysis_unit" %in% names(who_pathogens)) {
    who_pathogens %>%
      transmute(
        Family = family,
        Pathogens = analysis_unit,
        previous_name = source_previous_name,
        msl39_viral_name = source_msl39_viral_name,
        Disease_name = source_disease_name,
        in_gibb_etal = in_gibb_etal,
        in_empres_i = in_empres_i
      )
  } else {
    who_pathogens %>%
      mutate(
        in_gibb_etal = if ("in_gibb_etal" %in% names(.)) in_gibb_etal else FALSE,
        in_empres_i = if ("in_empres_i" %in% names(.)) in_empres_i else FALSE
      )
  }
}

run_esearch_diagnostic_http <- function(query, retmax = 20L, db = "nucleotide") {
  request_url <- build_esearch_url(query = query, retmax = retmax, db = db)

  response <- read_url_text(request_url)
  doc <- xml2::read_xml(response)
  error_nodes <- xml2::xml_find_all(doc, ".//ERROR | .//Error")
  error_text <- collapse_unique(xml2::xml_text(error_nodes))
  count_text <- xml2::xml_text(xml2::xml_find_first(doc, ".//Count"))
  ids <- clean_text(xml2::xml_text(xml2::xml_find_all(doc, ".//IdList/Id")))
  ids <- ids[!is.na(ids)]

  tibble(
    ok = is.na(error_text),
    records_found = suppressWarnings(as.integer(clean_text(count_text))),
    ids_returned = length(ids),
    sample_ids = collapse_unique(utils::head(ids, 10)),
    error_message = error_text,
    request_url = request_url,
    search_backend = paste0("http:", db)
  )
}

run_taxonomy_link_diagnostic <- function(taxid_query) {
  tax_ids <- extract_tax_ids_from_query(taxid_query)

  if (length(tax_ids) == 0) {
    return(tibble(
      ok = FALSE,
      records_found = NA_integer_,
      ids_returned = NA_integer_,
      sample_ids = NA_character_,
      error_message = "taxonomy_link_missing_taxids",
      request_url = NA_character_,
      search_backend = "taxonomy_link"
    ))
  }

  result <- tryCatch(
    {
      all_ids <- character(0)

      for (tax_id in tax_ids) {
        link_result <- rentrez::entrez_link(
          dbfrom = "taxonomy",
          db = "nuccore",
          id = tax_id
        )

        link_ids <- clean_text(link_result$links$taxonomy_nuccore)
        link_ids <- link_ids[!is.na(link_ids)]
        all_ids <- unique(c(all_ids, link_ids))
      }

      tibble(
        ok = length(all_ids) > 0,
        records_found = length(all_ids),
        ids_returned = length(all_ids),
        sample_ids = collapse_unique(utils::head(all_ids, 10)),
        error_message = if (length(all_ids) > 0) NA_character_ else "taxonomy_link_zero_records",
        request_url = NA_character_,
        search_backend = "taxonomy_link"
      )
    },
    error = function(e) {
      tibble(
        ok = FALSE,
        records_found = NA_integer_,
        ids_returned = NA_integer_,
        sample_ids = NA_character_,
        error_message = conditionMessage(e),
        request_url = NA_character_,
        search_backend = "taxonomy_link"
      )
    }
  )

  result
}

run_esearch_diagnostic <- function(query, retmax = 20L) {
  query <- clean_text(query)

  if (is.na(query)) {
    return(tibble(
      ok = FALSE,
      records_found = NA_integer_,
      ids_returned = NA_integer_,
      sample_ids = NA_character_,
      error_message = "query_missing",
      request_url = NA_character_,
      search_backend = NA_character_
    ))
  }

  attempts <- list(
    list(method = "http", db = "nucleotide"),
    list(method = "http", db = "nuccore")
  )
  last_error <- "NCBI esearch failed for all attempted backends."
  last_request_url <- NA_character_

  for (query_variant in search_query_variants(query)) {
    for (attempt in attempts) {
      parsed <- tryCatch(
        run_esearch_diagnostic_http(
          query = query_variant,
          retmax = retmax,
          db = attempt$db
        ),
        error = function(e) e
      )

      if (!inherits(parsed, "error") && isTRUE(parsed$ok[[1]])) {
        return(parsed)
      }

      if (!inherits(parsed, "error")) {
        last_error <- parsed$error_message[[1]]
        last_request_url <- parsed$request_url[[1]]
      } else {
        last_error <- conditionMessage(parsed)
      }
    }
  }

  tibble(
    ok = FALSE,
    records_found = NA_integer_,
    ids_returned = NA_integer_,
    sample_ids = NA_character_,
    error_message = last_error,
    request_url = last_request_url,
    search_backend = NA_character_
  )
}

# ------------------------------------------------------------------------------
# Parameters and paths ---------------------------------------------------------
# ------------------------------------------------------------------------------
max_pathogens <- parse_env_integer("GENBANK_MAX_PATHOGENS", default = NA_integer_)
diagnostic_retmax <- parse_env_integer("GENBANK_DIAGNOSTIC_RETMAX", default = 20L)
verbose <- parse_env_flag("GENBANK_VERBOSE", default = TRUE)
pathogen_filter <- clean_text(Sys.getenv("GENBANK_PATHOGEN_FILTER", unset = NA_character_))

if (diagnostic_retmax <= 0) {
  stop("GENBANK_DIAGNOSTIC_RETMAX must be greater than 0.")
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

who_dir <- here("pathogen_association_data", "WHO")
genbank_dir <- file.path(who_dir, "genbank")
dir.create(genbank_dir, recursive = TRUE, showWarnings = FALSE)

who_path <- who_working_pathogens_path()
virion_path <- file.path(who_dir, "virion", "who_pathogens_virion_taxid.csv")
clover_path <- file.path(who_dir, "clover", "who_bacteria_clover_taxid.csv")

diagnostic_manifest_path <- file.path(genbank_dir, "genbank_search_diagnostic_manifest.csv")
diagnostic_results_path <- file.path(genbank_dir, "genbank_search_diagnostic_results.csv")
diagnostic_summary_path <- file.path(genbank_dir, "genbank_search_diagnostic_summary.csv")

# ------------------------------------------------------------------------------
# Load and harmonize WHO pathogen inputs ---------------------------------------
# ------------------------------------------------------------------------------
who_pathogens <- read_csv(
  who_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  normalize_who_pathogen_input()

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
    previous_name = collapse_unique(previous_name),
    msl39_viral_name = collapse_unique(msl39_viral_name),
    in_gibb_etal = any(in_gibb_etal, na.rm = TRUE),
    in_empres_i = any(in_empres_i, na.rm = TRUE),
    .groups = "drop"
  )

virion_manifest <- virion_taxids %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    virion_tax_ids = collapse_unique(VirusTaxID),
    virion_names = collapse_unique(Virion_VirusName),
    .groups = "drop"
  )

clover_manifest <- clover_taxids %>%
  group_by(Pathogens, Disease_name) %>%
  summarise(
    clover_tax_ids = collapse_unique(PathogenTaxID),
    clover_names = collapse_unique(Clover_PathogenName),
    .groups = "drop"
  )

query_manifest <- who_manifest %>%
  left_join(virion_manifest, by = c("Pathogens", "Disease_name")) %>%
  left_join(clover_manifest, by = c("Pathogens", "Disease_name")) %>%
  mutate(
    manual_query_aliases = purrr::map2_chr(
      Pathogens,
      Disease_name,
      ~collapse_unique(manual_query_aliases(.x, .y))
    ),
    taxid_query = purrr::pmap_chr(
      list(virion_tax_ids, clover_tax_ids),
      ~build_taxid_query(..1, ..2)
    ),
    name_query = purrr::pmap_chr(
      list(Pathogens, previous_name, msl39_viral_name, virion_names, clover_names, manual_query_aliases),
      ~build_name_query(..1, ..2, ..3, ..4, ..5, ..6)
    ),
    preferred_metadata_source = purrr::pmap_chr(
      list(Pathogens, Family, virion_tax_ids, clover_tax_ids, msl39_viral_name),
      ~choose_metadata_source(..1, ..2, ..3, ..4, ..5)
    ),
    metadata_source_reason = purrr::pmap_chr(
      list(Pathogens, Family, virion_tax_ids, clover_tax_ids, msl39_viral_name),
      ~metadata_source_reason(..1, ..2, ..3, ..4, ..5)
    ),
    taxid_geo_query = purrr::map_chr(taxid_query, build_geo_query),
    name_geo_query = purrr::map_chr(name_query, build_geo_query),
    query_strategy = dplyr::case_when(
      purrr::map_lgl(Pathogens, prefer_name_query_pathogen) & !is.na(name_query) ~ "name_subtype",
      !is.na(taxid_query) ~ "taxid",
      !is.na(name_query) ~ "name",
      TRUE ~ "missing"
    )
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

write_csv(query_manifest, diagnostic_manifest_path, na = "")

if (nrow(query_manifest) == 0) {
  write_csv(tibble(), diagnostic_results_path, na = "")
  write_csv(tibble(), diagnostic_summary_path, na = "")
  cat("No pathogens matched the current filters; no NCBI requests were made.\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# Run diagnostics --------------------------------------------------------------
# ------------------------------------------------------------------------------
diagnostic_rows <- split(query_manifest, seq_len(nrow(query_manifest)))

diagnostic_results <- purrr::imap_dfr(
  diagnostic_rows,
  function(manifest_row, idx) {
    manifest_row <- tibble::as_tibble(manifest_row)
    progress_prefix <- paste0("[", idx, "/", length(diagnostic_rows), "] ")

    log_progress(
      progress_prefix,
      "Checking ",
      manifest_row$Pathogens,
      if (!is.na(manifest_row$Disease_name)) paste0(" | disease=", manifest_row$Disease_name) else "",
      if (!is.na(manifest_row$query_strategy)) paste0(" | strategy=", manifest_row$query_strategy) else "",
      verbose = verbose
    )

    manual_queries <- manual_query_candidates(
      pathogen = manifest_row$Pathogens,
      disease_name = manifest_row$Disease_name
    )
    manual_attempts <- if (length(manual_queries) == 0) {
      tibble(
        query_label = character(),
        query_text = character()
      )
    } else {
      tibble(
        query_label = paste0("manual_query_", seq_len(length(manual_queries))),
        query_text = manual_queries
      )
    }
    taxonomy_attempts <- if (is.na(manifest_row$taxid_query)) {
      tibble(
        query_label = character(),
        query_text = character()
      )
    } else {
      tibble(
        query_label = "taxonomy_link",
        query_text = manifest_row$taxid_query
      )
    }

    query_attempts <- bind_rows(
      manual_attempts,
      taxonomy_attempts,
      tibble(
        query_label = c("taxid_geo_query", "taxid_query", "name_geo_query", "name_query"),
        query_text = c(
          manifest_row$taxid_geo_query,
          manifest_row$taxid_query,
          manifest_row$name_geo_query,
          manifest_row$name_query
        )
      )
    ) %>%
      mutate(query_text = clean_text(query_text)) %>%
      filter(!is.na(query_text)) %>%
      distinct(query_text, .keep_all = TRUE)

    purrr::pmap_dfr(
      query_attempts,
      function(query_label, query_text) {
        result <- if (identical(query_label, "taxonomy_link")) {
          run_taxonomy_link_diagnostic(
            taxid_query = query_text
          )
        } else {
          run_esearch_diagnostic(
            query = query_text,
            retmax = diagnostic_retmax
          )
        }

        log_progress(
          progress_prefix,
          manifest_row$Pathogens,
          " | ",
          query_label,
          " | ok=",
          result$ok,
          " | records_found=",
          dplyr::coalesce(result$records_found, NA_integer_),
          if (!is.na(result$error_message)) paste0(" | error=", result$error_message) else "",
          verbose = verbose
        )

        result %>%
          mutate(
            Pathogens = manifest_row$Pathogens,
            Disease_name = manifest_row$Disease_name,
            Family = manifest_row$Family,
            preferred_metadata_source = manifest_row$preferred_metadata_source,
            metadata_source_reason = manifest_row$metadata_source_reason,
            query_strategy = manifest_row$query_strategy,
            query_label = query_label,
            query_text = query_text,
            .before = 1
          )
      }
    )
  }
)

diagnostic_summary <- diagnostic_results %>%
  group_by(
    Pathogens,
    Disease_name,
    Family,
    preferred_metadata_source,
    metadata_source_reason,
    query_strategy
  ) %>%
  summarise(
    any_success = any(ok, na.rm = TRUE),
    successful_queries = collapse_unique(query_label[ok %in% TRUE]),
    failed_queries = collapse_unique(query_label[ok %in% FALSE]),
    best_query_label = {
      successful_rows <- which(ok %in% TRUE)
      if (length(successful_rows) == 0) {
        NA_character_
      } else {
        query_label[successful_rows[which.max(dplyr::coalesce(records_found[successful_rows], -1L))]]
      }
    },
    best_records_found = suppressWarnings(max(records_found[ok %in% TRUE], na.rm = TRUE)),
    representative_error = collapse_unique(error_message[ok %in% FALSE]),
    .groups = "drop"
  ) %>%
  mutate(
    best_records_found = dplyr::if_else(
      is.infinite(best_records_found),
      NA_integer_,
      as.integer(best_records_found)
    )
  ) %>%
  arrange(desc(any_success), Pathogens)

write_csv(diagnostic_results, diagnostic_results_path, na = "")
write_csv(diagnostic_summary, diagnostic_summary_path, na = "")

cat("Pathogens checked:", nrow(query_manifest), "\n")
cat("Diagnostic result rows:", nrow(diagnostic_results), "\n")
cat("Diagnostic manifest path:", diagnostic_manifest_path, "\n")
cat("Diagnostic results path:", diagnostic_results_path, "\n")
cat("Diagnostic summary path:", diagnostic_summary_path, "\n")
