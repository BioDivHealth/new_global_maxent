# ------------------------------------------------------------------------------
# 1_4_NCBI_Broad_Taxa_Candidate_Metadata.R
# ------------------------------------------------------------------------------
# Purpose: Resolve accession.version values for the broad-taxa candidate strain
#          table, pull NCBI virus genome metadata via the local Datasets CLI,
#          and save analysis-ready metadata for downstream R work.
#
# Input  : pathogen_association_data/WHO/who_diseases/
#          who_broad_taxa_candidate_strains.csv
# Output : pathogen_association_data/WHO/who_diseases/
#          who_broad_taxa_candidate_strains_ncbi_metadata.csv
#          who_broad_taxa_candidate_strains_ncbi_enriched.csv
#          who_broad_taxa_candidate_strains_ncbi_raw.jsonl
# -------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, jsonlite, purrr, readr, stringr, tibble)

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_names <- function(x) {
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- clean_text(x)
  x <- unique(stats::na.omit(x))
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = "; ")
}

collapse_lineage <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  lineage_names <- purrr::map_chr(x, ~ clean_text(.x$name %||% NA_character_))
  lineage_names <- unique(stats::na.omit(lineage_names))

  if (length(lineage_names) == 0) {
    return(NA_character_)
  }

  paste(lineage_names, collapse = "; ")
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

who_dir <- here("pathogen_association_data", "WHO", "who_diseases")
candidate_path <- file.path(who_dir, "who_broad_taxa_candidate_strains.csv")
metadata_output_path <- file.path(
  who_dir,
  "who_broad_taxa_candidate_strains_ncbi_metadata.csv"
)
enriched_output_path <- file.path(
  who_dir,
  "who_broad_taxa_candidate_strains_ncbi_enriched.csv"
)
raw_output_path <- file.path(
  who_dir,
  "who_broad_taxa_candidate_strains_ncbi_raw.jsonl"
)

default_ncbi_bin_dir_candidates <- c(
  here::here("ncbi"),
  file.path(here::here(), "..", "ncbi")
)
default_ncbi_bin_dir <- default_ncbi_bin_dir_candidates[
  file.exists(file.path(default_ncbi_bin_dir_candidates, "datasets"))
][1]

if (is.na(default_ncbi_bin_dir) || length(default_ncbi_bin_dir) == 0) {
  default_ncbi_bin_dir <- default_ncbi_bin_dir_candidates[1]
}

default_ncbi_bin_dir <- normalizePath(
  default_ncbi_bin_dir,
  winslash = "/",
  mustWork = FALSE
)

ncbi_bin_dir <- Sys.getenv("NCBI_BIN_DIR", unset = default_ncbi_bin_dir)
datasets_bin <- file.path(ncbi_bin_dir, "datasets")

if (!file.exists(candidate_path)) {
  stop("Candidate strain table not found: ", candidate_path)
}

if (!file.exists(datasets_bin)) {
  stop(
    "NCBI datasets binary not found at ", datasets_bin,
    ". Set NCBI_BIN_DIR if your ncbi folder lives elsewhere."
  )
}

candidate_strains <- read_csv(
  candidate_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

accession_bases <- candidate_strains %>%
  transmute(accession_base = stringr::str_replace(accession, "\\.[0-9]+$", "")) %>%
  distinct() %>%
  filter(!is.na(accession_base))

fetch_accession_summary <- function(accession_base, max_version = 8, timeout_sec = 120) {
  accession_candidates <- c(
    accession_base,
    paste0(accession_base, ".", seq_len(max_version))
  ) %>%
    unique()

  for (accession_try in accession_candidates) {
    stdout_file <- tempfile(pattern = "ncbi_datasets_stdout_", fileext = ".txt")
    stderr_file <- tempfile(pattern = "ncbi_datasets_stderr_", fileext = ".txt")

    status <- tryCatch(
      system2(
        datasets_bin,
        args = c(
          "summary", "virus", "genome", "accession",
          accession_try,
          "--as-json-lines"
        ),
        stdout = stdout_file,
        stderr = stderr_file,
        timeout = timeout_sec
      ),
      error = function(e) e
    )

    stdout_text <- if (file.exists(stdout_file)) {
      paste(readLines(stdout_file, warn = FALSE), collapse = "\n")
    } else {
      ""
    }

    stderr_text <- if (file.exists(stderr_file)) {
      paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
    } else {
      ""
    }

    unlink(c(stdout_file, stderr_file), force = TRUE)

    if (inherits(status, "error")) {
      next
    }

    out <- paste(c(stdout_text, stderr_text), collapse = "\n")

    if (is.na(out) || out == "") {
      next
    }

    json_start <- regexpr("\\{", out)

    if (json_start[[1]] > 0) {
      out_json <- substr(out, json_start[[1]], nchar(out))
    } else {
      out_json <- out
    }

    out_json <- stringr::str_trim(out_json)

    parsed <- tryCatch(
      jsonlite::fromJSON(out_json, simplifyVector = FALSE),
      error = function(e) NULL
    )

    if (is.null(parsed)) {
      next
    }

    if ((parsed$total_count %||% 0) < 1) {
      next
    }

    report <- parsed$reports[[1]]

    return(list(
      accession_base = accession_base,
      accession_version = clean_text(report$accession %||% accession_try),
      raw_json = out_json,
      report = report
    ))
  }

  list(
    accession_base = accession_base,
    accession_version = NA_character_,
    raw_json = NA_character_,
    report = NULL
  )
}

extract_report_row <- function(x) {
  report <- x$report

  if (is.null(report)) {
    return(tibble(
      accession_base = x$accession_base,
      accession_version = NA_character_,
      ncbi_lookup_status = "not_found",
      completeness = NA_character_,
      is_annotated = NA,
      length = NA_real_,
      protein_count = NA_real_,
      source_database = NA_character_,
      release_date = NA_character_,
      update_date = NA_character_,
      isolate_name = NA_character_,
      geographic_location = NA_character_,
      geographic_region = NA_character_,
      virus_name_ncbi = NA_character_,
      virus_tax_id = NA_real_,
      virus_lineage = NA_character_,
      host_name_ncbi = NA_character_,
      host_tax_id = NA_real_,
      host_lineage = NA_character_,
      submitter_affiliation = NA_character_,
      submitter_country = NA_character_,
      submitter_names = NA_character_,
      raw_json = NA_character_
    ))
  }

  tibble(
    accession_base = x$accession_base,
    accession_version = clean_text(report$accession %||% NA_character_),
    ncbi_lookup_status = "ok",
    completeness = clean_text(report$completeness %||% NA_character_),
    is_annotated = report$is_annotated %||% NA,
    length = as.numeric(report$length %||% NA_real_),
    protein_count = as.numeric(report$protein_count %||% NA_real_),
    source_database = clean_text(report$source_database %||% NA_character_),
    release_date = clean_text(report$release_date %||% NA_character_),
    update_date = clean_text(report$update_date %||% NA_character_),
    isolate_name = clean_text((report$isolate %||% list())$name %||% NA_character_),
    geographic_location = clean_text((report$location %||% list())$geographic_location %||% NA_character_),
    geographic_region = clean_text((report$location %||% list())$geographic_region %||% NA_character_),
    virus_name_ncbi = clean_text((report$virus %||% list())$organism_name %||% NA_character_),
    virus_tax_id = as.numeric((report$virus %||% list())$tax_id %||% NA_real_),
    virus_lineage = collapse_lineage((report$virus %||% list())$lineage %||% list()),
    host_name_ncbi = clean_text((report$host %||% list())$organism_name %||% NA_character_),
    host_tax_id = as.numeric((report$host %||% list())$tax_id %||% NA_real_),
    host_lineage = collapse_lineage((report$host %||% list())$lineage %||% list()),
    submitter_affiliation = clean_text((report$submitter %||% list())$affiliation %||% NA_character_),
    submitter_country = clean_text((report$submitter %||% list())$country %||% NA_character_),
    submitter_names = collapse_names((report$submitter %||% list())$names %||% list()),
    raw_json = x$raw_json
  )
}

message("Resolving and querying ", nrow(accession_bases), " accession bases via NCBI Datasets...")

lookup_results <- purrr::map(accession_bases$accession_base, fetch_accession_summary)
metadata_tbl <- purrr::map_dfr(lookup_results, extract_report_row) %>%
  arrange(accession_base)

raw_json_lines <- metadata_tbl %>%
  filter(!is.na(raw_json)) %>%
  pull(raw_json)

if (length(raw_json_lines) > 0) {
  writeLines(raw_json_lines, raw_output_path, useBytes = TRUE)
} else {
  file.create(raw_output_path)
}

metadata_output <- metadata_tbl %>%
  select(-raw_json)

enriched_output <- candidate_strains %>%
  mutate(accession_base = stringr::str_replace(accession, "\\.[0-9]+$", "")) %>%
  left_join(metadata_output, by = "accession_base") %>%
  relocate(accession_base, accession_version, ncbi_lookup_status, .after = accession)

write_csv(metadata_output, metadata_output_path, na = "")
write_csv(enriched_output, enriched_output_path, na = "")

cat("Candidate rows read:", nrow(candidate_strains), "\n")
cat("Distinct accession bases queried:", nrow(accession_bases), "\n")
cat("Successful NCBI metadata lookups:", sum(metadata_output$ncbi_lookup_status == "ok"), "\n")
cat("Failed NCBI metadata lookups:", sum(metadata_output$ncbi_lookup_status != "ok"), "\n")
cat("Wrote metadata table to", metadata_output_path, "\n")
cat("Wrote enriched candidate table to", enriched_output_path, "\n")
cat("Wrote raw JSONL to", raw_output_path, "\n")
