# ------------------------------------------------------------------------------
# probe_vectraits_access.R
# ------------------------------------------------------------------------------
# Purpose: Optional exploratory probe of the VectorByte VecTraits API with a
#          small set of vector- and trait-oriented searches, then write a
#          compact manifest and schema summary for the datasets that come back.
#
# Inputs : VecTraits API search queries
# Outputs: pathogen_association_data/staged/vector_screening/vectraits/vectraits_probe/
#          - vectraits_probe_manifest.csv
#          - vectraits_probe_datasets.csv
#          - vectraits_probe_schema.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, httr, jsonlite, readr, stringr, tibble, purrr)

source(here("scripts", "associations", "working_inputs.R"))

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

first_non_empty <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

build_base_url <- function(endpoint, use_qa = FALSE) {
  paste0(
    "https://vectorbyte",
    if (use_qa) "-qa" else "",
    ".crc.nd.edu/portal/api/",
    endpoint
  )
}

fetch_json <- function(url) {
  response <- tryCatch(
    httr::GET(url),
    error = function(e) {
      return(list(ok = FALSE, status = NA_integer_, error = conditionMessage(e)))
    }
  )

  if (is.list(response) && identical(response$ok, FALSE)) {
    return(response)
  }

  status <- httr::status_code(response)
  if (status < 200 || status >= 300) {
    return(list(ok = FALSE, status = status, error = "HTTP request failed"))
  }

  payload <- tryCatch(
    jsonlite::fromJSON(
      httr::content(response, "text", encoding = "UTF-8"),
      flatten = TRUE
    ),
    error = function(e) {
      return(list(ok = FALSE, status = status, error = conditionMessage(e)))
    }
  )

  if (is.list(payload) && identical(payload$ok, FALSE)) {
    return(payload)
  }

  list(ok = TRUE, status = status, payload = payload)
}

search_vectraits <- function(keyword, use_qa = FALSE) {
  url <- paste0(
    build_base_url("vectraits-explorer/?format=json&keywords=", use_qa = use_qa),
    utils::URLencode(keyword, reserved = FALSE)
  )
  result <- fetch_json(url)

  if (!isTRUE(result$ok)) {
    return(tibble(
      query = keyword,
      dataset_id = NA_integer_,
      search_ok = FALSE,
      http_status = result$status %||% NA_integer_,
      title = NA_character_
    ))
  }

  ids <- result$payload$ids
  if (length(ids) == 0) {
    return(tibble(
      query = keyword,
      dataset_id = NA_integer_,
      search_ok = TRUE,
      http_status = result$status,
      title = NA_character_
    ))
  }

  tibble(
    query = keyword,
    dataset_id = as.integer(ids),
    search_ok = TRUE,
    http_status = result$status,
    title = NA_character_
  )
}

fetch_vectraits_dataset <- function(dataset_id, use_qa = FALSE) {
  url <- paste0(
    build_base_url(paste0("vectraits-dataset/", dataset_id, "/?format=json"), use_qa = use_qa)
  )
  result <- fetch_json(url)

  if (!isTRUE(result$ok)) {
    return(list(
      dataset_id = dataset_id,
      fetch_ok = FALSE,
      http_status = result$status %||% NA_integer_,
      error = result$error %||% "Unknown error",
      results = NULL,
      schema = tibble()
    ))
  }

  results <- NULL
  if (!is.null(result$payload$results) && is.data.frame(result$payload$results)) {
    results <- result$payload$results
  } else if (!is.null(result$payload$data) && is.data.frame(result$payload$data)) {
    results <- result$payload$data
  }

  if (is.null(results)) {
    results <- tibble()
  }

  schema <- if (ncol(results) > 0) {
    tibble(
      dataset_id = dataset_id,
      column_name = names(results),
      column_class = purrr::map_chr(results, ~ paste(class(.x), collapse = "|")),
      non_missing = purrr::map_int(results, ~ sum(!is.na(.x)))
    )
  } else {
    tibble(
      dataset_id = integer(),
      column_name = character(),
      column_class = character(),
      non_missing = integer()
    )
  }

  list(
    dataset_id = dataset_id,
    fetch_ok = TRUE,
    http_status = result$status,
    error = NA_character_,
    results = results,
    schema = schema
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(y)
  }
  x
}

# ------------------------------| Configuration |------------------------------
use_qa <- identical(tolower(Sys.getenv("VECTRAITS_USE_QA", "false")), "true")
max_datasets <- as.integer(Sys.getenv("VECTRAITS_MAX_DATASETS", "8"))

probe_keywords <- c(
  "mosquito",
  "aedes aegypti",
  "culex pipiens",
  "tick",
  "vector competence",
  "survival",
  "fecundity",
  "temperature"
)

output_dir <- file.path(vector_screening_vectraits_dir, "vectraits_probe")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(output_dir, "vectraits_probe_manifest.csv")
schema_path <- file.path(output_dir, "vectraits_probe_schema.csv")

# ------------------------------| Probe searches |------------------------------
search_results <- purrr::map_dfr(
  probe_keywords,
  ~ search_vectraits(.x, use_qa = use_qa)
)

search_results <- search_results %>%
  mutate(
    dataset_id = as.integer(dataset_id),
    query = clean_text(query)
  ) %>%
  arrange(query, dataset_id) %>%
  distinct(query, dataset_id, .keep_all = TRUE)

if (nrow(search_results) == 0 || all(is.na(search_results$dataset_id))) {
  write_csv(
    tibble(
      status = "no_search_results",
      use_qa = use_qa,
      probe_keywords = paste(probe_keywords, collapse = "; ")
    ),
    manifest_path
  )
  stop("VecTraits search returned no dataset IDs.")
}

candidate_ids <- search_results$dataset_id[!is.na(search_results$dataset_id)]
candidate_ids <- unique(candidate_ids)
candidate_ids <- candidate_ids[seq_len(min(length(candidate_ids), max_datasets))]

# ------------------------------| Fetch datasets |------------------------------
dataset_results <- purrr::map(candidate_ids, ~ fetch_vectraits_dataset(.x, use_qa = use_qa))

dataset_manifest <- purrr::map_dfr(dataset_results, function(x) {
  results <- x$results

  if (nrow(results) == 0) {
    return(tibble(
      dataset_id = x$dataset_id,
      fetch_ok = x$fetch_ok,
      http_status = x$http_status,
      error = x$error,
      n_rows = 0L,
      n_columns = 0L,
      result_fields = NA_character_,
      submitted_by = NA_character_,
      original_trait_name = NA_character_,
      interactor1_genus = NA_character_,
      interactor1_species = NA_character_,
      interactor2_genus = NA_character_,
      interactor2_species = NA_character_
    ))
  }

  tibble(
    dataset_id = x$dataset_id,
    fetch_ok = x$fetch_ok,
    http_status = x$http_status,
    error = x$error,
    n_rows = nrow(results),
    n_columns = ncol(results),
    result_fields = paste(names(results), collapse = "; "),
    submitted_by = first_non_empty(results$SubmittedBy),
    original_trait_name = first_non_empty(results$OriginalTraitName),
    interactor1_genus = first_non_empty(results$Interactor1Genus),
    interactor1_species = first_non_empty(results$Interactor1Species),
    interactor2_genus = first_non_empty(results$Interactor2Genus),
    interactor2_species = first_non_empty(results$Interactor2Species)
  )
})

schema_results <- purrr::map_dfr(dataset_results, ~ .x$schema)

write_csv(search_results, manifest_path, na = "")
write_csv(dataset_manifest, file.path(output_dir, "vectraits_probe_datasets.csv"), na = "")
write_csv(schema_results, schema_path, na = "")

cat("VecTraits probe complete.\n")
cat("Search rows:", nrow(search_results), "\n")
cat("Datasets fetched:", nrow(dataset_manifest), "\n")
cat("Manifest:", manifest_path, "\n")
cat("Dataset summary:", file.path(output_dir, "vectraits_probe_datasets.csv"), "\n")
cat("Schema summary:", schema_path, "\n")
