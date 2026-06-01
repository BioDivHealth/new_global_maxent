# ------------------------------------------------------------------------------
# summarize_vectraits_traits.R
# ------------------------------------------------------------------------------
# Purpose: Optional exploratory inspection of VecTraits datasets for a small set
#          of species and trait keywords, then summarise the trait names,
#          standardized names, and fields relevant to vector competence, blood
#          feeding, and transmission.
#
# Inputs : VecTraits API search queries
# Outputs: pathogen_association_data/staged/vector_screening/vectraits/vectraits_traits/
#          - vectraits_species_manifest.csv
#          - vectraits_trait_summary.csv
#          - vectraits_field_presence.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, httr, jsonlite, purrr, readr, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(y)
  }
  x
}

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
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
  fetch_json(url)
}

fetch_vectraits_dataset <- function(dataset_id, use_qa = FALSE) {
  url <- build_base_url(paste0("vectraits-dataset/", dataset_id, "/?format=json"), use_qa = use_qa)
  fetch_json(url)
}

first_non_empty <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

trait_signal <- function(x) {
  x <- clean_text(x)
  x_lower <- stringr::str_to_lower(paste(x, collapse = " | "))
  c(
    competence = str_detect(x_lower, "\\bcompetence\\b|vector competence|transmission potential|infective|infected saliva|saliva"),
    feeding = str_detect(x_lower, "blood feeding|\\bfeed\\b|host preference|host choice|bite|biting"),
    transmission = str_detect(x_lower, "\\btransmission\\b|dissemination|infection|vectorial capacity"),
    thermal = str_detect(x_lower, "temperature|thermal|degree|celsius|growth"),
    life_history = str_detect(x_lower, "fecund|survival|longevity|mortality|body size|development")
  )
}

# ------------------------------| Configuration |------------------------------
use_qa <- identical(tolower(Sys.getenv("VECTRAITS_USE_QA", "false")), "true")
max_datasets_per_query <- as.integer(Sys.getenv("VECTRAITS_MAX_DATASETS_PER_QUERY", "5"))

species_keywords <- c(
  "Aedes aegypti",
  "Aedes albopictus",
  "Culex pipiens",
  "Anopheles gambiae",
  "Culex quinquefasciatus"
)

trait_keywords <- c(
  "vector competence",
  "blood feeding",
  "transmission potential",
  "feeding",
  "temperature",
  "survival",
  "fecundity"
)

output_dir <- file.path(vector_screening_vectraits_dir, "vectraits_traits")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_manifest_path <- file.path(output_dir, "vectraits_species_manifest.csv")
trait_summary_path <- file.path(output_dir, "vectraits_trait_summary.csv")
field_presence_path <- file.path(output_dir, "vectraits_field_presence.csv")

queries <- tibble(
  query_type = c(rep("species", length(species_keywords)), rep("trait", length(trait_keywords))),
  query = c(species_keywords, trait_keywords)
)

query_hits <- purrr::map_dfr(seq_len(nrow(queries)), function(i) {
  q <- queries$query[[i]]
  result <- search_vectraits(q, use_qa = use_qa)

  if (!isTRUE(result$ok)) {
    return(tibble(
      query_type = queries$query_type[[i]],
      query = q,
      search_ok = FALSE,
      http_status = result$status %||% NA_integer_,
      dataset_id = NA_integer_
    ))
  }

  ids <- result$payload$ids %||% integer()
  if (length(ids) == 0) {
    return(tibble(
      query_type = queries$query_type[[i]],
      query = q,
      search_ok = TRUE,
      http_status = result$status,
      dataset_id = NA_integer_
    ))
  }

  tibble(
    query_type = queries$query_type[[i]],
    query = q,
    search_ok = TRUE,
    http_status = result$status,
    dataset_id = as.integer(ids)
  )
})

query_hits <- query_hits %>%
  mutate(query = clean_text(query)) %>%
  arrange(query_type, query, dataset_id) %>%
  distinct(query_type, query, dataset_id, .keep_all = TRUE)

candidate_ids <- query_hits %>%
  filter(!is.na(dataset_id)) %>%
  distinct(dataset_id) %>%
  pull(dataset_id)

candidate_ids <- head(candidate_ids, max_datasets_per_query * length(unique(queries$query)))

dataset_bundle <- purrr::map(candidate_ids, function(dataset_id) {
  result <- fetch_vectraits_dataset(dataset_id, use_qa = use_qa)
  if (!isTRUE(result$ok)) {
    return(list(
      dataset_id = dataset_id,
      fetch_ok = FALSE,
      http_status = result$status %||% NA_integer_,
      error = result$error %||% "Unknown error",
      data = tibble(),
      fields = tibble()
    ))
  }

  payload <- result$payload
  data <- payload$results %||% payload$data %||% tibble()
  if (!is.data.frame(data)) {
    data <- tibble()
  }

  fields <- tibble(
    dataset_id = dataset_id,
    field_name = names(data),
    field_class = purrr::map_chr(data, ~ paste(class(.x), collapse = "|")),
    non_missing = purrr::map_int(data, ~ sum(!is.na(.x)))
  )

  list(
    dataset_id = dataset_id,
    fetch_ok = TRUE,
    http_status = result$status,
    error = NA_character_,
    data = data,
    fields = fields
  )
})

species_manifest <- purrr::map_dfr(dataset_bundle, function(x) {
  data <- x$data
  if (nrow(data) == 0) {
    return(tibble(
      dataset_id = x$dataset_id,
      fetch_ok = x$fetch_ok,
      http_status = x$http_status,
      error = x$error,
      n_rows = 0L,
      n_columns = 0L,
      original_trait_names = NA_character_,
      standardised_trait_names = NA_character_,
      interactor1 = NA_character_,
      interactor2 = NA_character_
    ))
  }

  tibble(
    dataset_id = x$dataset_id,
    fetch_ok = x$fetch_ok,
    http_status = x$http_status,
    error = x$error,
    n_rows = nrow(data),
    n_columns = ncol(data),
    original_trait_names = paste(sort(unique(stats::na.omit(clean_text(data$OriginalTraitName)))), collapse = "; "),
    standardised_trait_names = paste(sort(unique(stats::na.omit(clean_text(data$StandardisedTraitName)))), collapse = "; "),
    interactor1 = paste(
      sort(unique(stats::na.omit(clean_text(paste(data$Interactor1Genus, data$Interactor1Species))))),
      collapse = "; "
    ),
    interactor2 = paste(
      sort(unique(stats::na.omit(clean_text(paste(data$Interactor2Genus, data$Interactor2Species))))),
      collapse = "; "
    )
  )
})

trait_summary <- purrr::map_dfr(dataset_bundle, function(x) {
  data <- x$data
  if (nrow(data) == 0) {
    return(tibble())
  }

  data %>%
    mutate(
      original_trait = clean_text(OriginalTraitName),
      standardized_trait = clean_text(StandardisedTraitName),
      original_def = clean_text(OriginalTraitDef),
      standardized_def = clean_text(StandardisedTraitDef),
      notes = clean_text(Notes),
      physical_process = clean_text(PhysicalProcess)
    ) %>%
    group_by(dataset_id = x$dataset_id, original_trait, standardized_trait) %>%
    summarise(
      n_rows = n(),
      original_def = first_non_empty(original_def),
      standardized_def = first_non_empty(standardized_def),
      notes = first_non_empty(notes),
      physical_process = first_non_empty(physical_process),
      competence_signal = any(trait_signal(c(original_trait, original_def, standardized_def, notes, physical_process))["competence"]),
      feeding_signal = any(trait_signal(c(original_trait, original_def, standardized_def, notes, physical_process))["feeding"]),
      transmission_signal = any(trait_signal(c(original_trait, original_def, standardized_def, notes, physical_process))["transmission"]),
      thermal_signal = any(trait_signal(c(original_trait, original_def, standardized_def, notes, physical_process))["thermal"]),
      life_history_signal = any(trait_signal(c(original_trait, original_def, standardized_def, notes, physical_process))["life_history"]),
      .groups = "drop"
    ) %>%
    arrange(desc(n_rows), original_trait)
})

field_presence <- purrr::map_dfr(dataset_bundle, ~ .x$fields) %>%
  mutate(
    field_name_lower = str_to_lower(field_name),
    competence_relevant = str_detect(field_name_lower, "trait|process|temp|feed|trans|compet|saliva|infection|blood|host|vector|fecund|survival")
  ) %>%
  arrange(dataset_id, desc(competence_relevant), field_name)

write_csv(query_hits, species_manifest_path, na = "")
write_csv(trait_summary, trait_summary_path, na = "")
write_csv(field_presence, field_presence_path, na = "")

cat("VecTraits trait summary complete.\n")
cat("Query hits:", nrow(query_hits), "\n")
cat("Datasets fetched:", length(dataset_bundle), "\n")
cat("Species manifest:", species_manifest_path, "\n")
cat("Trait summary:", trait_summary_path, "\n")
cat("Field presence:", field_presence_path, "\n")
