# ------------------------------------------------------------------------------
# build_vectraits_species_tables.R
# ------------------------------------------------------------------------------
# Purpose: Optional exploratory pull of VecTraits datasets for a few target
#          vector species, combining row-level trait evidence into one table per
#          species while preserving dataset provenance.
#
# Inputs : VecTraits API searches for the configured species names
# Outputs: pathogen_association_data/staged/vector_screening/vectraits/
#          vectraits_species_tables/
#          - vectraits_species_manifest.csv
#          - <species>_vectraits_combined.csv
#          - <species>_vectraits_datasets.csv
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

species_key <- function(genus, species) {
  paste0(
    stringr::str_to_lower(clean_text(genus)),
    " ",
    stringr::str_to_lower(clean_text(species))
  )
}

slugify <- function(x) {
  x <- stringr::str_to_lower(clean_text(x))
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
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

# ------------------------------| Configuration |------------------------------
use_qa <- identical(tolower(Sys.getenv("VECTRAITS_USE_QA", "false")), "true")

target_species <- tibble(
  species_name = c(
    "Aedes albopictus",
    "Aedes aegypti",
    "Aedes vittatus"
  )
)

output_dir <- file.path(vector_screening_vectraits_dir, "vectraits_species_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(output_dir, "vectraits_species_manifest.csv")

# ------------------------------| Search species |------------------------------
species_hits <- purrr::map_dfr(seq_len(nrow(target_species)), function(i) {
  species_name <- target_species$species_name[[i]]
  result <- search_vectraits(species_name, use_qa = use_qa)

  if (!isTRUE(result$ok)) {
    return(tibble(
      species_name = species_name,
      search_ok = FALSE,
      http_status = result$status %||% NA_integer_,
      dataset_id = NA_integer_
    ))
  }

  ids <- result$payload$ids %||% integer()
  if (length(ids) == 0) {
    return(tibble(
      species_name = species_name,
      search_ok = TRUE,
      http_status = result$status,
      dataset_id = NA_integer_
    ))
  }

  tibble(
    species_name = species_name,
    search_ok = TRUE,
    http_status = result$status,
    dataset_id = as.integer(ids)
  )
})

species_hits <- species_hits %>%
  mutate(
    species_name = clean_text(species_name)
  ) %>%
  arrange(species_name, dataset_id) %>%
  distinct(species_name, dataset_id, .keep_all = TRUE)

# ------------------------------| Fetch datasets |------------------------------
species_bundle <- purrr::map(target_species$species_name, function(species_name) {
  species_ids <- species_hits %>%
    filter(species_name == !!species_name, !is.na(dataset_id)) %>%
    pull(dataset_id) %>%
    unique()

  if (length(species_ids) == 0) {
    return(list(
      species_name = species_name,
      datasets = tibble(),
      rows = tibble()
    ))
  }

  dataset_bundle <- purrr::map(species_ids, function(dataset_id) {
    result <- fetch_vectraits_dataset(dataset_id, use_qa = use_qa)

    if (!isTRUE(result$ok)) {
      return(list(
        dataset_id = dataset_id,
        fetch_ok = FALSE,
        http_status = result$status %||% NA_integer_,
        error = result$error %||% "Unknown error",
        data = tibble()
      ))
    }

    payload <- result$payload
    data <- payload$results %||% payload$data %||% tibble()
    if (!is.data.frame(data)) {
      data <- tibble()
    }

    list(
      dataset_id = dataset_id,
      fetch_ok = TRUE,
      http_status = result$status,
      error = NA_character_,
      data = data
    )
  })

  dataset_summary <- purrr::map_dfr(dataset_bundle, function(x) {
    data <- x$data
    if (nrow(data) == 0) {
      return(tibble(
        species_name = species_name,
        dataset_id = x$dataset_id,
        fetch_ok = x$fetch_ok,
        http_status = x$http_status,
        error = x$error,
        n_rows = 0L,
        n_columns = 0L,
        original_trait_names = NA_character_,
        interactor1 = NA_character_,
        interactor2 = NA_character_
      ))
    }

    tibble(
      species_name = species_name,
      dataset_id = x$dataset_id,
      fetch_ok = x$fetch_ok,
      http_status = x$http_status,
      error = x$error,
      n_rows = nrow(data),
      n_columns = ncol(data),
      original_trait_names = paste(sort(unique(stats::na.omit(clean_text(data$OriginalTraitName)))), collapse = "; "),
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

  combined_rows <- purrr::map_dfr(dataset_bundle, function(x) {
    data <- x$data
    if (nrow(data) == 0) {
      return(tibble())
    }

    data %>%
      mutate(
        dataset_id = x$dataset_id,
        species_name = species_name,
        interactor1_key = species_key(Interactor1Genus, Interactor1Species),
        target_species_key = species_key(str_split_fixed(species_name, " ", 2)[, 1], str_split_fixed(species_name, " ", 2)[, 2])
      ) %>%
      filter(
        !is.na(interactor1_key),
        interactor1_key == target_species_key
      ) %>%
      mutate(
        across(everything(), ~ as.character(.x))
      ) %>%
      mutate(
        across(where(is.character), clean_text)
      )
  })

  list(
    species_name = species_name,
    datasets = dataset_summary,
    rows = combined_rows
  )
})

# ------------------------------| Write outputs |------------------------------
species_manifest <- purrr::map_dfr(species_bundle, function(x) {
  rows <- x$rows
  datasets <- x$datasets

  tibble(
    species_name = x$species_name,
    species_slug = slugify(x$species_name),
    n_datasets = nrow(datasets),
    n_rows = nrow(rows),
    n_unique_traits = if (nrow(rows) > 0 && "OriginalTraitName" %in% names(rows)) {
      n_distinct(clean_text(rows$OriginalTraitName))
    } else {
      0L
    }
  )
})

write_csv(species_manifest, manifest_path, na = "")

for (i in seq_along(species_bundle)) {
  species_name <- species_bundle[[i]]$species_name
  species_slug <- slugify(species_name)
  combined_path <- file.path(output_dir, paste0(species_slug, "_vectraits_combined.csv"))
  dataset_path <- file.path(output_dir, paste0(species_slug, "_vectraits_datasets.csv"))

  write_csv(species_bundle[[i]]$rows, combined_path, na = "")
  write_csv(species_bundle[[i]]$datasets, dataset_path, na = "")
}

cat("VecTraits species tables complete.\n")
cat("Species searched:", nrow(target_species), "\n")
cat("Manifest:", manifest_path, "\n")
