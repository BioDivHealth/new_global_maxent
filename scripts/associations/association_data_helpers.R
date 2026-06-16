################################################################################
# association_data_helpers.R
################################################################################
# Purpose: Shared low-level data helpers for association pipeline scripts.
#
# Notes  : Keep these helpers generic. Workflow-specific interpretation belongs
#          in the scoped role, readiness, vector, or network scripts.
################################################################################

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "NULL", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

clean_key <- function(x) {
  x %>%
    clean_text() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("&", " and ") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

stable_slug <- function(...) {
  paste(..., sep = "_") %>%
    clean_key() %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

stable_candidate_id <- function(prefix, ...) {
  slug <- stable_slug(...)
  paste0(prefix, "_", stringr::str_sub(slug, 1, 110))
}

is_true <- function(x) {
  as.character(x) %in% c("TRUE", "true", "True", "1", "yes", "Yes", "YES")
}

missing_as_false <- function(x) {
  out <- is_true(x)
  out[is.na(out)] <- FALSE
  out
}

read_csv_layer <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Required input is missing: ", path, call. = FALSE)
    }
    return(NULL)
  }

  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    dplyr::mutate(dplyr::across(where(is.character), clean_text))
}

first_non_empty <- function(...) {
  args <- list(...)
  if (length(args) == 0) {
    return(character())
  }

  out <- rep(NA_character_, length(args[[1]]))
  for (arg in args) {
    values <- clean_text(arg)
    fill <- is.na(out) & !is.na(values)
    out[fill] <- values[fill]
  }
  out
}

collapse_unique <- function(x, sep = "; ") {
  values <- sort(unique(clean_text(x)))
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(values, collapse = sep)
}

combine_reasons <- function(...) {
  values <- c(...)
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(unique(values), collapse = "; ")
}

max_tier <- function(strict, strong, supported, broad) {
  dplyr::case_when(
    strict ~ "strict",
    strong ~ "strong",
    supported ~ "supported",
    broad ~ "broad",
    TRUE ~ "excluded"
  )
}
