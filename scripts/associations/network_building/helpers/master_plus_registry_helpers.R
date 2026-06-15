# -----------------------------------------------------------------------------|
# master_plus_registry_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for master-plus registry stage scripts.
# -----------------------------------------------------------------------------|

registry_normalize_name <- function(x) {
  key <- stringr::str_to_lower(x)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[[:punct:]]+", " ")
  stringr::str_squish(key)
}

registry_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

registry_clean_key <- function(x) {
  key <- registry_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

registry_collapse_unique <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) {
    NA_character_
  } else {
    paste(x, collapse = "; ")
  }
}

registry_coalesce_chr <- function(...) {
  values <- purrr::map(list(...), as.character)
  dplyr::coalesce(!!!values)
}

registry_pick_preferred <- function(preferred_source, virion_value, clover_value) {
  dplyr::case_when(
    preferred_source == "virion" ~ as.character(virion_value),
    preferred_source == "clover" ~ as.character(clover_value),
    TRUE ~ NA_character_
  )
}
