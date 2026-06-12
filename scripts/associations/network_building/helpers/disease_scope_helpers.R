# -----------------------------------------------------------------------------|
# disease_scope_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for WHO disease-scope stage scripts.
# -----------------------------------------------------------------------------|

disease_scope_provenance_cols <- function() {
  c(
    "is_priority_pathogen",
    "is_prototype_pathogen",
    "in_gibb_etal",
    "in_empres_i",
    "priority_prototype_status",
    "region_africa",
    "region_americas",
    "region_europe",
    "region_mediterranean",
    "region_se_asia",
    "region_western_pacific"
  )
}

disease_scope_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

disease_scope_clean_key <- function(x) {
  key <- disease_scope_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

disease_scope_flag_from_mark <- function(x) {
  x <- disease_scope_clean_text(x)
  !is.na(x)
}

disease_scope_first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}
