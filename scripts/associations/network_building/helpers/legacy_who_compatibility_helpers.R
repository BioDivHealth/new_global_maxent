# -----------------------------------------------------------------------------|
# legacy_who_compatibility_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for retained legacy WHO compatibility stages.
# -----------------------------------------------------------------------------|

requireNamespace("dplyr", quietly = TRUE)
requireNamespace("stringr", quietly = TRUE)

legacy_who_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

legacy_who_collapse_vals <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x)])
  paste(x, collapse = sep)
}

legacy_who_clean_text <- function(x) {
  x <- ifelse(is.na(x), NA_character_, x)
  x <- dplyr::na_if(x, "")
  x <- dplyr::na_if(x, "NA")
  ifelse(is.na(x), NA_character_, stringr::str_squish(x))
}

legacy_who_first_non_missing <- function(x) {
  x <- legacy_who_clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

legacy_who_collapse_unique <- function(x) {
  x <- legacy_who_clean_text(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = "; ")
}

legacy_who_safe_lower <- function(x) {
  ifelse(is.na(x), NA_character_, stringr::str_to_lower(legacy_who_clean_text(x)))
}

legacy_who_normalize_pathogen <- function(x) {
  dplyr::case_when(
    legacy_who_safe_lower(x) == "subgenus sarbecovirus" ~ "Subgenus Sarbecovirus",
    legacy_who_safe_lower(x) == "subgenus merbecovirus" ~ "Subgenus Merbecovirus",
    TRUE ~ legacy_who_clean_text(x)
  )
}
