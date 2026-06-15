# -----------------------------------------------------------------------------|
# master_plus_host_network_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for master-plus host-network stage scripts.
# -----------------------------------------------------------------------------|

host_network_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

host_network_clean_key <- function(x) {
  key <- host_network_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

host_network_split_semicolon_values <- function(x) {
  x <- host_network_clean_text(x)
  if (is.na(x)) {
    return(character(0))
  }

  values <- stringr::str_split(x, ";", simplify = FALSE)[[1]]
  values <- stringr::str_squish(values)
  values <- purrr::discard(values, ~ .x == "")
  unique(values)
}

host_network_first_non_missing <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    NA_character_
  } else {
    x[[1]]
  }
}

host_network_collapse_reasons <- function(...) {
  reasons <- c(...)
  reasons <- reasons[!is.na(reasons) & reasons != ""]
  if (length(reasons) == 0) {
    NA_character_
  } else {
    paste(unique(reasons), collapse = "; ")
  }
}

host_network_collapse_unique <- function(x) {
  x <- host_network_clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

host_network_add_missing_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  for (col in missing) {
    data[[col]] <- NA
  }
  data
}

host_network_association_key <- function(data) {
  paste(
    host_network_clean_key(data$Disease_name),
    host_network_clean_text(data$PathogenTaxID),
    host_network_clean_key(data$Pathogen),
    host_network_clean_text(data$HostTaxID),
    host_network_clean_key(data$Host),
    sep = "|||"
  )
}

host_network_is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
}

host_network_collapse_true_flag <- function(x) {
  values <- x[!is.na(x)]

  if (length(values) == 0) {
    return(NA)
  }

  any(host_network_is_true(values))
}
