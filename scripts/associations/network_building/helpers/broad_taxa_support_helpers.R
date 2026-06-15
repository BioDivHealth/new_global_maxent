# -----------------------------------------------------------------------------|
# broad_taxa_support_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for broad-taxa support stage scripts.
# -----------------------------------------------------------------------------|

broad_taxa_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

broad_taxa_null_coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

`%||%` <- broad_taxa_null_coalesce

broad_taxa_collapse_names <- function(x) {
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- broad_taxa_clean_text(x)
  x <- unique(stats::na.omit(x))
  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

broad_taxa_collapse_lineage <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  lineage_names <- purrr::map_chr(
    x,
    ~ broad_taxa_clean_text(.x$name %||% NA_character_)
  )
  lineage_names <- unique(stats::na.omit(lineage_names))

  if (length(lineage_names) == 0) {
    return(NA_character_)
  }

  paste(lineage_names, collapse = "; ")
}
