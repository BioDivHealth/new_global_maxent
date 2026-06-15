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

legacy_who_manual_pathogen_map <- function() {
  tibble::tribble(
    ~Pathogen_raw, ~Pathogen_canonical, ~canonicalization_status,
    "Salmonella enterica", "Salmonella enterica non typhoidal serovars", "manual_specificity_map",
    "Shigella dysenteriae", "Shigella dysenteriae serotype 1", "manual_specificity_map",
    "Vibrio cholerae", "Vibrio cholerae serogroup 0139", "manual_specificity_map",
    "Betacoronavirus pandemicum", "Subgenus Sarbecovirus", "manual_synonym_map",
    "Severe acute respiratory syndrome-related coronavirus", "Subgenus Sarbecovirus", "manual_synonym_map",
    "Betacoronavirus cameli", "Subgenus Merbecovirus", "manual_synonym_map",
    "Zaire ebolavirus", "Orthoebolavirus zairense", "manual_synonym_map",
    "Enterovirus c", "Enterovirus coxsackiepol", "manual_synonym_map",
    "Human poliovirus", "Enterovirus coxsackiepol", "manual_synonym_map",
    "Enterovirus alphacoxsackie", "Enterovirus alphacoxsackie 71", "manual_synonym_map",
    "Enterovirus a", "Enterovirus alphacoxsackie 71", "manual_synonym_map",
    "Enterovirus deconjuncti", "Enterovirus deconjucti 68", "manual_synonym_map",
    "Alphainfluenzavirus influenzae", "Alphainfluenzavirus influenzae", "manual_group_retained",
    "Protoparvovirus carnivoran1", "Protoparvovirus carnivoran", "manual_group_map",
    "Protoparvovirus carnivoran3", "Protoparvovirus carnivoran", "manual_group_map",
    "Protoparvovirus carnivoran4", "Protoparvovirus carnivoran", "manual_group_map",
    "Protoparvovirus carnivoran5", "Protoparvovirus carnivoran", "manual_group_map"
  ) %>%
    dplyr::mutate(Pathogen_raw_key = legacy_who_safe_lower(Pathogen_raw))
}

legacy_who_zoonotic_override <- function() {
  tibble::tribble(
    ~Pathogen_canonical, ~is_zoonotic_override, ~zoonotic_status_override,
    "Alphainfluenzavirus influenzae", TRUE, "zoonotic_group_retained"
  ) %>%
    dplyr::mutate(Pathogen_canonical_key = legacy_who_safe_lower(Pathogen_canonical))
}
