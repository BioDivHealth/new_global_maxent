# ------------------------------------------------------------------------------
# 3_WHO_DON_Title_Country_Extraction.R
# ------------------------------------------------------------------------------
# Purpose: Recover country evidence from WHO DON titles using robust suffix
#          parsing, geography aliases, and WHO-region handling.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, purrr, readr, stringr, tibble, tidyr)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
records <- read_who_don_records("who_don_records_clean.csv")
gazetteer <- load_who_don_geography_gazetteer()
generic_labels <- generic_title_labels()

clean_title_suffix <- function(title) {
  title <- clean_scalar_text(title)

  if (is.na(title)) {
    return(NA_character_)
  }

  suffix <- stringr::str_match(
    title,
    ".*(?:\\s+[-–—:]\\s*|[-–—]\\s+|-(?=[A-Z]))(.+)$"
  )[, 2]

  clean_scalar_text(suffix)
}

normalize_title_suffix <- function(raw_suffix) {
  raw_suffix <- clean_scalar_text(raw_suffix)

  if (is.na(raw_suffix)) {
    return(NA_character_)
  }

  paren_match <- stringr::str_match(
    raw_suffix,
    "^(?:update|global update|situation report|regional update)\\s*\\((.+)\\)$"
  )

  if (!is.na(paren_match[, 2])) {
    return(clean_scalar_text(paren_match[, 2]))
  }

  in_match <- stringr::str_match(
    raw_suffix,
    "^(?:situation|outbreak|outbreaks|cases|case|confirmed|detection|detected|reported|infection|infections?)\\s+in\\s+(.+)$"
  )

  if (!is.na(in_match[, 2])) {
    return(clean_scalar_text(in_match[, 2]))
  }

  raw_suffix
}

match_single_geography <- function(raw_piece) {
  raw_piece <- clean_scalar_text(raw_piece)

  if (is.na(raw_piece)) {
    return(tibble())
  }

  piece_key <- normalize_match_text(raw_piece)

  generic_hit <- generic_labels %>%
    filter(piece_key == label_key | stringr::str_starts(piece_key, label_key))

  if (nrow(generic_hit) > 0) {
    return(
      tibble(
        geography_name = raw_piece,
        geography_type = generic_hit$label_type[[1]],
        country = NA_character_,
        country_standard = NA_character_,
        raw_match = raw_piece
      )
    )
  }

  direct_hit <- gazetteer %>%
    filter(alias_key == piece_key)

  if (nrow(direct_hit) > 0) {
    return(
      tibble(
        geography_name = raw_piece,
        geography_type = direct_hit$geography_type[[1]],
        country = ifelse(direct_hit$geography_type[[1]] == "country", raw_piece, NA_character_),
        country_standard = direct_hit$country_standard[[1]],
        raw_match = raw_piece
      )
    )
  }

  tibble(
    geography_name = raw_piece,
    geography_type = "unknown",
    country = NA_character_,
    country_standard = NA_character_,
    raw_match = raw_piece
  )
}

split_multi_geography <- function(raw_suffix) {
  raw_suffix <- clean_scalar_text(raw_suffix)

  if (is.na(raw_suffix)) {
    return(character(0))
  }

  direct_match <- match_single_geography(raw_suffix)

  if (nrow(direct_match) > 0 && direct_match$geography_type[[1]] != "unknown") {
    return(raw_suffix)
  }

  split_candidates <- stringr::str_split(
    raw_suffix,
    stringr::regex("\\s+(?:and|/)\\s+|\\s*;\\s*"),
    simplify = FALSE
  )[[1]]

  split_candidates <- clean_scalar_text(split_candidates)
  split_candidates <- split_candidates[!is.na(split_candidates)]

  if (length(split_candidates) <= 1) {
    comma_split <- stringr::str_split(raw_suffix, stringr::regex("\\s*,\\s*"), simplify = FALSE)[[1]]
    comma_split <- clean_scalar_text(comma_split)
    comma_split <- comma_split[!is.na(comma_split)]

    if (length(comma_split) > 1) {
      matched_parts <- purrr::map_lgl(comma_split, ~ nrow(match_single_geography(.x)) > 0)

      if (all(matched_parts)) {
        return(comma_split)
      }
    }

    return(raw_suffix)
  }

  matched_parts <- purrr::map_lgl(split_candidates, ~ {
    hit <- match_single_geography(.x)
    nrow(hit) > 0 && !all(hit$geography_type == "unknown")
  })

  if (all(matched_parts)) {
    return(split_candidates)
  }

  raw_suffix
}

message_if("Extracting title-based WHO DON country evidence...", verbose = verbose)

title_rows <- purrr::map_dfr(
  seq_len(nrow(records)),
  function(i) {
    suffix <- clean_title_suffix(records$Title[[i]])
    suffix <- normalize_title_suffix(suffix)

    if (is.na(suffix)) {
      return(tibble())
    }

    pieces <- split_multi_geography(suffix)
    hits <- purrr::map_dfr(pieces, match_single_geography)

    if (nrow(hits) == 0) {
      return(tibble())
    }

    hits %>%
      mutate(
        DonId = records$DonId[[i]],
        record_id = dplyr::coalesce(records$record_id[[i]], records$DonId[[i]], records$Id[[i]]),
        Title = records$Title[[i]],
        publication_datetime_utc = records$publication_datetime_utc[[i]],
        article_url = records$article_url[[i]],
        extraction_source = "title",
        extraction_basis = "title_suffix",
        confidence = "high"
      ) %>%
      select(
        DonId,
        record_id,
        Title,
        publication_datetime_utc,
        article_url,
        geography_name,
        geography_type,
        country,
        country_standard,
        extraction_source,
        extraction_basis,
        confidence,
        raw_match
      )
  }
)

title_rows <- title_rows %>%
  mutate(
    geography_type = dplyr::case_when(
      geography_type %in% c("generic_global", "generic_multicountry", "generic_regional", "generic_report") ~ "unknown",
      TRUE ~ geography_type
    )
  ) %>%
  distinct()

title_rows <- title_rows %>%
  filter(!raw_match %in% c("Global situation", "Global update", "Global Situation", "Global Update"))

title_path <- who_don_output_path("who_don_title_country_evidence.csv")
readr::write_csv(title_rows, title_path, na = "")

message_if(
  "Saved title-based WHO DON evidence to: ",
  title_path,
  " | rows=",
  nrow(title_rows),
  verbose = verbose
)
