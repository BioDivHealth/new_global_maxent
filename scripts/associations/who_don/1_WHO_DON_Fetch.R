# ------------------------------------------------------------------------------
# 1_WHO_DON_Fetch.R
# ------------------------------------------------------------------------------
# Purpose: Download the WHO Disease Outbreak News (DON) collection from the
#          public WHO API, save a raw JSON snapshot, and write flattened tables
#          that are easier to use in downstream country-evidence workflows.
#
# Inputs : WHO DON API endpoint
# Outputs: pathogen_association_data/WHO/disease_outbreak_news/who_don_raw.json
#          pathogen_association_data/WHO/disease_outbreak_news/who_don_records.csv
#
# Notes  : Optional environment variables:
#          `WHO_DON_START_DATE` (YYYY-MM-DD)
#          `WHO_DON_END_DATE`   (YYYY-MM-DD)
#          `WHO_DON_VERBOSE`    (TRUE/FALSE, default TRUE)
#          `WHO_DON_PAGE_SIZE`  (integer, default 50)
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, httr2, jsonlite, purrr, readr, stringr, tibble, tidyr)

# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
parse_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || value == "") {
    return(default)
  }

  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_env_date <- function(name, default = as.Date(NA)) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || value == "") {
    return(default)
  }

  parsed <- as.Date(value)

  if (is.na(parsed)) {
    stop(name, " must be a valid date in YYYY-MM-DD format.")
  }

  parsed
}

parse_env_integer <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || value == "") {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))

  if (is.na(parsed)) {
    stop(name, " must be an integer.")
  }

  parsed
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\r", "\n")
  x <- stringr::str_replace_all(x, "\n{2,}", "\n\n")
  x <- stringr::str_replace_all(x, "[ \t]+", " ")
  x <- stringr::str_replace_all(x, " ?\n ?", "\n")
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NULL")] <- NA_character_
  x
}

build_article_url <- function(item_default_url, url_name) {
  dplyr::case_when(
    !is.na(item_default_url) & item_default_url != "" ~ paste0(
      "https://www.who.int/emergencies/disease-outbreak-news/item",
      item_default_url
    ),
    !is.na(url_name) & url_name != "" ~ paste0(
      "https://www.who.int/emergencies/disease-outbreak-news/item/",
      url_name
    ),
    TRUE ~ NA_character_
  )
}

message_if <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

# ------------------------------------------------------------------------------
# Parameters and paths ---------------------------------------------------------
# ------------------------------------------------------------------------------
api_url <- "https://www.who.int/api/emergencies/diseaseoutbreaknews"
verbose <- parse_env_flag("WHO_DON_VERBOSE", default = TRUE)
start_date <- parse_env_date("WHO_DON_START_DATE")
end_date <- parse_env_date("WHO_DON_END_DATE")
page_size <- parse_env_integer("WHO_DON_PAGE_SIZE", default = 50L)

if (!is.na(start_date) && !is.na(end_date) && start_date > end_date) {
  stop("WHO_DON_START_DATE cannot be later than WHO_DON_END_DATE.")
}

if (is.na(page_size) || page_size <= 0) {
  stop("WHO_DON_PAGE_SIZE must be greater than 0.")
}

output_dir <- here("pathogen_association_data", "WHO", "disease_outbreak_news")
raw_json_path <- file.path(output_dir, "who_don_raw.json")
records_csv_path <- file.path(output_dir, "who_don_records.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Fetch DON collection ---------------------------------------------------------
# ------------------------------------------------------------------------------
message_if(
  "Fetching WHO DON collection from API with page size ", page_size, "...",
  verbose = verbose
)

fetch_don_page <- function(skip, top) {
  resp <- request(api_url) %>%
    req_url_query(`$skip` = skip, `$top` = top) %>%
    req_perform()

  raw_json <- resp_body_string(resp)
  parsed <- jsonlite::fromJSON(raw_json, simplifyVector = TRUE, flatten = TRUE)

  if (!"value" %in% names(parsed)) {
    stop("WHO DON API response did not contain a `value` field for skip=", skip, ".")
  }

  list(
    raw_json = raw_json,
    parsed = parsed,
    records = tibble::as_tibble(parsed$value)
  )
}

page_records <- list()
skip <- 0L
page_index <- 1L

repeat {
  page <- fetch_don_page(skip = skip, top = page_size)
  n_page <- nrow(page$records)

  message_if(
    "Fetched WHO DON page ", page_index,
    " | skip=", skip,
    " | records=", n_page,
    verbose = verbose
  )

  if (n_page == 0) {
    break
  }

  page_records[[page_index]] <- page$records

  if (n_page < page_size) {
    break
  }

  skip <- skip + page_size
  page_index <- page_index + 1L
}

if (length(page_records) == 0) {
  stop("WHO DON API returned zero records.")
}

records <- bind_rows(page_records)

raw_snapshot <- list(
  source_url = api_url,
  fetched_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  page_size = page_size,
  n_records = nrow(records),
  value = records
)

writeLines(
  jsonlite::toJSON(raw_snapshot, auto_unbox = TRUE, pretty = TRUE, na = "null"),
  raw_json_path
)

# ------------------------------------------------------------------------------
# Flatten and clean ------------------------------------------------------------
# ------------------------------------------------------------------------------
expected_cols <- c(
  "DonId", "Overview", "Assessment", "Advice", "OverrideTitle",
  "SystemSourceKey", "Title", "ItemDefaultUrl", "UseOverrideTitle",
  "TitleSuffix", "PublicationDateAndTime", "regionscountries", "Summary",
  "FurtherInformation", "Response", "UrlName", "Epidemiology",
  "IncludeInSitemap", "DateCreated", "PublicationDate", "LastModified",
  "Id", "FormattedDate"
)

missing_cols <- setdiff(expected_cols, names(records))

if (length(missing_cols) > 0) {
  for (col_name in missing_cols) {
    records[[col_name]] <- NA
  }
}

records <- records %>%
  select(all_of(expected_cols), everything()) %>%
  mutate(
    across(
      c(
        DonId, Overview, Assessment, Advice, OverrideTitle, SystemSourceKey,
        Title, ItemDefaultUrl, TitleSuffix, Summary, FurtherInformation,
        Response, UrlName, Epidemiology, FormattedDate, Id
      ),
      clean_text
    ),
    publication_datetime_utc = as.POSIXct(PublicationDateAndTime, tz = "UTC"),
    publication_date = as.Date(PublicationDate),
    date_created_utc = as.POSIXct(DateCreated, tz = "UTC"),
    last_modified_utc = as.POSIXct(LastModified, tz = "UTC"),
    record_id = dplyr::coalesce(DonId, Id),
    article_url = build_article_url(ItemDefaultUrl, UrlName)
  ) %>%
  arrange(desc(publication_datetime_utc), desc(last_modified_utc))

if (!is.na(start_date)) {
  records <- records %>%
    filter(is.na(publication_date) | publication_date >= start_date)
}

if (!is.na(end_date)) {
  records <- records %>%
    filter(is.na(publication_date) | publication_date <= end_date)
}

# ------------------------------------------------------------------------------
# Write outputs ----------------------------------------------------------------
# ------------------------------------------------------------------------------
readr::write_csv(records, records_csv_path, na = "")

# ------------------------------------------------------------------------------
# Console summary --------------------------------------------------------------
# ------------------------------------------------------------------------------
message_if(
  "Saved raw WHO DON JSON to: ", raw_json_path,
  verbose = verbose
)
message_if(
  "Saved flattened WHO DON records to: ", records_csv_path,
  verbose = verbose
)
message_if(
  "WHO DON records written: ", nrow(records),
  verbose = verbose
)
