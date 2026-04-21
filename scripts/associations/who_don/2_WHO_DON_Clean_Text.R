# ------------------------------------------------------------------------------
# 2_WHO_DON_Clean_Text.R
# ------------------------------------------------------------------------------
# Purpose: Convert the raw WHO DON record export into machine-readable plain
#          text and an inspection-friendly derivative for Excel/RStudio review.
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble, tidyr, xml2)

source(here("scripts", "associations", "who_don", "who_don_helpers.R"))

verbose <- TRUE
output_dir <- who_don_output_dir()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

records <- read_who_don_records("who_don_records.csv")

text_map <- c(
  Summary = "summary_text",
  Overview = "overview_text",
  Epidemiology = "epidemiology_text",
  Assessment = "assessment_text",
  Response = "response_text",
  Advice = "advice_text",
  FurtherInformation = "further_information_text"
)

message_if("Cleaning WHO DON HTML text fields...", verbose = verbose)

for (source_col in names(text_map)) {
  target_col <- text_map[[source_col]]
  records[[target_col]] <- html_to_text(records[[source_col]])
}

if (!"article_url" %in% names(records)) {
  records <- records %>%
    mutate(
      article_url = dplyr::case_when(
        !is.na(ItemDefaultUrl) & ItemDefaultUrl != "" ~ paste0(
          "https://www.who.int/emergencies/disease-outbreak-news/item",
          ItemDefaultUrl
        ),
        !is.na(UrlName) & UrlName != "" ~ paste0(
          "https://www.who.int/emergencies/disease-outbreak-news/item/",
          UrlName
        ),
        TRUE ~ NA_character_
      )
    )
}

if (!"publication_datetime_utc" %in% names(records)) {
  records <- records %>%
    mutate(
      publication_datetime_utc = as.POSIXct(PublicationDateAndTime, tz = "UTC")
    )
}

if (!"record_id" %in% names(records)) {
  records <- records %>%
    mutate(record_id = dplyr::coalesce(DonId, Id))
}

clean_path <- who_don_output_path("who_don_records_clean.csv")
inspection_path <- who_don_output_path("who_don_records_inspection.csv")

readr::write_csv(records, clean_path, na = "")

inspection <- records %>%
  transmute(
    DonId,
    Title,
    PublicationDateAndTime,
    publication_datetime_utc,
    FormattedDate,
    article_url,
    summary_text = truncate_text(summary_text),
    overview_text = truncate_text(overview_text),
    epidemiology_text = truncate_text(epidemiology_text),
    assessment_text = truncate_text(assessment_text),
    response_text = truncate_text(response_text),
    advice_text = truncate_text(advice_text),
    further_information_text = truncate_text(further_information_text)
  )

readr::write_csv(inspection, inspection_path, na = "")

message_if("Saved cleaned WHO DON records to: ", clean_path, verbose = verbose)
message_if(
  "Saved inspection-friendly WHO DON records to: ",
  inspection_path,
  verbose = verbose
)
