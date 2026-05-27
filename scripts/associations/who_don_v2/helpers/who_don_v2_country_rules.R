library(dplyr)
library(stringr)
library(tidyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_io.R"))

v2_country_alias_input_cols <- c(
  "country.name.en",
  "iso.name.en",
  "un.name.en",
  "cldr.name.en",
  "cldr.variant.en",
  "cow.name",
  "vdem.name"
)

v2_country_custom_aliases <- function() {
  tibble::tribble(
    ~alias, ~country_standard, ~alias_type, ~is_ambiguous, ~priority, ~notes,
    "Bolivia (Plurinational State of)", "Bolivia", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Brunei Darussalam", "Brunei", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Cabo Verde", "Cabo Verde", "who_variant", FALSE, 10L, "WHO standard spelling.",
    "Cape Verde", "Cabo Verde", "common_variant", FALSE, 30L, "Common historical spelling.",
    "Cap Verde", "Cabo Verde", "historical_variant", FALSE, 30L, "Historical spelling seen in accepted WHO DON evidence.",
    "Cote d'Ivoire", "Cote d'Ivoire", "ascii_variant", FALSE, 10L, "ASCII spelling used in current clean layer.",
    "Cote d`Ivoire", "Cote d'Ivoire", "ascii_variant", FALSE, 20L, "Backtick spelling seen in accepted WHO DON evidence.",
    "Cote d’Ivoire", "Cote d'Ivoire", "unicode_variant", FALSE, 20L, "ASCII Cote with curly apostrophe seen in accepted WHO DON evidence.",
    "Côte d’Ivoire", "Cote d'Ivoire", "unicode_variant", FALSE, 20L, "Unicode spelling with curly apostrophe.",
    "Côte d'Ivoire", "Cote d'Ivoire", "unicode_variant", FALSE, 20L, "Unicode spelling with straight apostrophe.",
    "Côte d'Ivoir", "Cote d'Ivoire", "accepted_typo_variant", FALSE, 30L, "Truncated spelling seen in accepted WHO DON evidence.",
    "Bosnia and Hezegovina", "Bosnia and Herzegovina", "accepted_typo_variant", FALSE, 30L, "Typo seen in accepted WHO DON evidence.",
    "Democratic Republic of Congo", "Democratic Republic of the Congo", "common_variant", FALSE, 20L, "Common shortened form.",
    "Democratic Republic of the Congo", "Democratic Republic of the Congo", "who_variant", FALSE, 10L, "WHO standard spelling.",
    "DRC", "Democratic Republic of the Congo", "abbreviation", TRUE, 80L, "Abbreviation can be ambiguous; keep reviewable.",
    "Congo-Kinshasa", "Democratic Republic of the Congo", "common_variant", FALSE, 40L, "Common disambiguating form.",
    "Zaire", "Democratic Republic of the Congo", "historical_variant", FALSE, 30L, "Historical country name seen in accepted WHO DON evidence.",
    "Republic of the Congo", "Republic of the Congo", "who_variant", FALSE, 10L, "WHO/UN style country name.",
    "Republic of Congo", "Republic of the Congo", "accepted_variant", FALSE, 20L, "Accepted WHO DON evidence sometimes omits 'the'.",
    "Congo-Brazzaville", "Republic of the Congo", "common_variant", FALSE, 40L, "Common disambiguating form.",
    "Gabonese", "Gabon", "demonym_event_context", FALSE, 45L, "Demonym seen in accepted Ministry of Health event evidence.",
    "Hong Kong", "China", "subnational_clean_compat", FALSE, 35L, "Accepted clean layer maps Hong Kong evidence to China.",
    "Hong Kong SAR", "China", "subnational_clean_compat", FALSE, 35L, "Accepted clean layer maps Hong Kong evidence to China.",
    "Islamic Republic of Iran", "Iran", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Iran (Islamic Republic of)", "Iran", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Lao People's Democratic Republic", "Laos", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Libyan Arab Jamahiriya", "Libya", "historical_variant", FALSE, 30L, "Historical WHO country name seen in accepted evidence.",
    "Luxemburg", "Luxembourg", "accepted_typo_variant", FALSE, 30L, "Historical spelling seen in accepted WHO DON evidence.",
    "Moldova, Republic of", "Moldova", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Burma", "Myanmar", "historical_variant", FALSE, 30L, "Historical country name.",
    "Russia", "Russia", "common_variant", FALSE, 20L, "Short common form.",
    "Russian Federation", "Russia", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "South Korea", "South Korea", "common_variant", FALSE, 20L, "Short common form.",
    "Republic of Korea", "South Korea", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Korea, Republic of", "South Korea", "who_variant", FALSE, 20L, "WHO/UN style country name seen in accepted evidence.",
    "North Korea", "North Korea", "common_variant", FALSE, 20L, "Short common form.",
    "Democratic People's Republic of Korea", "North Korea", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "DPR Korea", "North Korea", "who_variant", FALSE, 25L, "WHO short form seen in accepted evidence.",
    "Kosovo", "Kosovo", "who_variant", FALSE, 20L, "Current accepted layer country label.",
    "Micronesia", "Micronesia", "common_variant", FALSE, 20L, "Short common form used in current accepted layer.",
    "Micronesia (Federated States of)", "Micronesia", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Palestine", "Palestine", "common_variant", FALSE, 20L, "Short common form used in current accepted layer.",
    "State of Palestine", "Palestine", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Gaza Strip", "Palestine", "subnational_clean_compat", FALSE, 35L, "Accepted clean layer maps Gaza Strip event evidence to Palestine.",
    "West Bank and Gaza Strip", "Palestine", "subnational_clean_compat", FALSE, 35L, "Accepted clean layer maps West Bank and Gaza Strip evidence to Palestine.",
    "Syria", "Syria", "common_variant", FALSE, 20L, "Short common form.",
    "Syrian Arab Republic", "Syria", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Tanzania", "Tanzania", "common_variant", FALSE, 20L, "Short common form.",
    "United Republic of Tanzania", "Tanzania", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "United Kingdom", "United Kingdom", "common_variant", FALSE, 10L, "Common short form.",
    "United Kingdom of Great Britain and Northern Ireland", "United Kingdom", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "UK", "United Kingdom", "abbreviation", TRUE, 80L, "Abbreviation can appear outside country evidence.",
    "United States", "United States", "common_variant", FALSE, 10L, "Common short form used in current accepted layer.",
    "United States of America", "United States", "who_variant", FALSE, 10L, "WHO standard spelling mapped to current accepted layer.",
    "USA", "United States", "abbreviation", TRUE, 80L, "Abbreviation can appear outside country evidence.",
    "U.S.", "United States", "abbreviation", TRUE, 80L, "Abbreviation can appear outside country evidence.",
    "U.S.A.", "United States", "abbreviation", TRUE, 80L, "Abbreviation can appear outside country evidence.",
    "U S", "United States", "accepted_spaced_abbreviation", TRUE, 85L, "Spaced abbreviation appears in accepted legacy evidence; keep reviewable.",
    "US Centers for Disease Control", "United States", "institution_context", FALSE, 60L, "US CDC wording appears in accepted legacy evidence.",
    "Venezuela", "Venezuela", "common_variant", FALSE, 20L, "Short common form.",
    "Venezuela (Bolivarian Republic of)", "Venezuela", "who_variant", FALSE, 20L, "WHO/UN style country name.",
    "Viet Nam", "Vietnam", "who_variant", FALSE, 10L, "WHO spelling mapped to current accepted layer.",
    "Vietnam", "Vietnam", "common_variant", FALSE, 20L, "Common English spelling.",
    "Holland", "Netherlands", "historical_common_variant", FALSE, 30L, "Historical/common name seen in accepted WHO DON evidence.",
    "Tadjikistan", "Tajikistan", "accepted_typo_variant", FALSE, 30L, "Historical spelling seen in accepted WHO DON evidence.",
    "KSA", "Saudi Arabia", "abbreviation", TRUE, 80L, "Abbreviation seen in accepted evidence; keep reviewable outside strong context.",
    "Antigua/Barbuda", "Antigua and Barbuda", "slash_variant", FALSE, 30L, "Slash form seen in accepted H1N1 map/table evidence.",
    "St Kitts/Nevis", "Saint Kitts and Nevis", "slash_variant", FALSE, 30L, "Slash form seen in accepted H1N1 map/table evidence.",
    "St Vincent/Grenadines", "Saint Vincent and the Grenadines", "slash_variant", FALSE, 30L, "Slash form seen in accepted H1N1 map/table evidence.",
    "Trinidad/Tobago", "Trinidad and Tobago", "slash_variant", FALSE, 30L, "Slash form seen in accepted H1N1 map/table evidence."
  )
}

v2_default_country_aliases <- function() {
  if (!requireNamespace("countrycode", quietly = TRUE)) {
    stop("Package 'countrycode' is required to build native country aliases.", call. = FALSE)
  }

  accepted_standards <- v2_country_custom_aliases() %>%
    distinct(country_standard) %>%
    pull(country_standard)

  accepted_standards <- union(
    accepted_standards,
    countrycode::codelist %>%
      filter(!is.na(un), un != "", !is.na(country.name.en), country.name.en != "") %>%
      pull(country.name.en)
  )

  excluded_standards <- c(
    "Antarctica",
    "Congo - Brazzaville",
    "Congo - Kinshasa",
    "Global",
    "Hong Kong SAR China",
    "Unknown"
  )

  countrycode::codelist %>%
    filter(!is.na(un), un != "", !is.na(country.name.en), country.name.en != "") %>%
    select(any_of(v2_country_alias_input_cols)) %>%
    pivot_longer(everything(), names_to = "alias_type", values_to = "alias") %>%
    filter(!is.na(alias), alias != "") %>%
    mutate(
      country_standard = countrycode::countrycode(
        alias,
        origin = "country.name",
        destination = "country.name",
        warn = FALSE
      ),
      country_standard = coalesce(country_standard, alias),
      alias_type = paste0("countrycode_", alias_type),
      is_ambiguous = FALSE,
      priority = 50L,
      notes = "Generated from countrycode English country-name fields."
    ) %>%
    bind_rows(v2_country_custom_aliases()) %>%
    mutate(
      country_standard = recode(
        country_standard,
        "Congo - Brazzaville" = "Republic of the Congo",
        "Congo - Kinshasa" = "Democratic Republic of the Congo",
        "Antigua & Barbuda" = "Antigua and Barbuda",
        "Bosnia & Herzegovina" = "Bosnia and Herzegovina",
        "Cape Verde" = "Cabo Verde",
        "Côte d’Ivoire" = "Cote d'Ivoire",
        "Côte d'Ivoire" = "Cote d'Ivoire",
        "Micronesia (Federated States of)" = "Micronesia",
        "Myanmar (Burma)" = "Myanmar",
        "Palestinian Territories" = "Palestine",
        "State of Palestine" = "Palestine",
        "St. Lucia" = "Saint Lucia",
        "St. Kitts & Nevis" = "Saint Kitts and Nevis",
        "St. Vincent & Grenadines" = "Saint Vincent and the Grenadines",
        "São Tomé & Príncipe" = "Sao Tome and Principe",
        "Trinidad & Tobago" = "Trinidad and Tobago",
        "United States of America" = "United States",
        "Viet Nam" = "Vietnam"
      ),
      alias = str_squish(alias),
      country_standard = str_squish(country_standard),
      priority = if_else(alias == country_standard, pmin(priority, 10L), priority)
    ) %>%
    filter(
      alias != "",
      country_standard != "",
      country_standard %in% accepted_standards,
      !country_standard %in% excluded_standards,
      country_standard != "Viet Nam",
      !str_to_lower(alias) %in% c("global", "unknown")
    ) %>%
    arrange(priority, desc(nchar(alias)), alias) %>%
    distinct(alias, country_standard, .keep_all = TRUE)
}

v2_prepare_country_aliases <- function(path = who_don_v2_rules_dir("country_aliases.csv"), force = FALSE) {
  existing <- if (file.exists(path)) v2_read_csv(path) else tibble::tibble()
  if (nrow(existing) > 0 && !force) {
    return(existing)
  }

  aliases <- v2_default_country_aliases()
  v2_write_csv(aliases, path)
  aliases
}

v2_validate_country_aliases <- function(country_aliases) {
  required <- c("alias", "country_standard", "alias_type", "is_ambiguous", "priority", "notes")
  missing_cols <- setdiff(required, names(country_aliases))
  if (length(missing_cols) > 0) {
    return(tibble::tibble(
      severity = "blocking",
      issue = "missing_required_columns",
      detail = paste(missing_cols, collapse = ", ")
    ))
  }

  bind_rows(
    country_aliases %>%
      filter(is.na(alias) | str_squish(alias) == "") %>%
      transmute(severity = "blocking", issue = "blank_alias", detail = country_standard),
    country_aliases %>%
      filter(is.na(country_standard) | str_squish(country_standard) == "") %>%
      transmute(severity = "blocking", issue = "blank_country_standard", detail = alias),
    country_aliases %>%
      count(alias, name = "rows") %>%
      filter(rows > 1) %>%
      transmute(severity = "review", issue = "duplicate_alias", detail = alias)
  )
}

v2_default_country_policy_decisions <- function() {
  tibble::tribble(
    ~diff_category, ~policy_decision, ~policy_review_priority, ~policy_note,
    "exact_record_country_match", "accept_native", "low", "Record-country pair matches the accepted country layer and has native evidence.",
    "accepted_missing_native", "accept_legacy_exception", "medium", "Accepted country was not recovered by native extraction; retain as explicit legacy exception.",
    "native_new_country_candidate", "reject_native_unreviewed", "high", "Native-only country candidate is not adopted until reviewed or covered by policy."
  )
}

v2_prepare_country_policy_decisions <- function(
  path = who_don_v2_rules_dir("country_candidate_policy_decisions.csv"),
  force = FALSE
) {
  existing <- if (file.exists(path)) v2_read_csv(path) else tibble::tibble()
  if (nrow(existing) > 0 && !force) {
    return(existing)
  }

  policy <- v2_default_country_policy_decisions()
  v2_write_csv(policy, path)
  policy
}
