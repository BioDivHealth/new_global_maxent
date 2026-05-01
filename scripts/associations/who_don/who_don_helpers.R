# ------------------------------------------------------------------------------
# who_don_helpers.R
# ------------------------------------------------------------------------------
# Purpose: Shared helper functions for the WHO Disease Outbreak News (DON)
#          country-extraction workflow.
# ------------------------------------------------------------------------------

who_don_output_dir <- function() {
  here::here("pathogen_association_data", "WHO", "disease_outbreak_news")
}

who_don_output_path <- function(...) {
  file.path(who_don_output_dir(), ...)
}

message_if <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

clean_scalar_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\r", "\n")
  x <- stringr::str_replace_all(x, "\u00a0", " ")
  x <- stringr::str_replace_all(x, "[ \t]+", " ")
  x <- stringr::str_replace_all(x, "\n{2,}", "\n\n")
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NULL")] <- NA_character_
  x
}

normalize_whitespace <- function(x) {
  x <- clean_scalar_text(x)
  x <- stringr::str_replace_all(x, "\n+", " ")
  x <- stringr::str_replace_all(x, "[ ]{2,}", " ")
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NULL")] <- NA_character_
  x
}

strip_inline_media_payloads <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(
    x,
    stringr::regex("data:image/[^;]+;base64,[A-Za-z0-9+/=\\r\\n]+", ignore_case = TRUE),
    " "
  )
  x <- stringr::str_replace_all(
    x,
    stringr::regex("<img\\b[^>]*>", ignore_case = TRUE),
    " "
  )
  x <- stringr::str_replace_all(
    x,
    stringr::regex("<figure\\b[^>]*>.*?</figure>", ignore_case = TRUE, dotall = TRUE),
    " "
  )
  x <- stringr::str_replace_all(
    x,
    stringr::regex("<svg\\b[^>]*>.*?</svg>", ignore_case = TRUE, dotall = TRUE),
    " "
  )
  x
}

html_to_text <- function(x) {
  x <- clean_scalar_text(x)

  if (all(is.na(x))) {
    return(x)
  }

  vapply(
    x,
    FUN.VALUE = character(1),
    function(one_value) {
      if (is.na(one_value)) {
        return(NA_character_)
      }

      one_value <- strip_inline_media_payloads(one_value)
      one_value <- stringr::str_replace_all(
        one_value,
        stringr::regex("<br\\s*/?>", ignore_case = TRUE),
        "\n"
      )
      one_value <- stringr::str_replace_all(
        one_value,
        stringr::regex("</p>|</div>|</li>|</ul>|</ol>|</h[1-6]>", ignore_case = TRUE),
        "\n"
      )
      one_value <- stringr::str_replace_all(
        one_value,
        stringr::regex("<li\\b[^>]*>", ignore_case = TRUE),
        "\n- "
      )

      parsed <- tryCatch(
        xml2::read_html(paste0("<body>", one_value, "</body>")),
        error = function(e) NULL
      )

      plain <- if (!is.null(parsed)) {
        xml2::xml_text(xml2::xml_find_first(parsed, ".//body"), trim = FALSE)
      } else {
        stringr::str_replace_all(one_value, stringr::regex("<[^>]+>"), " ")
      }

      normalize_whitespace(plain)
    }
  )
}

truncate_text <- function(x, max_chars = 2000L) {
  x <- clean_scalar_text(x)

  vapply(
    x,
    FUN.VALUE = character(1),
    function(one_value) {
      if (is.na(one_value)) {
        return(NA_character_)
      }

      if (nchar(one_value) <= max_chars) {
        return(one_value)
      }

      paste0(substr(one_value, 1L, max_chars), "... [truncated]")
    }
  )
}

normalize_match_text <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = " ")
  x <- tolower(x)
  x <- stringr::str_replace_all(x, "[\u2012\u2013\u2014\u2212]", "-")
  x <- stringr::str_replace_all(x, "&", " and ")
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", " ")
  x <- stringr::str_replace_all(x, "\\s{2,}", " ")
  x <- stringr::str_trim(x)
  x[x %in% c("", "na", "null")] <- NA_character_
  x
}

split_semicolon_values <- function(x) {
  x <- clean_scalar_text(x)

  purrr::map(
    x,
    function(one_value) {
      if (is.na(one_value)) {
        return(character(0))
      }

      parts <- stringr::str_split(one_value, "\\s*;\\s*")[[1]]
      parts <- clean_scalar_text(parts)
      stats::na.omit(parts)
    }
  )
}

expand_semicolon_rows <- function(df, column) {
  column <- rlang::ensym(column)
  split_col <- paste0(rlang::as_string(column), "_split")

  df %>%
    dplyr::mutate(!!split_col := split_semicolon_values(!!column)) %>%
    tidyr::unnest(!!rlang::sym(split_col), keep_empty = FALSE) %>%
    dplyr::mutate(!!column := !!rlang::sym(split_col)) %>%
    dplyr::select(-!!rlang::sym(split_col))
}

build_phrase_boundary_pattern <- function(phrase_key) {
  if (is.na(phrase_key) || phrase_key == "") {
    return(NA_character_)
  }

  paste0("(?<![a-z0-9])", escape_regex(phrase_key), "(?![a-z0-9])")
}

detect_phrase_hits <- function(text, phrase_tbl, text_col = "alias_key") {
  text <- clean_scalar_text(text)

  if (is.na(text) || !nrow(phrase_tbl)) {
    return(phrase_tbl[0, , drop = FALSE])
  }

  text_key <- normalize_match_text(text)

  if (is.na(text_key) || text_key == "") {
    return(phrase_tbl[0, , drop = FALSE])
  }

  phrase_tbl %>%
    dplyr::filter(!is.na(.data[[text_col]]), .data[[text_col]] != "") %>%
    dplyr::mutate(
      pattern = vapply(.data[[text_col]], build_phrase_boundary_pattern, character(1)),
      is_match = stringr::str_detect(text_key, stringr::regex(pattern))
    ) %>%
    dplyr::filter(is_match) %>%
    dplyr::select(-pattern, -is_match)
}

who_don_base_countries <- function() {
  strsplit(
    paste(
      "Afghanistan|Albania|Algeria|Andorra|Angola|Antigua and Barbuda|Argentina|Armenia|Australia|Austria|Azerbaijan|Bahamas|Bahrain|Bangladesh|Barbados|Belarus|Belgium|Belize|Benin|Bhutan|Bolivia|Bosnia and Herzegovina|Botswana|Brazil|Brunei|Bulgaria|Burkina Faso|Burundi|Cabo Verde|Cambodia|Cameroon|Canada|Central African Republic|Chad|Chile|China|Colombia|Comoros|Costa Rica|Croatia|Cuba|Cyprus|Czechia|Denmark|Djibouti|Dominica|Dominican Republic|Ecuador|Egypt|El Salvador|Equatorial Guinea|Eritrea|Estonia|Eswatini|Ethiopia|Fiji|Finland|France|Gabon|Gambia|Georgia|Germany|Ghana|Greece|Grenada|Guatemala|Guinea|Guinea-Bissau|Guyana|Haiti|Honduras|Hungary|Iceland|India|Indonesia|Iran|Iraq|Ireland|Israel|Italy|Jamaica|Japan|Jordan|Kazakhstan|Kenya|Kiribati|Kuwait|Kyrgyzstan|Laos|Latvia|Lebanon|Lesotho|Liberia|Libya|Liechtenstein|Lithuania|Luxembourg|Madagascar|Malawi|Malaysia|Maldives|Mali|Malta|Marshall Islands|Mauritania|Mauritius|Mexico|Micronesia|Moldova|Monaco|Mongolia|Montenegro|Morocco|Mozambique|Myanmar|Namibia|Nauru|Nepal|Netherlands|New Zealand|Nicaragua|Niger|Nigeria|North Korea|North Macedonia|Norway|Oman|Pakistan|Palau|Panama|Papua New Guinea|Paraguay|Peru|Philippines|Poland|Portugal|Qatar|Romania|Russia|Rwanda|Saint Kitts and Nevis|Saint Lucia|Saint Vincent and the Grenadines|Samoa|San Marino|Sao Tome and Principe|Saudi Arabia|Senegal|Serbia|Seychelles|Sierra Leone|Singapore|Slovakia|Slovenia|Solomon Islands|Somalia|South Africa|South Korea|South Sudan|Spain|Sri Lanka|Sudan|Suriname|Sweden|Switzerland|Syria|Tajikistan|Tanzania|Thailand|Timor-Leste|Togo|Tonga|Trinidad and Tobago|Tunisia|Turkey|Turkmenistan|Tuvalu|Uganda|Ukraine|United Arab Emirates|United Kingdom|United States|Uruguay|Uzbekistan|Vanuatu|Venezuela|Vietnam|Yemen|Zambia|Zimbabwe|Palestine|Holy See|Kosovo",
        sep = ""
      ),
      "\\|"
    )[[1]]
}

load_who_don_geography_gazetteer <- function() {
  base_countries <- tibble::tibble(
    alias = who_don_base_countries(),
    country_standard = who_don_base_countries(),
    geography_type = "country",
    is_ambiguous = FALSE
  )

  special_path <- here::here(
    "scripts",
    "associations",
    "who_don",
    "who_don_country_aliases.csv"
  )

  special_aliases <- readr::read_csv(
    special_path,
    show_col_types = FALSE,
    na = c("", "NA")
  )

  dplyr::bind_rows(base_countries, special_aliases) %>%
    dplyr::mutate(
      alias = clean_scalar_text(alias),
      country_standard = clean_scalar_text(country_standard),
      geography_type = dplyr::coalesce(clean_scalar_text(geography_type), "country"),
      is_ambiguous = dplyr::coalesce(as.logical(is_ambiguous), FALSE),
      alias_key = normalize_match_text(alias),
      country_key = normalize_match_text(country_standard)
    ) %>%
    dplyr::filter(!is.na(alias_key)) %>%
    dplyr::distinct(alias_key, country_standard, geography_type, is_ambiguous, .keep_all = TRUE)
}

generic_title_labels <- function() {
  tibble::tribble(
    ~label, ~label_type,
    "Global situation", "generic_global",
    "Global Situation", "generic_global",
    "Global update", "generic_global",
    "Global Update", "generic_global",
    "Multi-country", "generic_multicountry",
    "Multi country", "generic_multicountry",
    "Regional situation report", "generic_regional",
    "Situation report", "generic_report",
    "Epidemic update", "generic_report",
    "Northern Hemisphere", "generic_region"
  ) %>%
    dplyr::mutate(label_key = normalize_match_text(label))
}

who_don_section_columns <- function() {
  c(
    "summary_text",
    "overview_text",
    "epidemiology_text",
    "response_text",
    "assessment_text",
    "advice_text",
    "further_information_text"
  )
}

who_don_record_key <- function(record_id, DonId = NULL, Id = NULL) {
  dplyr::coalesce(record_id, DonId, Id)
}

who_don_with_record_key <- function(df) {
  record_id_col <- if ("record_id" %in% names(df)) df$record_id else NULL
  DonId_col <- if ("DonId" %in% names(df)) df$DonId else NULL
  Id_col <- if ("Id" %in% names(df)) df$Id else NULL

  df %>%
    dplyr::mutate(record_key = who_don_record_key(record_id_col, DonId_col, Id_col))
}

who_don_confidence_rank <- function(confidence) {
  dplyr::case_when(
    confidence == "high" ~ 1L,
    confidence == "medium" ~ 2L,
    confidence == "low" ~ 3L,
    TRUE ~ 9L
  )
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x, perl = TRUE)
}

build_country_match_pattern <- function(gazetteer) {
  keys <- gazetteer %>%
    dplyr::filter(
      geography_type == "country",
      !is_ambiguous,
      !is.na(alias_key)
    ) %>%
    dplyr::pull(alias_key) %>%
    unique()

  keys <- keys[order(nchar(keys), decreasing = TRUE)]

  paste0(
    "(?<![a-z])(?:",
    paste(escape_regex(keys), collapse = "|"),
    ")(?![a-z])"
  )
}

who_don_event_patterns <- function() {
  c(
    "reported in",
    "reported from",
    "identified in",
    "detected in",
    "confirmed in",
    "outbreak in",
    "case in",
    "cases in",
    "notified by",
    "notified who of",
    "ministry of health of",
    "national ihr focal point for",
    "declared an outbreak in",
    "declared the outbreak in"
  )
}

who_don_first_sentences <- function(text, n_sentences = 2L) {
  text <- normalize_whitespace(text)

  vapply(
    text,
    FUN.VALUE = character(1),
    function(one_text) {
      if (is.na(one_text) || one_text == "") {
        return(NA_character_)
      }

      parts <- stringr::str_split(
        one_text,
        stringr::regex("(?<=[.!?])\\s+", multiline = TRUE)
      )[[1]]
      parts <- clean_scalar_text(parts)
      parts <- parts[!is.na(parts)]

      if (length(parts) == 0) {
        return(NA_character_)
      }

      paste(parts[seq_len(min(length(parts), n_sentences))], collapse = " ")
    }
  )
}

who_don_article_slug <- function(article_url) {
  article_url <- clean_scalar_text(article_url)

  vapply(
    article_url,
    FUN.VALUE = character(1),
    function(one_url) {
      if (is.na(one_url) || one_url == "") {
        return(NA_character_)
      }

      slug <- basename(one_url)
      slug <- stringr::str_replace(slug, "\\?.*$", "")
      clean_scalar_text(slug)
    }
  )
}

who_don_generic_region_patterns <- function() {
  c(
    "west africa",
    "central africa",
    "north africa",
    "southern africa",
    "east africa",
    "region of the americas",
    "african region",
    "european region",
    "south-east asia region",
    "eastern mediterranean region",
    "western pacific region",
    "northern hemisphere",
    "global situation",
    "global update",
    "multi-country",
    "multi country"
  )
}

who_don_region_hint <- function(title, title_suffix = NULL) {
  txt <- dplyr::coalesce(title_suffix, title)
  txt_norm <- normalize_match_text(txt)
  pats <- who_don_generic_region_patterns()

  vapply(
    txt_norm,
    FUN.VALUE = character(1),
    function(one_text) {
      if (is.na(one_text) || one_text == "") {
        return(NA_character_)
      }

      hit <- pats[stringr::str_detect(one_text, stringr::fixed(pats))]

      if (length(hit) == 0) {
        return(NA_character_)
      }

      hit[[1]]
    }
  )
}

who_don_title_suffix <- function(title) {
  title <- clean_scalar_text(title)

  suffix <- stringr::str_match(
    title,
    ".*(?:\\s+[-–—:]\\s*|[-–—]\\s+|-(?=[A-Z]))(.+)$"
  )[, 2]

  clean_scalar_text(suffix)
}

who_don_match_country_alias <- function(raw_text, gazetteer) {
  raw_text <- clean_scalar_text(raw_text)

  if (is.na(raw_text)) {
    return(tibble::tibble())
  }

  key <- normalize_match_text(raw_text)

  gazetteer %>%
    dplyr::filter(
      alias_key == key,
      geography_type == "country",
      !is_ambiguous
    ) %>%
    dplyr::slice(1) %>%
    dplyr::transmute(
      country = raw_text,
      country_standard,
      alias = raw_text,
      alias_key = key
    )
}

extract_context_window <- function(text, start, end, pad = 120L) {
  left <- max(1L, start - pad)
  right <- min(nchar(text), end + pad)
  substr(text, left, right)
}

read_who_don_records <- function(filename = "who_don_records.csv") {
  readr::read_csv(
    who_don_output_path(filename),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

safe_collapse_unique <- function(x, sep = "; ") {
  x <- unique(stats::na.omit(x))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = sep)
}
