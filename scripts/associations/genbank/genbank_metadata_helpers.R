# ------------------------------------------------------------------------------
# GenBank metadata helpers
# ------------------------------------------------------------------------------
# Shared helper functions for the GenBank search, fetch, and source-selection
# workflow used by `5_7_GenBank_Pathogen_Country_Metadata.R` and
# `5_7a_GenBank_Search_Diagnostics.R`.
# ------------------------------------------------------------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_unique <- function(x) {
  x <- clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

normalize_lookup_key <- function(x) {
  x <- clean_text(x)

  if (length(x) == 0) {
    return(character(0))
  }

  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "")
  x[x == ""] <- NA_character_
  x
}

pathogen_name_fuzzy_match <- function(candidate, target) {
  candidate <- clean_text(candidate)
  target <- clean_text(target)

  if (is.na(candidate) || is.na(target)) {
    return(FALSE)
  }

  candidate_key <- normalize_lookup_key(candidate)
  target_key <- normalize_lookup_key(target)

  if (is.na(candidate_key) || is.na(target_key)) {
    return(FALSE)
  }

  identical(candidate_key, target_key) ||
    stringr::str_detect(candidate_key, stringr::fixed(target_key)) ||
    stringr::str_detect(target_key, stringr::fixed(candidate_key))
}

split_multi_value <- function(x, pattern = ";") {
  x <- clean_text(x)

  if (length(x) == 0 || all(is.na(x))) {
    return(character(0))
  }

  x <- stats::na.omit(x)
  pieces <- stringr::str_split(x, pattern = pattern)
  pieces <- unlist(pieces, use.names = FALSE)
  pieces <- clean_text(pieces)
  pieces <- pieces[!is.na(pieces)]
  unique(pieces)
}

parse_env_integer <- function(name, default = NA_integer_) {
  value <- clean_text(Sys.getenv(name, unset = NA_character_))

  if (is.na(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))

  if (is.na(parsed)) {
    stop("Environment variable ", name, " must be an integer when set.")
  }

  parsed
}

parse_env_flag <- function(name, default = FALSE) {
  value <- clean_text(Sys.getenv(name, unset = NA_character_))

  if (is.na(value)) {
    return(default)
  }

  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

log_progress <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  cat(..., "\n", sep = "")
  invisible(NULL)
}

read_dotenv_value <- function(path, key) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  lines <- readLines(path, warn = FALSE)
  key_pattern <- paste0("^\\s*", key, "\\s*=\\s*")
  matches <- lines[stringr::str_detect(lines, key_pattern)]

  if (length(matches) == 0) {
    return(NA_character_)
  }

  value <- stringr::str_replace(matches[[1]], key_pattern, "")
  value <- stringr::str_replace(value, "\\s+#.*$", "")
  value <- clean_text(value)

  if (is.na(value)) {
    return(NA_character_)
  }

  if (
    stringr::str_length(value) >= 2 &&
      (
        (stringr::str_starts(value, "\"") && stringr::str_ends(value, "\"")) ||
          (stringr::str_starts(value, "'") && stringr::str_ends(value, "'"))
      )
  ) {
    value <- stringr::str_sub(value, 2, -2)
  }

  clean_text(value)
}

quote_organism_term <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(character(0))
  }

  paste0("\"", x, "\"[Organism]")
}

quote_all_fields_term <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(character(0))
  }

  paste0("\"", x, "\"[All Fields]")
}

looks_like_disease_label <- function(x) {
  x <- clean_text(x)

  if (is.na(x)) {
    return(FALSE)
  }

  stringr::str_detect(
    x,
    stringr::regex(
      "disease|fever|infection|syndrome|gastroenteritis|encephalitis|respiratory|hemorrhagic|haemorrhagic|aids|cholera|plague|monkeypox|vaccinia|poliomyelitis|hepatitis|influenza$",
      ignore_case = TRUE
    )
  )
}

token_count <- function(x) {
  x <- clean_text(x)

  if (is.na(x)) {
    return(0L)
  }

  tokens <- stringr::str_split(x, "\\s+")[[1]]
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  length(tokens)
}

filter_organism_aliases <- function(
  aliases,
  target_pathogen = NA_character_,
  disease_name = NA_character_,
  keep_broad_aliases = FALSE
) {
  aliases <- clean_text(aliases)
  aliases <- aliases[!is.na(aliases)]

  if (length(aliases) == 0) {
    return(character(0))
  }

  target_pathogen <- clean_text(target_pathogen)
  disease_name <- clean_text(disease_name)
  target_tokens <- token_count(target_pathogen)

  if (!is.na(disease_name)) {
    aliases <- aliases[aliases != disease_name]
  }
  aliases <- aliases[!purrr::map_lgl(aliases, looks_like_disease_label)]

  if (!isTRUE(keep_broad_aliases) && !is.na(target_pathogen)) {
    aliases <- aliases[
      !purrr::map_lgl(
        aliases,
        function(alias) {
          alias_tokens <- token_count(alias)

          alias_tokens > 0L &&
            target_tokens >= 3L &&
            alias_tokens <= 2L &&
            stringr::str_detect(
              alias,
              stringr::regex(
                "adenovirus|influenza|virus|vibrio|salmonella|shigella|enterovirus|poliovirus|picobirnavirus|parvovirus",
                ignore_case = TRUE
              )
            )
        }
      )
    ]
  }

  unique(aliases)
}

manual_query_aliases <- function(pathogen, disease_name = NA_character_) {
  pathogen <- clean_text(pathogen)

  aliases <- character(0)

  if (identical(pathogen, "Mamastrovirus 9 (GII.B-human)")) {
    aliases <- c("Mamastrovirus 9", "Mamastrovirus virginiaense")
  } else if (identical(pathogen, "Carivore protoparvoviruses (CPV)")) {
    aliases <- c("Carnivore protoparvoviruses (CPV)", "Protoparvovirus carnivoran")
  } else if (identical(pathogen, "Genus Vesiculovirus")) {
    aliases <- c("Vesiculovirus", "Vesicular stomatitis virus")
  } else if (identical(pathogen, "Genus Rotavirus")) {
    aliases <- c("Rotavirus", "Rotavirus A")
  }

  aliases <- clean_text(aliases)
  aliases[!is.na(aliases)]
}

is_influenza_network_target <- function(pathogen) {
  identical(clean_text(pathogen), "Alphainfluenzavirus influenzae")
}

is_sars_covid_network_target <- function(pathogen, disease_name) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)

  if (!is.na(disease_name) && identical(disease_name, "Severe Acute Respiratory Syndrome (SARS); COVID-19")) {
    return(TRUE)
  }

  pathogen %in% c(
    "Betacoronavirus pandemicum",
    "Severe acute respiratory syndrome-related coronavirus"
  )
}

extract_influenza_subtype_marker <- function(pathogen) {
  pathogen <- clean_text(pathogen)

  if (is.na(pathogen)) {
    return(NA_character_)
  }

  matched <- stringr::str_match(
    pathogen,
    "^Alphainfluenzavirus influenzae \\((H\\d+(?:N\\d+|Nx))\\)$"
  )

  clean_text(matched[, 2])
}

extract_influenza_subtype_token <- function(pathogen) {
  marker <- extract_influenza_subtype_marker(pathogen)

  if (is.na(marker)) {
    return(NA_character_)
  }

  if (stringr::str_detect(marker, "Nx$")) {
    return(clean_text(stringr::str_remove(marker, "Nx$")))
  }

  marker
}

is_influenza_family_fallback <- function(pathogen) {
  marker <- extract_influenza_subtype_marker(pathogen)
  !is.na(marker) && stringr::str_detect(marker, "Nx$")
}

split_influenza_candidate_sets <- function(who_pathogen_candidates) {
  candidates <- split_multi_value(who_pathogen_candidates)
  candidates <- clean_text(candidates)
  candidates <- candidates[!is.na(candidates)]

  subtype_candidates <- candidates[stringr::str_detect(
    candidates,
    "^Alphainfluenzavirus influenzae \\(H\\d+(?:N\\d+|Nx)\\)$"
  )]

  exact_order <- c(
    "Alphainfluenzavirus influenzae (H1N1)",
    "Alphainfluenzavirus influenzae (H2N1)",
    "Alphainfluenzavirus influenzae (H3N1)",
    "Alphainfluenzavirus influenzae (H3N2)",
    "Alphainfluenzavirus influenzae (H5N1)",
    "Alphainfluenzavirus influenzae (H6N1)",
    "Alphainfluenzavirus influenzae (H7N1)",
    "Alphainfluenzavirus influenzae (H10N1)"
  )

  fallback_order <- c(
    "Alphainfluenzavirus influenzae (H2Nx)",
    "Alphainfluenzavirus influenzae (H5Nx)",
    "Alphainfluenzavirus influenzae (H6Nx)",
    "Alphainfluenzavirus influenzae (H7Nx)",
    "Alphainfluenzavirus influenzae (H10Nx)"
  )

  exact_subtypes <- subtype_candidates[subtype_candidates %in% exact_order]
  exact_subtypes <- exact_order[exact_order %in% exact_subtypes]

  family_fallbacks <- subtype_candidates[subtype_candidates %in% fallback_order]
  family_fallbacks <- fallback_order[fallback_order %in% family_fallbacks]

  list(
    exact_subtypes = clean_text(exact_subtypes),
    family_fallbacks = clean_text(family_fallbacks)
  )
}

manual_genbank_queries <- function(pathogen, disease_name = NA_character_) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)

  queries <- character(0)

  if (identical(pathogen, "Human mastadenovirus B")) {
    queries <- c(
      "\"Human mastadenovirus B\"[All Fields]",
      "\"Human adenovirus B\"[All Fields]",
      "\"Human adenovirus\"[All Fields]"
    )
  } else if (identical(pathogen, "Mastadenovirus blackbeardi")) {
    queries <- c(
      "\"Mastadenovirus blackbeardi serotype 14\"[All Fields]",
      "\"Human adenovirus 14\"[All Fields]",
      "\"Adenovirus B14\"[All Fields]"
    )
  } else if (identical(pathogen, "Mastadenovirus blackbeardi serotype 14")) {
    queries <- c(
      "\"Human adenovirus 14\"[Organism]",
      "\"Human adenovirus 14\"[All Fields]",
      "\"Adenovirus B14\"[All Fields]"
    )
  } else if (identical(pathogen, "Mammarenavirus juninense")) {
    queries <- c(
      "\"Junin virus\"[Organism]",
      "\"Junin virus\"[All Fields]"
    )
  } else if (identical(pathogen, "Mammarenavirus lassaense")) {
    queries <- c(
      "\"Lassa virus\"[Organism]",
      "\"Lassa mammarenavirus\"[Organism]",
      "\"Lassa virus\"[All Fields]"
    )
  } else if (identical(pathogen, "Mammarenavirus lujoense")) {
    queries <- c(
      "\"Lujo virus\"[All Fields]",
      "\"Lujo mammarenavirus\"[Organism]"
    )
  } else if (identical(pathogen, "Mamastrovirus 9 (GII.B-human)")) {
    queries <- c(
      "\"Mamastrovirus virginiaense\"[Organism]",
      "\"Mamastrovirus 9\"[All Fields]"
    )
  } else if (identical(pathogen, "Salmonella enterica non typhoidal serovars")) {
    queries <- c(
      "\"Salmonella enterica\"[Organism]",
      "\"non-typhoidal salmonellosis\"[All Fields]"
    )
  } else if (identical(pathogen, "Vibrio cholerae") && identical(disease_name, "Cholera")) {
    queries <- c(
      "(\"Vibrio cholerae O1\"[All Fields] OR \"Vibrio cholerae O139\"[All Fields] OR \"Vibrio cholerae serogroup 0139\"[All Fields])",
      "((\"Vibrio cholerae\"[All Fields]) AND (O1[All Fields] OR O139[All Fields] OR \"serogroup 0139\"[All Fields]))"
    )
  } else if (identical(pathogen, "Orthoflavivirus encephalitidis")) {
    queries <- c(
      "\"Tick-borne encephalitis virus\"[Organism]",
      "\"tick-borne encephalitis virus\"[All Fields]",
      "\"TBEV\"[All Fields]"
    )
  } else if (identical(pathogen, "Orthoflavivirus nilense")) {
    queries <- c(
      "\"West Nile virus\"[Organism]",
      "\"west nile virus\"[All Fields]",
      "\"WNV\"[All Fields]"
    )
  } else if (identical(pathogen, "Henipavirus nipahense")) {
    queries <- c(
      "\"Nipah virus\"[Organism]",
      "\"Nipah virus\"[All Fields]"
    )
  } else if (identical(pathogen, "Metapneumovirus hominis")) {
    queries <- c(
      "\"Human metapneumovirus\"[Organism]",
      "\"Human metapneumovirus\"[All Fields]"
    )
  } else if (identical(pathogen, "Human immunodeficiency virus 1 (HIV-1)")) {
    queries <- c(
      "\"Human immunodeficiency virus 1\"[Organism]",
      "\"HIV-1\"[All Fields]"
    )
  } else if (identical(pathogen, "Alphavirus chikungunya")) {
    queries <- c(
      "\"Chikungunya virus\"[Organism]",
      "\"CHIKV\"[All Fields]"
    )
  } else if (identical(pathogen, "Alphavirus venezuelan")) {
    queries <- c(
      "\"Venezuelan equine encephalitis virus\"[Organism]",
      "\"VEEV\"[All Fields]"
    )
  } else if (identical(pathogen, "Betacoronavirus cameli")) {
    queries <- c(
      "\"Middle East respiratory syndrome-related coronavirus\"[Organism]",
      "\"MERS-CoV\"[All Fields]",
      "\"Middle East respiratory syndrome coronavirus\"[All Fields]"
    )
  } else if (identical(pathogen, "Zaire ebolavirus")) {
    queries <- c(
      "\"Zaire ebolavirus\"[Organism]",
      "\"Ebola virus\"[Organism]",
      "\"Ebola virus\"[All Fields]"
    )
  } else if (identical(pathogen, "Alphainfluenzavirus influenzae")) {
    queries <- c(
      paste0(
        "(\"Influenza A virus\"[Organism]) AND (",
        paste(
          c(
            "H1N1[All Fields]", "H2N1[All Fields]", "H3N1[All Fields]",
            "H3N2[All Fields]", "H5N1[All Fields]", "H6N1[All Fields]",
            "H7N1[All Fields]", "H10N1[All Fields]",
            "H2[All Fields]", "H5[All Fields]", "H6[All Fields]",
            "H7[All Fields]", "H10[All Fields]"
          ),
          collapse = " OR "
        ),
        ")"
      ),
      "\"Influenza A virus\"[Organism]"
    )
  } else if (!is.na(extract_influenza_subtype_marker(pathogen))) {
    subtype_token <- extract_influenza_subtype_token(pathogen)

    if (!is.na(subtype_token)) {
      queries <- c(
        paste0("(\"Influenza A virus\"[Organism]) AND (", subtype_token, "[All Fields])")
      )
    }
  }

  queries <- clean_text(queries)
  queries[!is.na(queries)]
}

build_geo_query <- function(base_query) {
  base_query <- clean_text(base_query)

  if (is.na(base_query)) {
    return(NA_character_)
  }

  paste0(
    "(", base_query, ")",
    " AND (\"country\"[All Fields] OR \"geo_loc_name\"[All Fields] OR \"lat_lon\"[All Fields])"
  )
}

manual_query_candidates <- function(pathogen, disease_name = NA_character_) {
  base_queries <- manual_genbank_queries(pathogen, disease_name)
  queries <- c(
    purrr::map_chr(base_queries, build_geo_query),
    base_queries
  )

  if (identical(clean_text(pathogen), "Orthoflavivirus encephalitidis")) {
    queries <- c(
      "(\"TBEV\"[All Fields]) AND (\"country\"[All Fields] OR \"geo_loc_name\"[All Fields] OR \"lat_lon\"[All Fields])",
      "\"tick-borne encephalitis virus\"[All Fields]",
      "(txid11084[Organism:exp]) AND (\"country\"[All Fields] OR \"geo_loc_name\"[All Fields] OR \"lat_lon\"[All Fields])",
      queries
    )
  }

  queries <- clean_text(queries)
  unique(queries[!is.na(queries)])
}

search_query_variants <- function(query) {
  query <- clean_text(query)

  if (is.na(query)) {
    return(NA_character_)
  }

  if (stringr::str_detect(query, "txid\\d+\\[Organism:exp\\]")) {
    return(query)
  }

  all_fields_query <- query %>%
    stringr::str_replace_all("\\[Organism:exp\\]", "[All Fields]") %>%
    stringr::str_replace_all("\\[Organism\\]", "[All Fields]")

  variants <- c(query, all_fields_query)
  variants <- clean_text(variants)
  unique(variants[!is.na(variants)])
}

prefer_taxonomy_link_fallback <- function(pathogen) {
  pathogen <- clean_text(pathogen)

  pathogen %in% c(
    "Human mastadenovirus B",
    "Mastadenovirus blackbeardi serotype 14",
    "Mammarenavirus juninense",
    "Mammarenavirus lassaense",
    "Mammarenavirus lujoense",
    "Mamastrovirus virginiaense"
  )
}

record_filter_pattern <- function(pathogen, disease_name = NA_character_) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)

  if (identical(pathogen, "Mastadenovirus blackbeardi")) {
    return("(blackbeardi|adenovirus\\s*14|adenovirus\\s*b14|b14)")
  }

  if (identical(pathogen, "Vibrio cholerae") && identical(disease_name, "Cholera")) {
    return("(vibrio cholerae|\\bo1\\b|\\bo139\\b|0139)")
  }

  if (identical(pathogen, "Alphainfluenzavirus influenzae (H2N1)")) {
    return("\\bH2N1\\b")
  }

  NA_character_
}

build_taxid_query <- function(...) {
  tax_ids <- c(...)
  tax_ids <- unlist(purrr::map(tax_ids, split_multi_value), use.names = FALSE)
  tax_ids <- clean_text(tax_ids)
  tax_ids <- tax_ids[!is.na(tax_ids)]

  if (length(tax_ids) == 0) {
    return(NA_character_)
  }

  paste0("txid", unique(tax_ids), "[Organism:exp]", collapse = " OR ")
}

extract_tax_ids_from_query <- function(taxid_query) {
  taxid_query <- clean_text(taxid_query)

  if (is.na(taxid_query)) {
    return(character(0))
  }

  matches <- stringr::str_match_all(taxid_query, "txid(\\d+)")[[1]]

  if (nrow(matches) == 0) {
    return(character(0))
  }

  clean_text(unique(matches[, 2]))
}

build_name_query <- function(...) {
  names_in <- c(...)
  names_in <- unlist(purrr::map(names_in, split_multi_value), use.names = FALSE)
  terms <- quote_organism_term(names_in)

  if (length(terms) == 0) {
    return(NA_character_)
  }

  paste(unique(terms), collapse = " OR ")
}

build_all_fields_query <- function(...) {
  names_in <- c(...)
  names_in <- unlist(purrr::map(names_in, split_multi_value), use.names = FALSE)
  terms <- quote_all_fields_term(names_in)

  if (length(terms) == 0) {
    return(NA_character_)
  }

  paste(unique(terms), collapse = " OR ")
}

resolve_target_resolution_type <- function(pathogen, disease_name, who_match_status, who_match_count) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)
  who_match_status <- clean_text(who_match_status)
  who_match_count <- suppressWarnings(as.integer(who_match_count))

  if (identical(pathogen, "Vibrio cholerae") && identical(disease_name, "Cholera")) {
    return("disease_profile")
  }

  if (is_influenza_network_target(pathogen)) {
    return("disease_profile")
  }

  if (is_sars_covid_network_target(pathogen, disease_name)) {
    return("disease_profile")
  }

  if (who_match_status == "exact_pair") {
    return("exact_candidate")
  }

  if (who_match_status == "disease_pathogen_fuzzy" && !is.na(who_match_count) && who_match_count == 1L) {
    return("narrow_fuzzy_candidate")
  }

  if (who_match_status == "disease_level") {
    return("disease_level")
  }

  if (who_match_status == "pathogen_only_fuzzy") {
    return("pathogen_only_fuzzy")
  }

  "network_only"
}

resolve_query_profile <- function(pathogen, disease_name, target_resolution_type) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)
  target_resolution_type <- clean_text(target_resolution_type)

  if (identical(pathogen, "Vibrio cholerae") && identical(disease_name, "Cholera")) {
    return("cholera_disease_profile")
  }

  if (is_influenza_network_target(pathogen)) {
    return("influenza_subtype_expansion")
  }

  if (is_sars_covid_network_target(pathogen, disease_name)) {
    return("sars_covid_profile")
  }

  if (identical(pathogen, "Mastadenovirus blackbeardi")) {
    return("adenovirus14_narrow")
  }

  dplyr::case_when(
    target_resolution_type == "exact_candidate" ~ "exact_candidate",
    target_resolution_type == "narrow_fuzzy_candidate" ~ "narrow_fuzzy_candidate",
    target_resolution_type == "disease_level" ~ "disease_level",
    target_resolution_type == "pathogen_only_fuzzy" ~ "pathogen_only_fuzzy",
    target_resolution_type == "network_only" ~ "network_only",
    TRUE ~ "default"
  )
}

resolve_organism_query <- function(
  pathogen,
  disease_name,
  who_pathogen_candidates,
  previous_name,
  msl39_viral_name,
  virion_names,
  clover_names,
  network_aliases,
  manual_query_aliases,
  query_profile
) {
  query_profile <- clean_text(query_profile)

  if (identical(query_profile, "sars_covid_profile")) {
    return(build_name_query(
      "Severe acute respiratory syndrome coronavirus 2",
      "Severe acute respiratory syndrome-related coronavirus",
      "SARS coronavirus"
    ))
  }

  candidate_pool <- dplyr::case_when(
    query_profile %in% c("adenovirus14_narrow", "narrow_fuzzy_candidate") ~ list(c(
      pathogen,
      who_pathogen_candidates,
      network_aliases,
      manual_query_aliases
    )),
    query_profile == "sars_covid_profile" ~ list(c(
      "Severe acute respiratory syndrome coronavirus 2",
      "Severe acute respiratory syndrome-related coronavirus",
      "SARS coronavirus"
    )),
    query_profile == "exact_candidate" ~ list(c(
      pathogen,
      who_pathogen_candidates,
      previous_name,
      msl39_viral_name,
      virion_names,
      clover_names,
      network_aliases,
      manual_query_aliases
    )),
    query_profile == "influenza_subtype_expansion" ~ list("Influenza A virus"),
    TRUE ~ list(c(
      pathogen,
      who_pathogen_candidates,
      network_aliases,
      manual_query_aliases
    ))
  )[[1]]

  keep_broad_aliases <- query_profile %in% c("cholera_disease_profile", "influenza_subtype_expansion")
  candidate_pool <- filter_organism_aliases(
    aliases = candidate_pool,
    target_pathogen = pathogen,
    disease_name = disease_name,
    keep_broad_aliases = keep_broad_aliases
  )

  build_name_query(candidate_pool)
}

resolve_all_fields_query <- function(pathogen, disease_name, who_pathogen_candidates, query_profile) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)
  query_profile <- clean_text(query_profile)

  if (query_profile == "cholera_disease_profile") {
    return(build_all_fields_query(
      "Vibrio cholerae O1",
      "Vibrio cholerae O139",
      "Vibrio cholerae serogroup 0139"
    ))
  }

  if (query_profile == "disease_level") {
    return(build_all_fields_query(who_pathogen_candidates, pathogen))
  }

  if (query_profile == "sars_covid_profile") {
    return(build_all_fields_query(
      "SARS-CoV-2",
      "COVID-19",
      "2019-nCoV",
      "SARS-CoV",
      "SARS coronavirus"
    ))
  }

  NA_character_
}

resolve_allow_taxonomy_link_fallback <- function(pathogen, who_match_status, taxid_query, query_profile) {
  pathogen <- clean_text(pathogen)
  who_match_status <- clean_text(who_match_status)
  taxid_query <- clean_text(taxid_query)
  query_profile <- clean_text(query_profile)

  !is.na(taxid_query) &&
    who_match_status == "exact_pair" &&
    query_profile == "exact_candidate" &&
    prefer_taxonomy_link_fallback(pathogen)
}

resolve_query_specificity_default <- function(target_resolution_type, query_profile) {
  target_resolution_type <- clean_text(target_resolution_type)
  query_profile <- clean_text(query_profile)

  dplyr::case_when(
    query_profile %in% c("exact_candidate", "narrow_fuzzy_candidate", "adenovirus14_narrow") ~ "specific",
    query_profile %in% c("cholera_disease_profile", "sars_covid_profile", "influenza_subtype_expansion", "influenza_exact_candidate", "influenza_family_fallback") ~ "disease_broad_but_acceptable",
    target_resolution_type %in% c("disease_level", "pathogen_only_fuzzy", "network_only") ~ "review_needed",
    TRUE ~ "review_needed"
  )
}

resolve_query_specificity_reason <- function(target_resolution_type, query_profile) {
  target_resolution_type <- clean_text(target_resolution_type)
  query_profile <- clean_text(query_profile)

  dplyr::case_when(
    query_profile == "exact_candidate" ~ "Exact WHO/network target with narrow organism-level search terms.",
    query_profile == "narrow_fuzzy_candidate" ~ "Single WHO candidate matched fuzzily; search narrowed to that candidate and close aliases only.",
    query_profile == "adenovirus14_narrow" ~ "Adenovirus row narrowed to serotype-14 style aliases to avoid broad adenovirus B retrieval.",
    query_profile == "cholera_disease_profile" ~ "Broad network target treated as a cholera-specific etiologic profile rather than literal Vibrio coverage.",
    query_profile == "sars_covid_profile" ~ "SARS/COVID network target is searched with curated SARS-CoV and SARS-CoV-2 terms for disease-level geographic coverage.",
    query_profile == "influenza_subtype_expansion" ~ "Broad influenza network target will be executed as a subtype-specific union with fallback family rows.",
    query_profile == "influenza_exact_candidate" ~ "Influenza network target is being retrieved via an exact subtype-specific search.",
    query_profile == "influenza_family_fallback" ~ "Influenza network target is being retrieved via a broader family fallback subtype search.",
    target_resolution_type %in% c("disease_level", "pathogen_only_fuzzy", "network_only") ~ "Target remains broad or weakly matched; results should be reviewed for specificity.",
    TRUE ~ "Specificity is uncertain and should be reviewed."
  )
}

resolve_network_target_candidates <- function(network_pathogen, network_disease_name, who_enriched) {
  network_pathogen <- clean_text(network_pathogen)
  network_disease_name <- clean_text(network_disease_name)

  exact_matches <- who_enriched %>%
    filter(
      Pathogens == network_pathogen,
      if (is.na(network_disease_name)) is.na(Disease_name) else Disease_name == network_disease_name
    )

  disease_matches <- if (is.na(network_disease_name)) {
    who_enriched %>% slice(0)
  } else {
    who_enriched %>%
      filter(Disease_name == network_disease_name)
  }

  fuzzy_matches <- if (nrow(disease_matches) == 0) {
    disease_matches
  } else {
    disease_matches %>%
      filter(purrr::map_lgl(Pathogens, pathogen_name_fuzzy_match, target = network_pathogen))
  }

  pathogen_only_matches <- if (is.na(network_pathogen)) {
    who_enriched %>% slice(0)
  } else {
    who_enriched %>%
      filter(purrr::map_lgl(Pathogens, pathogen_name_fuzzy_match, target = network_pathogen))
  }

  match_status <- dplyr::case_when(
    nrow(exact_matches) > 0 ~ "exact_pair",
    nrow(fuzzy_matches) > 0 ~ "disease_pathogen_fuzzy",
    nrow(disease_matches) > 0 ~ "disease_level",
    nrow(pathogen_only_matches) > 0 ~ "pathogen_only_fuzzy",
    TRUE ~ "network_only"
  )

  matched_rows <- dplyr::case_when(
    match_status == "exact_pair" ~ list(exact_matches),
    match_status == "disease_pathogen_fuzzy" ~ list(fuzzy_matches),
    match_status == "disease_level" ~ list(disease_matches),
    match_status == "pathogen_only_fuzzy" ~ list(pathogen_only_matches),
    TRUE ~ list(who_enriched %>% slice(0))
  )[[1]]

  tibble(
    Pathogens = network_pathogen,
    Disease_name = network_disease_name,
    network_pathogen = network_pathogen,
    network_disease_name = network_disease_name,
    who_match_status = match_status,
    who_match_count = nrow(matched_rows),
    who_pathogen_candidates = collapse_unique(matched_rows$Pathogens),
    Family = collapse_unique(matched_rows$Family),
    PHEIC_risk = collapse_unique(matched_rows$PHEIC_risk),
    previous_name = collapse_unique(matched_rows$previous_name),
    msl39_viral_name = collapse_unique(matched_rows$msl39_viral_name),
    virion_tax_ids = collapse_unique(matched_rows$virion_tax_ids),
    virion_names = collapse_unique(matched_rows$virion_names),
    matched_name_types = collapse_unique(matched_rows$matched_name_types),
    clover_tax_ids = collapse_unique(matched_rows$clover_tax_ids),
    clover_names = collapse_unique(matched_rows$clover_names),
    network_aliases = network_pathogen
  )
}

build_network_query_manifest <- function(
  network_data,
  who_manifest,
  virion_manifest,
  clover_manifest
) {
  who_enriched <- who_manifest %>%
    left_join(virion_manifest, by = c("Pathogens", "Disease_name")) %>%
    left_join(clover_manifest, by = c("Pathogens", "Disease_name"))

  network_targets <- network_data %>%
    transmute(
      network_pathogen = Pathogen,
      network_disease_name = Disease_name
    ) %>%
    distinct()

    network_targets %>%
    mutate(
      network_pathogen = clean_text(network_pathogen),
      network_disease_name = clean_text(network_disease_name)
    ) %>%
    arrange(network_disease_name, network_pathogen) %>%
    purrr::pmap_dfr(
      function(network_pathogen, network_disease_name) {
        resolve_network_target_candidates(
          network_pathogen = network_pathogen,
          network_disease_name = network_disease_name,
          who_enriched = who_enriched
        )
      }
    ) %>%
    mutate(
      influenza_candidate_sets = purrr::map(
        who_pathogen_candidates,
        split_influenza_candidate_sets
      ),
      influenza_search_mode = dplyr::if_else(
        purrr::map_lgl(Pathogens, is_influenza_network_target),
        "subtype_expansion",
        NA_character_
      ),
      influenza_exact_subtypes = purrr::map_chr(
        influenza_candidate_sets,
        ~collapse_unique(.x$exact_subtypes)
      ),
      influenza_family_fallbacks = purrr::map_chr(
        influenza_candidate_sets,
        ~collapse_unique(.x$family_fallbacks)
      ),
      influenza_family_policy = dplyr::if_else(
        !is.na(influenza_search_mode),
        "fallback_only",
        NA_character_
      ),
      manual_query_aliases = purrr::map2_chr(
        Pathogens,
        Disease_name,
        ~collapse_unique(c(manual_query_aliases(.x, .y), .x))
      ),
      target_resolution_type = purrr::pmap_chr(
        list(Pathogens, Disease_name, who_match_status, who_match_count),
        ~resolve_target_resolution_type(..1, ..2, ..3, ..4)
      ),
      query_profile = purrr::pmap_chr(
        list(Pathogens, Disease_name, target_resolution_type),
        ~resolve_query_profile(..1, ..2, ..3)
      ),
      taxid_query = purrr::pmap_chr(
        list(virion_tax_ids, clover_tax_ids),
        ~build_taxid_query(..1, ..2)
      ),
      organism_query = purrr::pmap_chr(
        list(Pathogens, Disease_name, who_pathogen_candidates, previous_name, msl39_viral_name, virion_names, clover_names, network_aliases, manual_query_aliases, query_profile),
        ~resolve_organism_query(..1, ..2, ..3, ..4, ..5, ..6, ..7, ..8, ..9, ..10)
      ),
      all_fields_query = purrr::pmap_chr(
        list(Pathogens, Disease_name, who_pathogen_candidates, query_profile),
        ~resolve_all_fields_query(..1, ..2, ..3, ..4)
      ),
      name_query = organism_query,
      allow_taxonomy_link_fallback = purrr::pmap_lgl(
        list(Pathogens, who_match_status, taxid_query, query_profile),
        ~resolve_allow_taxonomy_link_fallback(..1, ..2, ..3, ..4)
      ),
      query_specificity_status = purrr::pmap_chr(
        list(target_resolution_type, query_profile),
        ~resolve_query_specificity_default(..1, ..2)
      ),
      query_specificity_reason = purrr::pmap_chr(
        list(target_resolution_type, query_profile),
        ~resolve_query_specificity_reason(..1, ..2)
      ),
      preferred_metadata_source = purrr::pmap_chr(
        list(dplyr::coalesce(who_pathogen_candidates, Pathogens), Family, virion_tax_ids, clover_tax_ids, msl39_viral_name),
        ~choose_metadata_source(..1, ..2, ..3, ..4, ..5)
      ),
      metadata_source_reason = purrr::pmap_chr(
        list(dplyr::coalesce(who_pathogen_candidates, Pathogens), Family, virion_tax_ids, clover_tax_ids, msl39_viral_name),
        ~metadata_source_reason(..1, ..2, ..3, ..4, ..5)
      ),
      query_strategy = dplyr::case_when(
        query_profile == "influenza_subtype_expansion" ~ "name_subtype",
        query_profile == "sars_covid_profile" ~ "name",
        !is.na(taxid_query) ~ "taxid",
        !is.na(organism_query) | !is.na(all_fields_query) ~ "name",
        TRUE ~ "missing"
      ),
      search_query = dplyr::case_when(
        query_profile == "influenza_subtype_expansion" ~ dplyr::coalesce(organism_query, all_fields_query, taxid_query),
        query_profile == "sars_covid_profile" ~ dplyr::coalesce(all_fields_query, organism_query, taxid_query),
        query_profile == "cholera_disease_profile" ~ dplyr::coalesce(all_fields_query, organism_query, taxid_query),
        query_profile %in% c("adenovirus14_narrow", "narrow_fuzzy_candidate", "exact_candidate") ~ dplyr::coalesce(organism_query, all_fields_query, taxid_query),
        TRUE ~ dplyr::coalesce(taxid_query, organism_query, all_fields_query)
      ),
      geo_query = purrr::map_chr(search_query, build_geo_query)
    ) %>%
    select(-influenza_candidate_sets)
}

prefer_name_query_pathogen <- function(pathogen) {
  pathogen <- clean_text(pathogen)

  if (is.na(pathogen)) {
    return(FALSE)
  }

  stringr::str_detect(
    pathogen,
    stringr::regex("^Alphainfluenzavirus influenzae( \\(|$)", ignore_case = TRUE)
  )
}

is_strict_country_plateau_target <- function(pathogen, disease_name, query_profile = NA_character_) {
  pathogen <- clean_text(pathogen)
  disease_name <- clean_text(disease_name)
  query_profile <- clean_text(query_profile)

  broad_diseases <- c(
    "Severe Acute Respiratory Syndrome (SARS); COVID-19",
    "HIV infection / AIDS",
    "Influenza",
    "Non-typhoidal salmonellosis",
    "Dengue",
    "Cholera"
  )

  broad_pathogens <- c(
    "Betacoronavirus pandemicum",
    "Severe acute respiratory syndrome-related coronavirus",
    "Lentivirus humimdef1",
    "Alphainfluenzavirus influenzae",
    "Salmonella enterica",
    "Orthoflavivirus denguei",
    "Vibrio cholerae"
  )

  if (!is.na(disease_name) && disease_name %in% broad_diseases) {
    return(TRUE)
  }

  if (!is.na(pathogen) && pathogen %in% broad_pathogens) {
    return(TRUE)
  }

  identical(query_profile, "influenza_subtype_expansion")
}

resolve_plateau_controls <- function(
  pathogen,
  disease_name,
  query_profile,
  plateau_batches,
  min_batches_before_stop,
  min_records_before_stop,
  strict_plateau_batches,
  strict_min_batches_before_stop,
  strict_min_records_before_stop
) {
  if (is_strict_country_plateau_target(pathogen, disease_name, query_profile)) {
    return(list(
      plateau_batches = max(plateau_batches, strict_plateau_batches),
      min_batches_before_stop = max(min_batches_before_stop, strict_min_batches_before_stop),
      min_records_before_stop = max(min_records_before_stop, strict_min_records_before_stop)
    ))
  }

  list(
    plateau_batches = plateau_batches,
    min_batches_before_stop = min_batches_before_stop,
    min_records_before_stop = min_records_before_stop
  )
}

is_probable_viral_family <- function(family) {
  family <- clean_text(family)

  if (is.na(family)) {
    return(FALSE)
  }

  stringr::str_detect(
    family,
    stringr::regex("viridae|virinae|virales|virus|phenuiviridae|arenaviridae|flaviviridae|orthomyxoviridae|retroviridae|poxviridae|paramyxoviridae|picornaviridae|coronaviridae", ignore_case = TRUE)
  )
}

is_probable_viral_pathogen <- function(
  pathogen,
  family = NA_character_,
  virion_tax_ids = NA_character_,
  msl39_viral_name = NA_character_
) {
  pathogen <- clean_text(pathogen)

  if (!is.na(clean_text(virion_tax_ids)) || !is.na(clean_text(msl39_viral_name))) {
    return(TRUE)
  }

  if (is_probable_viral_family(family)) {
    return(TRUE)
  }

  if (is.na(pathogen)) {
    return(FALSE)
  }

  stringr::str_detect(
    pathogen,
    stringr::regex("virus|influenza|hiv|adenovirus|astrovirus|enterovirus|rotavirus|poxvirus|sarbecovirus|merbecovirus|ebolavirus|marburgvirus|flavivirus|alphavirus|hantavirus|hepadnavirus|henipavirus|mammarenavirus", ignore_case = TRUE)
  )
}

choose_metadata_source <- function(
  pathogen,
  family = NA_character_,
  virion_tax_ids = NA_character_,
  clover_tax_ids = NA_character_,
  msl39_viral_name = NA_character_
) {
  if (!is.na(clean_text(clover_tax_ids))) {
    return("biosample")
  }

  if (is_probable_viral_pathogen(
    pathogen = pathogen,
    family = family,
    virion_tax_ids = virion_tax_ids,
    msl39_viral_name = msl39_viral_name
  )) {
    return("nuccore")
  }

  "biosample"
}

metadata_source_reason <- function(
  pathogen,
  family = NA_character_,
  virion_tax_ids = NA_character_,
  clover_tax_ids = NA_character_,
  msl39_viral_name = NA_character_
) {
  source <- choose_metadata_source(
    pathogen = pathogen,
    family = family,
    virion_tax_ids = virion_tax_ids,
    clover_tax_ids = clover_tax_ids,
    msl39_viral_name = msl39_viral_name
  )

  if (identical(source, "biosample") && !is.na(clean_text(clover_tax_ids))) {
    return("CLOVER-linked bacterial pathogen; prefer sample-level geography over accession-heavy sequence records.")
  }

  if (identical(source, "nuccore")) {
    return("Virus or virus-like pathogen with VIRION/ICTV-style support; accession source qualifiers are a good first-pass geography source.")
  }

  "Non-viral or ambiguous pathogen without strong viral taxonomy support; prefer sample-level metadata first."
}

recommend_followup_metadata_source <- function(
  preferred_metadata_source,
  records_found = NA_integer_,
  countries_observed = NA_integer_,
  note = NA_character_
) {
  preferred_metadata_source <- clean_text(preferred_metadata_source)
  note <- clean_text(note)

  if (identical(preferred_metadata_source, "biosample")) {
    return("biosample")
  }

  if (
    !is.na(records_found) &&
      records_found >= 100000L &&
      !is.na(countries_observed) &&
      countries_observed <= 10L
  ) {
    return("biosample_review")
  }

  if (
    identical(note, "stopped_after_country_plateau") &&
      !is.na(records_found) &&
      records_found >= 1000000L &&
      !is.na(countries_observed) &&
      countries_observed <= 25L
  ) {
    return("biosample_review")
  }

  "nuccore"
}

followup_metadata_source_reason <- function(
  preferred_metadata_source,
  records_found = NA_integer_,
  countries_observed = NA_integer_,
  note = NA_character_
) {
  followup_source <- recommend_followup_metadata_source(
    preferred_metadata_source = preferred_metadata_source,
    records_found = records_found,
    countries_observed = countries_observed,
    note = note
  )

  if (identical(followup_source, "biosample_review")) {
    return("Very large accession count with weak country yield; sample metadata should be checked before trusting nuccore geography coverage.")
  }

  if (identical(followup_source, "biosample")) {
    return("Initial routing already points to BioSample.")
  }

  "Nuccore remains the preferred follow-up source for this pathogen."
}

evaluate_query_specificity <- function(manifest_row, fetched_records) {
  default_status <- clean_text(manifest_row$query_specificity_status)
  default_reason <- clean_text(manifest_row$query_specificity_reason)
  query_profile <- clean_text(manifest_row$query_profile)
  target_resolution_type <- clean_text(manifest_row$target_resolution_type)

  if (nrow(fetched_records) == 0) {
    return(list(
      status = dplyr::coalesce(default_status, "review_needed"),
      reason = dplyr::coalesce(default_reason, "No fetched records available to assess specificity."),
      organism_count = 0L
    ))
  }

  organisms <- fetched_records$organism %>%
    clean_text() %>%
    stats::na.omit() %>%
    unique()

  organism_count <- length(organisms)

  if (query_profile %in% c("exact_candidate", "narrow_fuzzy_candidate", "adenovirus14_narrow")) {
    if (organism_count <= 3L) {
      return(list(
        status = "specific",
        reason = "Fetched records stayed within a very small organism set for a narrow target.",
        organism_count = organism_count
      ))
    }

    return(list(
      status = "too_broad_review",
      reason = "Fetched records span many organism labels for a target that should be narrow.",
      organism_count = organism_count
    ))
  }

  if (query_profile %in% c(
    "cholera_disease_profile",
    "sars_covid_profile",
    "influenza_subtype_expansion",
    "influenza_exact_candidate",
    "influenza_family_fallback"
  )) {
    return(list(
      status = "disease_broad_but_acceptable",
      reason = "Broad disease profile query is acceptable, but should still be interpreted as disease-level rather than exact-pathogen coverage.",
      organism_count = organism_count
    ))
  }

  if (target_resolution_type %in% c("disease_level", "pathogen_only_fuzzy", "network_only")) {
    return(list(
      status = "too_broad_review",
      reason = "Target is broad or weakly matched, so results should be reviewed before being treated as pathogen-specific.",
      organism_count = organism_count
    ))
  }

  list(
    status = dplyr::coalesce(default_status, "review_needed"),
    reason = dplyr::coalesce(default_reason, "Specificity could not be assessed automatically."),
    organism_count = organism_count
  )
}

extract_gbqual_value <- function(source_feature, qualifier_name) {
  if (inherits(source_feature, "xml_missing")) {
    return(NA_character_)
  }

  qualifier_paths <- c(
    paste0(
      ".//GBQualifier[GBQualifier_name='",
      qualifier_name,
      "']/GBQualifier_value"
    ),
    paste0(
      ".//INSDQualifier[INSDQualifier_name='",
      qualifier_name,
      "']/INSDQualifier_value"
    )
  )

  node <- xml2::xml_find_first(source_feature, paste(qualifier_paths, collapse = " | "))

  value <- xml2::xml_text(node)

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

extract_seq_value <- function(seq_node, field_name) {
  field_paths <- c(
    paste0("./GBSeq_", field_name),
    paste0("./INSDSeq_", field_name)
  )

  node <- xml2::xml_find_first(seq_node, paste(field_paths, collapse = " | "))
  value <- xml2::xml_text(node)

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

standardize_country_name <- function(country_raw, geo_loc_name_raw) {
  location_raw <- dplyr::coalesce(clean_text(country_raw), clean_text(geo_loc_name_raw))

  if (is.na(location_raw)) {
    return(NA_character_)
  }

  country <- stringr::str_split_fixed(location_raw, ":", n = 2)[, 1]
  country <- stringr::str_split_fixed(country, ",", n = 2)[, 1]
  country <- clean_text(country)

  dplyr::case_when(
    country %in% c("USA", "U.S.A.", "United States of America") ~ "United States",
    country %in% c("UK", "U.K.") ~ "United Kingdom",
    country == "Viet Nam" ~ "Vietnam",
    country == "Russian Federation" ~ "Russia",
    country == "Czech Republic" ~ "Czechia",
    TRUE ~ country
  )
}

parse_gbseq_nodes <- function(xml_text, manifest_row) {
  if (is.na(xml_text) || !nzchar(xml_text)) {
    return(tibble())
  }

  doc <- xml2::read_xml(xml_text)
  gbseq_nodes <- xml2::xml_find_all(doc, ".//GBSeq | .//INSDSeq")

  if (length(gbseq_nodes) == 0) {
    return(tibble())
  }

  purrr::map_dfr(gbseq_nodes, function(seq_node) {
    source_feature <- xml2::xml_find_first(
      seq_node,
      ".//GBFeature[GBFeature_key='source'] | .//INSDFeature[INSDFeature_key='source']"
    )

    country_raw <- extract_gbqual_value(source_feature, "country")
    geo_loc_name_raw <- extract_gbqual_value(source_feature, "geo_loc_name")

    tibble(
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      network_pathogen = manifest_row$network_pathogen,
      network_disease_name = manifest_row$network_disease_name,
      Family = manifest_row$Family,
      who_match_status = manifest_row$who_match_status,
      who_match_count = manifest_row$who_match_count,
      who_pathogen_candidates = manifest_row$who_pathogen_candidates,
      target_resolution_type = manifest_row$target_resolution_type,
      query_profile = manifest_row$query_profile,
      influenza_search_mode = manifest_row$influenza_search_mode,
      influenza_query_label = manifest_row$influenza_query_label,
      influenza_query_class = manifest_row$influenza_query_class,
      preferred_metadata_source = manifest_row$preferred_metadata_source,
      metadata_source_reason = manifest_row$metadata_source_reason,
      query_strategy = manifest_row$query_strategy,
      query_used = manifest_row$query_used,
      taxid_query = manifest_row$taxid_query,
      organism_query = manifest_row$organism_query,
      all_fields_query = manifest_row$all_fields_query,
      name_query = manifest_row$name_query,
      accession_version = extract_seq_value(seq_node, "accession-version"),
      primary_accession = extract_seq_value(seq_node, "primary-accession"),
      definition = extract_seq_value(seq_node, "definition"),
      organism = extract_seq_value(seq_node, "organism"),
      taxonomy = extract_seq_value(seq_node, "taxonomy"),
      sequence_length = suppressWarnings(as.integer(extract_seq_value(seq_node, "length"))),
      country_raw = country_raw,
      geo_loc_name_raw = geo_loc_name_raw,
      country = standardize_country_name(country_raw, geo_loc_name_raw),
      lat_lon = extract_gbqual_value(source_feature, "lat_lon"),
      collection_date = extract_gbqual_value(source_feature, "collection_date"),
      host = extract_gbqual_value(source_feature, "host"),
      isolate = extract_gbqual_value(source_feature, "isolate"),
      strain = extract_gbqual_value(source_feature, "strain"),
      isolate_source = extract_gbqual_value(source_feature, "isolation_source"),
      db_xref = extract_gbqual_value(source_feature, "db_xref")
    )
  })
}

extract_biosample_attribute <- function(biosample_node, attribute_name = NA_character_, harmonized_name = NA_character_) {
  attribute_name <- clean_text(attribute_name)
  harmonized_name <- clean_text(harmonized_name)

  predicates <- character(0)

  if (!is.na(attribute_name)) {
    predicates <- c(predicates, paste0("@attribute_name='", attribute_name, "'"))
  }

  if (!is.na(harmonized_name)) {
    predicates <- c(predicates, paste0("@harmonized_name='", harmonized_name, "'"))
  }

  if (length(predicates) == 0) {
    return(NA_character_)
  }

  node <- xml2::xml_find_first(
    biosample_node,
    paste0(".//Attributes/Attribute[", paste(predicates, collapse = " or "), "]")
  )

  value <- xml2::xml_text(node)

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

extract_biosample_id <- function(biosample_node, db = "BioSample") {
  node <- xml2::xml_find_first(
    biosample_node,
    paste0(".//Ids/Id[@db='", db, "']")
  )

  value <- xml2::xml_text(node)

  if (identical(value, character(0)) || length(value) == 0) {
    return(NA_character_)
  }

  clean_text(value)
}

parse_biosample_nodes <- function(xml_text, manifest_row) {
  if (is.na(xml_text) || !nzchar(xml_text)) {
    return(tibble())
  }

  doc <- xml2::read_xml(xml_text)
  biosample_nodes <- xml2::xml_find_all(doc, ".//BioSample")

  if (length(biosample_nodes) == 0) {
    return(tibble())
  }

  purrr::map_dfr(biosample_nodes, function(biosample_node) {
    country_raw <- dplyr::coalesce(
      extract_biosample_attribute(biosample_node, attribute_name = "country", harmonized_name = "geo_loc_name"),
      extract_biosample_attribute(biosample_node, attribute_name = "geo_loc_name", harmonized_name = "geo_loc_name"),
      extract_biosample_attribute(biosample_node, attribute_name = "geographic location", harmonized_name = "geo_loc_name"),
      extract_biosample_attribute(biosample_node, attribute_name = "geographic location (country and/or sea)", harmonized_name = "geo_loc_name")
    )

    geo_loc_name_raw <- dplyr::coalesce(
      extract_biosample_attribute(biosample_node, attribute_name = "geo_loc_name", harmonized_name = "geo_loc_name"),
      extract_biosample_attribute(biosample_node, attribute_name = "geographic location", harmonized_name = "geo_loc_name"),
      extract_biosample_attribute(biosample_node, attribute_name = "geographic location (country and/or sea)", harmonized_name = "geo_loc_name")
    )

    biosample_accession <- clean_text(xml2::xml_attr(biosample_node, "accession"))
    biosample_id <- clean_text(xml2::xml_attr(biosample_node, "id"))

    tibble(
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      network_pathogen = manifest_row$network_pathogen,
      network_disease_name = manifest_row$network_disease_name,
      Family = manifest_row$Family,
      who_match_status = manifest_row$who_match_status,
      who_match_count = manifest_row$who_match_count,
      who_pathogen_candidates = manifest_row$who_pathogen_candidates,
      target_resolution_type = manifest_row$target_resolution_type,
      query_profile = manifest_row$query_profile,
      influenza_search_mode = manifest_row$influenza_search_mode,
      influenza_query_label = manifest_row$influenza_query_label,
      influenza_query_class = manifest_row$influenza_query_class,
      preferred_metadata_source = manifest_row$preferred_metadata_source,
      metadata_source_reason = manifest_row$metadata_source_reason,
      query_strategy = manifest_row$query_strategy,
      query_used = manifest_row$query_used,
      taxid_query = manifest_row$taxid_query,
      organism_query = manifest_row$organism_query,
      all_fields_query = manifest_row$all_fields_query,
      name_query = manifest_row$name_query,
      record_source_db = "biosample",
      biosample_accession = biosample_accession,
      biosample_id = biosample_id,
      accession_version = biosample_accession,
      primary_accession = biosample_accession,
      definition = clean_text(xml2::xml_text(xml2::xml_find_first(biosample_node, ".//Description/Title"))),
      organism = clean_text(xml2::xml_text(xml2::xml_find_first(biosample_node, ".//Description/Organism/OrganismName"))),
      taxonomy = clean_text(xml2::xml_attr(xml2::xml_find_first(biosample_node, ".//Description/Organism"), "taxonomy_id")),
      sequence_length = NA_integer_,
      country_raw = country_raw,
      geo_loc_name_raw = geo_loc_name_raw,
      country = standardize_country_name(country_raw, geo_loc_name_raw),
      lat_lon = extract_biosample_attribute(biosample_node, attribute_name = "lat_lon", harmonized_name = "lat_lon"),
      collection_date = extract_biosample_attribute(biosample_node, attribute_name = "collection_date", harmonized_name = "collection_date"),
      host = extract_biosample_attribute(biosample_node, attribute_name = "host", harmonized_name = "host"),
      isolate = dplyr::coalesce(
        extract_biosample_attribute(biosample_node, attribute_name = "isolate", harmonized_name = "isolate"),
        extract_biosample_attribute(biosample_node, attribute_name = "isolate_name_alias")
      ),
      strain = extract_biosample_attribute(biosample_node, attribute_name = "strain", harmonized_name = "strain"),
      isolate_source = extract_biosample_attribute(biosample_node, attribute_name = "isolation_source", harmonized_name = "isolation_source"),
      db_xref = extract_biosample_id(biosample_node, db = "SRA")
    )
  })
}

search_db <- "nucleotide"

build_esearch_url <- function(query, retmax = 0L, retstart = 0L, db = search_db) {
  query <- clean_text(query)

  if (is.na(query)) {
    stop("NCBI esearch query cannot be NA.")
  }

  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=", db,
    "&retmode=xml",
    "&usehistory=n",
    "&retmax=", as.integer(retmax),
    "&retstart=", as.integer(retstart),
    "&term=", utils::URLencode(query, reserved = FALSE)
  )

  if (exists("api_key", inherits = TRUE) && !is.na(api_key)) {
    url <- paste0(url, "&api_key=", utils::URLencode(api_key, reserved = FALSE))
  }

  url
}

read_url_text <- function(url) {
  con <- base::url(url, open = "rb")
  on.exit(close(con), add = TRUE)
  paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

parse_esearch_response <- function(xml_text) {
  if (is.na(xml_text) || !nzchar(xml_text)) {
    stop("NCBI esearch returned an empty response.")
  }

  doc <- xml2::read_xml(xml_text)
  error_nodes <- xml2::xml_find_all(doc, ".//ERROR | .//Error")

  if (length(error_nodes) > 0) {
    error_text <- collapse_unique(xml2::xml_text(error_nodes))
    stop("NCBI esearch error: ", error_text)
  }

  count_text <- xml2::xml_text(xml2::xml_find_first(doc, ".//Count"))
  ids <- xml2::xml_text(xml2::xml_find_all(doc, ".//IdList/Id"))

  list(
    count = suppressWarnings(as.integer(clean_text(count_text))),
    ids = clean_text(ids)
  )
}

search_nuccore_http <- function(query, retmax = 0L, retstart = 0L, db = search_db) {
  request_url <- build_esearch_url(
    query = query,
    retmax = retmax,
    retstart = retstart,
    db = db
  )

  log_progress(
    "NCBI ESearch request | method=http | db=",
    db,
    " | retmax=",
    retmax,
    " | retstart=",
    retstart,
    verbose = parse_env_flag("GENBANK_VERBOSE", default = TRUE)
  )

  response_xml <- read_url_text(request_url)
  result <- parse_esearch_response(response_xml)
  result$backend <- paste0("http:", db)
  result
}

search_nuccore_taxonomy_link <- function(taxid_query) {
  tax_ids <- extract_tax_ids_from_query(taxid_query)

  if (length(tax_ids) == 0) {
    stop("No usable taxonomy IDs found for taxonomy-link fallback.")
  }

  all_ids <- character(0)

  for (tax_id in tax_ids) {
    log_progress(
      "NCBI taxonomy link request | taxid=",
      tax_id,
      verbose = parse_env_flag("GENBANK_VERBOSE", default = TRUE)
    )

    link_result <- rentrez::entrez_link(
      dbfrom = "taxonomy",
      db = "nuccore",
      id = tax_id
    )

    link_ids <- clean_text(link_result$links$taxonomy_nuccore)
    link_ids <- link_ids[!is.na(link_ids)]
    all_ids <- unique(c(all_ids, link_ids))
  }

  list(
    count = length(all_ids),
    ids = all_ids,
    backend = "taxonomy_link"
  )
}

search_nuccore <- function(query, retmax = 0L, retstart = 0L) {
  query_variants <- search_query_variants(query)
  attempts <- list(
    list(method = "http", db = "nucleotide"),
    list(method = "http", db = "nuccore")
  )
  last_error <- "NCBI esearch failed for all attempted backends."

  for (query_variant in query_variants) {
    for (attempt in attempts) {
      result <- tryCatch(
        search_nuccore_http(
          query = query_variant,
          retmax = retmax,
          retstart = retstart,
          db = attempt$db
        ),
        error = function(e) e
      )

      if (!inherits(result, "error")) {
        return(result)
      }

      last_error <- conditionMessage(result)
    }
  }

  stop(last_error)
}

search_biosample <- function(query, retmax = 0L, retstart = 0L) {
  query_variants <- search_query_variants(query)
  last_error <- "NCBI BioSample esearch failed."

  for (query_variant in query_variants) {
    result <- tryCatch(
      search_nuccore_http(
        query = query_variant,
        retmax = retmax,
        retstart = retstart,
        db = "biosample"
      ),
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- conditionMessage(result)
  }

  stop(last_error)
}

fetch_nuccore_batch <- function(ids) {
  rentrez::entrez_fetch(
    db = "nuccore",
    id = ids,
    rettype = "gbc",
    retmode = "xml"
  )
}

fetch_biosample_batch <- function(ids) {
  rentrez::entrez_fetch(
    db = "biosample",
    id = ids,
    rettype = "full",
    retmode = "xml"
  )
}

is_request_too_large_error <- function(x) {
  if (!inherits(x, "error")) {
    return(FALSE)
  }

  stringr::str_detect(
    conditionMessage(x),
    stringr::regex("HTTP failure 414|request is too large", ignore_case = TRUE)
  )
}

fetch_nuccore_batch_with_retry <- function(
  ids,
  max_attempts = 3L,
  retry_wait_seconds = 1
) {
  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      fetch_nuccore_batch(ids = ids),
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- result

    if (attempt < max_attempts) {
      Sys.sleep(retry_wait_seconds * attempt)
    }
  }

  if (is_request_too_large_error(last_error) && length(ids) > 1) {
    split_index <- floor(length(ids) / 2)

    left_result <- fetch_nuccore_batch_with_retry(
      ids = ids[seq_len(split_index)],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )

    if (inherits(left_result, "error")) {
      return(left_result)
    }

    right_result <- fetch_nuccore_batch_with_retry(
      ids = ids[seq.int(split_index + 1L, length(ids))],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )

    if (inherits(right_result, "error")) {
      return(right_result)
    }

    return(c(left_result, right_result))
  }

  last_error
}

fetch_biosample_batch_with_retry <- function(
  ids,
  max_attempts = 3L,
  retry_wait_seconds = 1
) {
  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      fetch_biosample_batch(ids = ids),
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- result

    if (attempt < max_attempts) {
      Sys.sleep(retry_wait_seconds * attempt)
    }
  }

  if (is_request_too_large_error(last_error) && length(ids) > 1) {
    split_index <- floor(length(ids) / 2)

    left_result <- fetch_biosample_batch_with_retry(
      ids = ids[seq_len(split_index)],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )

    if (inherits(left_result, "error")) {
      return(left_result)
    }

    right_result <- fetch_biosample_batch_with_retry(
      ids = ids[seq.int(split_index + 1L, length(ids))],
      max_attempts = max_attempts,
      retry_wait_seconds = retry_wait_seconds
    )

    if (inherits(right_result, "error")) {
      return(right_result)
    }

    return(c(left_result, right_result))
  }

  last_error
}

collect_interval_sample_ids <- function(
  query,
  total_count,
  target_n,
  batch_size,
  search_fun,
  progress_prefix = "",
  pathogen = NA_character_,
  verbose = TRUE
) {
  if (is.na(total_count) || total_count <= 0 || target_n <= 0) {
    return(character(0))
  }

  n_batches <- ceiling(target_n / batch_size)
  retstarts <- floor(seq(0, n_batches - 1L) * total_count / n_batches)
  retstarts <- unique(pmax(0L, pmin(as.integer(total_count - 1L), as.integer(retstarts))))

  sampled_ids <- character(0)

  for (i in seq_along(retstarts)) {
    current_retstart <- retstarts[[i]]
    current_batch_size <- min(batch_size, target_n - length(sampled_ids))

    if (current_batch_size <= 0) {
      break
    }

    result <- tryCatch(
      search_fun(
        query = query,
        retmax = current_batch_size,
        retstart = current_retstart
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      log_progress(
        progress_prefix,
        "Interval ID sample failed for ",
        pathogen,
        " | sample_batch=",
        i,
        " | retstart=",
        current_retstart,
        " | error=",
        conditionMessage(result),
        verbose = verbose
      )
      next
    }

    batch_ids <- clean_text(result$ids)
    batch_ids <- batch_ids[!is.na(batch_ids)]
    sampled_ids <- unique(c(sampled_ids, batch_ids))

    log_progress(
      progress_prefix,
      "Interval ID sample ",
      i,
      "/",
      length(retstarts),
      " for ",
      pathogen,
      " | retstart=",
      current_retstart,
      " | ids_collected=",
      length(sampled_ids),
      verbose = verbose
    )

    if (length(sampled_ids) >= target_n) {
      break
    }
  }

  utils::head(sampled_ids, target_n)
}

build_influenza_subquery_row <- function(manifest_row, subtype_pathogen, query_class) {
  subtype_pathogen <- clean_text(subtype_pathogen)
  query_class <- clean_text(query_class)
  subtype_token <- extract_influenza_subtype_token(subtype_pathogen)

  sub_row <- manifest_row
  sub_row$Pathogens <- subtype_pathogen
  sub_row$who_pathogen_candidates <- subtype_pathogen
  sub_row$target_resolution_type <- "influenza_subtype_candidate"
  sub_row$query_profile <- dplyr::case_when(
    query_class == "exact_subtype" ~ "influenza_exact_candidate",
    query_class == "family_fallback" ~ "influenza_family_fallback",
    TRUE ~ "influenza_exact_candidate"
  )
  sub_row$taxid_query <- NA_character_
  sub_row$allow_taxonomy_link_fallback <- FALSE
  sub_row$query_strategy <- "name_subtype"
  sub_row$organism_query <- "\"Influenza A virus\"[Organism]"
  sub_row$all_fields_query <- if (!is.na(subtype_token)) {
    build_all_fields_query(subtype_token)
  } else {
    NA_character_
  }
  sub_row$name_query <- sub_row$organism_query
  sub_row$search_query <- dplyr::coalesce(sub_row$all_fields_query, sub_row$organism_query)
  sub_row$geo_query <- build_geo_query(sub_row$search_query)
  sub_row$query_specificity_status <- "disease_broad_but_acceptable"
  sub_row$query_specificity_reason <- dplyr::case_when(
    query_class == "exact_subtype" ~ "Influenza network target is being retrieved via an exact subtype-specific search.",
    query_class == "family_fallback" ~ "Influenza network target is being retrieved via a broader family fallback subtype search.",
    TRUE ~ "Influenza network target is being retrieved via a subtype-specific search."
  )
  sub_row$influenza_search_mode <- "subtype_expansion"
  sub_row$influenza_query_label <- subtype_pathogen
  sub_row$influenza_query_class <- query_class

  sub_row
}

summarize_influenza_subquery_logs <- function(manifest_row, sub_logs, fetched_records) {
  sub_logs <- sub_logs %>%
    mutate(
      records_found = suppressWarnings(as.integer(records_found)),
      records_fetched = suppressWarnings(as.integer(records_fetched)),
      batches_fetched = suppressWarnings(as.integer(batches_fetched)),
      countries_observed = suppressWarnings(as.integer(countries_observed)),
      organism_count_observed = suppressWarnings(as.integer(organism_count_observed)),
      last_batch_new_countries = suppressWarnings(as.integer(last_batch_new_countries))
    )

  specificity_assessment <- evaluate_query_specificity(
    manifest_row = manifest_row,
    fetched_records = fetched_records
  )
  countries_observed <- if ("country" %in% names(fetched_records)) {
    dplyr::n_distinct(fetched_records$country, na.rm = TRUE)
  } else {
    0L
  }

  total_records_found <- sum(sub_logs$records_found, na.rm = TRUE)
  if (all(is.na(sub_logs$records_found))) {
    total_records_found <- NA_integer_
  }

  stop_reason <- dplyr::case_when(
    any(sub_logs$note == "truncated_at_retmax", na.rm = TRUE) ~ "truncated_at_retmax",
    any(sub_logs$note == "stopped_after_country_plateau", na.rm = TRUE) ~ "stopped_after_country_plateau",
    any(sub_logs$note == "fetch_error", na.rm = TRUE) ~ "fetch_error",
    any(sub_logs$note == "search_failed", na.rm = TRUE) ~ "search_failed",
    TRUE ~ "fetched_all_available_records"
  )

  tibble(
    Pathogens = manifest_row$Pathogens,
    Disease_name = manifest_row$Disease_name,
    network_pathogen = manifest_row$network_pathogen,
    network_disease_name = manifest_row$network_disease_name,
    Family = manifest_row$Family,
    who_match_status = manifest_row$who_match_status,
    who_match_count = manifest_row$who_match_count,
    who_pathogen_candidates = manifest_row$who_pathogen_candidates,
    target_resolution_type = manifest_row$target_resolution_type,
    query_profile = manifest_row$query_profile,
    influenza_search_mode = manifest_row$influenza_search_mode,
    preferred_metadata_source = manifest_row$preferred_metadata_source,
    metadata_source_reason = manifest_row$metadata_source_reason,
    query_strategy = manifest_row$query_strategy,
    query_used = collapse_unique(sub_logs$query_used),
    records_found = total_records_found,
    records_fetched = nrow(fetched_records),
    batches_fetched = sum(sub_logs$batches_fetched, na.rm = TRUE),
    countries_observed = countries_observed,
    organism_count_observed = specificity_assessment$organism_count,
    last_batch_new_countries = dplyr::last(sub_logs$last_batch_new_countries, order_by = seq_len(nrow(sub_logs))),
    query_specificity_status = specificity_assessment$status,
    query_specificity_reason = specificity_assessment$reason,
    status = if (all(sub_logs$status == "ok", na.rm = TRUE)) "ok" else collapse_unique(sub_logs$status),
    note = stop_reason
  )
}

run_influenza_subtype_queries <- function(
  manifest_row,
  retmax,
  batch_size,
  plateau_batches,
  min_new_countries,
  min_batches_before_stop,
  min_records_before_stop,
  strict_plateau_batches,
  strict_min_batches_before_stop,
  strict_min_records_before_stop,
  interval_sample_threshold,
  manifest_index = NA_integer_,
  manifest_total = NA_integer_,
  verbose = TRUE
) {
  exact_subtypes <- split_multi_value(manifest_row$influenza_exact_subtypes)
  exact_subtypes <- clean_text(exact_subtypes)
  exact_subtypes <- exact_subtypes[!is.na(exact_subtypes)]

  family_fallbacks <- split_multi_value(manifest_row$influenza_family_fallbacks)
  family_fallbacks <- clean_text(family_fallbacks)
  family_fallbacks <- family_fallbacks[!is.na(family_fallbacks)]

  if (length(exact_subtypes) == 0 && length(family_fallbacks) == 0) {
    return(list(
      log = tibble(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        Family = manifest_row$Family,
        who_match_status = manifest_row$who_match_status,
        who_match_count = manifest_row$who_match_count,
        who_pathogen_candidates = manifest_row$who_pathogen_candidates,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        preferred_metadata_source = manifest_row$preferred_metadata_source,
        metadata_source_reason = manifest_row$metadata_source_reason,
        query_strategy = manifest_row$query_strategy,
        query_used = NA_character_,
        records_found = 0L,
        records_fetched = 0L,
        batches_fetched = 0L,
        countries_observed = 0L,
        organism_count_observed = 0L,
        last_batch_new_countries = 0L,
        query_specificity_status = "review_needed",
        query_specificity_reason = "Influenza subtype expansion was requested, but no subtype candidates were available.",
        status = "no_records",
        note = "no_influenza_subtypes"
      ),
      records = tibble()
    ))
  }

  sub_logs <- list()
  retained_records <- list()
  seen_countries <- character(0)

  run_one_subquery <- function(subtype_pathogen, query_class) {
    sub_row <- build_influenza_subquery_row(
      manifest_row = manifest_row,
      subtype_pathogen = subtype_pathogen,
      query_class = query_class
    )

    sub_result <- run_manifest_query(
      manifest_row = sub_row,
      retmax = retmax,
      batch_size = batch_size,
      plateau_batches = plateau_batches,
      min_new_countries = min_new_countries,
      min_batches_before_stop = min_batches_before_stop,
      min_records_before_stop = min_records_before_stop,
      strict_plateau_batches = strict_plateau_batches,
      strict_min_batches_before_stop = strict_min_batches_before_stop,
      strict_min_records_before_stop = strict_min_records_before_stop,
      interval_sample_threshold = interval_sample_threshold,
      manifest_index = manifest_index,
      manifest_total = manifest_total,
      verbose = verbose
    )

    sub_log <- sub_result$log %>%
      mutate(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode
      )

    sub_records <- sub_result$records %>%
      mutate(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        influenza_query_label = subtype_pathogen,
        influenza_query_class = query_class
      )

    list(log = sub_log, records = sub_records)
  }

  for (subtype_pathogen in exact_subtypes) {
    sub_result <- run_one_subquery(subtype_pathogen, "exact_subtype")
    sub_logs[[length(sub_logs) + 1L]] <- sub_result$log

    if (nrow(sub_result$records) > 0) {
      retained_records[[length(retained_records) + 1L]] <- sub_result$records
      seen_countries <- unique(c(seen_countries, clean_text(sub_result$records$country)))
    }
  }

  exact_record_count <- sum(vapply(retained_records, nrow, integer(1)))

  for (subtype_pathogen in family_fallbacks) {
    sub_result <- run_one_subquery(subtype_pathogen, "family_fallback")
    sub_logs[[length(sub_logs) + 1L]] <- sub_result$log

    if (nrow(sub_result$records) == 0) {
      next
    }

    fallback_countries <- clean_text(sub_result$records$country)
    fallback_countries <- unique(fallback_countries[!is.na(fallback_countries)])
    adds_new_countries <- length(setdiff(fallback_countries, seen_countries)) > 0

    if (exact_record_count == 0 || adds_new_countries) {
      retained_records[[length(retained_records) + 1L]] <- sub_result$records
      seen_countries <- unique(c(seen_countries, fallback_countries))
    }
  }

  combined_logs <- dplyr::bind_rows(sub_logs)
  combined_records <- dplyr::bind_rows(retained_records)

  list(
    log = summarize_influenza_subquery_logs(
      manifest_row = manifest_row,
      sub_logs = combined_logs,
      fetched_records = combined_records
    ),
    records = combined_records
  )
}

run_manifest_query <- function(
  manifest_row,
  retmax,
  batch_size,
  plateau_batches,
  min_new_countries,
  min_batches_before_stop,
  min_records_before_stop,
  strict_plateau_batches,
  strict_min_batches_before_stop,
  strict_min_records_before_stop,
  interval_sample_threshold,
  manifest_index = NA_integer_,
  manifest_total = NA_integer_,
  verbose = TRUE
) {
  if (identical(clean_text(manifest_row$query_profile), "influenza_subtype_expansion")) {
    return(run_influenza_subtype_queries(
      manifest_row = manifest_row,
      retmax = retmax,
      batch_size = batch_size,
      plateau_batches = plateau_batches,
      min_new_countries = min_new_countries,
      min_batches_before_stop = min_batches_before_stop,
      min_records_before_stop = min_records_before_stop,
      strict_plateau_batches = strict_plateau_batches,
      strict_min_batches_before_stop = strict_min_batches_before_stop,
      strict_min_records_before_stop = strict_min_records_before_stop,
      interval_sample_threshold = interval_sample_threshold,
      manifest_index = manifest_index,
      manifest_total = manifest_total,
      verbose = verbose
    ))
  }

  plateau_controls <- resolve_plateau_controls(
    pathogen = manifest_row$Pathogens,
    disease_name = manifest_row$Disease_name,
    query_profile = manifest_row$query_profile,
    plateau_batches = plateau_batches,
    min_batches_before_stop = min_batches_before_stop,
    min_records_before_stop = min_records_before_stop,
    strict_plateau_batches = strict_plateau_batches,
    strict_min_batches_before_stop = strict_min_batches_before_stop,
    strict_min_records_before_stop = strict_min_records_before_stop
  )

  plateau_batches <- plateau_controls$plateau_batches
  min_batches_before_stop <- plateau_controls$min_batches_before_stop
  min_records_before_stop <- plateau_controls$min_records_before_stop

  geo_query <- clean_text(manifest_row$geo_query)
  base_query <- clean_text(manifest_row$search_query)
  organism_query <- clean_text(manifest_row$organism_query)
  organism_geo_query <- build_geo_query(organism_query)
  all_fields_query <- clean_text(manifest_row$all_fields_query)
  all_fields_geo_query <- build_geo_query(all_fields_query)
  manual_queries <- manual_query_candidates(
    pathogen = manifest_row$Pathogens,
    disease_name = manifest_row$Disease_name
  )
  use_name_queries_only <- isTRUE(clean_text(manifest_row$query_profile) %in% c(
    "cholera_disease_profile",
    "influenza_exact_candidate",
    "influenza_family_fallback"
  ))
  filter_pattern <- record_filter_pattern(
    pathogen = manifest_row$Pathogens,
    disease_name = manifest_row$Disease_name
  )
  allow_taxonomy_link_fallback <- isTRUE(manifest_row$allow_taxonomy_link_fallback)
  source_db <- dplyr::coalesce(clean_text(manifest_row$preferred_metadata_source), "nuccore")
  progress_prefix <- if (
    !is.na(manifest_index) &&
      !is.na(manifest_total)
  ) {
    paste0("[", manifest_index, "/", manifest_total, "] ")
  } else {
    ""
  }

  log_progress(
    progress_prefix,
    "Starting ",
    manifest_row$Pathogens,
    if (!is.na(manifest_row$Disease_name)) paste0(" | disease=", manifest_row$Disease_name) else "",
    if (!is.na(manifest_row$query_strategy)) paste0(" | strategy=", manifest_row$query_strategy) else "",
    verbose = verbose
  )

  queries_to_try <- if (use_name_queries_only) {
    c(
      manual_queries,
      all_fields_geo_query,
      all_fields_query,
      organism_geo_query,
      organism_query
    )
  } else {
    c(
      manual_queries,
      organism_geo_query,
      organism_query,
      all_fields_geo_query,
      all_fields_query,
      geo_query,
      base_query
    )
  }
  queries_to_try <- queries_to_try[!is.na(queries_to_try)]
  queries_to_try <- unique(queries_to_try)

  search_result <- NULL
  query_used <- NA_character_
  search_error <- NA_character_

  if (
    identical(source_db, "nuccore") &&
    !use_name_queries_only &&
    allow_taxonomy_link_fallback &&
      !is.na(manifest_row$taxid_query)
  ) {
    taxonomy_attempt <- tryCatch(
      search_nuccore_taxonomy_link(
        taxid_query = manifest_row$taxid_query
      ),
      error = function(e) e
    )

    if (!inherits(taxonomy_attempt, "error") && taxonomy_attempt$count > 0) {
      search_result <- taxonomy_attempt
      query_used <- paste0("taxonomy_link:", manifest_row$taxid_query)
    } else if (inherits(taxonomy_attempt, "error")) {
      search_error <- conditionMessage(taxonomy_attempt)
    }
  }

  for (query in queries_to_try) {
    if (!is.null(search_result) && search_result$count > 0) {
      break
    }

    attempt <- tryCatch(
      if (identical(source_db, "biosample")) {
        search_biosample(query = query, retmax = retmax, retstart = 0L)
      } else {
        search_nuccore(query = query, retmax = retmax, retstart = 0L)
      },
      error = function(e) e
    )

    if (inherits(attempt, "error")) {
      search_error <- conditionMessage(attempt)
      next
    }

    search_result <- attempt
    query_used <- query

    if (attempt$count > 0) {
      break
    }
  }

  if (
    identical(source_db, "nuccore") &&
    !use_name_queries_only &&
    (is.null(search_result) || search_result$count == 0) &&
      allow_taxonomy_link_fallback &&
      !is.na(manifest_row$taxid_query)
  ) {
    taxonomy_attempt <- tryCatch(
      search_nuccore_taxonomy_link(
        taxid_query = manifest_row$taxid_query
      ),
      error = function(e) e
    )

    if (!inherits(taxonomy_attempt, "error") && taxonomy_attempt$count > 0) {
      search_result <- taxonomy_attempt
      query_used <- paste0("taxonomy_link:", manifest_row$taxid_query)
    } else if (inherits(taxonomy_attempt, "error")) {
      search_error <- conditionMessage(taxonomy_attempt)
    }
  }

  if (is.null(search_result)) {
    log_progress(
      progress_prefix,
      "Search failed for ",
      manifest_row$Pathogens,
      if (!is.na(search_error)) paste0(" | error=", search_error) else "",
      verbose = verbose
    )
    return(list(
      log = tibble(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        Family = manifest_row$Family,
        who_match_status = manifest_row$who_match_status,
        who_match_count = manifest_row$who_match_count,
        who_pathogen_candidates = manifest_row$who_pathogen_candidates,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        preferred_metadata_source = manifest_row$preferred_metadata_source,
        metadata_source_reason = manifest_row$metadata_source_reason,
        query_strategy = manifest_row$query_strategy,
        query_used = query_used,
        records_found = NA_integer_,
        records_fetched = 0L,
        batches_fetched = 0L,
        countries_observed = 0L,
        organism_count_observed = 0L,
        last_batch_new_countries = 0L,
        query_specificity_status = manifest_row$query_specificity_status,
        query_specificity_reason = dplyr::coalesce(search_error, manifest_row$query_specificity_reason),
        status = "search_failed",
        note = search_error
      ),
      records = tibble()
    ))
  }

  if (search_result$count == 0) {
    log_progress(
      progress_prefix,
      "No records found for ",
      manifest_row$Pathogens,
      verbose = verbose
    )
    return(list(
      log = tibble(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        Family = manifest_row$Family,
        who_match_status = manifest_row$who_match_status,
        who_match_count = manifest_row$who_match_count,
        who_pathogen_candidates = manifest_row$who_pathogen_candidates,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        preferred_metadata_source = manifest_row$preferred_metadata_source,
        metadata_source_reason = manifest_row$metadata_source_reason,
        query_strategy = manifest_row$query_strategy,
        query_used = query_used,
        records_found = 0L,
        records_fetched = 0L,
        batches_fetched = 0L,
        countries_observed = 0L,
        organism_count_observed = 0L,
        last_batch_new_countries = 0L,
        query_specificity_status = manifest_row$query_specificity_status,
        query_specificity_reason = manifest_row$query_specificity_reason,
        status = "no_records",
        note = NA_character_
      ),
      records = tibble()
    ))
  }

  records_to_fetch <- min(search_result$count, retmax)
  use_interval_sampling <- !is.na(search_result$count) &&
    search_result$count >= interval_sample_threshold

  if (use_interval_sampling) {
    search_fun <- if (identical(source_db, "biosample")) {
      search_biosample
    } else {
      search_nuccore
    }

    result_ids <- collect_interval_sample_ids(
      query = query_used,
      total_count = search_result$count,
      target_n = records_to_fetch,
      batch_size = batch_size,
      search_fun = search_fun,
      progress_prefix = progress_prefix,
      pathogen = manifest_row$Pathogens,
      verbose = verbose
    )
  } else {
    result_ids <- clean_text(search_result$ids)
    result_ids <- result_ids[!is.na(result_ids)]
  }

  if (records_to_fetch > 0 && length(result_ids) == 0) {
    log_progress(
      progress_prefix,
      "Search found records but returned no IDs for ",
      manifest_row$Pathogens,
      " | query=",
      query_used,
      verbose = verbose
    )
    return(list(
      log = tibble(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        Family = manifest_row$Family,
        who_match_status = manifest_row$who_match_status,
        who_match_count = manifest_row$who_match_count,
        who_pathogen_candidates = manifest_row$who_pathogen_candidates,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        preferred_metadata_source = manifest_row$preferred_metadata_source,
        metadata_source_reason = manifest_row$metadata_source_reason,
        query_strategy = manifest_row$query_strategy,
        query_used = query_used,
        records_found = as.integer(search_result$count),
        records_fetched = 0L,
        batches_fetched = 0L,
        countries_observed = 0L,
        organism_count_observed = 0L,
        last_batch_new_countries = 0L,
        query_specificity_status = manifest_row$query_specificity_status,
        query_specificity_reason = "Search returned matches but no IDs were available for fetch.",
        status = "search_failed",
        note = "search_returned_no_ids"
      ),
      records = tibble()
    ))
  }

  records_to_fetch <- min(records_to_fetch, length(result_ids))
  fetched_batches <- list()
  seen_countries <- character(0)
  plateau_counter <- 0L
  batch_counter <- 0L
  last_batch_new_countries <- 0L
  fetched_total <- 0L
  retstart <- 0L
  stop_reason <- if (search_result$count > retmax) "truncated_at_retmax" else "fetched_all_available_records"

  log_progress(
    progress_prefix,
    "Search returned ",
    search_result$count,
    " records for ",
    manifest_row$Pathogens,
    "; fetching up to ",
    records_to_fetch,
    " in batches of ",
    batch_size,
    " | source_db=",
    source_db,
    if (use_interval_sampling) paste0(" | sampling=interval_across_", search_result$count) else "",
    verbose = verbose
  )

  while (retstart < records_to_fetch) {
    current_batch_size <- min(batch_size, records_to_fetch - retstart)
    batch_ids <- result_ids[seq.int(retstart + 1L, retstart + current_batch_size)]

    if (length(batch_ids) == 0 || all(is.na(batch_ids))) {
      stop_reason <- "search_returned_no_ids"
      log_progress(
        progress_prefix,
        "Local ID slice returned no IDs for ",
        manifest_row$Pathogens,
        " | retstart=",
        retstart,
        verbose = verbose
      )
      break
    }

    batch_xml <- if (identical(source_db, "biosample")) {
      fetch_biosample_batch_with_retry(ids = batch_ids)
    } else {
      fetch_nuccore_batch_with_retry(ids = batch_ids)
    }
    batch_counter <- batch_counter + 1L

    if (inherits(batch_xml, "error")) {
      fetched_batches[[length(fetched_batches) + 1L]] <- tibble(
        Pathogens = manifest_row$Pathogens,
        Disease_name = manifest_row$Disease_name,
        network_pathogen = manifest_row$network_pathogen,
        network_disease_name = manifest_row$network_disease_name,
        Family = manifest_row$Family,
        who_match_status = manifest_row$who_match_status,
        who_match_count = manifest_row$who_match_count,
        who_pathogen_candidates = manifest_row$who_pathogen_candidates,
        target_resolution_type = manifest_row$target_resolution_type,
        query_profile = manifest_row$query_profile,
        influenza_search_mode = manifest_row$influenza_search_mode,
        influenza_query_label = manifest_row$influenza_query_label,
        influenza_query_class = manifest_row$influenza_query_class,
        preferred_metadata_source = manifest_row$preferred_metadata_source,
        metadata_source_reason = manifest_row$metadata_source_reason,
        query_strategy = manifest_row$query_strategy,
        query_used = query_used,
        record_source_db = source_db,
        biosample_accession = NA_character_,
        biosample_id = NA_character_,
        taxid_query = manifest_row$taxid_query,
        organism_query = manifest_row$organism_query,
        all_fields_query = manifest_row$all_fields_query,
        name_query = manifest_row$name_query,
        accession_version = NA_character_,
        primary_accession = NA_character_,
        definition = NA_character_,
        organism = NA_character_,
        taxonomy = NA_character_,
        sequence_length = NA_integer_,
        country_raw = NA_character_,
        geo_loc_name_raw = NA_character_,
        country = NA_character_,
        lat_lon = NA_character_,
        collection_date = NA_character_,
        host = NA_character_,
        isolate = NA_character_,
        strain = NA_character_,
        isolate_source = NA_character_,
        db_xref = NA_character_,
        fetch_error = conditionMessage(batch_xml)
      )
      stop_reason <- "fetch_error"
      log_progress(
        progress_prefix,
        "Batch ",
        batch_counter,
        " failed for ",
        manifest_row$Pathogens,
        " | id_start=",
        retstart,
        " | error=",
        conditionMessage(batch_xml),
        verbose = verbose
      )
      break
    }

    parsed <- purrr::map_dfr(
      batch_xml,
      ~if (identical(source_db, "biosample")) {
        parse_biosample_nodes(
          .x,
          manifest_row = dplyr::mutate(manifest_row, query_used = query_used)
        )
      } else {
        parse_gbseq_nodes(
          .x,
          manifest_row = dplyr::mutate(manifest_row, query_used = query_used)
        )
      }
    )

    if (nrow(parsed) > 0) {
      if (!is.na(filter_pattern)) {
        parsed <- parsed %>%
          dplyr::filter(
            stringr::str_detect(
              string = stringr::str_c(
                dplyr::coalesce(definition, ""),
                dplyr::coalesce(organism, ""),
                dplyr::coalesce(taxonomy, ""),
                dplyr::coalesce(isolate, ""),
                dplyr::coalesce(strain, ""),
                sep = " | "
              ),
              pattern = stringr::regex(filter_pattern, ignore_case = TRUE)
            )
          )
      }

      parsed <- parsed %>%
        dplyr::mutate(fetch_error = NA_character_)

      fetched_batches[[length(fetched_batches) + 1L]] <- parsed

      batch_countries <- parsed %>%
        dplyr::filter(!is.na(country)) %>%
        dplyr::pull(country) %>%
        unique()

      new_country_count <- length(setdiff(batch_countries, seen_countries))
      seen_countries <- unique(c(seen_countries, batch_countries))
    } else {
      new_country_count <- 0L
    }

    last_batch_new_countries <- new_country_count
    fetched_total <- sum(vapply(fetched_batches, nrow, integer(1)))

    log_progress(
      progress_prefix,
      "Batch ",
      batch_counter,
      " complete for ",
      manifest_row$Pathogens,
      " | fetched_this_batch=",
      nrow(parsed),
      " | fetched_total=",
      fetched_total,
      " | new_countries=",
      new_country_count,
      " | countries_observed=",
      length(seen_countries),
      verbose = verbose
    )

    if (
      fetched_total >= min_records_before_stop &&
        batch_counter >= min_batches_before_stop &&
        new_country_count <= min_new_countries
    ) {
      plateau_counter <- plateau_counter + 1L
    } else {
      plateau_counter <- 0L
    }

    retstart <- retstart + current_batch_size

    if (plateau_counter >= plateau_batches && retstart < records_to_fetch) {
      stop_reason <- "stopped_after_country_plateau"
      log_progress(
        progress_prefix,
        "Stopping early for ",
        manifest_row$Pathogens,
        " after country plateau | plateau_batches=",
        plateau_batches,
        " | min_records_before_stop=",
        min_records_before_stop,
        " | countries_observed=",
        length(seen_countries),
        verbose = verbose
      )
      break
    }
  }

  fetched_records <- dplyr::bind_rows(fetched_batches)
  specificity_assessment <- evaluate_query_specificity(
    manifest_row = manifest_row,
    fetched_records = fetched_records
  )

  log_progress(
    progress_prefix,
    "Finished ",
    manifest_row$Pathogens,
    " | records_found=",
    search_result$count,
    " | records_fetched=",
    nrow(fetched_records),
    " | batches_fetched=",
    batch_counter,
    " | countries_observed=",
    length(seen_countries),
    " | note=",
    stop_reason,
    verbose = verbose
  )

  list(
    log = tibble(
      Pathogens = manifest_row$Pathogens,
      Disease_name = manifest_row$Disease_name,
      network_pathogen = manifest_row$network_pathogen,
      network_disease_name = manifest_row$network_disease_name,
      Family = manifest_row$Family,
      who_match_status = manifest_row$who_match_status,
      who_match_count = manifest_row$who_match_count,
      who_pathogen_candidates = manifest_row$who_pathogen_candidates,
      target_resolution_type = manifest_row$target_resolution_type,
      query_profile = manifest_row$query_profile,
      influenza_search_mode = manifest_row$influenza_search_mode,
      preferred_metadata_source = manifest_row$preferred_metadata_source,
      metadata_source_reason = manifest_row$metadata_source_reason,
      query_strategy = manifest_row$query_strategy,
      query_used = query_used,
      records_found = as.integer(search_result$count),
      records_fetched = nrow(fetched_records),
      batches_fetched = batch_counter,
      countries_observed = length(seen_countries),
      organism_count_observed = specificity_assessment$organism_count,
      last_batch_new_countries = last_batch_new_countries,
      query_specificity_status = specificity_assessment$status,
      query_specificity_reason = specificity_assessment$reason,
      status = "ok",
      note = stop_reason
    ),
    records = fetched_records
  )
}
