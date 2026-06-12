# ------------------------------------------------------------------------------
# Filter VectorMap host-vector links to hosts present in the WHO network
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

# Clean text while keeping raw source values as intact as possible.
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

# Build a normalized key for exact and fuzzy host matching.
normalize_host_key <- function(x) {
  x <- clean_text(x)

  transliterated <- suppressWarnings(iconv(x, from = "", to = "ASCII//TRANSLIT"))
  transliterated[is.na(transliterated)] <- x[is.na(transliterated)]

  transliterated %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9 ]+", " ") %>%
    str_squish() %>%
    na_if("")
}

# Count word-like tokens in a host label after normalization.
count_host_words <- function(x) {
  key <- normalize_host_key(x)
  if_else(is.na(key), 0L, str_count(key, " ") + 1L)
}

# Flag broad or non-actionable labels that should not pass directly into the
# final species-focused host-vector table.
is_non_actionable_host <- function(x) {
  key <- normalize_host_key(x)
  word_count <- count_host_words(x)

  is.na(key) |
    key == "" |
    word_count < 2L |
    str_detect(
      key,
      "\\b(sp|sp\\.|spp|spp\\.|species complex|complex|group|subgroup|unknown|unidentified|organism|organisms|mixed|multiple|multi|other)\\b"
    ) |
    str_detect(key, "\\band\\b|/|,")
}

# Parse scientific-name strings down to a binomial matching key while keeping
# stripped infraspecific and authorship detail in separate columns.
canonicalize_scientific_name <- function(x) {
  cleaned <- clean_text(x)
  cleaned <- str_replace_all(cleaned, "_", " ")
  cleaned <- str_squish(cleaned)
  cleaned[cleaned == ""] <- NA_character_

  qualifier_tokens <- c(
    "sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.",
    "nr", "nr.", "group", "complex", "subgroup"
  )

  parsed <- lapply(cleaned, function(value) {
    if (is.na(value)) {
      return(list(
        host_scientific_name_clean = NA_character_,
        host_binomial = NA_character_,
        host_infraspecific_epithet = NA_character_,
        host_authorship_year_suffix = NA_character_,
        host_canonicalization_flag = "missing"
      ))
    }

    tokens <- str_split(value, "\\s+", simplify = TRUE)
    tokens <- tokens[tokens != ""]

    if (length(tokens) < 2) {
      return(list(
        host_scientific_name_clean = value,
        host_binomial = NA_character_,
        host_infraspecific_epithet = NA_character_,
        host_authorship_year_suffix = NA_character_,
        host_canonicalization_flag = "non_actionable"
      ))
    }

    genus <- tokens[[1]]
    species <- tokens[[2]]
    species_key <- str_to_lower(species)

    if (species_key %in% qualifier_tokens) {
      return(list(
        host_scientific_name_clean = value,
        host_binomial = NA_character_,
        host_infraspecific_epithet = NA_character_,
        host_authorship_year_suffix = NA_character_,
        host_canonicalization_flag = "non_actionable"
      ))
    }

    host_binomial <- paste(genus, species)

    infraspecific <- NA_character_
    suffix_start <- 3L

    if (length(tokens) >= 3) {
      third_token <- tokens[[3]]
      third_key <- str_to_lower(third_token)

      if (
        str_detect(third_token, "^[a-z-]+$") &&
        !third_key %in% qualifier_tokens
      ) {
        infraspecific <- third_token
        suffix_start <- 4L
      }
    }

    suffix_tokens <- if (length(tokens) >= suffix_start) {
      tokens[suffix_start:length(tokens)]
    } else {
      character(0)
    }

    suffix <- if (length(suffix_tokens) > 0) {
      paste(suffix_tokens, collapse = " ")
    } else {
      NA_character_
    }

    flag <- case_when(
      !is.na(infraspecific) & !is.na(suffix) ~ "trinomial_and_suffix",
      !is.na(infraspecific) ~ "trinomial",
      !is.na(suffix) ~ "authorship_or_suffix",
      TRUE ~ "already_binomial"
    )

    list(
      host_scientific_name_clean = value,
      host_binomial = host_binomial,
      host_infraspecific_epithet = infraspecific,
      host_authorship_year_suffix = suffix,
      host_canonicalization_flag = flag
    )
  })

  tibble(
    host_scientific_name_clean = vapply(parsed, `[[`, character(1), "host_scientific_name_clean"),
    host_binomial = vapply(parsed, `[[`, character(1), "host_binomial"),
    host_infraspecific_epithet = vapply(parsed, `[[`, character(1), "host_infraspecific_epithet"),
    host_authorship_year_suffix = vapply(parsed, `[[`, character(1), "host_authorship_year_suffix"),
    host_canonicalization_flag = vapply(parsed, `[[`, character(1), "host_canonicalization_flag")
  ) %>%
    mutate(across(everything(), ~ na_if(.x, "")))
}

# Restrict verbatim rescue to labels that already look like clean Latin binomials.
looks_like_latin_binomial <- function(x) {
  cleaned <- clean_text(x)
  word_count <- count_host_words(cleaned)
  tokens <- str_split(cleaned, "\\s+", simplify = TRUE)
  key_tokens <- str_split(normalize_host_key(cleaned), " ", simplify = TRUE)
  qualifier_tokens <- c("sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.")

  !is.na(cleaned) &
    word_count == 2L &
    str_detect(tokens[, 1], "^[A-Z][a-z-]+$") &
    str_detect(tokens[, 2], "^[a-z-]+$") &
    !key_tokens[, 2] %in% qualifier_tokens
}

# Return the first non-missing value from a vector.
first_non_missing <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

# Parse manual crosswalk flags consistently from empty, logical, or text values.
parse_logical_flag <- function(x) {
  x <- clean_text(as.character(x))
  x <- str_to_lower(x)

  case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ FALSE
  )
}

# Suggest one fuzzy WHO host candidate per unresolved raw label.
suggest_fuzzy_match <- function(raw_key, who_lookup, min_score = 0.75) {
  if (is.na(raw_key) || raw_key == "" || nchar(raw_key) < 5) {
    return(list(suggested_who_host = NA_character_, fuzzy_score = NA_real_))
  }

  raw_tokens <- str_split(raw_key, " ", simplify = TRUE)
  raw_tokens <- raw_tokens[raw_tokens != ""]
  raw_first <- raw_tokens[[1]]
  raw_token_count <- length(raw_tokens)

  candidate_pool <- who_lookup %>%
    filter(abs(token_count - raw_token_count) <= 1)

  if (nrow(candidate_pool) == 0) {
    candidate_pool <- who_lookup
  }

  genus_pool <- candidate_pool %>%
    filter(word(who_host_key, 1) == raw_first)

  if (nrow(genus_pool) > 0) {
    candidate_pool <- genus_pool
  }

  distances <- utils::adist(raw_key, candidate_pool$who_host_key, ignore.case = TRUE)
  distances <- as.numeric(distances[1, ])
  scores <- 1 - distances / pmax(nchar(raw_key), nchar(candidate_pool$who_host_key))

  best_index <- which.max(scores)
  best_score <- scores[[best_index]]

  if (length(best_score) == 0 || is.na(best_score) || best_score < min_score) {
    return(list(suggested_who_host = NA_character_, fuzzy_score = best_score))
  }

  list(
    suggested_who_host = candidate_pool$matched_who_host[[best_index]],
    fuzzy_score = round(best_score, 3)
  )
}

outputs_dir <- vectormap_outputs_dir
manual_dir <- vectormap_manual_dir
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)

raw_links_path <- file.path(outputs_dir, "vectormap_vector_host_links_raw.csv")
manual_crosswalk_path <- file.path(manual_dir, "vectormap_host_manual_crosswalk.csv")

exact_output_path <- file.path(outputs_dir, "vectormap_vector_host_links_who_exact.csv")
review_output_path <- file.path(outputs_dir, "vectormap_host_crosswalk_review.csv")
filtered_output_path <- file.path(outputs_dir, "vectormap_vector_host_links_who_filtered.csv")
package_candidates_path <- file.path(outputs_dir, "vectormap_host_package_candidates.csv")
manual_candidates_path <- file.path(outputs_dir, "vectormap_host_manual_crosswalk_candidates.csv")

if (!file.exists(manual_crosswalk_path)) {
  manual_crosswalk_template <- tibble(
    raw_host_label = character(),
    source_field = character(),
    matched_who_host = character(),
    include_in_final = logical(),
    note = character()
  )

  write_csv(manual_crosswalk_template, manual_crosswalk_path, na = "")
}

# Load the authoritative WHO host universe and keep one row per host.
who_hosts <- read_legacy_compatible_master_plus_network() %>%
  transmute(
    matched_who_host = clean_text(Host),
    matched_who_host_tax_id = clean_text(HostTaxID),
    matched_who_host_class = clean_text(HostClass),
    matched_who_host_order = clean_text(HostOrder),
    matched_who_host_family = clean_text(HostFamily)
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  group_by(matched_who_host) %>%
  summarise(
    matched_who_host_tax_id = first_non_missing(matched_who_host_tax_id),
    matched_who_host_class = first_non_missing(matched_who_host_class),
    matched_who_host_order = first_non_missing(matched_who_host_order),
    matched_who_host_family = first_non_missing(matched_who_host_family),
    .groups = "drop"
  ) %>%
  mutate(
    who_host_key = normalize_host_key(matched_who_host),
    token_count = str_count(who_host_key, " ") + 1L
  )

duplicate_who_keys <- who_hosts %>%
  count(who_host_key, sort = TRUE) %>%
  filter(!is.na(who_host_key), n > 1)

if (nrow(duplicate_who_keys) > 0) {
  duplicate_examples <- who_hosts %>%
    semi_join(duplicate_who_keys, by = "who_host_key") %>%
    group_by(who_host_key) %>%
    summarise(
      host_options = paste(sort(unique(matched_who_host)), collapse = "; "),
      taxid_options = paste(
        sort(unique(stats::na.omit(matched_who_host_tax_id))),
        collapse = "; "
      ),
      .groups = "drop"
    )

  cat(
    "Collapsing",
    nrow(duplicate_who_keys),
    "duplicated normalized WHO host keys in the VectorMap lookup.\n"
  )
  print(duplicate_examples, n = Inf)

  who_hosts <- who_hosts %>%
    mutate(
      .taxonomy_field_count =
        as.integer(!is.na(matched_who_host_tax_id)) +
        as.integer(!is.na(matched_who_host_class)) +
        as.integer(!is.na(matched_who_host_order)) +
        as.integer(!is.na(matched_who_host_family)),
      .title_case_label = str_detect(matched_who_host, "^[[:upper:]]")
    ) %>%
    arrange(
      who_host_key,
      desc(.taxonomy_field_count),
      desc(.title_case_label),
      matched_who_host,
      matched_who_host_tax_id
    ) %>%
    group_by(who_host_key) %>%
    slice(1) %>%
    ungroup() %>%
    select(-.taxonomy_field_count, -.title_case_label)
}

# Load the raw VectorMap host-vector evidence and preserve one row id per record.
vectormap_raw <- read_csv(
  raw_links_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  bind_cols(canonicalize_scientific_name(.$host_scientific_name)) %>%
  mutate(
    vectormap_row_id = row_number(),
    host_scientific_name = clean_text(host_scientific_name),
    verbatim_host_name = clean_text(verbatim_host_name),
    host_scientific_clean_key = normalize_host_key(host_scientific_name_clean),
    host_binomial_key = normalize_host_key(host_binomial),
    verbatim_host_key = normalize_host_key(verbatim_host_name),
    scientific_non_actionable = is.na(host_binomial) | is_non_actionable_host(host_scientific_name_clean),
    scientific_missing = is.na(host_scientific_name_clean),
    verbatim_non_actionable = is_non_actionable_host(verbatim_host_name),
    verbatim_binomial = if_else(
      looks_like_latin_binomial(verbatim_host_name),
      clean_text(verbatim_host_name),
      NA_character_
    ),
    verbatim_binomial_key = normalize_host_key(verbatim_binomial)
  )

raw_row_count <- nrow(vectormap_raw)

# Apply exact scientific-name matching on the canonicalized binomial first.
exact_scientific_binomial <- vectormap_raw %>%
  filter(!is.na(host_binomial_key)) %>%
  left_join(
    who_hosts,
    by = c("host_binomial_key" = "who_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  mutate(
    match_method = "exact_scientific_binomial",
    match_source_field = "host_binomial"
  )

remaining_after_scientific_binomial <- vectormap_raw %>%
  filter(!vectormap_row_id %in% exact_scientific_binomial$vectormap_row_id)

# Then try the cleaned raw scientific label for any rows not caught by the
# binomial key.
exact_scientific_raw <- remaining_after_scientific_binomial %>%
  filter(!scientific_non_actionable, !is.na(host_scientific_clean_key)) %>%
  left_join(
    who_hosts,
    by = c("host_scientific_clean_key" = "who_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  mutate(
    match_method = "exact_scientific_raw",
    match_source_field = "host_scientific_name_clean"
  )

exact_matches <- bind_rows(
  exact_scientific_binomial,
  exact_scientific_raw
) %>%
  arrange(vectormap_row_id)

# Load any reviewed manual crosswalk entries and apply them after exact matching.
manual_crosswalk <- read_csv(
  manual_crosswalk_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    raw_host_label = clean_text(raw_host_label),
    raw_host_key = normalize_host_key(raw_host_label),
    source_field = clean_text(source_field),
    matched_who_host_raw = clean_text(matched_who_host),
    matched_who_host_key = normalize_host_key(matched_who_host),
    include_in_final = parse_logical_flag(include_in_final),
    note = clean_text(note)
  ) %>%
  filter(
    !is.na(raw_host_key),
    !is.na(source_field),
    !is.na(matched_who_host_key),
    include_in_final
  ) %>%
  left_join(
    who_hosts %>%
      select(
        matched_who_host_key = who_host_key,
        matched_who_host_resolved = matched_who_host
      ),
    by = "matched_who_host_key"
  )

invalid_manual_sources <- manual_crosswalk %>%
  filter(!source_field %in% c("host_binomial", "host_scientific_name_clean", "verbatim_binomial"))

if (nrow(invalid_manual_sources) > 0) {
  stop("Manual crosswalk source_field must be host_binomial, host_scientific_name_clean, or verbatim_binomial.")
}

invalid_manual_targets <- manual_crosswalk %>%
  filter(is.na(matched_who_host_resolved))

if (nrow(invalid_manual_targets) > 0) {
  cat(
    "Skipping",
    nrow(invalid_manual_targets),
    "manual crosswalk rows with targets outside the current WHO host universe.\n"
  )
  print(
    invalid_manual_targets %>%
      select(raw_host_label, source_field, matched_who_host_raw, note),
    n = Inf
  )
}

manual_crosswalk <- manual_crosswalk %>%
  filter(!is.na(matched_who_host_resolved)) %>%
  mutate(matched_who_host = matched_who_host_resolved) %>%
  select(-matched_who_host_resolved)

remaining_after_exact <- vectormap_raw %>%
  filter(!vectormap_row_id %in% exact_matches$vectormap_row_id)

manual_scientific_binomial <- remaining_after_exact %>%
  filter(!is.na(host_binomial_key)) %>%
  left_join(
    manual_crosswalk %>%
      filter(source_field == "host_binomial") %>%
      select(raw_host_key, matched_who_host),
    by = c("host_binomial_key" = "raw_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  left_join(who_hosts, by = "matched_who_host") %>%
  mutate(
    match_method = "manual_crosswalk",
    match_source_field = "host_binomial"
  )

remaining_after_manual_scientific_binomial <- remaining_after_exact %>%
  filter(!vectormap_row_id %in% manual_scientific_binomial$vectormap_row_id)

manual_scientific_raw <- remaining_after_manual_scientific_binomial %>%
  filter(!is.na(host_scientific_clean_key)) %>%
  left_join(
    manual_crosswalk %>%
      filter(source_field == "host_scientific_name_clean") %>%
      select(raw_host_key, matched_who_host),
    by = c("host_scientific_clean_key" = "raw_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  left_join(who_hosts, by = "matched_who_host") %>%
  mutate(
    match_method = "manual_crosswalk",
    match_source_field = "host_scientific_name_clean"
  )

remaining_after_manual_scientific <- remaining_after_manual_scientific_binomial %>%
  filter(!vectormap_row_id %in% manual_scientific_raw$vectormap_row_id)

manual_verbatim_binomial <- remaining_after_manual_scientific %>%
  filter(!is.na(verbatim_binomial_key)) %>%
  left_join(
    manual_crosswalk %>%
      filter(source_field == "verbatim_binomial") %>%
      select(raw_host_key, matched_who_host),
    by = c("verbatim_binomial_key" = "raw_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  left_join(who_hosts, by = "matched_who_host") %>%
  mutate(
    match_method = "manual_crosswalk",
    match_source_field = "verbatim_binomial"
  )

remaining_after_manual_verbatim <- remaining_after_manual_scientific %>%
  filter(!vectormap_row_id %in% manual_verbatim_binomial$vectormap_row_id)

# Only use verbatim names as an exact rescue when the scientific label is
# missing or non-actionable and the verbatim value is already a clean binomial.
exact_verbatim_binomial_rescue <- remaining_after_manual_verbatim %>%
  filter(
    (scientific_missing | scientific_non_actionable),
    !is.na(verbatim_binomial_key)
  ) %>%
  left_join(
    who_hosts,
    by = c("verbatim_binomial_key" = "who_host_key")
  ) %>%
  filter(!is.na(matched_who_host)) %>%
  mutate(
    match_method = "exact_verbatim_binomial_rescue",
    match_source_field = "verbatim_binomial"
  )

filtered_matches <- bind_rows(
  exact_matches,
  manual_scientific_binomial,
  manual_scientific_raw,
  manual_verbatim_binomial,
  exact_verbatim_binomial_rescue
) %>%
  arrange(vectormap_row_id)

# Build a review table from unresolved labels after exact and manual matching.
remaining_after_manual <- vectormap_raw %>%
  filter(!vectormap_row_id %in% filtered_matches$vectormap_row_id)

review_labels <- bind_rows(
  remaining_after_manual %>%
    transmute(
      raw_host_label = host_scientific_name,
      source_field = "host_scientific_name",
      scientific_name_clean = host_scientific_name_clean,
      host_binomial = host_binomial,
      host_infraspecific_epithet = host_infraspecific_epithet,
      host_authorship_year_suffix = host_authorship_year_suffix,
      host_canonicalization_flag = host_canonicalization_flag,
      review_match_key = coalesce(host_binomial_key, host_scientific_clean_key),
      word_count = count_host_words(host_scientific_name),
      binomial_like = !is.na(host_binomial)
    ),
  remaining_after_manual %>%
    transmute(
      raw_host_label = verbatim_host_name,
      source_field = "verbatim_host_name",
      scientific_name_clean = NA_character_,
      host_binomial = verbatim_binomial,
      host_infraspecific_epithet = NA_character_,
      host_authorship_year_suffix = NA_character_,
      host_canonicalization_flag = NA_character_,
      review_match_key = verbatim_binomial_key,
      word_count = count_host_words(verbatim_host_name),
      binomial_like = !is.na(verbatim_binomial)
    )
) %>%
  mutate(
    raw_host_label = clean_text(raw_host_label),
    raw_host_key = normalize_host_key(raw_host_label),
    non_actionable = case_when(
      source_field == "host_scientific_name" ~ is.na(host_binomial) | is_non_actionable_host(scientific_name_clean),
      source_field == "verbatim_host_name" ~ is.na(host_binomial),
      TRUE ~ is_non_actionable_host(raw_host_label)
    )
  ) %>%
  group_by(
    raw_host_label,
    source_field,
    raw_host_key,
    scientific_name_clean,
    host_binomial,
    host_infraspecific_epithet,
    host_authorship_year_suffix,
    host_canonicalization_flag,
    review_match_key,
    word_count,
    binomial_like,
    non_actionable
  ) %>%
  summarise(row_count = n(), .groups = "drop")

review_missing <- review_labels %>%
  filter(is.na(raw_host_label)) %>%
  transmute(
    raw_host_label = NA_character_,
    source_field,
    scientific_name_clean,
    host_binomial,
    host_infraspecific_epithet,
    host_authorship_year_suffix,
    host_canonicalization_flag,
    review_bucket = "drop_from_species_level",
    package_check_target = FALSE,
    suggested_who_host = NA_character_,
    fuzzy_score = NA_real_,
    reason_flag = "missing_host_label",
    row_count
  )

review_non_missing <- review_labels %>%
  filter(!is.na(raw_host_label)) %>%
  mutate(
    fuzzy_result = lapply(
      review_match_key,
      suggest_fuzzy_match,
      who_lookup = who_hosts
    ),
    suggested_who_host = vapply(
      fuzzy_result,
      function(x) x$suggested_who_host,
      character(1)
    ),
    fuzzy_score = vapply(
      fuzzy_result,
      function(x) {
        if (is.null(x$fuzzy_score)) {
          return(NA_real_)
        }

        as.numeric(x$fuzzy_score)
      },
      numeric(1)
    ),
    review_bucket = case_when(
      non_actionable ~ "drop_from_species_level",
      source_field == "host_scientific_name" & !is.na(host_binomial) ~ "manual_crosswalk_candidate",
      source_field == "verbatim_host_name" & binomial_like ~ "package_candidate_check",
      TRUE ~ "drop_from_species_level"
    ),
    package_check_target = review_bucket %in% c(
      "manual_crosswalk_candidate",
      "package_candidate_check"
    ),
    reason_flag = case_when(
      review_bucket == "drop_from_species_level" ~ "drop_from_species_level",
      source_field == "host_scientific_name" &
        !is.na(host_canonicalization_flag) &
        host_canonicalization_flag != "already_binomial" ~ "needs_binomialization",
      source_field == "host_scientific_name" ~ "needs_manual_synonym_map",
      source_field == "verbatim_host_name" & binomial_like ~ "verbatim_binomial_candidate",
      TRUE ~ "drop_from_species_level"
    )
  ) %>%
  transmute(
    raw_host_label,
    source_field,
    scientific_name_clean,
    host_binomial,
    host_infraspecific_epithet,
    host_authorship_year_suffix,
    host_canonicalization_flag,
    review_bucket,
    package_check_target,
    suggested_who_host,
    fuzzy_score,
    reason_flag,
    row_count
  )

crosswalk_review <- bind_rows(review_missing, review_non_missing) %>%
  arrange(desc(row_count), source_field, raw_host_label)

package_candidates <- crosswalk_review %>%
  filter(package_check_target) %>%
  arrange(desc(row_count), source_field, raw_host_label)

manual_crosswalk_candidates <- crosswalk_review %>%
  filter(review_bucket == "manual_crosswalk_candidate") %>%
  arrange(desc(row_count), raw_host_label)

unresolved_row_count <- nrow(remaining_after_manual)

output_match_cols <- c(
  "vectormap_row_id",
  "matched_who_host",
  "matched_who_host_tax_id",
  "matched_who_host_class",
  "matched_who_host_order",
  "matched_who_host_family",
  "host_scientific_name_clean",
  "host_binomial",
  "host_infraspecific_epithet",
  "host_authorship_year_suffix",
  "host_canonicalization_flag",
  "verbatim_binomial",
  "match_method",
  "match_source_field"
)

exact_output <- exact_matches %>%
  select(all_of(output_match_cols), everything())

filtered_output <- filtered_matches %>%
  select(all_of(output_match_cols), everything())

if (!all(filtered_output$matched_who_host %in% who_hosts$matched_who_host)) {
  stop("Filtered output contains matched hosts outside the WHO host universe.")
}

if (any(filtered_output$match_method == "fuzzy_candidate", na.rm = TRUE)) {
  stop("Fuzzy candidates should not appear in the final filtered output.")
}

write_csv(exact_output, exact_output_path, na = "")
write_csv(crosswalk_review, review_output_path, na = "")
write_csv(filtered_output, filtered_output_path, na = "")
write_csv(package_candidates, package_candidates_path, na = "")
write_csv(manual_crosswalk_candidates, manual_candidates_path, na = "")

cat("Raw VectorMap rows:", raw_row_count, "\n")
cat("Unique WHO hosts:", nrow(who_hosts), "\n")
cat(
  "Exact scientific binomial matches:",
  nrow(exact_scientific_binomial),
  "rows across",
  n_distinct(exact_scientific_binomial$matched_who_host),
  "hosts\n"
)
cat(
  "Exact scientific raw matches:",
  nrow(exact_scientific_raw),
  "rows across",
  n_distinct(exact_scientific_raw$matched_who_host),
  "hosts\n"
)
cat(
  "Manual scientific crosswalk matches:",
  nrow(bind_rows(manual_scientific_binomial, manual_scientific_raw)),
  "rows across",
  n_distinct(bind_rows(manual_scientific_binomial, manual_scientific_raw)$matched_who_host),
  "hosts\n"
)
cat(
  "Exact verbatim binomial rescue matches:",
  nrow(exact_verbatim_binomial_rescue),
  "rows across",
  n_distinct(exact_verbatim_binomial_rescue$matched_who_host),
  "hosts\n"
)
cat("Unresolved rows sent to crosswalk review:", unresolved_row_count, "\n")
cat("Crosswalk review entries:", nrow(crosswalk_review), "\n")
cat("Package-check candidate entries:", nrow(package_candidates), "\n")
cat("Manual crosswalk candidate entries:", nrow(manual_crosswalk_candidates), "\n")
cat("Final filtered rows:", nrow(filtered_output), "\n")
cat("Wrote exact matches to", exact_output_path, "\n")
cat("Wrote review table to", review_output_path, "\n")
cat("Wrote filtered output to", filtered_output_path, "\n")
cat("Wrote package candidate table to", package_candidates_path, "\n")
cat("Wrote manual crosswalk candidate table to", manual_candidates_path, "\n")
