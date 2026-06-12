suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

# ------------------------------------------------------------------------------
# Clean EFSA appendices A and G, then crosswalk EFSA pathogens to the combined
# WHO host-pathogen network using conservative, auditable matching rules.
# ------------------------------------------------------------------------------

snake_case_names <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "([a-z0-9])([A-Z])", "\\1_\\2")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_trim(x)
  x <- str_to_lower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  x <- ifelse(is.na(x) | x == "", "x", x)
  make.unique(x, sep = "_")
}

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

normalize_yes_no <- function(x) {
  x <- str_to_lower(clean_text(x))
  case_when(
    x %in% c("yes", "y", "true", "1") ~ "yes",
    x %in% c("no", "n", "false", "0") ~ "no",
    TRUE ~ x
  )
}

normalize_name_for_match <- function(x) {
  x <- clean_text(x)
  x <- str_to_lower(x)
  x <- str_replace_all(x, "haemorrh", "hemorrh")
  x <- str_replace_all(x, "&", " and ")
  x <- str_replace_all(x, "[/]", " ")
  x <- str_replace_all(x, "[-–—]", " ")
  x <- str_replace_all(x, "[()\\[\\],.;:*'`\"]", " ")
  x <- str_replace_all(x, "\\bviruses\\b", "virus")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  x
}

remove_virus_suffix <- function(x) {
  x <- clean_text(x)
  x <- str_replace(x, regex("\\s+virus(?:es)?$", ignore_case = TRUE), "")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

normalize_status <- function(x) {
  x <- str_to_lower(clean_text(x))
  case_when(
    x %in% c("highly likely", "highly_likely") ~ "highly_likely",
    x %in% c("potential") ~ "potential",
    TRUE ~ x
  )
}

normalize_vector_group <- function(x) {
  x <- str_to_upper(clean_text(x))
  case_when(
    x == "BM" ~ "biting_midge",
    x == "M" ~ "mosquito",
    x == "T" ~ "tick",
    x == "SF" ~ "sand_fly",
    TRUE ~ str_to_lower(x)
  )
}

normalize_positive_field <- function(x) {
  x_raw <- clean_text(x)
  x_clean <- str_to_lower(x_raw)
  case_when(
    is.na(x_clean) ~ NA_character_,
    str_detect(x_clean, "positive") ~ "positive",
    TRUE ~ x_clean
  )
}

drop_blank_rows_and_cols <- function(df) {
  blank_cols <- vapply(
    df,
    function(col) all(is.na(clean_text(col))),
    logical(1)
  )
  df <- df[, !blank_cols, drop = FALSE]

  if (nrow(df) == 0) {
    return(df)
  }

  blank_rows <- apply(
    df,
    1,
    function(row) all(is.na(clean_text(row)))
  )

  df[!blank_rows, , drop = FALSE]
}

detect_header_row <- function(path, expected_tokens, sheet = 1, preview_rows = 10) {
  preview <- read_excel(
    path,
    sheet = sheet,
    col_names = FALSE,
    n_max = preview_rows
  )

  expected_tokens <- normalize_name_for_match(expected_tokens)

  row_scores <- apply(
    preview,
    1,
    function(row_values) {
      row_values <- normalize_name_for_match(row_values)
      sum(vapply(
        expected_tokens,
        function(token) any(!is.na(row_values) & str_detect(row_values, fixed(token))),
        logical(1)
      ))
    }
  )

  which.max(row_scores)
}

read_excel_detect_header <- function(path, expected_tokens, sheet = 1) {
  header_row <- detect_header_row(
    path = path,
    expected_tokens = expected_tokens,
    sheet = sheet
  )

  dat <- read_excel(
    path,
    sheet = sheet,
    skip = header_row - 1
  )

  list(data = dat, header_row = header_row)
}

safe_left_join_unique <- function(x, y, by, prefix) {
  y_match <- y %>%
    add_count(across(all_of(by)), name = paste0(prefix, "_candidate_count"))

  x %>%
    left_join(y_match, by = by)
}

load_manual_map <- function(path) {
  if (!file.exists(path)) {
    stop("Manual mapping file not found: ", path)
  }

  read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(
      across(everything(), clean_text),
      source_name_clean = normalize_name_for_match(source_name),
      canonical_name_clean = normalize_name_for_match(canonical_name)
    )
}

apply_manual_name_override <- function(name_vector, manual_map) {
  cleaned <- normalize_name_for_match(name_vector)
  override_lookup <- manual_map %>%
    filter(applies_to %in% c("efsa_internal", "both")) %>%
    filter(!is.na(source_name_clean), !is.na(canonical_name)) %>%
    distinct(source_name_clean, canonical_name)

  matched <- override_lookup$canonical_name[
    match(cleaned, override_lookup$source_name_clean)
  ]

  ifelse(!is.na(matched), matched, clean_text(name_vector))
}

build_best_fuzzy_candidate <- function(efsa_clean, combined_reference) {
  if (is.na(efsa_clean) || efsa_clean == "") {
    return(tibble())
  }

  pathogen_pool <- combined_reference %>%
    filter(!is.na(pathogen_clean)) %>%
    transmute(
      candidate_type = "pathogen",
      combined_pathogen_raw,
      combined_pathogen_clean = pathogen_clean,
      combined_disease_name_raw,
      candidate_clean = pathogen_clean
    )

  disease_pool <- combined_reference %>%
    filter(!is.na(disease_name_clean)) %>%
    transmute(
      candidate_type = "disease",
      combined_pathogen_raw,
      combined_pathogen_clean = pathogen_clean,
      combined_disease_name_raw,
      candidate_clean = disease_name_clean
    )

  candidates <- bind_rows(pathogen_pool, disease_pool) %>%
    distinct()

  if (nrow(candidates) == 0) {
    return(tibble())
  }

  distances <- adist(efsa_clean, candidates$candidate_clean, ignore.case = TRUE)

  candidates %>%
    mutate(
      distance = as.numeric(distances[1, ]),
      relative_distance = distance / pmax(nchar(efsa_clean), nchar(candidate_clean)),
      shared_tokens = vapply(
        candidate_clean,
        function(candidate) {
          efsa_tokens <- unique(str_split(efsa_clean, " ", simplify = TRUE))
          candidate_tokens <- unique(str_split(candidate, " ", simplify = TRUE))
          sum(efsa_tokens != "" & efsa_tokens %in% candidate_tokens)
        },
        numeric(1)
      )
    ) %>%
    arrange(relative_distance, distance, desc(shared_tokens)) %>%
    slice(1) %>%
    filter(relative_distance <= 0.12, distance <= 4, shared_tokens >= 2)
}

output_dir <- vector_screening_efsa_outputs_dir

appendix_a_path <- vector_screening_efsa_source_path("efsa_report_appendix_a.xlsx")
appendix_g_path <- vector_screening_efsa_source_path("efsa_report_appendix_g.xlsx")
screening_path <- vector_screening_manual_path("disease_vector_screening.csv")
manual_map_path <- vector_screening_efsa_manual_path("efsa_name_manual_map.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manual_map <- load_manual_map(manual_map_path)

appendix_a_raw <- read_excel_detect_header(
  appendix_a_path,
  expected_tokens = c(
    "abbreviation",
    "species",
    "isolate virus name",
    "family",
    "genus",
    "oie notifiable"
  )
)

appendix_g_raw <- read_excel_detect_header(
  appendix_g_path,
  expected_tokens = c(
    "vector borne pathogen",
    "arthropod species",
    "vector group",
    "vector status",
    "field results",
    "laboratory results"
  )
)

appendix_a_clean <- appendix_a_raw$data %>%
  setNames(snake_case_names(names(.))) %>%
  mutate(across(everything(), clean_text)) %>%
  drop_blank_rows_and_cols() %>%
  rename(
    abbreviation = abbreviation,
    species_name = species,
    isolate_virus_name = isolate_virus_name,
    pathogen_family = family,
    pathogen_genus = genus,
    oie_notifiable_raw = oie_notifiable
  ) %>%
  tidyr::fill(pathogen_family, pathogen_genus, oie_notifiable_raw, .direction = "down") %>%
  mutate(
    species_name = apply_manual_name_override(species_name, manual_map),
    isolate_virus_name = apply_manual_name_override(isolate_virus_name, manual_map),
    abbreviation = clean_text(abbreviation),
    oie_notifiable = normalize_yes_no(oie_notifiable_raw),
    species_name_clean = normalize_name_for_match(species_name),
    isolate_virus_name_clean = normalize_name_for_match(isolate_virus_name),
    abbreviation_clean = normalize_name_for_match(abbreviation),
    efsa_pathogen_raw = coalesce(isolate_virus_name, species_name),
    efsa_pathogen_clean = normalize_name_for_match(efsa_pathogen_raw),
    efsa_disease_candidate = remove_virus_suffix(efsa_pathogen_raw),
    efsa_disease_candidate_clean = normalize_name_for_match(efsa_disease_candidate),
    efsa_pathogen_id = if_else(
      !is.na(abbreviation_clean),
      paste0("efsa_", str_replace_all(abbreviation_clean, " ", "_")),
      paste0("efsa_", str_replace_all(efsa_pathogen_clean, " ", "_"))
    )
  ) %>%
  select(
    efsa_pathogen_id,
    abbreviation,
    abbreviation_clean,
    species_name,
    species_name_clean,
    isolate_virus_name,
    isolate_virus_name_clean,
    efsa_pathogen_raw,
    efsa_pathogen_clean,
    efsa_disease_candidate,
    efsa_disease_candidate_clean,
    pathogen_family,
    pathogen_genus,
    oie_notifiable,
    oie_notifiable_raw
  ) %>%
  distinct()

appendix_a_aliases <- bind_rows(
  appendix_a_clean %>%
    transmute(
      efsa_pathogen_id,
      alias_source = "species_name",
      alias_raw = species_name,
      alias_clean = species_name_clean
    ),
  appendix_a_clean %>%
    transmute(
      efsa_pathogen_id,
      alias_source = "isolate_virus_name",
      alias_raw = isolate_virus_name,
      alias_clean = isolate_virus_name_clean
    ),
  appendix_a_clean %>%
    transmute(
      efsa_pathogen_id,
      alias_source = "abbreviation",
      alias_raw = abbreviation,
      alias_clean = abbreviation_clean
    )
) %>%
  filter(!is.na(alias_raw), !is.na(alias_clean)) %>%
  mutate(
    alias_disease_candidate = remove_virus_suffix(alias_raw),
    alias_disease_candidate_clean = normalize_name_for_match(alias_disease_candidate)
  ) %>%
  distinct(efsa_pathogen_id, alias_source, alias_clean, .keep_all = TRUE)

appendix_g_clean <- appendix_g_raw$data %>%
  setNames(snake_case_names(names(.))) %>%
  mutate(across(everything(), clean_text)) %>%
  drop_blank_rows_and_cols() %>%
  rename(
    efsa_pathogen_raw = vector_borne_pathogen,
    vector_species_raw = arthropod_species,
    vector_group_raw = vector_group,
    vector_status_raw = vector_status,
    field_results_raw = field_results,
    laboratory_results_raw = laboratory_results,
    vn_database_species_raw = vn_database_species,
    vn_priority_species_raw = vn_priority_species
  ) %>%
  mutate(
    efsa_pathogen_standardized = apply_manual_name_override(efsa_pathogen_raw, manual_map),
    efsa_pathogen_clean = normalize_name_for_match(efsa_pathogen_standardized),
    vector_species = str_replace(clean_text(vector_species_raw), "\\*+$", ""),
    vector_species_clean = normalize_name_for_match(vector_species),
    vector_group = normalize_vector_group(vector_group_raw),
    vector_status = normalize_status(vector_status_raw),
    field_results = normalize_positive_field(field_results_raw),
    laboratory_results = normalize_positive_field(laboratory_results_raw),
    vn_database_species = normalize_yes_no(vn_database_species_raw),
    vn_priority_species = normalize_yes_no(vn_priority_species_raw)
  ) %>%
  filter(!is.na(efsa_pathogen_clean), !is.na(vector_species_clean)) %>%
  distinct()

appendix_g_clean <- appendix_g_clean %>%
  left_join(
    appendix_a_aliases %>%
      distinct(alias_clean, efsa_pathogen_id, .keep_all = TRUE) %>%
      add_count(alias_clean, name = "appendix_a_alias_match_count"),
    by = c("efsa_pathogen_clean" = "alias_clean")
  ) %>%
  left_join(
    appendix_a_clean %>%
      select(
        efsa_pathogen_id,
        abbreviation,
        species_name,
        isolate_virus_name,
        pathogen_family,
        pathogen_genus,
        oie_notifiable
      ),
    by = "efsa_pathogen_id"
  ) %>%
  mutate(
    appendix_a_match_status = case_when(
      is.na(efsa_pathogen_id) ~ "unmatched_to_appendix_a",
      appendix_a_alias_match_count > 1 ~ "ambiguous_appendix_a_match",
      TRUE ~ "matched_to_appendix_a"
    )
  )

combined_network_raw <- read_legacy_compatible_master_plus_network()

combined_network_clean <- combined_network_raw %>%
  setNames(snake_case_names(names(.))) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    pathogen_raw = pathogen,
    disease_name_raw = disease_name,
    pathogen_family_raw = pathogen_family,
    pathogen_genus_raw = pathogen_genus,
    pathogen_type_raw = pathogen_type,
    pathogen_clean = normalize_name_for_match(pathogen_raw),
    disease_name_clean = normalize_name_for_match(disease_name_raw),
    pathogen_family_clean = normalize_name_for_match(pathogen_family_raw),
    pathogen_genus_clean = normalize_name_for_match(pathogen_genus_raw),
    pathogen_type_clean = normalize_name_for_match(pathogen_type_raw)
  )

combined_reference <- combined_network_clean %>%
  distinct(
    combined_pathogen_raw = pathogen_raw,
    pathogen_clean,
    combined_disease_name_raw = disease_name_raw,
    disease_name_clean,
    pathogen_family_raw,
    pathogen_genus_raw,
    pathogen_type_raw
  ) %>%
  filter(!is.na(combined_pathogen_raw) | !is.na(combined_disease_name_raw))

exact_pathogen_matches <- combined_reference %>%
  filter(!is.na(pathogen_clean)) %>%
  distinct(pathogen_clean, .keep_all = TRUE)

exact_disease_matches <- combined_reference %>%
  filter(!is.na(disease_name_clean)) %>%
  distinct(disease_name_clean, .keep_all = TRUE)

exact_alias_pathogen_matches <- appendix_a_aliases %>%
  inner_join(
    exact_pathogen_matches %>%
      transmute(
        pathogen_clean,
        alias_pathogen_combined_pathogen = combined_pathogen_raw,
        alias_pathogen_combined_disease = combined_disease_name_raw
      ),
    by = c("alias_clean" = "pathogen_clean")
  ) %>%
  group_by(efsa_pathogen_id) %>%
  mutate(alias_pathogen_candidate_count = n_distinct(paste(alias_pathogen_combined_pathogen, alias_pathogen_combined_disease, sep = "|||"))) %>%
  filter(alias_pathogen_candidate_count == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    efsa_pathogen_id,
    alias_pathogen_combined_pathogen,
    alias_pathogen_combined_disease
  )

exact_alias_disease_matches <- appendix_a_aliases %>%
  filter(!is.na(alias_disease_candidate_clean)) %>%
  inner_join(
    exact_disease_matches %>%
      transmute(
        disease_name_clean,
        alias_disease_combined_pathogen = combined_pathogen_raw,
        alias_disease_combined_disease = combined_disease_name_raw
      ),
    by = c("alias_disease_candidate_clean" = "disease_name_clean")
  ) %>%
  group_by(efsa_pathogen_id) %>%
  mutate(alias_disease_candidate_count = n_distinct(paste(alias_disease_combined_pathogen, alias_disease_combined_disease, sep = "|||"))) %>%
  filter(alias_disease_candidate_count == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    efsa_pathogen_id,
    alias_disease_combined_pathogen,
    alias_disease_combined_disease
  )

manual_combined_map <- manual_map %>%
  filter(applies_to %in% c("combined_mapping", "both")) %>%
  filter(!is.na(source_name_clean))

crosswalk <- appendix_a_clean %>%
  left_join(
    appendix_g_clean %>%
      distinct(efsa_pathogen_id, efsa_pathogen_standardized) %>%
      count(efsa_pathogen_id, name = "appendix_g_vector_rows"),
    by = "efsa_pathogen_id"
  ) %>%
  mutate(
    appendix_g_vector_rows = replace_na(appendix_g_vector_rows, 0L),
    present_in_appendix_g = appendix_g_vector_rows > 0
  ) %>%
  left_join(
    exact_pathogen_matches %>%
      transmute(
        pathogen_clean,
        exact_pathogen_combined_pathogen = combined_pathogen_raw,
        exact_pathogen_combined_disease = combined_disease_name_raw
      ),
    by = c("efsa_pathogen_clean" = "pathogen_clean")
  ) %>%
  left_join(
    exact_disease_matches %>%
      transmute(
        disease_name_clean,
        exact_disease_combined_pathogen = combined_pathogen_raw,
        exact_disease_combined_disease = combined_disease_name_raw
      ),
    by = c("efsa_disease_candidate_clean" = "disease_name_clean")
  ) %>%
  left_join(
    exact_alias_pathogen_matches,
    by = "efsa_pathogen_id"
  ) %>%
  left_join(
    exact_alias_disease_matches,
    by = "efsa_pathogen_id"
  ) %>%
  left_join(
    manual_combined_map %>%
      transmute(
        source_name_clean,
        manual_combined_pathogen = combined_pathogen,
        manual_combined_disease = combined_disease_name,
        manual_notes = notes
      ),
    by = c("efsa_pathogen_clean" = "source_name_clean")
  ) %>%
  mutate(
    combined_pathogen_raw = case_when(
      !is.na(exact_pathogen_combined_pathogen) ~ exact_pathogen_combined_pathogen,
      !is.na(exact_disease_combined_pathogen) ~ exact_disease_combined_pathogen,
      !is.na(alias_pathogen_combined_pathogen) ~ alias_pathogen_combined_pathogen,
      !is.na(alias_disease_combined_pathogen) ~ alias_disease_combined_pathogen,
      !is.na(manual_combined_pathogen) ~ manual_combined_pathogen,
      TRUE ~ NA_character_
    ),
    combined_disease_name = case_when(
      !is.na(exact_pathogen_combined_disease) ~ exact_pathogen_combined_disease,
      !is.na(exact_disease_combined_disease) ~ exact_disease_combined_disease,
      !is.na(alias_pathogen_combined_disease) ~ alias_pathogen_combined_disease,
      !is.na(alias_disease_combined_disease) ~ alias_disease_combined_disease,
      !is.na(manual_combined_disease) ~ manual_combined_disease,
      TRUE ~ NA_character_
    ),
    match_method = case_when(
      !is.na(exact_pathogen_combined_pathogen) ~ "exact_pathogen",
      !is.na(exact_disease_combined_disease) ~ "exact_disease",
      !is.na(alias_pathogen_combined_pathogen) ~ "exact_alias_pathogen",
      !is.na(alias_disease_combined_disease) ~ "exact_alias_disease",
      !is.na(manual_combined_pathogen) | !is.na(manual_combined_disease) ~ "manual_synonym",
      TRUE ~ NA_character_
    ),
    match_confidence = case_when(
      match_method %in% c("exact_pathogen", "exact_disease", "exact_alias_pathogen", "exact_alias_disease") ~ "exact",
      match_method == "manual_synonym" ~ "manual",
      TRUE ~ NA_character_
    ),
    review_needed = FALSE,
    notes = manual_notes
  )

unresolved_crosswalk <- crosswalk %>%
  filter(is.na(match_method)) %>%
  mutate(fuzzy_candidate = map(efsa_pathogen_clean, build_best_fuzzy_candidate, combined_reference = combined_reference))

fuzzy_updates <- unresolved_crosswalk %>%
  mutate(has_fuzzy_candidate = map_int(fuzzy_candidate, nrow) > 0) %>%
  filter(has_fuzzy_candidate) %>%
  transmute(
    efsa_pathogen_id,
    fuzzy_combined_pathogen = map_chr(fuzzy_candidate, ~ .x$combined_pathogen_raw[[1]]),
    fuzzy_combined_disease = map_chr(fuzzy_candidate, ~ .x$combined_disease_name_raw[[1]]),
    fuzzy_candidate_type = map_chr(fuzzy_candidate, ~ .x$candidate_type[[1]]),
    fuzzy_distance = map_dbl(fuzzy_candidate, ~ .x$distance[[1]]),
    fuzzy_relative_distance = map_dbl(fuzzy_candidate, ~ .x$relative_distance[[1]])
  )

crosswalk <- crosswalk %>%
  left_join(fuzzy_updates, by = "efsa_pathogen_id") %>%
  mutate(
    combined_pathogen_raw = if_else(
      is.na(combined_pathogen_raw) & !is.na(fuzzy_combined_pathogen),
      fuzzy_combined_pathogen,
      combined_pathogen_raw
    ),
    combined_disease_name = if_else(
      is.na(combined_disease_name) & !is.na(fuzzy_combined_disease),
      fuzzy_combined_disease,
      combined_disease_name
    ),
    match_method = case_when(
      !is.na(match_method) ~ match_method,
      !is.na(fuzzy_combined_pathogen) ~ paste0("fuzzy_", fuzzy_candidate_type),
      TRUE ~ "unmatched"
    ),
    match_confidence = case_when(
      !is.na(match_confidence) ~ match_confidence,
      !is.na(fuzzy_combined_pathogen) ~ "fuzzy_review",
      TRUE ~ "unmatched"
    ),
    review_needed = case_when(
      match_confidence %in% c("manual", "fuzzy_review", "unmatched") ~ TRUE,
      TRUE ~ FALSE
    ),
    notes = case_when(
      !is.na(notes) ~ notes,
      !is.na(fuzzy_combined_pathogen) ~ paste0(
        "Fuzzy candidate surfaced for review; distance=",
        fuzzy_distance,
        ", relative_distance=",
        round(fuzzy_relative_distance, 3)
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    combined_pathogen_clean = normalize_name_for_match(combined_pathogen_raw)
  ) %>%
  select(
    efsa_pathogen_id,
    efsa_pathogen_raw,
    efsa_pathogen_clean,
    efsa_abbreviation = abbreviation,
    appendix_a_species_name = species_name,
    appendix_a_isolate_virus_name = isolate_virus_name,
    pathogen_family,
    pathogen_genus,
    oie_notifiable,
    present_in_appendix_g,
    appendix_g_vector_rows,
    combined_pathogen_raw,
    combined_pathogen_clean,
    combined_disease_name,
    match_method,
    match_confidence,
    review_needed,
    notes
  )

pathogen_vector_links_efsa <- appendix_g_clean %>%
  left_join(
    crosswalk %>%
      select(
        efsa_pathogen_id,
        combined_pathogen = combined_pathogen_raw,
        combined_disease_name,
        match_method,
        match_confidence,
        review_needed
      ),
    by = "efsa_pathogen_id"
  ) %>%
  mutate(
    efsa_pathogen = coalesce(efsa_pathogen_standardized, efsa_pathogen_raw),
    source = "EFSA Appendix G",
    review_needed = if_else(
      is.na(review_needed),
      TRUE,
      review_needed | appendix_a_match_status != "matched_to_appendix_a"
    ),
    match_method = case_when(
      appendix_a_match_status != "matched_to_appendix_a" ~ paste0("appendix_a_", appendix_a_match_status),
      TRUE ~ match_method
    )
  ) %>%
  transmute(
    combined_pathogen,
    combined_disease_name,
    efsa_pathogen,
    efsa_pathogen_clean,
    efsa_abbreviation = abbreviation,
    vector_species = vector_species,
    vector_group,
    vector_status,
    field_results,
    laboratory_results,
    vn_database_species,
    vn_priority_species,
    source,
    match_method,
    match_confidence,
    review_needed,
    vector_species_raw,
    vector_group_raw,
    vector_status_raw,
    field_results_raw,
    laboratory_results_raw
  )

review_table <- crosswalk %>%
  filter(review_needed) %>%
  arrange(desc(present_in_appendix_g), match_confidence, efsa_pathogen_raw)

screening <- read_csv(screening_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(everything(), clean_text))

matched_combined_diseases <- crosswalk %>%
  filter(match_confidence %in% c("exact", "manual")) %>%
  pull(combined_disease_name) %>%
  unique() %>%
  discard(is.na)

vector_screening_not_covered <- screening %>%
  filter(screen_status %in% c("clear", "review")) %>%
  filter(!disease_name %in% matched_combined_diseases) %>%
  arrange(screen_status, disease_name)

summary_lines <- c(
  "# EFSA to Combined WHO Mapping Summary",
  "",
  paste0("- Appendix A header row detected: ", appendix_a_raw$header_row),
  paste0("- Appendix A pathogens read: ", nrow(appendix_a_clean)),
  paste0("- Appendix G header row detected: ", appendix_g_raw$header_row),
  paste0("- Appendix G rows read: ", nrow(appendix_g_clean)),
  paste0("- Unique EFSA pathogens represented in Appendix G: ", n_distinct(appendix_g_clean$efsa_pathogen_clean)),
  paste0("- Unique pathogens in combined WHO network: ", n_distinct(combined_reference$combined_pathogen_raw)),
  paste0("- Unique diseases in combined WHO network: ", n_distinct(combined_reference$combined_disease_name_raw)),
  paste0("- EFSA Appendix A pathogens matched to combined WHO network: ", sum(crosswalk$match_confidence %in% c("exact", "manual"))),
  paste0("- Exact matches: ", sum(crosswalk$match_confidence == "exact")),
  paste0("- Manual matches: ", sum(crosswalk$match_confidence == "manual")),
  paste0("- Fuzzy review candidates: ", sum(crosswalk$match_confidence == "fuzzy_review")),
  paste0("- Unmatched EFSA pathogens: ", sum(crosswalk$match_confidence == "unmatched")),
  "",
  "## Unmatched EFSA pathogens",
  if (sum(crosswalk$match_confidence == "unmatched") == 0) {
    "- None"
  } else {
    paste0("- ", crosswalk$efsa_pathogen_raw[crosswalk$match_confidence == "unmatched"])
  },
  "",
  "## Combined WHO vector-screened diseases not covered by EFSA matches",
  if (nrow(vector_screening_not_covered) == 0) {
    "- None"
  } else {
    paste0("- ", vector_screening_not_covered$disease_name, " [", vector_screening_not_covered$screen_status, "]")
  },
  "",
  "## Rows needing manual review",
  if (nrow(review_table) == 0) {
    "- None"
  } else {
    paste0("- ", review_table$efsa_pathogen_raw, " (", review_table$match_confidence, "; ", review_table$match_method, ")")
  }
)

write_csv(appendix_a_clean, file.path(output_dir, "efsa_appendix_a_clean.csv"), na = "")
write_csv(appendix_g_clean, file.path(output_dir, "efsa_appendix_g_clean.csv"), na = "")
write_csv(
  combined_network_clean,
  file.path(output_dir, "combined_who_network_clean_for_efsa_mapping.csv"),
  na = ""
)
write_csv(crosswalk, file.path(output_dir, "efsa_to_combinedwho_crosswalk.csv"), na = "")
write_csv(pathogen_vector_links_efsa, file.path(output_dir, "pathogen_vector_links_efsa.csv"), na = "")
write_csv(review_table, file.path(output_dir, "efsa_crosswalk_review_needed.csv"), na = "")
write_lines(summary_lines, file.path(output_dir, "mapping_summary.md"))

cat("Wrote cleaned Appendix A to", file.path(output_dir, "efsa_appendix_a_clean.csv"), "\n")
cat("Wrote cleaned Appendix G to", file.path(output_dir, "efsa_appendix_g_clean.csv"), "\n")
cat("Wrote cleaned combined WHO network to", file.path(output_dir, "combined_who_network_clean_for_efsa_mapping.csv"), "\n")
cat("Wrote EFSA crosswalk to", file.path(output_dir, "efsa_to_combinedwho_crosswalk.csv"), "\n")
cat("Wrote EFSA pathogen-vector links to", file.path(output_dir, "pathogen_vector_links_efsa.csv"), "\n")
cat("Wrote review table to", file.path(output_dir, "efsa_crosswalk_review_needed.csv"), "\n")
cat("Wrote mapping summary to", file.path(output_dir, "mapping_summary.md"), "\n")
