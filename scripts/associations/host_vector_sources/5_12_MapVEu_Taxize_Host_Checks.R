# ------------------------------------------------------------------------------
# 5_12_MapVEu_Taxize_Host_Checks.R
# ------------------------------------------------------------------------------
# Purpose: Use taxize to review specific unmatched MapVEu host names against the
#          WHO host list, mirroring the VectorMap host-taxize stage but scoped
#          only to the small unresolved MapVEu host set.
#
# Inputs : pathogen_association_data/staged/mapveu/outputs/
#          mapveu_vector_host_links_raw.csv
#          pathogen_association_data/manual/mapveu/
#          mapveu_host_manual_crosswalk.csv
#          WHO network helper path for combined_who_network_canonical_zoonotic.csv
# Outputs: pathogen_association_data/staged/mapveu/outputs/
#          mapveu_host_taxize_candidates.csv
#          pathogen_association_data/staged/mapveu/outputs/
#          mapveu_host_taxize_review.csv
#          pathogen_association_data/staged/mapveu/outputs/
#          mapveu_host_taxize_who_hits.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, taxize, tibble)

source(here("scripts", "associations", "working_inputs.R"))

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

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

count_host_words <- function(x) {
  key <- normalize_host_key(x)
  if_else(is.na(key), 0L, str_count(key, " ") + 1L)
}

is_non_actionable_host <- function(x) {
  key <- normalize_host_key(x)
  word_count <- count_host_words(x)

  is.na(key) |
    key == "" |
    word_count < 2L |
    str_detect(
      key,
      "\\b(sp|sp\\.|spp|spp\\.|species complex|complex|group|subgroup|unknown|unidentified|organism|organisms|multi species|multi species collection|multi|mixed|other|vertebrates)\\b"
    ) |
    str_detect(key, "\\band\\b|/")
}

pull_first <- function(df, candidates) {
  cols <- intersect(candidates, names(df))

  if (length(cols) == 0 || nrow(df) == 0) {
    return(NA_character_)
  }

  value <- df[[cols[[1]]]][[1]]
  clean_text(value)
}

resolve_with_taxize <- function(query_name) {
  gna_result <- tryCatch(
    taxize::gna_verifier(
      names = query_name,
      all_matches = FALSE
    ),
    error = function(e) e
  )

  if (inherits(gna_result, "error")) {
    return(tibble(
      taxize_status = "error",
      taxize_note = conditionMessage(gna_result),
      matched_name = NA_character_,
      canonical_form = NA_character_,
      current_name_string = NA_character_,
      score = NA_real_,
      data_source_title = NA_character_,
      classification_path = NA_character_
    ))
  }

  if (nrow(gna_result) == 0) {
    return(tibble(
      taxize_status = "no_match",
      taxize_note = "No taxize match returned.",
      matched_name = NA_character_,
      canonical_form = NA_character_,
      current_name_string = NA_character_,
      score = NA_real_,
      data_source_title = NA_character_,
      classification_path = NA_character_
    ))
  }

  tibble(
    taxize_status = "matched",
    taxize_note = NA_character_,
    matched_name = pull_first(gna_result, c("matchedName", "submittedName")),
    canonical_form = pull_first(gna_result, c("matchedCanonicalSimple", "matchedCanonicalFull", "matchedName")),
    current_name_string = pull_first(gna_result, c("currentCanonicalSimple", "currentCanonicalFull", "currentName", "matchedCanonicalSimple")),
    score = suppressWarnings(as.numeric(pull_first(gna_result, c("sortScore", "acceptedNameScore", "editDistance")))),
    data_source_title = pull_first(gna_result, c("dataSourceTitleShort")),
    classification_path = paste(
      na.omit(c(
        pull_first(gna_result, c("taxonomicStatus")),
        ifelse(isTRUE(pull_first(gna_result, c("isSynonym")) == "TRUE"), "synonym", NA_character_),
        pull_first(gna_result, c("matchType"))
      )),
      collapse = "; "
    )
  )
}

outputs_dir <- mapveu_outputs_dir
manual_dir <- mapveu_manual_dir

raw_links_path <- file.path(outputs_dir, "mapveu_vector_host_links_raw.csv")
host_crosswalk_path <- file.path(manual_dir, "mapveu_host_manual_crosswalk.csv")
combined_network_path <- who_working_network_path()

candidate_path <- file.path(outputs_dir, "mapveu_host_taxize_candidates.csv")
taxize_review_path <- file.path(outputs_dir, "mapveu_host_taxize_review.csv")
taxize_who_hits_path <- file.path(outputs_dir, "mapveu_host_taxize_who_hits.csv")

max_queries_env <- Sys.getenv("TAXIZE_MAX_QUERIES", unset = "")
max_queries <- if (max_queries_env == "") {
  Inf
} else {
  suppressWarnings(as.numeric(max_queries_env))
}

sleep_seconds <- suppressWarnings(as.numeric(Sys.getenv("TAXIZE_SLEEP_SECONDS", unset = "0")))
if (is.na(sleep_seconds)) {
  sleep_seconds <- 0
}

who_hosts <- read_csv(
  combined_network_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    matched_who_host = clean_text(Host),
    matched_who_host_key = normalize_host_key(Host)
  ) %>%
  distinct()

host_crosswalk <- if (file.exists(host_crosswalk_path)) {
  read_csv(
    host_crosswalk_path,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    transmute(
      raw_host_label = clean_text(raw_host_label),
      raw_host_key = normalize_host_key(raw_host_label),
      source_field = clean_text(source_field),
      matched_who_host = clean_text(matched_who_host),
      include_in_final = tolower(clean_text(include_in_final)) %in% c("true", "t", "1", "yes", "y")
    ) %>%
    filter(
      source_field == "host_organism_raw",
      include_in_final
    )
} else {
  tibble(
    raw_host_label = character(),
    raw_host_key = character(),
    source_field = character(),
    matched_who_host = character(),
    include_in_final = logical()
  )
}

candidates <- read_csv(
  raw_links_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(
    host_organism_raw = clean_text(host_organism_raw),
    host_key = normalize_host_key(host_organism_raw),
    host_non_actionable = is_non_actionable_host(host_organism_raw)
  ) %>%
  filter(!is.na(host_organism_raw)) %>%
  left_join(
    who_hosts,
    by = c("host_key" = "matched_who_host_key")
  ) %>%
  left_join(
    host_crosswalk %>%
      select(raw_host_key, crosswalk_who_host = matched_who_host),
    by = c("host_key" = "raw_host_key")
  ) %>%
  mutate(
    final_who_host = coalesce(matched_who_host, crosswalk_who_host)
  ) %>%
  filter(
    is.na(final_who_host),
    !host_non_actionable
  ) %>%
  group_by(host_organism_raw, host_key) %>%
  summarise(
    total_row_count = n(),
    source_fields = "host_organism_raw",
    example_reason_flags = "specific_unmatched",
    .groups = "drop"
  ) %>%
  rename(candidate_query = host_organism_raw, candidate_query_key = host_key) %>%
  arrange(desc(total_row_count), candidate_query)

if (is.finite(max_queries)) {
  candidates <- candidates %>%
    slice_head(n = max_queries)
}

write_csv(candidates, candidate_path, na = "")

cat("Taxize candidate queries:", nrow(candidates), "\n")

resolved_results <- vector("list", length = nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  query_name <- candidates$candidate_query[[i]]

  if (sleep_seconds > 0 && i > 1) {
    Sys.sleep(sleep_seconds)
  }

  cat("Resolving", i, "of", nrow(candidates), ":", query_name, "\n")
  resolved_results[[i]] <- resolve_with_taxize(query_name)
}

resolved_tbl <- if (length(resolved_results) == 0) {
  tibble(
    taxize_status = character(),
    taxize_note = character(),
    matched_name = character(),
    canonical_form = character(),
    current_name_string = character(),
    score = numeric(),
    data_source_title = character(),
    classification_path = character()
  )
} else {
  bind_rows(resolved_results)
}

taxize_review <- resolved_tbl %>%
  bind_cols(candidates) %>%
  mutate(
    resolved_name = coalesce(current_name_string, canonical_form, matched_name),
    resolved_name_key = normalize_host_key(resolved_name),
    resolved_matches_who = resolved_name_key %in% who_hosts$matched_who_host_key,
    resolved_who_host = who_hosts$matched_who_host[match(resolved_name_key, who_hosts$matched_who_host_key)]
  ) %>%
  relocate(
    candidate_query,
    total_row_count,
    source_fields,
    example_reason_flags,
    taxize_status,
    taxize_note,
    matched_name,
    canonical_form,
    current_name_string,
    resolved_name,
    score,
    data_source_title,
    classification_path,
    resolved_matches_who,
    resolved_who_host
  )

taxize_who_hits <- taxize_review %>%
  filter(resolved_matches_who) %>%
  arrange(desc(total_row_count), candidate_query)

write_csv(taxize_review, taxize_review_path, na = "")
write_csv(taxize_who_hits, taxize_who_hits_path, na = "")

cat("Wrote taxize candidates to", candidate_path, "\n")
cat("Wrote taxize review to", taxize_review_path, "\n")
cat("Wrote taxize WHO hits to", taxize_who_hits_path, "\n")
cat("Resolved WHO-host hits:", nrow(taxize_who_hits), "\n")
