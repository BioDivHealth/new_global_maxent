# ------------------------------------------------------------------------------
# Use taxize to review unresolved VectorMap host names against the WHO host list
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, taxize, tibble)

source(here("scripts", "associations", "working_inputs.R"))

# Clean text fields while preserving the original review tables on disk.
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x[x == ""] <- NA_character_
  x
}

# Normalize host strings for exact post-resolution matching back to WHO hosts.
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

# Pull a single value from a result data frame without assuming every field exists.
pull_first <- function(df, candidates) {
  cols <- intersect(candidates, names(df))

  if (length(cols) == 0 || nrow(df) == 0) {
    return(NA_character_)
  }

  value <- df[[cols[[1]]]][[1]]
  clean_text(value)
}

# Resolve one query through taxize and keep the best available match only.
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

outputs_dir <- vectormap_outputs_dir

candidate_path <- file.path(outputs_dir, "vectormap_host_package_candidates.csv")
combined_network_path <- who_working_network_path()

taxize_review_path <- file.path(outputs_dir, "vectormap_host_taxize_review.csv")
taxize_who_hits_path <- file.path(outputs_dir, "vectormap_host_taxize_who_hits.csv")

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

# Load WHO hosts and keep one authoritative label per name for exact back-matching.
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

# Restrict taxize work to unresolved labels that already look taxonomic.
candidates <- read_csv(
  candidate_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(
    raw_host_label = clean_text(raw_host_label),
    source_field = clean_text(source_field),
    host_binomial = clean_text(host_binomial),
    scientific_name_clean = clean_text(scientific_name_clean),
    candidate_query = coalesce(host_binomial, scientific_name_clean),
    candidate_query_key = normalize_host_key(candidate_query)
  ) %>%
  filter(
    package_check_target,
    !is.na(candidate_query),
    candidate_query != ""
  ) %>%
  group_by(candidate_query, candidate_query_key) %>%
  summarise(
    total_row_count = sum(row_count, na.rm = TRUE),
    source_fields = paste(sort(unique(source_field)), collapse = "; "),
    example_raw_labels = paste(head(sort(unique(raw_host_label)), 5), collapse = " | "),
    example_reason_flags = paste(sort(unique(reason_flag)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_row_count), candidate_query)

if (is.finite(max_queries)) {
  candidates <- candidates %>%
    slice_head(n = max_queries)
}

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

taxize_review <- bind_rows(resolved_results) %>%
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
    example_raw_labels,
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

cat("Wrote taxize review to", taxize_review_path, "\n")
cat("Wrote taxize WHO hits to", taxize_who_hits_path, "\n")
cat("Resolved WHO-host hits:", nrow(taxize_who_hits), "\n")
