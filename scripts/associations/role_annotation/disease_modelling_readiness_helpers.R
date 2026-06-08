################################################################################
# disease_modelling_readiness_helpers.R
################################################################################
# Purpose: Pure helpers used by 6_10_Build_Disease_Modelling_Readiness.R.
#
# Notes  : Keep path declarations, input reads, table assembly, validation, and
#          writes in the numbered script. Helpers here should not write files or
#          change upstream evidence interpretation.
################################################################################

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "NULL", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

clean_key <- function(x) {
  x %>%
    clean_text() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("&", " and ") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

is_true <- function(x) {
  as.character(x) %in% c("TRUE", "true", "True", "1", "yes", "Yes", "YES")
}

read_csv_layer <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Required input is missing: ", path, call. = FALSE)
    }
    return(NULL)
  }

  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_text))
}

prefer_existing_path <- function(primary, fallback) {
  if (file.exists(primary)) {
    return(primary)
  }

  fallback
}

empty_disease_summary <- function() {
  tibble(disease_name = character())
}

first_non_empty <- function(...) {
  args <- list(...)
  if (length(args) == 0) {
    return(character())
  }

  out <- rep(NA_character_, length(args[[1]]))
  for (arg in args) {
    values <- clean_text(arg)
    fill <- is.na(out) & !is.na(values)
    out[fill] <- values[fill]
  }
  out
}

is_species_like_name <- function(x) {
  x <- clean_text(x)
  !is.na(x) &
    stringr::str_count(x, "\\S+") >= 2 &
    !stringr::str_detect(
      stringr::str_to_lower(x),
      "^(genus|subgenus|family|order|group|complex)\\b"
    )
}

format_taxon_name <- function(x) {
  x <- clean_text(x)
  ok <- !is.na(x)
  x[ok] <- paste0(
    stringr::str_to_upper(stringr::str_sub(x[ok], 1, 1)),
    stringr::str_sub(x[ok], 2)
  )
  x
}

coalesce_column <- function(data, left_col, right_col) {
  left <- if (left_col %in% names(data)) data[[left_col]] else rep(NA_character_, nrow(data))
  right <- if (right_col %in% names(data)) data[[right_col]] else rep(NA_character_, nrow(data))
  first_non_empty(left, right)
}

add_missing_count_cols <- function(data, cols) {
  for (col in cols) {
    if (!col %in% names(data)) {
      data[[col]] <- 0
    }
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
    data[[col]][is.na(data[[col]])] <- 0
  }
  data
}

add_missing_flag_cols <- function(data, cols) {
  for (col in cols) {
    if (!col %in% names(data)) {
      data[[col]] <- FALSE
    }
    data[[col]] <- is_true(data[[col]])
    data[[col]][is.na(data[[col]])] <- FALSE
  }
  data
}

build_taxid_lookup <- function(data, taxid_col, name_cols, source_name) {
  empty_lookup <- tibble(
    .lookup_key = character(),
    pathogen_taxid_lookup = character(),
    pathogen_taxid_source_lookup = character()
  )

  if (is.null(data) || !taxid_col %in% names(data)) {
    return(empty_lookup)
  }

  name_cols <- intersect(name_cols, names(data))
  if (length(name_cols) == 0) {
    return(empty_lookup)
  }

  purrr::map_dfr(name_cols, function(name_col) {
    data %>%
      transmute(
        .lookup_key = clean_key(.data[[name_col]]),
        pathogen_taxid_lookup = clean_text(.data[[taxid_col]]),
        pathogen_taxid_source_lookup = paste(source_name, name_col, sep = ":")
      )
  }) %>%
    filter(
      !is.na(.lookup_key),
      .lookup_key != "",
      !is.na(pathogen_taxid_lookup),
      pathogen_taxid_lookup != ""
    ) %>%
    group_by(.lookup_key) %>%
    summarise(
      pathogen_taxid_lookup = paste(sort(unique(pathogen_taxid_lookup)), collapse = "; "),
      pathogen_taxid_source_lookup = paste(sort(unique(pathogen_taxid_source_lookup)), collapse = "; "),
      .groups = "drop"
    )
}

lookup_pathogen_taxid <- function(lookup, ..., return_source = FALSE) {
  if (nrow(lookup) == 0) {
    return(NA_character_)
  }

  keys <- clean_key(c(...))
  keys <- keys[!is.na(keys) & keys != ""]
  if (length(keys) == 0) {
    return(NA_character_)
  }

  matched_index <- match(keys, lookup$.lookup_key)
  matched_index <- matched_index[!is.na(matched_index)]
  if (length(matched_index) == 0) {
    return(NA_character_)
  }

  if (return_source) {
    lookup$pathogen_taxid_source_lookup[[matched_index[[1]]]]
  } else {
    lookup$pathogen_taxid_lookup[[matched_index[[1]]]]
  }
}

report_unmatched <- function(label, evidence_names, matched_names) {
  evidence_names <- sort(unique(clean_text(evidence_names)))
  matched_names <- sort(unique(clean_text(matched_names)))
  unmatched <- setdiff(evidence_names[!is.na(evidence_names)], matched_names[!is.na(matched_names)])

  message(label, " unmatched disease labels: ", length(unmatched))
  if (length(unmatched) > 0) {
    message("  ", paste(unmatched, collapse = "; "))
  }
}

collapse_unique <- function(x) {
  values <- sort(unique(clean_text(x)))
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(values, collapse = "; ")
}

collapse_latest <- function(x) {
  values <- sort(unique(clean_text(x)))
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }

  values[[length(values)]]
}

join_evidence_by_names <- function(base, evidence, prefix, candidate_fields) {
  status_col <- paste0(prefix, "_join_status")
  field_col <- paste0(prefix, "_match_field")
  name_col <- paste0(prefix, "_matched_disease_name")

  empty_matches <- tibble(
    analysis_unit_id = base$analysis_unit_id,
    !!status_col := "unmatched_no_evidence",
    !!field_col := NA_character_,
    !!name_col := NA_character_,
    .evidence_key = NA_character_
  )

  if (is.null(evidence) || nrow(evidence) == 0 || !"disease_name" %in% names(evidence)) {
    return(base %>% left_join(empty_matches %>% select(-.evidence_key), by = "analysis_unit_id"))
  }

  evidence_lookup <- evidence %>%
    rename(.evidence_disease_name = disease_name) %>%
    mutate(.evidence_key = clean_key(.evidence_disease_name)) %>%
    filter(!is.na(.evidence_key), .evidence_key != "")

  if (nrow(evidence_lookup) == 0) {
    return(base %>% left_join(empty_matches %>% select(-.evidence_key), by = "analysis_unit_id"))
  }

  key_status <- evidence_lookup %>%
    group_by(.evidence_key) %>%
    summarise(
      .evidence_names = paste(sort(unique(.evidence_disease_name)), collapse = "; "),
      .n_evidence_names = dplyr::n_distinct(.evidence_disease_name),
      .groups = "drop"
    )

  unique_lookup <- evidence_lookup %>%
    inner_join(key_status %>% filter(.n_evidence_names == 1), by = ".evidence_key") %>%
    arrange(.evidence_disease_name) %>%
    distinct(.evidence_key, .keep_all = TRUE) %>%
    select(-.evidence_names, -.n_evidence_names)

  initial_matches <- purrr::map_dfr(seq_len(nrow(base)), function(row_index) {
    for (field in candidate_fields) {
      value <- if (field %in% names(base)) base[[field]][[row_index]] else NA_character_
      key <- clean_key(value)
      if (is.na(key) || key == "") {
        next
      }

      key_info <- key_status %>% filter(.evidence_key == key)
      if (nrow(key_info) == 0) {
        next
      }

      if (key_info$.n_evidence_names[[1]] > 1) {
        return(tibble(
          analysis_unit_id = base$analysis_unit_id[[row_index]],
          !!status_col := "ambiguous_not_joined",
          !!field_col := field,
          !!name_col := key_info$.evidence_names[[1]],
          .evidence_key = NA_character_
        ))
      }

      matched_name <- unique_lookup$.evidence_disease_name[unique_lookup$.evidence_key == key][[1]]
      join_status <- if (field == candidate_fields[[1]]) {
        "matched_primary_name"
      } else {
        "matched_alternate_name"
      }

      return(tibble(
        analysis_unit_id = base$analysis_unit_id[[row_index]],
        !!status_col := join_status,
        !!field_col := field,
        !!name_col := matched_name,
        .evidence_key = key
      ))
    }

    tibble(
      analysis_unit_id = base$analysis_unit_id[[row_index]],
      !!status_col := "unmatched_no_evidence",
      !!field_col := NA_character_,
      !!name_col := NA_character_,
      .evidence_key = NA_character_
    )
  })

  duplicate_selected_keys <- initial_matches %>%
    filter(!is.na(.evidence_key)) %>%
    count(.evidence_key, name = "matched_analysis_units") %>%
    filter(matched_analysis_units > 1) %>%
    pull(.evidence_key)

  matches <- initial_matches %>%
    mutate(
      "{status_col}" := if_else(
        .evidence_key %in% duplicate_selected_keys,
        "ambiguous_not_joined",
        .data[[status_col]]
      ),
      "{field_col}" := if_else(.evidence_key %in% duplicate_selected_keys, NA_character_, .data[[field_col]]),
      "{name_col}" := if_else(.evidence_key %in% duplicate_selected_keys, NA_character_, .data[[name_col]]),
      .evidence_key = if_else(.evidence_key %in% duplicate_selected_keys, NA_character_, .evidence_key)
    )

  evidence_values <- unique_lookup %>%
    select(-.evidence_disease_name)

  base %>%
    left_join(matches, by = "analysis_unit_id") %>%
    left_join(evidence_values, by = ".evidence_key") %>%
    select(-.evidence_key)
}

available_species <- function(data) {
  if (is.null(data) || !"species" %in% names(data)) {
    return(character())
  }

  if ("status" %in% names(data)) {
    data <- data %>% filter(status == "success")
  }

  unique(clean_key(data$species))
}

build_combined_pathogen_taxid_lookup <- function(virion_taxid, clover_taxid) {
  bind_rows(
    build_taxid_lookup(
      virion_taxid,
      taxid_col = "VirusTaxID",
      name_cols = c("Pathogens", "previous_name", "msl39_viral_name", "Virion_VirusName", "matched_virus_name"),
      source_name = "virion_taxid"
    ),
    build_taxid_lookup(
      clover_taxid,
      taxid_col = "PathogenTaxID",
      name_cols = c("Pathogens", "previous_name", "msl39_viral_name", "Clover_PathogenName"),
      source_name = "clover_taxid"
    )
  ) %>%
    group_by(.lookup_key) %>%
    summarise(
      pathogen_taxid_lookup = paste(sort(unique(pathogen_taxid_lookup)), collapse = "; "),
      pathogen_taxid_source_lookup = paste(sort(unique(pathogen_taxid_source_lookup)), collapse = "; "),
      .groups = "drop"
    )
}

readiness_output_specs <- function() {
  core_cols <- c(
    "analysis_unit_id",
    "readiness_disease_name",
    "source_disease_name",
    "disease_master_name",
    "analysis_unit",
    "analysis_unit_label",
    "row_type",
    "family",
    "priority_prototype_status",
    "master_tier",
    "master_guild",
    "master_livestock_amplified",
    "master_key_host_vector",
    "include_as_analysis_unit",
    "modelling_scope_status",
    "modelling_scope_reason"
  )

  readiness_cols <- c(
    "vectored_status",
    "generalist_status",
    "transmission_complexity",
    "guild",
    "host_sdm_needed",
    "vector_sdm_needed",
    "host_range_rule",
    "vector_range_rule",
    "range_limiting_layer",
    "transmission_rule_review_status",
    "transmission_rule_notes"
  )

  evidence_cols <- c(
    "host_candidate_rows",
    "roster_host_rows",
    "roster_vector_rows",
    "host_role_evidence_rows",
    "host_role_assignment_rows",
    "direct_disease_vector_rows",
    "direct_vector_confirmed_rows",
    "direct_vector_probable_rows",
    "direct_vector_candidate_rows",
    "direct_vector_competence_joined_rows",
    "direct_vector_transmission_yes_rows",
    "direct_vector_natural_infection_yes_rows",
    "roster_vectors_with_host_vector_evidence",
    "roster_vectors_with_human_biting",
    "genbank_country_rows",
    "genbank_distinct_countries_or_territories",
    "genbank_records_with_country",
    "who_don_focal_country_rows",
    "who_don_distinct_countries",
    "who_don_distinct_records",
    "has_substantive_host_role_evidence",
    "has_substantive_vector_role_evidence",
    "role_assignment_status"
  )

  sdm_cols <- c(
    "host_sdm_species_available",
    "vector_sdm_species_available",
    "host_sdm_projection_species_available",
    "vector_sdm_projection_species_available",
    "host_sdm_comparison_species_available",
    "vector_sdm_comparison_species_available"
  )

  decision_cols <- c(
    "has_direct_vector_evidence",
    "direct_vector_evidence_status",
    "country_evidence_status",
    "sdm_availability_status",
    "readiness_blocker",
    "recommended_next_action",
    "evidence_join_status"
  )

  match_cols <- c(
    "evidence_qa_join_status",
    "evidence_qa_match_field",
    "evidence_qa_matched_disease_name",
    "direct_vector_join_status",
    "direct_vector_match_field",
    "direct_vector_matched_disease_name",
    "genbank_join_status",
    "genbank_match_field",
    "genbank_matched_disease_name",
    "who_don_join_status",
    "who_don_match_field",
    "who_don_matched_disease_name",
    "sdm_join_status",
    "sdm_match_field",
    "sdm_matched_disease_name"
  )

  slim_cols <- c(
    "analysis_unit_id",
    "readiness_disease_name",
    "pathogen_species_name",
    "pathogen_taxid",
    "analysis_unit_label",
    "family",
    "in_master_who",
    "priority_prototype_status",
    "master_tier",
    "modelling_scope_status",
    "transmission_rule_review_status",
    "recommended_next_action",
    "readiness_blocker",
    "vectored_status",
    "generalist_status",
    "transmission_complexity",
    "guild",
    "range_limiting_layer",
    "host_sdm_needed",
    "vector_sdm_needed",
    "direct_vector_evidence_status",
    "country_evidence_status",
    "sdm_availability_status",
    "role_assignment_status",
    "evidence_join_status",
    "has_direct_vector_evidence",
    "host_sdm_species_available",
    "vector_sdm_species_available",
    "genbank_distinct_countries_or_territories",
    "who_don_distinct_countries"
  )

  list(
    core = core_cols,
    readiness = readiness_cols,
    evidence = evidence_cols,
    sdm = sdm_cols,
    decision = decision_cols,
    match = match_cols,
    slim = slim_cols
  )
}

validate_readiness_outputs <- function(
  readiness_full,
  readiness_slim,
  readiness_pilot,
  master,
  slim_cols,
  pilot_cols,
  held_analysis_unit_ids
) {
  if (nrow(readiness_full) != nrow(master)) {
    stop("Full readiness output row count does not match master table row count.", call. = FALSE)
  }

  if (any(readiness_slim$analysis_unit_id %in% held_analysis_unit_ids)) {
    stop("Slim readiness output still contains held analysis units.", call. = FALSE)
  }

  if (any(readiness_pilot$analysis_unit_id %in% held_analysis_unit_ids)) {
    stop("Pilot readiness output still contains held analysis units.", call. = FALSE)
  }

  if (
    any(is.na(readiness_full$analysis_unit_id)) || any(readiness_full$analysis_unit_id == "") ||
      any(is.na(readiness_slim$analysis_unit_id)) || any(readiness_slim$analysis_unit_id == "") ||
      any(is.na(readiness_pilot$analysis_unit_id)) || any(readiness_pilot$analysis_unit_id == "")
  ) {
    stop("Readiness output has missing `analysis_unit_id` values.", call. = FALSE)
  }

  if (
    any(duplicated(readiness_full$analysis_unit_id)) ||
      any(duplicated(readiness_slim$analysis_unit_id)) ||
      any(duplicated(readiness_pilot$analysis_unit_id))
  ) {
    stop("Readiness output has duplicate `analysis_unit_id` values.", call. = FALSE)
  }

  if (ncol(readiness_slim) != length(slim_cols)) {
    stop("Slim readiness output does not have the expected column count.", call. = FALSE)
  }

  if (ncol(readiness_pilot) != length(pilot_cols)) {
    stop("Pilot readiness output does not have the expected column count.", call. = FALSE)
  }

  has_cols <- names(readiness_full)[stringr::str_starts(names(readiness_full), "has_")]
  has_na_counts <- purrr::map_int(readiness_full[has_cols], ~ sum(is.na(.x)))
  if (any(has_na_counts > 0)) {
    stop("Derived `has_*` flags contain NA values.", call. = FALSE)
  }

  invisible(TRUE)
}

pilot_package_source_descriptions <- function() {
  c(
    disease_modelling_pilot = "Generated from disease_modelling_readiness_full.csv rows in the non-held WHO pilot subset.",
    pilot_hosts = "Filtered host rows from species_host_vector_roster.csv using readiness-script disease matches.",
    pilot_vectors = "Filtered vector rows from species_host_vector_roster.csv using readiness-script disease matches.",
    pilot_countries = "Filtered and summarised GenBank disease-country and WHO DON country evidence using readiness-script disease matches.",
    pilot_sdm_species = "Filtered host/vector roster species joined to the accessible SDM species inventory, with projection and comparison status where available.",
    pilot_evidence_summary = "Disease-level evidence, role, country, SDM, and join-count fields from the full readiness audit table."
  )
}

build_pilot_package_manifest <- function(data_tables, generated_at_utc, sources) {
  purrr::imap_dfr(data_tables, function(table, table_name) {
    tibble(
      generated_at_utc = generated_at_utc,
      table_name = table_name,
      file_name = paste0(table_name, ".csv"),
      rows = nrow(table),
      columns = ncol(table),
      source_description = sources[[table_name]]
    )
  })
}

validate_pilot_package_tables <- function(data_tables, pilot_ids) {
  for (table_name in names(data_tables)) {
    table <- data_tables[[table_name]]
    if ("analysis_unit_id" %in% names(table) && any(!table$analysis_unit_id %in% pilot_ids)) {
      stop("Pilot package table contains analysis units outside the pilot subset: ", table_name, call. = FALSE)
    }
  }

  invisible(TRUE)
}

build_pilot_match_lookup <- function(pilot_context, layers) {
  match_cols <- paste0(layers, "_matched_disease_name")
  match_cols <- intersect(match_cols, names(pilot_context))

  if (length(match_cols) == 0) {
    return(tibble(
      analysis_unit_id = character(),
      match_layer = character(),
      package_matched_disease_name = character(),
      matched_key = character()
    ))
  }

  layer_priority <- tibble(
    match_layer = layers,
    match_priority = seq_along(layers)
  )

  pilot_context %>%
    select(analysis_unit_id, all_of(match_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(match_cols),
      names_to = "match_layer",
      values_to = "package_matched_disease_name"
    ) %>%
    mutate(
      match_layer = stringr::str_remove(match_layer, "_matched_disease_name$"),
      package_matched_disease_name = clean_text(package_matched_disease_name),
      matched_key = clean_key(package_matched_disease_name)
    ) %>%
    filter(!is.na(matched_key), matched_key != "") %>%
    inner_join(layer_priority, by = "match_layer") %>%
    arrange(match_priority, match_layer, package_matched_disease_name) %>%
    distinct(analysis_unit_id, matched_key, .keep_all = TRUE) %>%
    select(analysis_unit_id, match_layer, package_matched_disease_name, matched_key)
}

join_pilot_layer <- function(data, disease_col, layers, pilot_context, context_cols) {
  empty_layer <- tibble()

  if (is.null(data) || nrow(data) == 0 || !disease_col %in% names(data)) {
    return(empty_layer)
  }

  match_lookup <- build_pilot_match_lookup(pilot_context, layers)
  if (nrow(match_lookup) == 0) {
    return(empty_layer)
  }

  data %>%
    mutate(.matched_key = clean_key(.data[[disease_col]])) %>%
    inner_join(match_lookup, by = c(".matched_key" = "matched_key")) %>%
    left_join(pilot_context %>% select(any_of(context_cols)), by = "analysis_unit_id") %>%
    mutate(package_match_layer = match_layer) %>%
    select(
      any_of(context_cols),
      package_match_layer,
      package_matched_disease_name,
      everything(),
      -any_of(c(".matched_key", "match_layer"))
    )
}

pilot_package_readme_lines <- function() {
  c(
    "# Disease Modelling Pilot Package",
    "",
    "Generated by `scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R`.",
    "",
    "This package is a modelling handoff bundle keyed by `analysis_unit_id`. The",
    "disease-level pilot table is the spine; the host, vector, country, SDM, and",
    "evidence-summary tables are long companion tables for one-to-many layers.",
    "Companion tables keep only a small repeated disease context block:",
    "`analysis_unit_id` and `readiness_disease_name`. Disease-level readiness",
    "status fields live in `disease_modelling_pilot.csv` and",
    "`pilot_evidence_summary.csv`.",
    "",
    "## Tables",
    "",
    "- `manifest.csv`: row counts, column counts, source descriptions, and generation time.",
    "- `disease_modelling_pilot.csv`: one row per non-held WHO-origin pilot analysis unit.",
    "- `pilot_hosts.csv`: host species rows for pilot diseases from `species_host_vector_roster.csv`.",
    "- `pilot_vectors.csv`: vector taxon rows for pilot diseases from `species_host_vector_roster.csv`.",
    "- `pilot_countries.csv`: GenBank and WHO DON disease-country evidence summaries.",
    "- `pilot_sdm_species.csv`: host/vector species with accessible SDM model status, plus projection and comparison status where available.",
    "- `pilot_evidence_summary.csv`: disease-level evidence, country, role, and SDM counts.",
    "",
    "The package is also written as `../disease_modelling_pilot_package.rds` and,",
    "when `writexl` is available, `../disease_modelling_pilot_package.xlsx`."
  )
}
