#!/usr/bin/env Rscript
################################################################################
# 6_11_Derive_Modelling_Evidence_Tiers.R
################################################################################
# Purpose: Build a read-only prototype evidence-tier surface for modelling.
#
# Output : pathogen_association_data/readiness/evidence_tiers/
#
# Notes  : This script derives a compact modelling-readiness handoff surface
#          from the current readiness package. It does not edit canonical
#          evidence, source-check ledgers, role assignments, or SDM manifests.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tibble, tidyr)

source(here::here("scripts", "associations", "working_inputs.R"))
source(here::here(
  "scripts",
  "associations",
  "role_annotation",
  "disease_modelling_readiness_helpers.R"
))

# ------------------------------------------------------------------------------|
#      Paths -------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_package_dir <- file.path(readiness_dir, "disease_modelling_pilot_package")
output_dir <- file.path(readiness_dir, "evidence_tiers")
staged_output_dir <- file.path(
  staged_data_dir,
  "role_annotation",
  "modelling_evidence_tiers"
)

chikungunya_delivery_dir <- Sys.getenv(
  "CHIKUNGUNYA_DELIVERY_DIR",
  unset = "/Volumes/LaCie/new_global_maxent/sdms/delivery/chikungunya_vector_sdm_delivery_20260609"
)
delivery_package_dir <- file.path(chikungunya_delivery_dir, "readiness", "disease_modelling_pilot_package")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staged_output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  repo_hosts = file.path(repo_package_dir, "pilot_hosts.csv"),
  repo_vectors = file.path(repo_package_dir, "pilot_vectors.csv"),
  repo_sdm_species = file.path(repo_package_dir, "pilot_sdm_species.csv"),
  delivery_hosts = file.path(delivery_package_dir, "pilot_hosts.csv"),
  delivery_vectors = file.path(delivery_package_dir, "pilot_vectors.csv"),
  delivery_sdm_species = file.path(delivery_package_dir, "pilot_sdm_species.csv"),
  delivery_model_qc = file.path(chikungunya_delivery_dir, "model_qc_summary.csv")
)

# ------------------------------------------------------------------------------|
#      Small Helpers -----------------------------------------------------------|
# ------------------------------------------------------------------------------|
read_optional_csv <- function(path) {
  read_csv_layer(path, required = FALSE)
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(x))
}

missing_as_false <- function(x) {
  out <- is_true(x)
  out[is.na(out)] <- FALSE
  out
}

combine_reasons <- function(...) {
  values <- c(...)
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(unique(values), collapse = "; ")
}

max_tier <- function(strict, strong, supported, broad) {
  case_when(
    strict ~ "strict",
    strong ~ "strong",
    supported ~ "supported",
    broad ~ "broad",
    TRUE ~ "excluded"
  )
}

is_avian_influenza_generalist <- function(x) {
  x %in% c(
    "Influenza (H5N1 avian influenza)",
    "Influenza (H7N9 avian influenza)"
  )
}

is_west_nile <- function(x) {
  x == "West Nile fever"
}

confidence_weight <- function(confidence, bucket) {
  case_when(
    bucket == "unknown_or_unreviewed" ~ 0,
    bucket == "host_presence_only" & confidence == "low" ~ 0.1,
    confidence == "high" ~ 1,
    confidence == "medium" ~ 0.75,
    confidence == "medium_low" ~ 0.5,
    confidence == "low" ~ 0.25,
    TRUE ~ 0.1
  )
}

delivery_available <- function() {
  all(file.exists(c(paths$delivery_hosts, paths$delivery_vectors, paths$delivery_sdm_species)))
}

build_sdm_lookup <- function(sdm_species) {
  if (is.null(sdm_species)) {
    return(tibble(
      analysis_unit_id = character(),
      species_name = character(),
      species_role = character(),
      sdm_needed_for_disease = character(),
      sdm_available = logical(),
      sdm_species = character()
    ))
  }

  sdm_species %>%
    mutate(
      species_role = clean_text(species_role),
      sdm_available = missing_as_false(sdm_available)
    ) %>%
    select(any_of(c(
      "analysis_unit_id",
      "species_name",
      "species_role",
      "sdm_needed_for_disease",
      "sdm_available",
      "sdm_species"
    ))) %>%
    distinct(analysis_unit_id, species_name, species_role, .keep_all = TRUE)
}

add_sdm_fields <- function(data, sdm_lookup) {
  joined <- data %>%
    left_join(
      sdm_lookup,
      by = c("analysis_unit_id", "species_name", "species_role"),
      suffix = c("", "_sdm")
    )

  sdm_cols <- c(
    "sdm_needed_for_disease",
    "sdm_needed_for_disease_sdm",
    "sdm_available",
    "sdm_available_sdm",
    "sdm_species",
    "sdm_species_sdm"
  )
  for (col in sdm_cols) {
    if (!col %in% names(joined)) {
      joined[[col]] <- NA
    }
  }

  joined %>%
    mutate(
      sdm_needed_for_disease = first_non_empty(sdm_needed_for_disease, sdm_needed_for_disease_sdm),
      sdm_available = missing_as_false(first_non_empty(sdm_available, sdm_available_sdm)),
      sdm_species = first_non_empty(sdm_species, sdm_species_sdm)
    ) %>%
    select(-any_of(c("sdm_needed_for_disease_sdm", "sdm_available_sdm", "sdm_species_sdm")))
}

add_model_quality <- function(data, model_qc = NULL) {
  if (!is.null(model_qc)) {
    qc <- model_qc %>%
      select(any_of(c(
        "species_name",
        "species_role",
        "model_quality",
        "retained_models",
        "min_boyce",
        "max_boyce",
        "min_test_auc",
        "max_test_auc",
        "min_max_tss",
        "max_max_tss",
        "prediction_tif_path"
      ))) %>%
      distinct(species_name, species_role, .keep_all = TRUE)

    data <- data %>%
      left_join(qc, by = c("species_name", "species_role"))
  }

  model_cols <- c(
    "model_quality",
    "retained_models",
    "min_boyce",
    "max_boyce",
    "min_test_auc",
    "max_test_auc",
    "min_max_tss",
    "max_max_tss",
    "prediction_tif_path"
  )
  for (col in model_cols) {
    if (!col %in% names(data)) {
      data[[col]] <- NA
    }
  }

  data %>%
    mutate(
      retained_models = as_numeric_clean(retained_models),
      min_boyce = as_numeric_clean(min_boyce),
      model_quality_raw = clean_text(model_quality),
      model_quality_tier = case_when(
        !sdm_available ~ "no_sdm",
        str_detect(model_quality_raw, regex("^diagnostic", ignore_case = TRUE)) ~ "diagnostic",
        str_detect(model_quality_raw, regex("production", ignore_case = TRUE)) &
          !is.na(min_boyce) & min_boyce >= 0.5 &
          !is.na(retained_models) & retained_models >= 5 ~ "usable_sdm",
        str_detect(model_quality_raw, regex("existing_host_model", ignore_case = TRUE)) ~ "existing_host_model",
        sdm_available ~ "available_unscored",
        TRUE ~ "no_sdm"
      ),
      model_quality_weight = case_when(
        model_quality_tier == "usable_sdm" ~ 1,
        model_quality_tier == "existing_host_model" ~ 0.75,
        model_quality_tier == "available_unscored" ~ 0.5,
        model_quality_tier == "diagnostic" ~ 0.25,
        TRUE ~ 0
      )
    )
}

add_host_modelling_proxy <- function(data) {
  data %>%
    mutate(
      .species_key = clean_key(species_name),
      .tax_id_key = clean_text(tax_id),
      .avian_influenza_generalist = is_avian_influenza_generalist(readiness_disease_name),
      .west_nile = is_west_nile(readiness_disease_name),
      .source_backed_specific_role = host_role_source_backed & host_role_specific,
      .reviewed_specific_role = host_role_specific &
        host_role_assignment_status %in% c("draft_source_backed", "draft_needs_review"),
      .wild_aquatic_bird_group =
        .avian_influenza_generalist &
          host_class == "aves" &
          host_order %in% c("anseriformes", "charadriiformes"),
      .avian_host_group =
        .avian_influenza_generalist &
          host_class == "aves",
      .human_spillover_host =
        .avian_influenza_generalist &
          (.species_key == "homo sapiens" | .tax_id_key == "9606"),
      .mammalian_spillover_host =
        .avian_influenza_generalist &
          host_class == "mammalia",
      .west_nile_avian_group =
        .west_nile &
          host_class == "aves",
      .west_nile_corvidae_group =
        .west_nile_avian_group &
          host_family == "corvidae",
      .west_nile_passeriform_group =
        .west_nile_avian_group &
          host_order == "passeriformes",
      .west_nile_charadriiform_group =
        .west_nile_avian_group &
          host_order == "charadriiformes",
      modelling_role_proxy = case_when(
        .reviewed_specific_role ~ host_role_assignment,
        .west_nile_corvidae_group ~ "elevated_avian_amplifying_competence_proxy",
        .west_nile_passeriform_group ~ "elevated_avian_amplifying_competence_proxy",
        .west_nile_charadriiform_group ~ "elevated_avian_competence_proxy",
        .west_nile_avian_group ~ "avian_reservoir_amplifying_group_proxy",
        .wild_aquatic_bird_group ~ "wild_aquatic_bird_reservoir_group_proxy",
        .avian_host_group & host_order == "galliformes" ~ "galliform_avian_host_group_proxy",
        .avian_host_group ~ "avian_host_presence_group_proxy",
        .human_spillover_host ~ "incidental_spillover_host_proxy",
        .mammalian_spillover_host ~ "mammalian_susceptible_or_spillover_host_proxy",
        .avian_influenza_generalist ~ "host_presence_group_proxy",
        TRUE ~ host_role_assignment
      ),
      modelling_role_proxy_basis = case_when(
        .source_backed_specific_role ~ "source_backed_role_assignment",
        .reviewed_specific_role ~ "reviewed_role_assignment_needs_review",
        .west_nile_corvidae_group ~ "west_nile_weighted_rule_corvidae",
        .west_nile_passeriform_group ~ "west_nile_weighted_rule_passeriformes",
        .west_nile_charadriiform_group ~ "west_nile_weighted_rule_charadriiformes",
        .west_nile_avian_group ~ "west_nile_group_rule_aves",
        .wild_aquatic_bird_group ~ "avian_influenza_group_rule_wild_aquatic_birds",
        .avian_host_group & host_order == "galliformes" ~ "avian_influenza_group_rule_galliform_birds",
        .avian_host_group ~ "avian_influenza_group_rule_avian_host_presence",
        .human_spillover_host ~ "avian_influenza_group_rule_human_spillover",
        .mammalian_spillover_host ~ "avian_influenza_group_rule_mammalian_spillover_or_susceptible_host",
        .avian_influenza_generalist ~ "avian_influenza_group_rule_taxonomy_missing_or_other_host",
        TRUE ~ "candidate_role_assignment"
      ),
      modelling_role_proxy_confidence = case_when(
        .reviewed_specific_role ~ host_role_confidence,
        .west_nile_corvidae_group ~ "medium",
        .west_nile_passeriform_group ~ "medium_low",
        .west_nile_charadriiform_group ~ "medium_low",
        .west_nile_avian_group ~ "low",
        .wild_aquatic_bird_group ~ "medium",
        .human_spillover_host ~ "medium",
        .avian_influenza_generalist ~ "low",
        TRUE ~ host_role_confidence
      ),
      modelling_role_proxy_rule_id = case_when(
        .source_backed_specific_role ~ "source_backed_role_v0_1",
        .reviewed_specific_role ~ "reviewed_role_needs_review_v0_1",
        .west_nile_corvidae_group ~ "wnv_corvidae_weighted_proxy_v0_1",
        .west_nile_passeriform_group ~ "wnv_passeriformes_weighted_proxy_v0_1",
        .west_nile_charadriiform_group ~ "wnv_charadriiformes_weighted_proxy_v0_1",
        .west_nile_avian_group ~ "wnv_aves_group_proxy_v0_1",
        .avian_influenza_generalist ~ "avian_influenza_group_proxy_v0_1",
        TRUE ~ "candidate_role_v0_1"
      ),
      modelling_role_proxy_needs_review = case_when(
        .reviewed_specific_role ~ host_role_needs_manual_review,
        .west_nile_avian_group ~ TRUE,
        .avian_influenza_generalist ~ TRUE,
        TRUE ~ host_role_needs_manual_review
      ),
      host_role_bucket = case_when(
        str_detect(modelling_role_proxy, regex("reservoir|amplifying|competence", ignore_case = TRUE)) ~
          "reservoir_or_amplifying_host",
        str_detect(modelling_role_proxy, regex("dead_end|incidental", ignore_case = TRUE)) ~
          "dead_end_or_incidental_host",
        modelling_role_proxy == "galliform_avian_host_group_proxy" ~
          "susceptible_or_spillover_host",
        str_detect(modelling_role_proxy, regex("susceptible|spillover", ignore_case = TRUE)) ~
          "susceptible_or_spillover_host",
        modelling_role_proxy == "host_presence_only" |
          str_detect(modelling_role_proxy, regex("host_presence", ignore_case = TRUE)) ~
          "host_presence_only",
        TRUE ~ "unknown_or_unreviewed"
      ),
      host_role_evidence_basis = case_when(
        modelling_role_proxy_basis == "source_backed_role_assignment" ~ "exact_source_backed",
        modelling_role_proxy_basis == "reviewed_role_assignment_needs_review" ~ "exact_reviewed_needs_review",
        str_detect(modelling_role_proxy_basis, regex("^west_nile_weighted_rule", ignore_case = TRUE)) ~
          "weighted_taxonomic_proxy",
        str_detect(modelling_role_proxy_basis, regex("_group_rule_", ignore_case = TRUE)) ~
          "disease_group_proxy",
        modelling_role_proxy_basis == "candidate_role_assignment" ~ "candidate_only",
        TRUE ~ "candidate_only"
      ),
      host_role_weight = confidence_weight(
        modelling_role_proxy_confidence,
        host_role_bucket
      ),
      role_proxy_applied = !is.na(modelling_role_proxy) &
        modelling_role_proxy != "host_presence_only",
      group_proxy_applied = str_detect(
        modelling_role_proxy_rule_id,
        regex("avian_influenza_group_proxy|wnv_.*proxy", ignore_case = TRUE)
      ),
      profile_group_proxy = group_proxy_applied & taxonomy_ok
    ) %>%
    select(-starts_with("."))
}

derive_host_tiers <- function(hosts, sdm_lookup, source_dataset, model_qc = NULL) {
  hosts %>%
    mutate(
      source_dataset = source_dataset,
      species_role = "host"
    ) %>%
    add_sdm_fields(sdm_lookup) %>%
    add_model_quality(model_qc) %>%
    mutate(
      taxonomy_ok = !missing_as_false(taxonomy_caution),
      host_detection_method = clean_text(host_detection_method),
      host_detection_tier = case_when(
        str_detect(host_detection_method, regex("PCR|Sequencing", ignore_case = TRUE)) ~ "pcr_or_sequencing",
        str_detect(host_detection_method, regex("Isolation|Observation", ignore_case = TRUE)) ~ "isolation_or_observation",
        str_detect(host_detection_method, regex("Antibod", ignore_case = TRUE)) ~ "serology",
        is.na(host_detection_method) | str_detect(host_detection_method, regex("not specified", ignore_case = TRUE)) ~ "not_specified",
        TRUE ~ "other_detection"
      ),
      host_direct_detection_supported = host_detection_tier %in% c(
        "pcr_or_sequencing",
        "isolation_or_observation"
      ),
      host_role_assignment = coalesce(clean_text(host_role_assignment), "host_presence_only"),
      host_role_confidence = coalesce(clean_text(host_role_confidence), "low"),
      host_role_assignment_status = coalesce(clean_text(host_role_assignment_status), "candidate_only"),
      host_role_needs_manual_review = missing_as_false(host_role_needs_manual_review),
      host_role_specific = host_role_assignment != "host_presence_only",
      host_role_source_backed = host_role_assignment_status == "draft_source_backed",
      host_role_medium_high = host_role_confidence %in% c("medium", "high"),
      profile_broad = taxonomy_ok,
      profile_supported = taxonomy_ok & (
        host_direct_detection_supported |
          (host_role_source_backed & host_role_specific)
      ),
      profile_strong = taxonomy_ok &
        host_role_source_backed &
        host_role_specific &
        host_role_medium_high,
      profile_strict = profile_strong & !host_role_needs_manual_review,
      no_spatial_layer = !sdm_available,
      biological_evidence_tier = max_tier(
        profile_strict,
        profile_strong,
        profile_supported,
        profile_broad
      ),
      tier_rule_id = "host_v0_1",
      missingness_reason = pmap_chr(
        list(
          if_else(taxonomy_ok, NA_character_, "taxonomy_caution"),
          if_else(sdm_available, NA_character_, "no_sdm"),
          if_else(host_role_specific, NA_character_, "host_role_presence_only"),
          if_else(host_role_needs_manual_review, "host_role_review_needed", NA_character_)
        ),
        combine_reasons
      )
    ) %>%
    add_host_modelling_proxy()
}

derive_vector_tiers <- function(vectors, sdm_lookup, source_dataset, model_qc = NULL) {
  vectors %>%
    mutate(
      source_dataset = source_dataset,
      species_role = "vector"
    ) %>%
    add_sdm_fields(sdm_lookup) %>%
    add_model_quality(model_qc) %>%
    mutate(
      taxonomy_ok = !missing_as_false(taxonomy_caution),
      has_disease_vector_evidence = missing_as_false(has_disease_vector_evidence),
      has_host_vector_evidence = missing_as_false(has_host_vector_evidence),
      has_competence_evidence = missing_as_false(has_competence_evidence),
      best_evidence_level = clean_text(best_evidence_level),
      vector_competence_status = clean_text(vector_competence_status),
      transmission_demonstrated = clean_text(transmission_demonstrated),
      natural_infection_reported = clean_text(natural_infection_reported),
      bites_humans_known = !is.na(clean_text(bites_humans)),
      bites_humans_true = missing_as_false(bites_humans),
      evidence_level_supported = best_evidence_level %in% c("probable", "confirmed"),
      competence_or_transmission_supported =
        vector_competence_status %in% c("competent", "mixed") |
          transmission_demonstrated %in% c("yes", "mixed"),
      profile_broad = taxonomy_ok & has_disease_vector_evidence,
      profile_supported = taxonomy_ok & (
        evidence_level_supported |
          has_competence_evidence |
          has_host_vector_evidence
      ),
      profile_strong = taxonomy_ok & competence_or_transmission_supported,
      profile_strict = profile_strong &
        evidence_level_supported &
        bites_humans_true,
      no_spatial_layer = !sdm_available,
      biological_evidence_tier = max_tier(
        profile_strict,
        profile_strong,
        profile_supported,
        profile_broad
      ),
      tier_rule_id = "vector_v0_1",
      missingness_reason = pmap_chr(
        list(
          if_else(taxonomy_ok, NA_character_, "taxonomy_caution"),
          if_else(is.na(clean_text(tax_id)), "tax_id_missing", NA_character_),
          if_else(sdm_available, NA_character_, "no_sdm"),
          if_else(bites_humans_known, NA_character_, "bites_humans_unknown"),
          if_else(!is.na(clean_text(vector_role_hint)), NA_character_, "vector_role_hint_blank"),
          if_else(!is.na(transmission_demonstrated), NA_character_, "transmission_demonstrated_unknown"),
          if_else(!is.na(natural_infection_reported), NA_character_, "natural_infection_unknown")
        ),
        combine_reasons
      )
    )
}

build_host_role_bucket_counts <- function(tiered_rows) {
  tiered_rows %>%
    filter(species_role == "host") %>%
    group_by(
      source_dataset,
      analysis_unit_id,
      readiness_disease_name,
      host_role_bucket,
      host_role_evidence_basis,
      modelling_role_proxy_confidence,
      modelling_role_proxy_needs_review
    ) %>%
    summarise(rows = n(), .groups = "drop") %>%
    arrange(source_dataset, readiness_disease_name, desc(rows), host_role_bucket)
}

build_review_queues <- function(tiered_rows) {
  vector_rows <- tiered_rows %>%
    filter(species_role == "vector")

  host_rows <- tiered_rows %>%
    filter(species_role == "host")

  list(
    vector_taxonomy = vector_rows %>%
      filter(is.na(clean_text(tax_id)) | !taxonomy_ok | vector_taxon_rank != "species") %>%
      transmute(
        source_dataset,
        readiness_disease_name,
        species_name,
        tax_id,
        vector_taxon_rank,
        vector_join_key,
        taxonomy_ok,
        sdm_available,
        review_reason = missingness_reason
      ) %>%
      arrange(source_dataset, readiness_disease_name, species_name),

    vector_human_biting = vector_rows %>%
      filter(!bites_humans_known, profile_supported | profile_strong) %>%
      transmute(
        source_dataset,
        readiness_disease_name,
        species_name,
        vector_group,
        best_evidence_level,
        vector_competence_status,
        transmission_demonstrated,
        has_host_vector_evidence,
        has_competence_evidence,
        sdm_available,
        review_reason = "bites_humans_unknown"
      ) %>%
      arrange(source_dataset, readiness_disease_name, species_name),

    vector_role_hints = vector_rows %>%
      filter(is.na(clean_text(vector_role_hint)), profile_strong | (bites_humans_true & profile_supported)) %>%
      transmute(
        source_dataset,
        readiness_disease_name,
        species_name,
        vector_group,
        best_evidence_level,
        vector_competence_status,
        transmission_demonstrated,
        bites_humans,
        has_host_vector_evidence,
        sdm_available,
        review_reason = "strong_vector_without_role_hint"
      ) %>%
      arrange(source_dataset, readiness_disease_name, species_name),

    host_roles = host_rows %>%
      mutate(
        review_priority_score =
          4 * as.integer(host_direct_detection_supported) +
          3 * as.integer(readiness_disease_name %in% c(
            "Chikungunya fever",
            "Dengue",
            "West Nile fever",
            "Rift Valley fever",
            "Influenza (H5N1 avian influenza)"
          )) +
          2 * as.integer(host_class %in% c("mammalia", "aves")) +
          2 * as.integer(host_order %in% c("primates", "rodentia"))
      ) %>%
      filter(host_role_assignment == "host_presence_only", review_priority_score > 0) %>%
      transmute(
        source_dataset,
        readiness_disease_name,
        species_name,
        tax_id,
        host_class,
        host_order,
        host_family,
        host_detection_method,
        sdm_available,
        review_priority_score,
        review_reason = missingness_reason
      ) %>%
      arrange(desc(review_priority_score), source_dataset, readiness_disease_name, species_name)
  )
}

write_table <- function(data, filename, dir = output_dir) {
  path <- file.path(dir, filename)
  write_csv(data, path, na = "")
  path
}

csv_row_count <- function(path) {
  max(length(readLines(path, warn = FALSE)) - 1, 0)
}

csv_column_count <- function(path) {
  length(names(readr::read_csv(path, n_max = 0, show_col_types = FALSE)))
}

# ------------------------------------------------------------------------------|
#      Inputs ------------------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_hosts <- read_optional_csv(paths$repo_hosts)
repo_vectors <- read_optional_csv(paths$repo_vectors)
repo_sdm <- read_optional_csv(paths$repo_sdm_species)

if (is.null(repo_hosts) || is.null(repo_vectors) || is.null(repo_sdm)) {
  stop("Repository pilot package inputs are missing under: ", repo_package_dir, call. = FALSE)
}

delivery_hosts <- if (delivery_available()) read_optional_csv(paths$delivery_hosts) else NULL
delivery_vectors <- if (delivery_available()) read_optional_csv(paths$delivery_vectors) else NULL
delivery_sdm <- if (delivery_available()) read_optional_csv(paths$delivery_sdm_species) else NULL
delivery_model_qc <- read_optional_csv(paths$delivery_model_qc)

# ------------------------------------------------------------------------------|
#      Tier Derivation ---------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_sdm_lookup <- build_sdm_lookup(repo_sdm)

repo_tiered <- bind_rows(
  derive_host_tiers(repo_hosts, repo_sdm_lookup, "repo_pilot"),
  derive_vector_tiers(repo_vectors, repo_sdm_lookup, "repo_pilot")
)

delivery_tiered <- tibble()
if (!is.null(delivery_hosts) && !is.null(delivery_vectors) && !is.null(delivery_sdm)) {
  delivery_sdm_lookup <- build_sdm_lookup(delivery_sdm)
  delivery_tiered <- bind_rows(
    derive_host_tiers(delivery_hosts, delivery_sdm_lookup, "chikungunya_delivery", delivery_model_qc),
    derive_vector_tiers(delivery_vectors, delivery_sdm_lookup, "chikungunya_delivery", delivery_model_qc)
  )
}

tiered_rows <- bind_rows(repo_tiered, delivery_tiered) %>%
  select(any_of(c(
    "source_dataset",
    "analysis_unit_id",
    "readiness_disease_name",
    "disease_name",
    "species_role",
    "species_name",
    "tax_id",
    "host_class",
    "host_order",
    "host_family",
    "host_detection_method",
    "host_detection_tier",
    "host_direct_detection_supported",
    "host_role_assignment",
    "host_role_confidence",
    "host_role_assignment_status",
    "host_role_needs_manual_review",
    "modelling_role_proxy",
    "modelling_role_proxy_basis",
    "modelling_role_proxy_confidence",
    "modelling_role_proxy_rule_id",
    "modelling_role_proxy_needs_review",
    "host_role_bucket",
    "host_role_evidence_basis",
    "host_role_weight",
    "role_proxy_applied",
    "group_proxy_applied",
    "vector_group",
    "vector_taxon_rank",
    "vector_join_key",
    "best_evidence_level",
    "best_evidence_basis",
    "has_disease_vector_evidence",
    "has_host_vector_evidence",
    "has_competence_evidence",
    "bites_humans",
    "bites_humans_known",
    "bites_humans_true",
    "vector_competence_status",
    "transmission_demonstrated",
    "natural_infection_reported",
    "vector_role_hint",
    "taxonomy_ok",
    "sdm_needed_for_disease",
    "sdm_available",
    "sdm_species",
    "model_quality_raw",
    "model_quality_tier",
    "model_quality_weight",
    "retained_models",
    "min_boyce",
    "max_boyce",
    "min_test_auc",
    "max_test_auc",
    "min_max_tss",
    "max_max_tss",
    "profile_broad",
    "profile_supported",
    "profile_strong",
    "profile_strict",
    "profile_group_proxy",
    "biological_evidence_tier",
    "no_spatial_layer",
    "missingness_reason",
    "tier_rule_id"
  ))) %>%
  arrange(desc(species_role == "vector"), source_dataset, readiness_disease_name, species_name)

review_queues <- build_review_queues(tiered_rows)
host_role_bucket_counts <- build_host_role_bucket_counts(tiered_rows)

# ------------------------------------------------------------------------------|
#      Minimal Writes ----------------------------------------------------------|
# ------------------------------------------------------------------------------|
outputs <- list(
  tiered_species = write_table(tiered_rows, "tiered_species.csv"),
  host_role_bucket_counts = write_table(host_role_bucket_counts, "host_role_bucket_counts.csv")
)

host_review_queue_path <- write_table(
  review_queues$host_roles,
  "review_queue_host_roles.csv",
  staged_output_dir
)

manifest <- enframe(outputs, name = "table_name", value = "path") %>%
  mutate(
    relative_path = stringr::str_remove(path, paste0("^", stringr::fixed(here::here()), "/?")),
    row_count = map_int(path, csv_row_count),
    column_count = map_int(path, csv_column_count),
    description = case_when(
      table_name == "tiered_species" ~ "One row per readiness species/taxon with derived tier/profile fields.",
      table_name == "host_role_bucket_counts" ~ "Counts of broad modelling-facing host-role buckets and evidence bases.",
      TRUE ~ "Derived prototype output."
    )
  ) %>%
  select(table_name, relative_path, row_count, column_count, description)

invisible(write_table(manifest, "manifest.csv"))

readme <- c(
  "# Modelling Evidence Tiers",
  "",
  "Generated by `scripts/associations/role_annotation/6_11_Derive_Modelling_Evidence_Tiers.R`.",
  "",
  "These files are prototype modelling/readiness outputs. They do not replace",
  "canonical role evidence, vector evidence, source-check decisions, or SDM",
  "manifests.",
  "",
  "Current scope: evidence tiers and review queues. SDM/model-quality fields are",
  "retained as contextual metadata for later, but they do not drive the current",
  "review-queue priority while SDM runs are still in progress.",
  "",
  "Host rows include both a small modelling-facing `host_role_bucket` and a",
  "more detailed `modelling_role_proxy`/rule layer. The bucket keeps modelling",
  "categories compact, while the proxy fields preserve exact-vs-group evidence,",
  "taxonomic weighting, confidence, and review status.",
  "",
  "H5N1/H7N9 avian-influenza host rows and West Nile fever bird rows include",
  "derived modelling-role proxies. These proxy roles are for sensitivity",
  "modelling and review triage; they are not source-backed species-level role",
  "assignments.",
  "",
  "## Storage Rules",
  "",
  "- Derived modelling-readiness handoff outputs live in this folder.",
  "- Generated role-curation queues live under `pathogen_association_data/staged/role_annotation/modelling_evidence_tiers/`.",
  "- Human review decisions should live under `pathogen_association_data/manual/`.",
  "- Canonical accepted role evidence remains under `pathogen_association_data/evidence/role_annotation/`.",
  "- Spatial admin extracts should later live under `pathogen_association_data/readiness/admin_modelling/`.",
  "",
  "## Main Files",
  "",
  "- `tiered_species.csv`: derived rows for the repository pilot package and, when available, the Chikungunya delivery bundle.",
  "- `host_role_bucket_counts.csv`: counts for compact host-role buckets by evidence basis.",
  "- `manifest.csv`: row/column counts for the readiness files in this folder.",
  "",
  "The generated host-role review queue is written separately to",
  "`pathogen_association_data/staged/role_annotation/modelling_evidence_tiers/review_queue_host_roles.csv`.",
  "",
  "## Inline Contract",
  "",
  "- `biological_evidence_tier`: `excluded`, `broad`, `supported`, `strong`, or `strict`.",
  "- `model_quality_tier`: `no_sdm`, `diagnostic`, `available_unscored`, `existing_host_model`, or `usable_sdm`.",
  "- `host_role_bucket`: `reservoir_or_amplifying_host`, `dead_end_or_incidental_host`, `susceptible_or_spillover_host`, `host_presence_only`, or `unknown_or_unreviewed`.",
  "- `host_role_evidence_basis`: `exact_source_backed`, `exact_reviewed_needs_review`, `disease_group_proxy`, `weighted_taxonomic_proxy`, or `candidate_only`.",
  "- `host_role_weight`: draft role weight derived from host-role confidence and bucket; it is separate from SDM/model-quality weight.",
  "",
  "Other exploratory counts can be regenerated from `tiered_species.csv` when",
  "needed; they are intentionally not written as separate files in this minimal",
  "handoff surface."
)
writeLines(readme, file.path(output_dir, "README.md"))

message("Wrote modelling evidence-tier prototype outputs to: ", output_dir)
message("Wrote generated host-role review queue to: ", host_review_queue_path)
message("Rows in tiered_species.csv: ", nrow(tiered_rows))
message("Review queue rows: host=", nrow(review_queues$host_roles))
