# -----------------------------------------------------------------------------|
# master_plus_registry_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for master-plus registry stage scripts.
# -----------------------------------------------------------------------------|

registry_normalize_name <- function(x) {
  key <- stringr::str_to_lower(x)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[[:punct:]]+", " ")
  stringr::str_squish(key)
}

registry_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

registry_clean_key <- function(x) {
  key <- registry_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

registry_collapse_unique <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) {
    NA_character_
  } else {
    paste(x, collapse = "; ")
  }
}

registry_coalesce_chr <- function(...) {
  values <- purrr::map(list(...), as.character)
  dplyr::coalesce(!!!values)
}

registry_pick_preferred <- function(preferred_source, virion_value, clover_value) {
  dplyr::case_when(
    preferred_source == "virion" ~ as.character(virion_value),
    preferred_source == "clover" ~ as.character(clover_value),
    TRUE ~ NA_character_
  )
}

registry_build_master_bridge <- function(master_units, manual_units, master_matches) {
  master_units %>%
    dplyr::left_join(manual_units, by = "master_row") %>%
    dplyr::left_join(master_matches, by = "master_row") %>%
    dplyr::mutate(
      analysis_unit_id = dplyr::coalesce(analysis_unit_id, manual_analysis_unit_id, paste0("master_", master_row)),
      bridge_source = "disease_master_list",
      bridge_row_type = dplyr::case_when(
        combined_row_type == "existing_who_analysis_unit" & !is.na(manual_resolved_pathogen_name) ~ "master_and_who_manually_resolved",
        combined_row_type == "existing_who_analysis_unit" ~ "master_and_who_existing_resolution",
        !is.na(manual_resolved_pathogen_name) ~ "master_manual_resolution",
        TRUE ~ "master_unresolved"
      ),
      bridge_to_existing_who = combined_row_type == "existing_who_analysis_unit",
      active_master_analysis_unit = include_as_analysis_unit == "yes",
      resolved_disease_name_final = registry_coalesce_chr(
        match_resolved_disease_name,
        manual_resolved_disease_name,
        source_disease_name,
        disease_master_name
      ),
      resolved_pathogen_name_final = registry_coalesce_chr(
        match_resolved_pathogen_name,
        manual_resolved_pathogen_name,
        analysis_unit,
        source_pathogen
      ),
      resolved_pathogen_rank_final = registry_coalesce_chr(
        match_resolved_pathogen_rank,
        manual_resolved_pathogen_rank,
        analysis_unit_rank
      ),
      include_status_final = dplyr::case_when(
        !is.na(include_as_analysis_unit) ~ include_as_analysis_unit,
        analysis_decision == "keep" ~ "yes_existing_who",
        analysis_decision == "review_name_resolution" ~ "review",
        TRUE ~ "not_reviewed"
      ),
      host_query_source = preferred_match_source,
      host_query_pathogen_names = registry_pick_preferred(
        preferred_match_source,
        virion_matched_pathogen_names,
        clover_matched_pathogen_names
      ),
      host_query_taxids = registry_pick_preferred(
        preferred_match_source,
        virion_matched_taxids,
        clover_matched_taxids
      ),
      host_query_include_default = active_master_analysis_unit &
        preferred_source_match_status == "preferred_source_matched" &
        !dplyr::coalesce(match_review_flag, FALSE) &
        !dplyr::coalesce(shared_species_proxy_flag, FALSE),
      host_query_bucket = dplyr::case_when(
        include_status_final == "hold" ~ "hold",
        include_status_final == "review" ~ "manual_review",
        !active_master_analysis_unit ~ "inactive_or_not_reviewed",
        dplyr::coalesce(shared_species_proxy_flag, FALSE) ~ "shared_species_proxy_review",
        dplyr::coalesce(match_review_flag, FALSE) ~ "match_review",
        preferred_source_match_status == "preferred_source_matched" ~ "default_clean",
        overall_match_status == "matched_or_candidate" ~ "fallback_source_review",
        TRUE ~ "unmatched"
      )
    ) %>%
    dplyr::select(
      bridge_source,
      bridge_row_type,
      bridge_to_existing_who,
      analysis_unit_id,
      master_row,
      disease_master_name,
      resolved_disease_name_final,
      resolved_pathogen_name_final,
      resolved_pathogen_rank_final,
      include_status_final,
      active_master_analysis_unit,
      split_group,
      pathogen_aliases,
      pathogen_family_master,
      master_tier,
      master_guild,
      master_livestock_amplified,
      master_key_host_vector,
      in_master_who,
      in_master_gibb,
      in_master_empres_i,
      in_master_atlas,
      master_gbif_checked,
      master_notes,
      name_resolution_status,
      existing_lookup_name,
      match_field,
      row_type,
      family,
      pheic_risk,
      source_pathogen,
      source_previous_name,
      source_msl39_viral_name,
      source_disease_name,
      is_priority_pathogen,
      is_prototype_pathogen,
      in_gibb_etal,
      in_empres_i,
      priority_prototype_status,
      region_africa,
      region_americas,
      region_europe,
      region_mediterranean,
      region_se_asia,
      region_western_pacific,
      source_unit_scope,
      analysis_unit,
      analysis_unit_label,
      analysis_unit_rank,
      analysis_decision,
      decision_rule_trigger,
      transmission_context,
      human_infection_status,
      host_link_status,
      vector_data_status,
      amplifier_data_status,
      example_members,
      rationale,
      notes,
      resolution_source,
      resolution_notes,
      clover_matched_pathogen_names,
      virion_matched_pathogen_names,
      clover_matched_taxids,
      virion_matched_taxids,
      clover_matched_families,
      virion_matched_families,
      clover_matched_source_types,
      virion_matched_source_types,
      clover_best_match_type,
      virion_best_match_type,
      clover_match_status,
      virion_match_status,
      preferred_match_source,
      overall_match_status,
      preferred_source_match_status,
      match_review_flag,
      shared_species_proxy_flag,
      match_review_notes,
      host_query_include_default,
      host_query_bucket,
      host_query_source,
      host_query_pathogen_names,
      host_query_taxids
    ) %>%
    dplyr::arrange(master_row)
}

registry_build_who_only_units <- function(who_units, bridge) {
  master_keys_for_who_overlap <- bridge %>%
    dplyr::transmute(
      analysis_unit_key = registry_clean_key(analysis_unit),
      analysis_unit_label_key = registry_clean_key(analysis_unit_label),
      source_disease_key = registry_clean_key(source_disease_name),
      source_pathogen_key = registry_clean_key(source_pathogen)
    )

  who_units %>%
    dplyr::filter(
      !who_analysis_unit_key %in% master_keys_for_who_overlap$analysis_unit_key,
      !who_analysis_unit_label_key %in% master_keys_for_who_overlap$analysis_unit_label_key,
      !who_source_disease_key %in% master_keys_for_who_overlap$source_disease_key,
      !who_source_pathogen_key %in% master_keys_for_who_overlap$source_pathogen_key
    ) %>%
    dplyr::transmute(
      bridge_source = "who_pathogen_analysis_units",
      bridge_row_type = "who_only_existing_analysis_unit",
      bridge_to_existing_who = TRUE,
      analysis_unit_id = paste0("who_", who_unit_row),
      master_row = NA_integer_,
      disease_master_name = NA_character_,
      resolved_disease_name_final = registry_coalesce_chr(source_disease_name, analysis_unit_label),
      resolved_pathogen_name_final = registry_coalesce_chr(analysis_unit, source_pathogen),
      resolved_pathogen_rank_final = analysis_unit_rank,
      include_status_final = dplyr::if_else(analysis_decision == "keep", "yes_existing_who", "review"),
      active_master_analysis_unit = FALSE,
      split_group = NA_character_,
      pathogen_aliases = NA_character_,
      pathogen_family_master = NA_character_,
      master_tier = NA_character_,
      master_guild = NA_character_,
      master_livestock_amplified = NA,
      master_key_host_vector = NA_character_,
      in_master_who = FALSE,
      in_master_gibb = FALSE,
      in_master_empres_i = FALSE,
      in_master_atlas = FALSE,
      master_gbif_checked = FALSE,
      master_notes = NA_character_,
      name_resolution_status = "who_only_existing_unit",
      existing_lookup_name = NA_character_,
      match_field = NA_character_,
      row_type,
      family,
      pheic_risk,
      source_pathogen,
      source_previous_name,
      source_msl39_viral_name,
      source_disease_name,
      is_priority_pathogen,
      is_prototype_pathogen,
      in_gibb_etal,
      in_empres_i,
      priority_prototype_status,
      region_africa,
      region_americas,
      region_europe,
      region_mediterranean,
      region_se_asia,
      region_western_pacific,
      source_unit_scope,
      analysis_unit,
      analysis_unit_label,
      analysis_unit_rank,
      analysis_decision,
      decision_rule_trigger,
      transmission_context,
      human_infection_status,
      host_link_status,
      vector_data_status,
      amplifier_data_status,
      example_members,
      rationale,
      notes,
      resolution_source = "existing_who_analysis_unit",
      resolution_notes = NA_character_,
      clover_matched_pathogen_names = NA_character_,
      virion_matched_pathogen_names = NA_character_,
      clover_matched_taxids = NA_character_,
      virion_matched_taxids = NA_character_,
      clover_matched_families = NA_character_,
      virion_matched_families = NA_character_,
      clover_matched_source_types = NA_character_,
      virion_matched_source_types = NA_character_,
      clover_best_match_type = NA_character_,
      virion_best_match_type = NA_character_,
      clover_match_status = NA_character_,
      virion_match_status = NA_character_,
      preferred_match_source = NA_character_,
      overall_match_status = NA_character_,
      preferred_source_match_status = NA_character_,
      match_review_flag = NA,
      shared_species_proxy_flag = NA,
      match_review_notes = NA_character_,
      host_query_include_default = FALSE,
      host_query_bucket = "who_only_needs_master_match",
      host_query_source = NA_character_,
      host_query_pathogen_names = NA_character_,
      host_query_taxids = NA_character_
    )
}

registry_compact_analysis_units <- function(combined_units) {
  combined_units %>%
    dplyr::transmute(
      row_type,
      family,
      pheic_risk,
      source_pathogen,
      source_previous_name,
      source_msl39_viral_name,
      source_disease_name,
      is_priority_pathogen,
      is_prototype_pathogen,
      in_gibb_etal,
      in_empres_i,
      priority_prototype_status,
      region_africa,
      region_americas,
      region_europe,
      region_mediterranean,
      region_se_asia,
      region_western_pacific,
      source_unit_scope,
      analysis_unit,
      analysis_unit_label,
      analysis_unit_rank,
      analysis_decision,
      decision_rule_trigger,
      transmission_context,
      human_infection_status,
      host_link_status,
      vector_data_status,
      amplifier_data_status,
      example_members,
      rationale,
      notes,
      vectored_status,
      generalist_status,
      transmission_complexity,
      guild,
      host_sdm_needed,
      vector_sdm_needed,
      host_range_rule,
      vector_range_rule,
      range_limiting_layer,
      transmission_rule_notes,
      transmission_rule_review_status,
      modelling_scope_status,
      modelling_scope_reason,
      analysis_unit_id,
      master_row,
      disease_master_name,
      master_tier,
      master_guild,
      master_livestock_amplified,
      master_key_host_vector,
      include_as_analysis_unit = include_status_final,
      preferred_match_source,
      matched_pathogen_names = host_query_pathogen_names,
      matched_taxids = host_query_taxids,
      match_review_flag = dplyr::coalesce(match_review_flag, FALSE),
      shared_species_proxy_flag = dplyr::coalesce(shared_species_proxy_flag, FALSE),
      match_review_notes,
      host_query_bucket
    )
}

registry_build_host_query_units <- function(combined_units) {
  combined_units %>%
    dplyr::filter(
      bridge_source == "disease_master_list",
      include_status_final == "yes"
    ) %>%
    dplyr::transmute(
      analysis_unit_id,
      master_row,
      disease_master_name,
      resolved_disease_name = resolved_disease_name_final,
      resolved_pathogen_name = resolved_pathogen_name_final,
      resolved_pathogen_rank = resolved_pathogen_rank_final,
      preferred_match_source,
      host_query_include_default,
      host_query_bucket,
      host_query_source,
      host_query_pathogen_names,
      host_query_taxids,
      match_review_flag = dplyr::coalesce(match_review_flag, FALSE),
      shared_species_proxy_flag = dplyr::coalesce(shared_species_proxy_flag, FALSE),
      match_review_notes,
      split_group,
      pathogen_family_master,
      master_tier,
      master_guild,
      master_livestock_amplified,
      master_key_host_vector
    ) %>%
    dplyr::arrange(host_query_bucket, master_row)
}

registry_make_source_matches <- function(query, source_table, source_name, manual_aliases, max_dist = 0.08) {
  source_proc <- source_table %>%
    dplyr::filter(!is.na(source_pathogen_name), source_pathogen_name != "") %>%
    dplyr::mutate(source_key = registry_normalize_name(source_pathogen_name)) %>%
    dplyr::distinct(source, source_pathogen_name, source_taxid, source_family, source_type, source_key)

  exact_matches <- query %>%
    dplyr::inner_join(source_proc, by = c("query_key" = "source_key")) %>%
    dplyr::mutate(match_type = "exact", match_distance = 0)

  alias_matches <- query %>%
    dplyr::inner_join(
      manual_aliases %>%
        dplyr::filter(source == .env$source_name) %>%
        dplyr::select(query_key, source_key, alias_source_name = source_name, alias_type, alias_notes, alias_review_flag),
      by = "query_key"
    ) %>%
    dplyr::inner_join(source_proc, by = "source_key", relationship = "many-to-many") %>%
    dplyr::mutate(match_type = "manual_alias", match_distance = 0)

  unmatched <- query %>%
    dplyr::filter(!analysis_unit_id %in% c(exact_matches$analysis_unit_id, alias_matches$analysis_unit_id))

  fuzzy_matches <- tibble::tibble()
  if (nrow(unmatched) > 0 && nrow(source_proc) > 0) {
    fuzzy_matches <- purrr::map_dfr(seq_len(nrow(unmatched)), function(i) {
      query_row <- unmatched[i, ]
      distances <- stringdist::stringdist(
        query_row$query_key,
        source_proc$source_key,
        method = "jw"
      )
      keep <- which(distances <= max_dist)
      if (length(keep) == 0) {
        return(tibble::tibble())
      }
      keep <- keep[order(distances[keep], source_proc$source_pathogen_name[keep])]
      keep <- head(keep, 5)

      dplyr::bind_cols(
        query_row[rep(1, length(keep)), ],
        source_proc[keep, ] %>% dplyr::select(-source)
      ) %>%
        dplyr::mutate(
          match_distance = distances[keep],
          match_type = "fuzzy_candidate"
        )
    })
  }

  dplyr::bind_rows(exact_matches, alias_matches, fuzzy_matches) %>%
    dplyr::mutate(source = source_name) %>%
    dplyr::select(
      analysis_unit_id, master_row, disease_master_name, resolved_disease_name,
      resolved_pathogen_name, resolved_pathogen_rank, include_as_analysis_unit,
      split_group, source, source_pathogen_name, source_taxid, source_family,
      source_type, match_type, match_distance, alias_type, alias_notes, alias_review_flag
    ) %>%
    dplyr::distinct()
}
