################################################################################
# host_proxy_rules.R
################################################################################
# Purpose: Derive modelling-facing host-role proxy fields from reviewed role
#          assignments and explicit disease/taxonomic proxy rules.
#
# Rule contract:
# - exact source-backed assignments win;
# - reviewed assignments win before broad proxies;
# - group proxies stay review-needed unless explicitly accepted upstream;
# - role weights are separate from SDM or model-quality weights.
################################################################################

host_proxy_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "NULL", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

host_proxy_clean_key <- function(x) {
  x %>%
    host_proxy_clean_text() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("&", " and ") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

host_proxy_missing_as_false <- function(x) {
  out <- as.character(x) %in% c("TRUE", "true", "True", "1", "yes", "Yes", "YES")
  out[is.na(out)] <- FALSE
  out
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

host_role_confidence_weight <- function(confidence, bucket) {
  dplyr::case_when(
    bucket == "unknown_or_unreviewed" ~ 0,
    bucket == "host_presence_only" & confidence == "low" ~ 0.1,
    confidence == "high" ~ 1,
    confidence == "medium" ~ 0.75,
    confidence == "medium_low" ~ 0.5,
    confidence == "low" ~ 0.25,
    TRUE ~ 0.1
  )
}

add_host_modelling_proxy <- function(data) {
  data %>%
    dplyr::mutate(
      .species_key = host_proxy_clean_key(species_name),
      .tax_id_key = host_proxy_clean_text(tax_id),
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
      modelling_role_proxy = dplyr::case_when(
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
      modelling_role_proxy_basis = dplyr::case_when(
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
      modelling_role_proxy_confidence = dplyr::case_when(
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
      modelling_role_proxy_rule_id = dplyr::case_when(
        .source_backed_specific_role ~ "source_backed_role_v0_1",
        .reviewed_specific_role ~ "reviewed_role_needs_review_v0_1",
        .west_nile_corvidae_group ~ "wnv_corvidae_weighted_proxy_v0_1",
        .west_nile_passeriform_group ~ "wnv_passeriformes_weighted_proxy_v0_1",
        .west_nile_charadriiform_group ~ "wnv_charadriiformes_weighted_proxy_v0_1",
        .west_nile_avian_group ~ "wnv_aves_group_proxy_v0_1",
        .avian_influenza_generalist ~ "avian_influenza_group_proxy_v0_1",
        TRUE ~ "candidate_role_v0_1"
      ),
      modelling_role_proxy_needs_review = dplyr::case_when(
        .reviewed_specific_role ~ host_proxy_missing_as_false(host_role_needs_manual_review),
        .west_nile_avian_group ~ TRUE,
        .avian_influenza_generalist ~ TRUE,
        TRUE ~ host_proxy_missing_as_false(host_role_needs_manual_review)
      ),
      host_role_bucket = dplyr::case_when(
        stringr::str_detect(
          modelling_role_proxy,
          stringr::regex("reservoir|amplifying|competence", ignore_case = TRUE)
        ) ~ "reservoir_or_amplifying_host",
        stringr::str_detect(
          modelling_role_proxy,
          stringr::regex("dead_end|incidental", ignore_case = TRUE)
        ) ~ "dead_end_or_incidental_host",
        modelling_role_proxy == "galliform_avian_host_group_proxy" ~
          "susceptible_or_spillover_host",
        stringr::str_detect(
          modelling_role_proxy,
          stringr::regex("susceptible|spillover", ignore_case = TRUE)
        ) ~ "susceptible_or_spillover_host",
        modelling_role_proxy == "host_presence_only" |
          stringr::str_detect(
            modelling_role_proxy,
            stringr::regex("host_presence", ignore_case = TRUE)
          ) ~ "host_presence_only",
        TRUE ~ "unknown_or_unreviewed"
      ),
      host_role_evidence_basis = dplyr::case_when(
        modelling_role_proxy_basis == "source_backed_role_assignment" ~ "exact_source_backed",
        modelling_role_proxy_basis == "reviewed_role_assignment_needs_review" ~
          "exact_reviewed_needs_review",
        stringr::str_detect(
          modelling_role_proxy_basis,
          stringr::regex("^west_nile_weighted_rule", ignore_case = TRUE)
        ) ~ "weighted_taxonomic_proxy",
        stringr::str_detect(
          modelling_role_proxy_basis,
          stringr::regex("_group_rule_", ignore_case = TRUE)
        ) ~ "disease_group_proxy",
        modelling_role_proxy_basis == "candidate_role_assignment" ~ "candidate_only",
        TRUE ~ "candidate_only"
      ),
      host_role_weight = host_role_confidence_weight(
        modelling_role_proxy_confidence,
        host_role_bucket
      ),
      role_proxy_applied = !is.na(modelling_role_proxy) &
        modelling_role_proxy != "host_presence_only",
      group_proxy_applied = stringr::str_detect(
        modelling_role_proxy_rule_id,
        stringr::regex("avian_influenza_group_proxy|wnv_.*proxy", ignore_case = TRUE)
      ),
      profile_group_proxy = group_proxy_applied & taxonomy_ok
    ) %>%
    dplyr::select(-dplyr::starts_with("."))
}
