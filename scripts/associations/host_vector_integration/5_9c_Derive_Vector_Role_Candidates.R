# ------------------------------------------------------------------------------
# 5_9c_Derive_Vector_Role_Candidates.R
# ------------------------------------------------------------------------------
# Purpose: Derive conservative vector-role candidate flags from the integrated
#          disease-host-vector table plus the joined competence annotations.
#
# Inputs : WHO host-vector helper path for
#          disease_host_vector_links_competence_annotated.csv
# Outputs: role-annotation helper paths for WHO vector role candidates
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

collapse_flag_values <- function(x) {
  x <- clean_text(x)
  x <- unlist(stringr::str_split(stats::na.omit(x), "\\|"), use.names = FALSE)
  x <- sort(unique(stats::na.omit(clean_text(x))))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "|")
}

summarise_yes_no_mixed <- function(x) {
  values <- unique(stats::na.omit(clean_text(x)))

  if (length(values) == 0) {
    return(NA_character_)
  }

  has_yes <- "yes" %in% values
  has_no <- "no" %in% values
  has_mixed <- "mixed" %in% values

  dplyr::case_when(
    has_mixed | (has_yes & has_no) ~ "mixed",
    has_yes ~ "yes",
    has_no ~ "no",
    TRUE ~ NA_character_
  )
}

summarise_competence_status <- function(x) {
  values <- unique(stats::na.omit(clean_text(x)))

  if (length(values) == 0) {
    return(NA_character_)
  }

  has_competent <- "competent" %in% values
  has_mixed <- "mixed" %in% values
  has_not_competent <- "not_competent" %in% values

  dplyr::case_when(
    has_competent & (has_mixed | has_not_competent) ~ "mixed",
    has_competent ~ "competent",
    has_mixed ~ "mixed",
    has_not_competent ~ "not_competent",
    TRUE ~ "unclear"
  )
}

input_path <- who_network_host_vector_path(
  "disease_host_vector_links_competence_annotated.csv"
)
output_path <- role_vector_candidate_path("who")
summary_path <- role_vector_candidate_summary_path("who")

livestock_species <- c(
  "Bos taurus",
  "Bubalus bubalis",
  "Camelus bactrianus",
  "Camelus dromedarius",
  "Capra hircus",
  "Equus ferus",
  "Ovis aries",
  "Sus scrofa"
)

role_input <- read_csv(
  input_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    host_type = case_when(
      host_tax_id == 9606 ~ "human",
      !is.na(host_tax_id) ~ "nonhuman",
      TRUE ~ "unknown"
    ),
    is_human_host = host_type == "human",
    is_nonhuman_host = host_type == "nonhuman",
    is_livestock_host = host %in% livestock_species,
    is_competence_supported = !is.na(vector_competence_status) &
      vector_competence_status != "not_competent"
  ) %>%
  filter(is_competence_supported)

vector_role_candidates <- role_input %>%
  group_by(disease_name, vector_join_key) %>%
  summarise(
    vector_species = collapse_unique(dplyr::coalesce(vector_species, host_vector_species)),
    vector_group = collapse_unique(dplyr::coalesce(vector_group, vector_competence_group)),
    best_evidence_level = collapse_unique(best_evidence_level),
    best_evidence_basis = collapse_unique(best_evidence_basis),
    record_sources = collapse_unique(record_sources),
    vector_competence_status = summarise_competence_status(vector_competence_status),
    vector_competence_evidence_types = collapse_unique(vector_competence_evidence_types),
    transmission_demonstrated = summarise_yes_no_mixed(transmission_demonstrated),
    natural_infection_reported = summarise_yes_no_mixed(natural_infection_reported),
    vector_role_hint = collapse_flag_values(vector_role_hint),
    uncertainty_reason = collapse_flag_values(uncertainty_reason),
    host_species_count = n_distinct(host),
    human_host_count = n_distinct(host[is_human_host]),
    nonhuman_host_count = n_distinct(host[is_nonhuman_host]),
    livestock_host_count = n_distinct(host[is_livestock_host]),
    human_host_examples = collapse_unique(host[is_human_host]),
    nonhuman_host_examples = collapse_unique(host[is_nonhuman_host]),
    livestock_host_examples = collapse_unique(host[is_livestock_host]),
    interaction_type_examples = collapse_unique(interaction_type_examples),
    source_platform_examples = collapse_unique(source_platform_examples),
    source_dataset_examples = collapse_unique(source_dataset_examples),
    country_examples = collapse_unique(country_examples),
    competence_source_examples = collapse_unique(competence_source_examples),
    competence_note_examples = collapse_unique(competence_note_examples),
    competence_row_count = max(competence_row_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bridge_vector_candidate = human_host_count > 0 & nonhuman_host_count > 0,
    enzootic_maintenance_candidate = nonhuman_host_count > 0,
    human_amplification_candidate = human_host_count > 0,
    role_confidence = case_when(
      vector_competence_status == "competent" &
        transmission_demonstrated == "yes" &
        bridge_vector_candidate ~ "moderate",
      vector_competence_status %in% c("competent", "mixed") &
        (bridge_vector_candidate |
           enzootic_maintenance_candidate |
           human_amplification_candidate) ~ "moderate",
      vector_competence_status == "unclear" &
        (bridge_vector_candidate |
           enzootic_maintenance_candidate |
           human_amplification_candidate) ~ "low",
      TRUE ~ "low"
    ),
    role_assignment_basis = dplyr::case_when(
      bridge_vector_candidate ~
        "vector in curated disease-vector table with competence support, observed on humans and non-human hosts",
      enzootic_maintenance_candidate & human_amplification_candidate ~
        "vector in curated disease-vector table with competence support, observed on both humans and non-human hosts",
      enzootic_maintenance_candidate ~
        "vector in curated disease-vector table with competence support, observed on non-human hosts",
      human_amplification_candidate ~
        "vector in curated disease-vector table with competence support, observed on humans",
      TRUE ~ "no candidate role assigned"
    ),
    needs_manual_review = TRUE
  ) %>%
  arrange(disease_name, desc(bridge_vector_candidate), desc(nonhuman_host_count), vector_species)

summary_table <- bind_rows(
  vector_role_candidates %>%
    transmute(disease_name, metric = "bridge_vector_candidate", flag_value = bridge_vector_candidate),
  vector_role_candidates %>%
    transmute(disease_name, metric = "enzootic_maintenance_candidate", flag_value = enzootic_maintenance_candidate),
  vector_role_candidates %>%
    transmute(disease_name, metric = "human_amplification_candidate", flag_value = human_amplification_candidate)
) %>%
  count(disease_name, metric, flag_value, name = "row_count") %>%
  arrange(disease_name, metric, desc(flag_value))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(vector_role_candidates, output_path, na = "")
write_csv(summary_table, summary_path, na = "")

message("Wrote vector role candidates: ", output_path)
message("Wrote vector role summary: ", summary_path)
