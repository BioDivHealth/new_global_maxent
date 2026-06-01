# ------------------------------------------------------------------------------
# 5_6c_Join_Vector_Competence_Evidence.R
# ------------------------------------------------------------------------------
# Purpose: Collapse the extracted vector competence evidence to disease-vector
#          grain and left-join it onto the canonical disease-vector and
#          disease-host-vector outputs without changing their row inclusion logic.
#
# Inputs : diseases/vector_competence.csv
#          disease_vector_links_taxonomy_cleaned.csv
#          disease_host_vector_links.csv
#          disease_host_vector_links_expanded.csv
# Outputs: vector_competence_collapsed.csv
#          vector_competence_join_unmatched.csv
#          disease_vector_links_taxonomy_cleaned_competence_annotated.csv
#          disease_host_vector_links_competence_annotated.csv
#          disease_host_vector_links_expanded_competence_annotated.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr)

source(here("scripts", "associations", "working_inputs.R"))
source(here("scripts", "associations", "association_text_helpers.R"))

normalize_vector_key <- function(x) {
  x <- clean_text(x)
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "[/]", " ")
  x <- stringr::str_replace_all(x, "[-–—]", " ")
  x <- stringr::str_replace_all(x, "[()\\[\\],.;:*'`\"]", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_flag_values <- function(x) {
  x <- clean_text(x)
  x <- unlist(stringr::str_split(stats::na.omit(x), "\\|"), use.names = FALSE)
  x <- sort(unique(stats::na.omit(clean_text(x))))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "|")
}

summarise_competence_status <- function(x) {
  statuses <- unique(stats::na.omit(clean_text(x)))

  if (length(statuses) == 0) {
    return(NA_character_)
  }

  has_competent <- "competent" %in% statuses
  has_mixed <- "mixed" %in% statuses
  has_not_competent <- "not_competent" %in% statuses

  dplyr::case_when(
    has_competent & (has_mixed | has_not_competent) ~ "mixed",
    has_competent ~ "competent",
    has_mixed ~ "mixed",
    has_not_competent ~ "not_competent",
    TRUE ~ "unclear"
  )
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

apply_vector_name_cleanup <- function(x, manual_map) {
  cleaned <- normalize_vector_key(x)
  method <- rep("normalized_name", length(cleaned))

  abbreviations <- c(
    "ae\\." = "aedes",
    "cx\\." = "culex",
    "oc\\." = "ochlerotatus",
    "an\\." = "anopheles"
  )

  for (pattern in names(abbreviations)) {
    matched <- !is.na(cleaned) & stringr::str_detect(cleaned, paste0("^", pattern, "\\s+"))
    cleaned[matched] <- stringr::str_replace(
      cleaned[matched],
      paste0("^", pattern),
      abbreviations[[pattern]]
    )
    method[matched] <- "rule_expand_genus_abbreviation"
  }

  abbreviations_no_period <- c(
    "ae" = "aedes",
    "cx" = "culex",
    "oc" = "ochlerotatus",
    "an" = "anopheles"
  )

  for (token in names(abbreviations_no_period)) {
    matched <- !is.na(cleaned) & stringr::str_detect(cleaned, paste0("^", token, "\\s+"))
    cleaned[matched] <- stringr::str_replace(
      cleaned[matched],
      paste0("^", token, "\\b"),
      abbreviations_no_period[[token]]
    )
    method[matched] <- "rule_expand_genus_abbreviation"
  }

  parenthetical_subgenus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^[a-z]+ \\([a-z]+\\)\\s+")
  cleaned[parenthetical_subgenus] <- stringr::str_replace(
    cleaned[parenthetical_subgenus],
    "^([a-z]+) \\([a-z]+\\)\\s+",
    "\\1 "
  )
  method[parenthetical_subgenus] <- "rule_drop_parenthetical_subgenus"

  repeated_genus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^([a-z]+) \\1\\b")
  cleaned[repeated_genus] <- stringr::str_replace(
    cleaned[repeated_genus],
    "^([a-z]+) \\1\\s+",
    "\\1 "
  )
  method[repeated_genus] <- "rule_repeated_genus"

  aedes_subgenus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^aedes (ochlerotatus|neomelaniconion)\\b")
  cleaned[aedes_subgenus] <- stringr::str_replace(
    cleaned[aedes_subgenus],
    "^aedes (ochlerotatus|neomelaniconion)\\s+",
    "aedes "
  )
  method[aedes_subgenus] <- "rule_drop_subgenus_token"

  culex_subgenus <- !is.na(cleaned) &
    stringr::str_detect(cleaned, "^culex (melanoconion|culex)\\b")
  cleaned[culex_subgenus] <- stringr::str_replace(
    cleaned[culex_subgenus],
    "^culex (melanoconion|culex)\\s+",
    "culex "
  )
  method[culex_subgenus] <- "rule_drop_subgenus_token"

  map_key <- normalize_vector_key(cleaned)
  manual_source <- normalize_vector_key(manual_map$source_name)
  mapped <- manual_map$canonical_name[match(map_key, manual_source)]
  has_manual_map <- !is.na(mapped)

  cleaned[has_manual_map] <- mapped[has_manual_map]
  method[has_manual_map] <- "manual_map"

  tibble::tibble(
    vector_join_key = normalize_vector_key(cleaned),
    vector_competence_name_cleanup_method = method
  )
}

standardise_competence_disease_name <- function(x) {
  x_clean <- clean_text(x)
  x_join <- normalize_name_for_match(x_clean)

  dplyr::case_when(
    x_join == "west nile virus" ~ "West Nile fever",
    x_join == "zika" ~ "Zika virus disease",
    x_join == "sftsv" ~ "Severe fever with thrombocytopenia syndrome (SFTS)",
    x_join == "oropouche virus" ~ "Oropouche fever",
    x_join == "venezuelan equine encephalitis virus" ~ "Venezuelan equine encephalitis",
    x_join == "crimean congo hemorrhagic fever" ~ "Crimean-Congo hemorrhagic fever",
    x_join == "dengue" ~ "Dengue",
    x_join == "plague" ~ "Plague",
    x_join == "yellow fever" ~ "Yellow fever",
    x_join == "tick borne encephalitis" ~ "Tick-borne encephalitis",
    x_join == "rift valley fever" ~ "Rift Valley fever",
    x_join == "chikungunya" ~ "Chikungunya fever",
    x_join == "chikungunya virus" ~ "Chikungunya fever",
    TRUE ~ x_clean
  )
}

annotate_with_competence <- function(input_path, output_path, competence_collapsed) {
  input_table <- read_csv(
    input_path,
    show_col_types = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(across(where(is.character), clean_text))

  if ("vector_join_key" %in% names(input_table)) {
    join_vector_key <- normalize_vector_key(input_table$vector_join_key)
  } else if ("vector_species_taxonomy_cleaned" %in% names(input_table)) {
    join_vector_key <- normalize_vector_key(input_table$vector_species_taxonomy_cleaned)
  } else {
    stop("No vector join column found in ", input_path)
  }

  joinable <- input_table %>%
    mutate(
      disease_name_join = normalize_name_for_match(disease_name),
      vector_join_key = join_vector_key
    )

  annotated <- joinable %>%
    left_join(
      competence_collapsed,
      by = c("disease_name_join", "vector_join_key")
    ) %>%
    select(-disease_name_join)

  if (nrow(annotated) != nrow(input_table)) {
    stop("Row count changed while annotating ", input_path)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(annotated, output_path, na = "")

  tibble::tibble(
    output_path = output_path,
    input_rows = nrow(input_table),
    annotated_rows = sum(!is.na(clean_text(annotated$vector_competence_status)))
  )
}

competence_path <- here("diseases", "vector_competence.csv")

disease_vector_path <- vector_screening_evidence_path(
  "disease_vector_links_taxonomy_cleaned.csv"
)
disease_vector_annotated_path <- file.path(
  vector_screening_evidence_dir,
  "disease_vector_links_taxonomy_cleaned_competence_annotated.csv"
)
dhv_path <- who_network_host_vector_path("disease_host_vector_links.csv")
dhv_annotated_path <- who_network_host_vector_path(
  "disease_host_vector_links_competence_annotated.csv"
)
dhv_expanded_path <- who_network_host_vector_path("disease_host_vector_links_expanded.csv")
dhv_expanded_annotated_path <- who_network_host_vector_path(
  "disease_host_vector_links_expanded_competence_annotated.csv"
)
collapsed_path <- file.path(vector_screening_evidence_dir, "vector_competence_collapsed.csv")
unmatched_path <- file.path(vector_screening_qa_dir, "vector_competence_join_unmatched.csv")
manual_map_path <- vector_screening_taxonomy_manual_path("vector_taxonomy_manual_map.csv")
dir.create(vector_screening_evidence_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_screening_qa_dir, recursive = TRUE, showWarnings = FALSE)

manual_map <- read_csv(
  manual_map_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(source_name = stringr::str_to_lower(source_name)) %>%
  distinct(source_name, .keep_all = TRUE)

competence <- read_csv(
  competence_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

required_cols <- c(
  "disease",
  "v_species",
  "v_group",
  "competence_status",
  "evidence_type",
  "transmission_demonstrated",
  "natural_infection_reported",
  "location",
  "source",
  "notes",
  "vector_role_hint",
  "uncertainty_reason"
)

missing_cols <- setdiff(required_cols, names(competence))
if (length(missing_cols) > 0) {
  stop(
    "vector_competence.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

vector_cleanup <- apply_vector_name_cleanup(competence$v_species, manual_map)

competence_joinable <- competence %>%
  bind_cols(vector_cleanup) %>%
  mutate(
    disease_name = standardise_competence_disease_name(disease),
    disease_name_join = normalize_name_for_match(disease_name)
  )

competence_collapsed <- competence_joinable %>%
  filter(!is.na(disease_name_join), !is.na(vector_join_key)) %>%
  group_by(disease_name_join, vector_join_key) %>%
  summarise(
    competence_disease_name = collapse_unique(disease_name),
    vector_competence_species = collapse_unique(v_species),
    vector_competence_group = collapse_unique(v_group),
    vector_competence_status = summarise_competence_status(competence_status),
    competence_statuses = collapse_unique(competence_status),
    vector_competence_evidence_types = collapse_unique(evidence_type),
    transmission_demonstrated = summarise_yes_no_mixed(transmission_demonstrated),
    natural_infection_reported = summarise_yes_no_mixed(natural_infection_reported),
    vector_role_hint = collapse_flag_values(vector_role_hint),
    uncertainty_reason = collapse_flag_values(uncertainty_reason),
    competence_locations_summary = collapse_unique(location),
    competence_source_examples = collapse_unique(source),
    competence_note_examples = collapse_unique(notes),
    competence_row_count = dplyr::n(),
    vector_competence_name_cleanup_methods = collapse_unique(vector_competence_name_cleanup_method),
    .groups = "drop"
  ) %>%
  arrange(competence_disease_name, vector_join_key)

write_csv(competence_collapsed, collapsed_path, na = "")

disease_vector_join_keys <- read_csv(
  disease_vector_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text)) %>%
  transmute(
    disease_name_join = normalize_name_for_match(disease_name),
    vector_join_key = normalize_vector_key(vector_species_taxonomy_cleaned)
  ) %>%
  distinct()

unmatched_competence <- competence_collapsed %>%
  anti_join(
    disease_vector_join_keys,
    by = c("disease_name_join", "vector_join_key")
  ) %>%
  select(
    competence_disease_name,
    vector_join_key,
    vector_competence_species,
    vector_competence_status,
    vector_competence_evidence_types,
    transmission_demonstrated,
    natural_infection_reported,
    competence_source_examples,
    competence_note_examples
  ) %>%
  arrange(competence_disease_name, vector_join_key)

write_csv(unmatched_competence, unmatched_path, na = "")

annotation_summaries <- bind_rows(
  annotate_with_competence(
    disease_vector_path,
    disease_vector_annotated_path,
    competence_collapsed
  ),
  annotate_with_competence(
    dhv_path,
    dhv_annotated_path,
    competence_collapsed
  ),
  annotate_with_competence(
    dhv_expanded_path,
    dhv_expanded_annotated_path,
    competence_collapsed
  )
)

cat("Vector competence rows read:", nrow(competence), "\n")
cat("Collapsed disease-vector competence rows:", nrow(competence_collapsed), "\n")
cat("Competence rows not matching disease-vector table:", nrow(unmatched_competence), "\n")
print(annotation_summaries)
cat("Wrote collapsed competence table to", collapsed_path, "\n")
cat("Wrote unmatched competence QA table to", unmatched_path, "\n")
