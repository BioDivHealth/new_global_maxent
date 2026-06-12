# ------------------------------------------------------------------------------
# Build the first pathogen_vector_links scaffold from disease-level screening
# ------------------------------------------------------------------------------

library(pacman)
p_load(here, readr)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "master_plus_compatibility_helpers.R"
))

# Normalize disease labels before joining to the screening table.
clean_disease_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x
}

# Preserve the first populated taxonomy field when duplicate network rows disagree.
first_non_empty <- function(x) {
  x <- x[!is.na(x) & trimws(x) != ""]

  if (length(x) == 0) {
    return("")
  }

  x[1]
}

screening_path <- vector_screening_manual_path("disease_vector_screening.csv")
output_dir <- vector_screening_staged_outputs_dir
output_path <- file.path(output_dir, "pathogen_vector_links.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the combined host-pathogen network and the disease-level vector screen.
combined_network <- read_legacy_compatible_master_plus_network()

screening <- read_csv(
  screening_path,
  show_col_types = FALSE
)

# Standardize disease labels so small whitespace differences do not break joins.
combined_network$disease_name_clean <- clean_disease_name(combined_network$Disease_name)
screening$disease_name_clean <- clean_disease_name(screening$disease_name_clean)

combined_network <- combined_network[
  !is.na(combined_network$disease_name_clean) &
    combined_network$disease_name_clean != "",
]

screening_lookup <- screening[
  !duplicated(screening$disease_name_clean),
]

if (any(duplicated(screening$disease_name_clean))) {
  stop("Duplicate disease_name_clean values found in disease_vector_screening.csv")
}

missing_in_screening <- setdiff(
  sort(unique(combined_network$disease_name_clean)),
  sort(unique(screening_lookup$disease_name_clean))
)

if (length(missing_in_screening) > 0) {
  stop(
    "The screening table is missing diseases from the canonical zoonotic WHO network: ",
    paste(missing_in_screening, collapse = ", ")
  )
}

# Use disease + pathogen + taxid as the seed unit for vector curation.
combined_keys <- paste(
  combined_network$disease_name_clean,
  combined_network$Pathogen,
  combined_network$PathogenTaxID,
  sep = "|||"
)

pair_seed <- aggregate(
  combined_network[
    ,
    c(
      "Disease_name",
      "PathogenType",
      "PHEIC risk",
      "PathogenFamily",
      "PathogenGenus"
    )
  ],
  by = list(
    disease_name_clean = combined_network$disease_name_clean,
    pathogen = combined_network$Pathogen,
    pathogen_tax_id = combined_network$PathogenTaxID
  ),
  FUN = first_non_empty
)

names(pair_seed)[names(pair_seed) == "Disease_name"] <- "disease_name"
names(pair_seed)[names(pair_seed) == "PathogenType"] <- "pathogen_type"
names(pair_seed)[names(pair_seed) == "PHEIC risk"] <- "pheic_risk"
names(pair_seed)[names(pair_seed) == "PathogenFamily"] <- "pathogen_family"
names(pair_seed)[names(pair_seed) == "PathogenGenus"] <- "pathogen_genus"

# Add context from the existing network so high-yield rows can be prioritised first.
pair_seed$key <- paste(
  pair_seed$disease_name_clean,
  pair_seed$pathogen,
  pair_seed$pathogen_tax_id,
  sep = "|||"
)

row_count_lookup <- as.data.frame(table(combined_keys), stringsAsFactors = FALSE)
names(row_count_lookup) <- c("key", "row_count_in_network")

host_count_lookup <- aggregate(
  combined_network$Host,
  by = list(key = combined_keys),
  FUN = function(x) length(unique(x))
)
names(host_count_lookup)[2] <- "host_count_in_network"

pair_seed <- merge(pair_seed, row_count_lookup, by = "key", all.x = TRUE)
pair_seed <- merge(pair_seed, host_count_lookup, by = "key", all.x = TRUE)
pair_seed <- merge(
  pair_seed,
  screening_lookup,
  by = "disease_name_clean",
  all.x = TRUE,
  suffixes = c("", "_screen")
)

# Only carry forward diseases that are clear vector-borne candidates or review cases.
eligible <- pair_seed[
  pair_seed$screen_status %in% c("clear", "review"),
]

eligible <- eligible[
  order(eligible$priority_tier, eligible$disease_name_clean, eligible$pathogen),
]

# Create the first pathogen-level template with blank fields ready for manual curation.
pathogen_vector_links <- data.frame(
  disease_name = eligible$disease_name,
  disease_name_clean = eligible$disease_name_clean,
  pathogen = eligible$pathogen,
  pathogen_tax_id = eligible$pathogen_tax_id,
  pathogen_type = eligible$pathogen_type,
  pheic_risk = eligible$pheic_risk,
  pathogen_family = eligible$pathogen_family,
  pathogen_genus = eligible$pathogen_genus,
  screen_status = eligible$screen_status,
  priority_tier = eligible$priority_tier,
  likely_vector_group = eligible$likely_vector_group,
  screening_source_key = eligible$screening_source_key,
  screening_basis = eligible$screening_basis,
  disease_rows_in_network = eligible$row_count_in_network,
  host_count_in_network = eligible$host_count_in_network,
  candidate_vector_species = "",
  candidate_vector_genus = "",
  candidate_vector_family = "",
  candidate_vector_group = "",
  vector_status = "",
  evidence_type = "",
  evidence_strength = "",
  source_org = "",
  source_title = "",
  source_url = "",
  source_accessed = "",
  notes = eligible$notes,
  next_action = eligible$next_action,
  stringsAsFactors = FALSE
)

# Write the scaffold CSV that will be populated with vector evidence next.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(pathogen_vector_links, output_path, na = "")

cat("Wrote", nrow(pathogen_vector_links), "pathogen-level rows to", output_path, "\n")
cat("Diseases in scope:", length(unique(pathogen_vector_links$disease_name_clean)), "\n")
