# ------------------------------------------------------------------------------
# 5_3_Combine_LitReview_EFSA_Vector_Table.R
# ------------------------------------------------------------------------------
# Purpose: Merge EFSA pathogen-vector links into the curated `vector_table.xlsx`
#          while keeping the same table structure used for lit-review rows.
#
# Inputs : vector_table.xlsx
#          pathogen_vector_links_efsa.csv (from 5_2_EFSA_Vector_Crosswalk.R)
#
# Output : vector_table_with_efsa.csv
#          vector_table_with_efsa.xlsx
# ------------------------------------------------------------------------------

# ------------------------------| Load libraries |------------------------------
library(pacman)
p_load(dplyr, here, readr, readxl, stringr, writexl)

source(here("scripts", "associations", "working_inputs.R"))

# ------------------------------| Helper functions |----------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "NaN")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

map_efsa_status_to_evidence_level <- function(x) {
  x <- stringr::str_to_lower(clean_text(x))
  dplyr::case_when(
    x == "highly_likely" ~ "probable",
    x == "potential" ~ "candidate",
    TRUE ~ "candidate"
  )
}

to_vector_group <- function(x) {
  x <- stringr::str_to_lower(clean_text(x))
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  x
}

# ------------------------------| Define paths |--------------------------------
root_dir <- here()
vector_table_path <- prefer_existing_path(
  file.path(root_dir, "diseases", "vector_table.xlsx"),
  file.path(root_dir, "vector_table.xlsx")
)
efsa_links_path <- vector_screening_efsa_staged_path("pathogen_vector_links_efsa.csv")

output_dir <- vector_screening_efsa_outputs_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_csv <- file.path(output_dir, "vector_table_with_efsa.csv")
output_xlsx <- file.path(output_dir, "vector_table_with_efsa.xlsx")

# ------------------------------| Expected schema |-----------------------------
required_vector_cols <- c(
  "disease",
  "v_species",
  "v_group",
  "evidence_level",
  "evidence_basis",
  "location",
  "source",
  "notes",
  "record_source"
)
required_input_cols <- setdiff(required_vector_cols, "record_source")

# ------------------------------| Load lit-review table |-----------------------
vector_table <- read_excel(vector_table_path) %>%
  mutate(across(everything(), clean_text))

missing_vector_cols <- setdiff(required_input_cols, names(vector_table))
if (length(missing_vector_cols) > 0) {
  stop(
    "vector_table.xlsx is missing required columns: ",
    paste(missing_vector_cols, collapse = ", ")
  )
}

vector_table <- vector_table %>%
  select(all_of(required_input_cols)) %>%
  mutate(record_source = "lit_review")

# ------------------------------| Load EFSA links |-----------------------------
efsa_links <- read_csv(
  efsa_links_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

# ------------------------------| Map EFSA -> vector schema |-------------------
efsa_rows <- efsa_links %>%
  transmute(
    disease = coalesce(combined_disease_name, efsa_pathogen),
    v_species = vector_species,
    v_group = to_vector_group(vector_group),
    evidence_level = map_efsa_status_to_evidence_level(vector_status),
    evidence_basis = "efsa_appendix_g",
    location = "not specified",
    source = "EFSA Appendix G",
    notes = case_when(
      !is.na(efsa_pathogen) & !is.na(match_confidence) ~ paste0(
        "efsa_pathogen=", efsa_pathogen,
        "; match_confidence=", match_confidence,
        "; match_method=", match_method
      ),
      !is.na(efsa_pathogen) ~ paste0("efsa_pathogen=", efsa_pathogen),
      TRUE ~ NA_character_
    ),
    record_source = "EFSA"
  ) %>%
  filter(!is.na(disease), !is.na(v_species)) %>%
  select(all_of(required_vector_cols))

# ------------------------------| Combine and save |----------------------------
# Requested behavior: append all EFSA rows without deduplicating against
# existing lit-review rows in vector_table.xlsx.
vector_table_with_efsa <- bind_rows(vector_table, efsa_rows) %>%
  select(all_of(required_vector_cols))

write_csv(vector_table_with_efsa, output_csv, na = "")
write_xlsx(vector_table_with_efsa, output_xlsx)

# ------------------------------| Console summary |-----------------------------
cat("Original vector_table rows:", nrow(vector_table), "\n")
cat("EFSA rows candidate after cleaning:", nrow(efsa_rows), "\n")
cat("EFSA rows appended:", nrow(efsa_rows), "\n")
cat("Combined rows written:", nrow(vector_table_with_efsa), "\n")
cat("Wrote combined CSV to", output_csv, "\n")
cat("Wrote combined XLSX to", output_xlsx, "\n")
