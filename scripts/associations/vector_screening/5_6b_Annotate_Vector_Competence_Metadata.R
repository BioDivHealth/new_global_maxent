# ------------------------------------------------------------------------------
# 5_6b_Annotate_Vector_Competence_Metadata.R
# ------------------------------------------------------------------------------
# Purpose: Add conservative derived metadata columns to the combined
#          vector_competence.csv artifact using source-language hints captured
#          in the notes field.
#
# Input  : diseases/vector_competence.csv
# Output : diseases/vector_competence.csv
# Adds   : vector_role_hint
#          uncertainty_reason
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, purrr)

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

collapse_flags <- function(flags) {
  flags <- flags[!is.na(flags) & flags != ""]

  if (length(flags) == 0) {
    return(NA_character_)
  }

  paste(unique(flags), collapse = "|")
}

extract_vector_role_hint <- function(note) {
  note <- clean_text(note)

  if (is.na(note)) {
    return(NA_character_)
  }

  note_lower <- stringr::str_to_lower(note)
  flags <- character()

  if (stringr::str_detect(note_lower, "\\burban\\b")) {
    flags <- c(flags, "urban_vector")
  }

  if (stringr::str_detect(note_lower, "\\bsylvatic\\b")) {
    flags <- c(flags, "sylvatic_vector")
  }

  if (stringr::str_detect(note_lower, "\\bbridge\\b")) {
    flags <- c(flags, "bridge_vector")
  }

  if (stringr::str_detect(note_lower, "\\benzootic\\b")) {
    flags <- c(flags, "enzootic_vector")
  }

  if (stringr::str_detect(note_lower, "\\bepizootic\\b")) {
    flags <- c(flags, "epizootic_vector")
  }

  if (stringr::str_detect(note_lower, "semi[- ]domestic")) {
    flags <- c(flags, "semidomestic_vector")
  }

  if (stringr::str_detect(note_lower, "\\bmaintenance\\b")) {
    flags <- c(flags, "maintenance_vector")
  }

  if (stringr::str_detect(
    note_lower,
    "\\b(primary|principal|main|major|dominant|most important|strongest)\\b"
  )) {
    flags <- c(flags, "primary_vector")
  }

  if (stringr::str_detect(
    note_lower,
    "\\b(secondary|accessory|minor|less important)\\b"
  )) {
    flags <- c(flags, "secondary_vector")
  }

  if (stringr::str_detect(
    note_lower,
    "\\b(candidate|putative vector|potential vector|regional candidate)\\b"
  )) {
    flags <- c(flags, "candidate_vector")
  }

  collapse_flags(flags)
}

extract_uncertainty_reason <- function(note) {
  note <- clean_text(note)

  if (is.na(note)) {
    return(NA_character_)
  }

  note_lower <- stringr::str_to_lower(note)
  flags <- character()

  if (stringr::str_detect(
    note_lower,
    "temperature[- ]dependent|higher temperature|tested temperatures|\\b18 c\\b|\\b27 c\\b"
  )) {
    flags <- c(flags, "temperature_dependent")
  }

  if (stringr::str_detect(
    note_lower,
    "strain[- ]dependent|genotype|lineage|subtype"
  )) {
    flags <- c(flags, "strain_dependent")
  }

  if (stringr::str_detect(
    note_lower,
    "population|populations|geographic variation|variable results across experiments|variable results across studies"
  )) {
    flags <- c(flags, "population_dependent")
  }

  if (stringr::str_detect(
    note_lower,
    "mixed evidence|one positive and one negative|controversial|debated|inconsistent|varies across studies|transmission varies|mixed lab evidence|debated role"
  )) {
    flags <- c(flags, "mixed_or_debated_evidence")
  }

  if (stringr::str_detect(
    note_lower,
    "field detection only|field infection only|field association; no transmission confirmed|positive only in reviewed studies|wnv[- ]positive pools|\\bpcr only\\b|rt[- ]pcr positive only|rt[- ]qpcr positive|field prevalence only|field isolate|review cites isolation$|review cites isolation and field detection|virus isolated from ticks|detected in tick|found infected|veev isolated|confirmed cchfv[- ]positive tick|positive in questing ticks|field positive only|field positive role unclear|associated with cchfv|rna detected"
  )) {
    flags <- c(flags, "field_detection_only")
  }

  if (stringr::str_detect(
    note_lower,
    "no transmission|did not transmit|not able to transmit|no infectious saliva|salivary rna only|no confirmed transmission|no bite transmission|transmission potential not infectious virus|weak competence; no infectious saliva|infected but did not transmit|no saliva|no dissemination or transmission"
  )) {
    flags <- c(flags, "no_transmission_demonstrated")
  }

  if (stringr::str_detect(
    note_lower,
    "unable to be infected|could not be infected|no infection"
  )) {
    flags <- c(flags, "no_infection_detected")
  }

  if (stringr::str_detect(
    note_lower,
    "poor vector|low vector competence|low competence|low efficiency|limited competence|low or no infection|very low|near-zero|moderately competent|poor laboratory vector|low to moderate efficiency|poor biological vector|weak competence|unlikely major vector|significance unclear"
  )) {
    flags <- c(flags, "low_competence_or_efficiency")
  }

  if (stringr::str_detect(note_lower, "gbm predicted|predicted")) {
    flags <- c(flags, "predicted_only")
  }

  if (stringr::str_detect(note_lower, "genus[- ]level|family[- ]level")) {
    flags <- c(flags, "taxon_level_not_specific")
  }

  if (stringr::str_detect(
    note_lower,
    "limited support|limited sample|no clear role|no conclusive studies|not ascertained|role uncertain|role still unclear|competence not assessed|competence not established|not a known vector|only in the cited passage|remains unknown"
  )) {
    flags <- c(flags, "limited_or_indirect_support")
  }

  if (stringr::str_detect(note_lower, "spelling|spelled|ocr")) {
    flags <- c(flags, "spelling_or_ocr_issue")
  }

  collapse_flags(flags)
}

input_path <- here("diseases", "vector_competence.csv")

vector_competence <- read_csv(
  input_path,
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
  "notes"
)

missing_cols <- setdiff(required_cols, names(vector_competence))
if (length(missing_cols) > 0) {
  stop(
    "vector_competence.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

vector_competence_annotated <- vector_competence %>%
  mutate(
    vector_role_hint = purrr::map_chr(notes, extract_vector_role_hint),
    uncertainty_reason = purrr::map_chr(notes, extract_uncertainty_reason)
  ) %>%
  relocate(vector_role_hint, uncertainty_reason, .after = notes)

write_csv(vector_competence_annotated, input_path, na = "")

cat("Rows written:", nrow(vector_competence_annotated), "\n")
cat(
  "Rows with vector_role_hint:",
  sum(!is.na(vector_competence_annotated$vector_role_hint)),
  "\n"
)
cat(
  "Rows with uncertainty_reason:",
  sum(!is.na(vector_competence_annotated$uncertainty_reason)),
  "\n"
)
cat("Wrote annotated competence table to", input_path, "\n")
