#!/usr/bin/env Rscript
################################################################################
# 6_4_Prepare_Deep_Research_Batch_Inputs.R
################################################################################
# Purpose: Prepare prompt-specific attachment folders for Deep Research
#          role-evidence review batches.
#
# Inputs : ROLE_EVIDENCE_FULL_CURATION_PLAN.md
#          species_host_vector_roster.csv
#          host_role_candidates.csv
#          disease_vector_links_taxonomy_cleaned_competence_annotated.csv
#          existing role review markdowns
#
# Outputs: pathogen_association_data/WHO/role_annotation/deep_research_inputs/
#            README.md
#            batch_manifest.csv
#            <batch_id>/PROMPT.md
#            <batch_id>/ATTACHMENT_LIST.md
#            <batch_id>/ROLE_EVIDENCE_FULL_CURATION_PLAN.md
#            <batch_id>/species_host_vector_roster_batch.csv
#            <batch_id>/host_role_candidates_batch.csv
#            <batch_id>/vector_candidates_batch.csv, for Phase V only
#            <batch_id>/existing_reviews_combined.md
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, readr, stringr, tibble)

source(here::here("scripts", "associations", "working_inputs.R"))

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

collapse_prompt <- function(lines) {
  paste(lines, collapse = "\n")
}

slug_review_file <- function(slug) {
  file.path(reviews_dir, paste0(slug, "_role_review.md"))
}

# Preserve missing-review markers in the batch attachment so absent markdowns are
# visible to the reviewer instead of silently dropping disease context.
combine_reviews <- function(review_paths, output_path) {
  blocks <- lapply(review_paths, function(path) {
    if (!file.exists(path)) {
      return(c(
        paste0("# Missing Review: ", basename(path)),
        "",
        "This expected review markdown was not present when the batch input pack was generated."
      ))
    }

    c(
      paste0("<!-- Source file: ", path, " -->"),
      readLines(path, warn = FALSE)
    )
  })

  write_text(unlist(Map(c, blocks, list(c("", "\n---\n"))), use.names = FALSE), output_path)
}

# Build phase-specific instructions while keeping vector and non-vector batches
# under the same table schema expectations for downstream parsing.
make_prompt <- function(batch_title, diseases, phase, has_vector_csv, cautions) {
  disease_lines <- paste0(seq_along(diseases), ". ", diseases)
  attached <- c(
    "- ROLE_EVIDENCE_FULL_CURATION_PLAN.md",
    "- Batch-filtered species_host_vector_roster_batch.csv",
    "- Batch-filtered host_role_candidates_batch.csv",
    if (has_vector_csv) "- Batch-filtered vector_candidates_batch.csv" else NULL,
    "- existing_reviews_combined.md"
  )

  phase_text <- if (phase == "Phase V") {
    c(
      "I am curating source-backed host/vector ecological role evidence for a disease modelling repo.",
      "",
      "These are vectored diseases in this curation phase."
    )
  } else {
    c(
      "I am curating source-backed host ecological role evidence for a disease modelling repo.",
      "",
      "These are non-vectored diseases in this curation phase. Do not create vector role rows unless strong source evidence shows an arthropod vector role, which would be unexpected and should be flagged as an upstream scope issue."
    )
  }

  vector_warning <- if (has_vector_csv) {
    "Use the attached vector candidate/competence CSV only as candidate/context, not as final vector-role proof."
  } else {
    "For this non-vectored batch, record vector status as not applicable unless strong source evidence indicates an arthropod vector role; flag any such finding as an upstream scope issue."
  }

  vector_section <- if (phase == "Phase V") {
    c(
      "3. Vector role evidence table:",
      "   disease_name, vector_or_group, taxonomic_grain, role_claim, evidence_direction, confidence, source_id, source_title, source_type, source_url, doi_pmid_pmc, source_access, short_evidence_span, caveat, assignment_recommendation.",
      "4. Assignment/staging recommendations table:",
      "   disease_name, entity_type, entity_name, role_assignment, assignment_confidence, action, evidence_source_ids, evidence_basis, needs_manual_review, review_reason.",
      "5. Sources used table:",
      "   source_id, disease_name, source_title, authors_or_organization, year, source_type, source_url, doi, pmid, pmcid, source_access, rows_supported, reliability_note.",
      "6. Deferred candidates / unresolved issues."
    )
  } else {
    c(
      "3. Vector non-applicability note:",
      "   disease_name, vector_status, evidence_basis, caveat.",
      "4. Assignment/staging recommendations table:",
      "   disease_name, entity_name, role_assignment, assignment_confidence, action, evidence_source_ids, evidence_basis, needs_manual_review, review_reason.",
      "5. Sources used table:",
      "   source_id, disease_name, source_title, authors_or_organization, year, source_type, source_url, doi, pmid, pmcid, source_access, rows_supported, reliability_note.",
      "6. Deferred candidates / unresolved issues."
    )
  }

  source_policy <- if (phase == "Phase V") {
    c(
      "Source policy:",
      "Prioritize peer-reviewed reviews, primary ecological/vector/reservoir studies, WHO/WOAH/FAO/ECDC technical reports, and official CDC/WHO/ECDC pages.",
      "Do not use news, generic health websites, pest-control websites, Wikipedia-style summaries, blogs, or uncited secondary pages as evidence for role assignments.",
      "Official public-health pages may support broad transmission-cycle claims, but do not use them to make species-level assignments unless the page explicitly names the species and role.",
      "For disputed or species-level reservoir/vector assignments, prefer peer-reviewed review or primary research evidence."
    )
  } else {
    c(
      "Source policy:",
      "Prioritize peer-reviewed reviews, primary reservoir/ecology studies, WHO/WOAH/FAO/ECDC technical reports, and official CDC/WHO pages.",
      "Do not use news, generic health websites, Wikipedia-style summaries, blogs, or uncited secondary pages as evidence for role assignments.",
      "Official pages may support broad transmission-cycle claims, but do not use them for species-level reservoir assignments unless they explicitly name the species and role.",
      "For disputed reservoirs, prefer peer-reviewed reviews and primary studies, and keep confidence conservative."
    )
  }

  role_distinctions <- if (phase == "Phase V") {
    "Distinguish reservoir host, amplifying host, incidental host, dead-end host, spillover host, primary vector, bridge vector, sylvatic vector, enzootic vector, epidemic vector, candidate vector, and field-detection-only evidence."
  } else {
    "Distinguish suspected reservoir, confirmed reservoir, spillover host, amplifying host, incidental host, susceptible host only, and human-to-human transmission contexts."
  }

  collapse_prompt(c(
    phase_text,
    "",
    "Batch:",
    batch_title,
    "",
    "Batch diseases:",
    disease_lines,
    "",
    "Attached files:",
    attached,
    "",
    "Use the attached CSVs only to understand which host/vector taxa are present in my repo.",
    "Do not treat candidate presence, host presence, vector presence, GenBank geography, outbreak country evidence, natural infection, pathogen detection, or vector competence alone as ecological-role proof.",
    vector_warning,
    "",
    "Goal:",
    "Produce auditable evidence rows suitable for later conversion into CSV. Do not produce generic disease summaries.",
    "I will later convert your output into staging CSVs and then into repo role evidence tables, so favor explicit rows and stable source metadata over prose.",
    "",
    source_policy,
    "",
    "Important rules:",
    "- Separate evidence rows from assignment recommendations.",
    "- Do not assign species-level roles from group-level source language.",
    "- Mark group-level, regional, cycle-specific, outbreak-specific, or disputed evidence clearly.",
    paste0("- ", role_distinctions),
    "- Treat attached CSV rows as candidates/context only, even if a species appears many times.",
    "- If a claim is based only on vector competence, natural infection, PCR/pathogen detection, host presence, outbreak geography, or GenBank geography, label it as evidence-only or deferred rather than a role assignment.",
    "- Do not use Deep Research citation handles alone as provenance. Include stable source identifiers: title, authors or organization, year, URL, and DOI/PMID/PMCID where available.",
    "- Use plain URLs/DOIs/PMIDs/PMCIDs in the tables, not only hidden links or turn/file citation handles.",
    "- For `source_access`, use one of: full_text, abstract_only, official_page, report_pdf, dataset_page, unclear.",
    "- For `action`, use one of: already_covered, candidate_add_after_source_check, evidence_only_group, defer_vocabulary_or_taxonomy, defer_insufficient_evidence, reject_or_negative_evidence_only.",
    "- Use `already_covered` when the existing review markdown indicates the same row is already represented in my repo.",
    "- Use `candidate_add_after_source_check` only when the source support looks strong enough that a human should verify and potentially add the row.",
    "- Use `evidence_only_group` for class/order/genus/group rows that should not be propagated to species.",
    "- Use `defer_vocabulary_or_taxonomy` for species pairs, complexes, uncertain synonyms, missing candidate taxa, or useful roles not represented by the current vocabulary.",
    "- Use `reject_or_negative_evidence_only` for negative, poor competence, unsupported, or contradicted role claims.",
    "",
    "Return these sections for each disease:",
    "1. Brief source-search summary.",
    "2. Host role evidence table:",
    "   disease_name, host_or_group, taxonomic_grain, role_claim, evidence_direction, confidence, source_id, source_title, source_type, source_url, doi_pmid_pmc, source_access, short_evidence_span, caveat, assignment_recommendation.",
    vector_section,
    "",
    "At the end of the whole report, add a cross-batch summary table:",
    "disease_name, entity_type, entity_name, proposed_role_or_issue, action, priority_for_manual_review, reason.",
    "",
    "Do not omit low-confidence/deferred findings if they are important for preventing overclaiming; include them with the correct action and caveat.",
    "",
    "Batch-specific cautions:",
    paste0("- ", cautions)
  ))
}

# ------------------------------------------------------------------------------|
#      Define paths ------------------------------------------------------------|
# ------------------------------------------------------------------------------|
repo_root <- here::here()
repo_relative_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  sub(paste0("^", root, "/?"), "", path)
}

role_dir <- role_annotation_dir
reviews_dir <- role_reviews_dir
output_root <- role_deep_research_dir
plan_path <- here::here("ROLE_EVIDENCE_FULL_CURATION_PLAN.md")
roster_path <- file.path(role_roster_dir, "species_host_vector_roster.csv")
host_candidates_path <- file.path(role_candidates_dir, "host_role_candidates.csv")
vector_candidates_path <- here::here(
  "pathogen_association_data", "WHO", "vector_screening", "outputs",
  "disease_vector_links_taxonomy_cleaned_competence_annotated.csv"
)

# ------------------------------------------------------------------------------|
#      Load input layers -------------------------------------------------------|
# ------------------------------------------------------------------------------|
if (!file.exists(plan_path)) {
  stop("Missing plan file: ", plan_path)
}

roster <- read_csv(roster_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

host_candidates <- read_csv(host_candidates_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

vector_candidates <- read_csv(vector_candidates_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(across(where(is.character), clean_text))

# ------------------------------------------------------------------------------|
#      Define review batches ---------------------------------------------------|
# ------------------------------------------------------------------------------|
# Batches are intentionally small enough for manual source review, with cautions
# carrying disease-specific overclaim guardrails into the prompt.
batches <- list(
  list(
    id = "v_batch_1_west_nile_yellow_fever_dengue_plague",
    phase = "Phase V",
    title = "V Batch 1: West Nile, Yellow fever, Dengue, Plague",
    diseases = c("West Nile fever", "Yellow fever", "Dengue", "Plague"),
    review_slugs = c("west_nile_fever", "yellow_fever", "dengue", "plague"),
    cautions = c(
      "West Nile: keep bird reservoir evidence group-level unless species-specific support exists; distinguish human/horse dead-end roles.",
      "Yellow fever: distinguish sylvatic, intermediate, and urban cycles.",
      "Dengue: distinguish Aedes aegypti primary vector from Aedes albopictus secondary/candidate vector; humans are urban-cycle amplifying hosts only.",
      "Plague: distinguish maintenance hosts, amplifying epizootic hosts, incidental hosts, and flea vectors."
    )
  ),
  list(
    id = "v_batch_2_rvf_chikungunya_zika_vee",
    phase = "Phase V",
    title = "V Batch 2: Rift Valley fever, Chikungunya, Zika, VEE",
    diseases = c(
      "Rift Valley fever",
      "Chikungunya fever",
      "Zika virus disease",
      "Venezuelan equine encephalitis"
    ),
    review_slugs = c("rift_valley_fever", "chikungunya_fever", "zika_virus_disease", "venezuelan_equine_encephalitis"),
    cautions = c(
      "Rift Valley fever: separate livestock amplifying hosts from wildlife reservoir uncertainty; be careful with vertical transmission and mechanical-vector wording.",
      "Chikungunya: distinguish urban Aedes vectors from sylvatic/enzootic cycles and competence-only taxa.",
      "Zika: distinguish urban Aedes vectors, sylvatic non-human primate evidence, human amplification, and incidental hosts.",
      "Venezuelan equine encephalitis: distinguish enzootic rodent/mosquito cycles, equine amplification, epidemic vectors, and human roles."
    )
  ),
  list(
    id = "v_batch_3_cchf_tbe_sfts_oropouche",
    phase = "Phase V",
    title = "V Batch 3: CCHF, TBE, SFTS, Oropouche",
    diseases = c(
      "Crimean-Congo hemorrhagic fever",
      "Tick-borne encephalitis",
      "Severe fever with thrombocytopenia syndrome (SFTS)",
      "Oropouche fever"
    ),
    review_slugs = c("crimean_congo_hemorrhagic_fever", "tick_borne_encephalitis", "sfts", "oropouche_fever"),
    cautions = c(
      "CCHF: distinguish Hyalomma tick vector/reservoir roles from livestock amplifying/asymptomatic host roles and human spillover/nosocomial transmission.",
      "TBE: distinguish small mammal reservoir/amplifying hosts, tick vectors/reservoirs, indicator hosts, and human incidental roles.",
      "SFTS: distinguish tick vectors, animal hosts, possible amplifying hosts, and human-to-human transmission.",
      "Oropouche: carefully assess Culicoides paraensis primary-vector evidence, mosquito candidate evidence, and uncertain sylvatic vertebrate hosts."
    )
  ),
  list(
    id = "n_batch_1_ebola_sudan_marburg_nipah",
    phase = "Phase N",
    title = "N Batch 1: Ebola, Sudan virus disease, Marburg, Nipah",
    diseases = c(
      "Ebola virus disease",
      "Sudan virus disease (Ebola virus disease)",
      "Marburg virus disease",
      "Nipah virus disease"
    ),
    review_slugs = c("ebola_virus_disease", "sudan_virus_disease", "marburg_virus_disease", "nipah_virus_disease"),
    cautions = c(
      "Ebola virus disease: avoid species-level bat assignments unless source-supported; distinguish non-human primate spillover/source animals from reservoirs.",
      "Sudan virus disease: do a targeted check for Sudan-virus-specific reservoir evidence; do not rely only on generic ebolavirus statements.",
      "Marburg: assess Rousettus aegyptiacus evidence carefully and distinguish confirmed reservoir from broader bat group language.",
      "Nipah: separate Pteropus reservoir evidence, pig amplifying host evidence, human infection, and context-dependent human-to-human transmission."
    )
  ),
  list(
    id = "n_batch_2_lassa_hantaan_argentine_hf",
    phase = "Phase N",
    title = "N Batch 2: Lassa, Hantaan/HFRS, Argentine hemorrhagic fever",
    diseases = c(
      "Lassa fever",
      "Hemorrhagic fever with renal syndrome (Hantaan virus)",
      "Argentine hemorrhagic fever"
    ),
    review_slugs = c("lassa_fever", "hfrs_hantaan", "argentine_hemorrhagic_fever"),
    cautions = c(
      "Lassa fever: distinguish Mastomys natalensis from broader Mastomys/multimammate mouse/rodent language; note human-to-human transmission context separately.",
      "Hantaan/HFRS: keep Hantaan virus reservoir evidence separate from broader hantavirus/HFRS evidence; check Apodemus agrarius support.",
      "Argentine hemorrhagic fever: keep Junin virus / Calomys musculinus evidence distinct from other New World arenaviruses."
    )
  ),
  list(
    id = "n_batch_3_h5n1_mpox",
    phase = "Phase N",
    title = "N Batch 3: H5N1, Mpox",
    diseases = c(
      "Influenza (H5N1 avian influenza)",
      "Mpox (Monkeypox)"
    ),
    review_slugs = c("influenza_h5n1_avian_influenza", "mpox_monkeypox"),
    cautions = c(
      "H5N1: separate wild aquatic bird reservoir/natural-host evidence, domestic poultry amplification, dairy cattle susceptibility/outbreak context, mammalian spillover, and human incidental/spillover infection. Do not assign all birds identically.",
      "Mpox: the animal reservoir is uncertain. Distinguish suspected reservoir groups, susceptible animals, outbreak-associated animals, non-human primates, rodents/small mammals, and human amplification/person-to-person transmission."
    )
  )
)

# ------------------------------------------------------------------------------|
#      Write batch packs -------------------------------------------------------|
# ------------------------------------------------------------------------------|
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

manifest_rows <- list()

for (batch in batches) {
  batch_dir <- file.path(output_root, batch$id)
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

  batch_roster <- roster %>%
    filter(disease_name %in% batch$diseases) %>%
    arrange(disease_name, species_role, species_name)

  batch_hosts <- host_candidates %>%
    filter(disease_name %in% batch$diseases) %>%
    arrange(disease_name, host)

  write_csv(batch_roster, file.path(batch_dir, "species_host_vector_roster_batch.csv"), na = "")
  write_csv(batch_hosts, file.path(batch_dir, "host_role_candidates_batch.csv"), na = "")

  vector_rows <- NA_integer_
  if (batch$phase == "Phase V") {
    batch_vectors <- vector_candidates %>%
      filter(disease_name %in% batch$diseases) %>%
      arrange(disease_name, vector_species_taxonomy_cleaned)

    write_csv(batch_vectors, file.path(batch_dir, "vector_candidates_batch.csv"), na = "")
    vector_rows <- nrow(batch_vectors)
  }

  # The full curation plan and existing reviews give the external review the
  # current repo guardrails and already-reviewed role decisions.
  file.copy(plan_path, file.path(batch_dir, "ROLE_EVIDENCE_FULL_CURATION_PLAN.md"), overwrite = TRUE)

  review_paths <- vapply(batch$review_slugs, slug_review_file, character(1))
  combine_reviews(review_paths, file.path(batch_dir, "existing_reviews_combined.md"))

  prompt <- make_prompt(
    batch_title = batch$title,
    diseases = batch$diseases,
    phase = batch$phase,
    has_vector_csv = batch$phase == "Phase V",
    cautions = batch$cautions
  )
  write_text(prompt, file.path(batch_dir, "PROMPT.md"))

  attachment_lines <- c(
    paste0("# Attachment List: ", batch$title),
    "",
    "Paste the contents of `PROMPT.md` into Deep Research and attach these files:",
    "",
    "- `ROLE_EVIDENCE_FULL_CURATION_PLAN.md`",
    "- `species_host_vector_roster_batch.csv`",
    "- `host_role_candidates_batch.csv`",
    if (batch$phase == "Phase V") "- `vector_candidates_batch.csv`" else NULL,
    "- `existing_reviews_combined.md`",
    "",
    "Do not attach GenBank country tables or WHO DON outputs for this task.",
    "Those are geography/outbreak evidence layers, not host/vector role proof.",
    "",
    "Batch diseases:",
    paste0("- ", batch$diseases)
  )
  write_text(attachment_lines, file.path(batch_dir, "ATTACHMENT_LIST.md"))

  manifest_rows[[batch$id]] <- tibble(
    batch_id = batch$id,
    phase = batch$phase,
    batch_title = batch$title,
    disease_count = length(batch$diseases),
    roster_rows = nrow(batch_roster),
    host_candidate_rows = nrow(batch_hosts),
    vector_candidate_rows = vector_rows,
    folder = repo_relative_path(file.path(output_root, batch$id))
  )
}

# ------------------------------------------------------------------------------|
#      Write manifest and README ----------------------------------------------|
# ------------------------------------------------------------------------------|
manifest <- bind_rows(manifest_rows)
write_csv(manifest, file.path(output_root, "batch_manifest.csv"), na = "")

readme_lines <- c(
  "# Deep Research Input Packs",
  "",
  "These folders contain prompt-specific files to upload with each Deep Research batch.",
  "",
  "For each batch:",
  "",
  "1. Open the batch folder.",
  "2. Paste `PROMPT.md` into Deep Research.",
  "3. Attach the files listed in `ATTACHMENT_LIST.md`.",
  "4. After Deep Research finishes, bring the output back to the repo for schema conversion and overclaim review.",
  "",
  "The batch CSVs are filtered candidate/context files only. Candidate presence, host presence, vector presence, GenBank geography, WHO DON outbreak-country evidence, natural infection, pathogen detection, and vector competence are not role proof by themselves.",
  "",
  "Batch folders:",
  paste0("- `", manifest$batch_id, "`: ", manifest$batch_title),
  "",
  "Generated by:",
  "",
  "`Rscript scripts/associations/role_annotation/6_4_Prepare_Deep_Research_Batch_Inputs.R`"
)
write_text(readme_lines, file.path(output_root, "README.md"))

message("Wrote Deep Research input packs to: ", output_root)
print(manifest)
