#!/usr/bin/env Rscript
################################################################################
# 6_5_Reformat_Deep_Research_Reports.R
################################################################################
# Purpose: Reformat Deep Research markdown reports into cleaned,
#          reviewable staging artifacts without changing official role tables.
#
# Inputs : Deep Research markdown reports from ~/Downloads
# Outputs: pathogen_association_data/WHO/role_annotation/deep_research_inputs/
#            <batch_id>/reformatted/
#              DEEP_RESEARCH_REFORMATTED.md
#              extracted_<table_type>.csv
#              extraction_quality_issues.csv
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    stop("Package `pacman` is required.", call. = FALSE)
  }
})

pacman::p_load(dplyr, purrr, readr, stringr, tibble)

# ------------------------------------------------------------------------------|
#      Define paths and report map --------------------------------------------|
# ------------------------------------------------------------------------------|
role_dir <- here::here("pathogen_association_data", "WHO", "role_annotation")
input_root <- file.path(role_dir, "deep_research_inputs")

report_specs <- tribble(
  ~source_path, ~batch_id, ~batch_title, ~phase,
  "~/Downloads/deep-research-report(6).md",
  "v_batch_1_west_nile_yellow_fever_dengue_plague",
  "V Batch 1: West Nile fever, Yellow fever, Dengue, Plague",
  "Phase V",
  "~/Downloads/deep-research-report(1).md",
  "v_batch_2_rvf_chikungunya_zika_vee",
  "V Batch 2: Rift Valley fever, Chikungunya, Zika, VEE",
  "Phase V",
  "~/Downloads/deep-research-report(2).md",
  "v_batch_3_cchf_tbe_sfts_oropouche",
  "V Batch 3: CCHF, TBE, SFTS, Oropouche",
  "Phase V",
  "~/Downloads/deep-research-report(3).md",
  "n_batch_1_ebola_sudan_marburg_nipah",
  "N Batch 1: Ebola, Sudan virus disease, Marburg, Nipah",
  "Phase N",
  "~/Downloads/deep-research-report(5).md",
  "n_batch_2_lassa_hantaan_argentine_hf",
  "N Batch 2: Lassa fever, HFRS/Hantaan, Argentine hemorrhagic fever",
  "Phase N",
  "~/Downloads/deep-research-report(4).md",
  "n_batch_3_h5n1_mpox",
  "N Batch 3: H5N1 avian influenza, Mpox",
  "Phase N"
) %>%
  mutate(source_path = path.expand(source_path))

# ------------------------------------------------------------------------------|
#      Define allowed vocabulary ----------------------------------------------|
# ------------------------------------------------------------------------------|
# These allowed sets flag drift from the repo's current role/action vocabulary;
# flagged rows remain staging artifacts until manually source-checked.
allowed_actions <- c(
  "already_covered",
  "candidate_add_after_source_check",
  "evidence_only_group",
  "defer_vocabulary_or_taxonomy",
  "defer_insufficient_evidence",
  "reject_or_negative_evidence_only"
)

allowed_host_roles <- c(
  "reservoir_host",
  "reservoir_host_group",
  "amplifying_host",
  "maintenance_host",
  "incidental_host",
  "dead_end_host",
  "dead_end_incidental_host",
  "spillover_host",
  "susceptible_host_only",
  "host_presence_only",
  "negative_or_unsupported_role",
  "unclear_role"
)

allowed_vector_roles <- c(
  "primary_vector",
  "principal_vector_genus",
  "main_vector",
  "epidemic_vector",
  "enzootic_maintenance_vector",
  "bridge_vector",
  "amplificatory_vector",
  "sylvatic_vector",
  "candidate_vector",
  "competent_vector_only",
  "field_detection_only",
  "mechanical_vector",
  "not_competent",
  "negative_or_unsupported_role",
  "unclear_role"
)

allowed_source_access <- c(
  "full_text",
  "abstract_only",
  "official_page",
  "report_pdf",
  "dataset_page",
  "unclear"
)

allowed_confidence <- c("high", "medium", "low")

# ------------------------------------------------------------------------------|
#      Helpers -----------------------------------------------------------------|
# ------------------------------------------------------------------------------|
write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

# External reports may contain citation wrappers or encoded URL markers. Strip
# those wrappers so downstream CSVs keep plain, reviewable source metadata.
clean_deep_research_markup <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(
    x,
    "\ue200entity\ue202\\[\"disease\",\"([^\"]+)\"[^\\]]*\\]\ue201",
    "\\1"
  )
  x <- str_replace_all(
    x,
    "\ue200url\ue202([^\ue202\ue201]+)(?:\ue202[^\ue201]*)?\ue201",
    "\\1"
  )
  x <- str_replace_all(x, "\ue200(?:cite|filecite)\ue202[^\ue201]+\ue201", "")
  x <- str_replace_all(x, "\ue200[^\ue201]*\ue201", "")
  x <- str_replace_all(x, "\\s+([.;,])", "\\1")
  x <- str_replace_all(x, "[ \t]+", " ")
  str_trim(x)
}

split_markdown_row <- function(line) {
  placeholder <- "\ue010"
  line <- str_trim(line)
  line <- str_remove(line, "^\\|")
  line <- str_remove(line, "\\|$")
  line <- str_replace_all(line, "\\\\[|]", placeholder)
  cells <- str_split(line, "\\|", simplify = FALSE)[[1]]
  cells <- str_replace_all(cells, fixed(placeholder), "|")
  clean_deep_research_markup(cells)
}

# Markdown table separator rows are used to identify tables without relying on
# section names, which vary across reports.
is_table_separator <- function(line) {
  str_detect(str_trim(line), "^\\|?[[:space:]:\\-\\|]+\\|?[[:space:]]*$")
}

normalize_header <- function(x) {
  x %>%
    clean_deep_research_markup() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

infer_table_type <- function(header) {
  if ("host_or_group" %in% header) return("host_role_evidence")
  if ("vector_or_group" %in% header) return("vector_role_evidence")
  if ("vector_status" %in% header) return("vector_non_applicability")
  if ("proposed_role_or_issue" %in% header) return("cross_batch_summary")
  if ("candidate_or_issue" %in% header) return("deferred_unresolved")
  if ("source_id" %in% header && "rows_supported" %in% header) return("sources_used")
  if ("role_assignment" %in% header && "action" %in% header) return("assignment_staging")
  "other_table"
}

drop_empty_columns <- function(df) {
  keep <- vapply(df, function(column) {
    if (!is.character(column)) {
      return(!all(is.na(column)))
    }
    !all(is.na(column) | column == "")
  }, logical(1))
  df[, keep, drop = FALSE]
}

extract_disease_heading <- function(line) {
  line <- clean_deep_research_markup(line)
  if (!str_detect(line, "^## ")) return(NA_character_)
  line <- str_remove(line, "^## +")
  if (str_detect(str_to_lower(line), "^cross-batch")) return(NA_character_)
  str_trim(line)
}

# ------------------------------------------------------------------------------|
#      Extract markdown tables -------------------------------------------------|
# ------------------------------------------------------------------------------|
extract_tables <- function(lines, spec) {
  disease <- NA_character_
  tables <- list()
  i <- 1L

  while (i <= length(lines)) {
    disease_hit <- extract_disease_heading(lines[[i]])
    if (!is.na(disease_hit)) {
      disease <- disease_hit
    }

    is_table_start <- str_detect(str_trim(lines[[i]]), "^\\|") &&
      i < length(lines) &&
      is_table_separator(lines[[i + 1L]])

    if (!is_table_start) {
      i <- i + 1L
      next
    }

    header <- normalize_header(split_markdown_row(lines[[i]]))
    table_type <- infer_table_type(header)
    j <- i + 2L
    row_lines <- character()

    while (j <= length(lines) && str_detect(str_trim(lines[[j]]), "^\\|")) {
      row_lines <- c(row_lines, lines[[j]])
      j <- j + 1L
    }

    if (length(row_lines) > 0) {
      rows <- lapply(row_lines, split_markdown_row)
      width <- length(header)
      rows <- lapply(rows, function(row) {
        if (length(row) < width) {
          row <- c(row, rep(NA_character_, width - length(row)))
        }
        if (length(row) > width) {
          row <- c(row[seq_len(width - 1L)], paste(row[width:length(row)], collapse = " | "))
        }
        row
      })

      table <- as_tibble(do.call(rbind, rows), .name_repair = "minimal")
      names(table) <- header
      section_disease <- if (table_type == "cross_batch_summary") {
        "Cross-batch summary"
      } else {
        disease
      }
      table <- table %>%
        mutate(
          batch_id = spec$batch_id,
          batch_title = spec$batch_title,
          phase = spec$phase,
          source_report = spec$source_path,
          section_disease = section_disease,
          table_type = table_type,
          .before = 1
        ) %>%
        mutate(across(where(is.character), clean_deep_research_markup))

      tables[[length(tables) + 1L]] <- table
    }

    i <- j
  }

  bind_rows(tables)
}

# ------------------------------------------------------------------------------|
#      Build parser and vocabulary QA -----------------------------------------|
# ------------------------------------------------------------------------------|
issue_table <- function(tables) {
  issues <- list()

  if (nrow(tables) == 0) {
    return(tibble(issue_type = "no_tables_extracted", detail = "No markdown tables were extracted."))
  }

  optional_columns <- c(
    "action",
    "source_access",
    "role_claim",
    "role_assignment",
    "confidence",
    "assignment_confidence",
    "disease_name",
    "entity_name",
    "host_or_group",
    "vector_or_group"
  )
  for (column in optional_columns) {
    if (!column %in% names(tables)) {
      tables[[column]] <- NA_character_
    }
  }

  if ("action" %in% names(tables)) {
    issues[[length(issues) + 1L]] <- tables %>%
      filter(!is.na(action), action != "", !action %in% allowed_actions) %>%
      transmute(
        batch_id,
        table_type,
        disease_name = coalesce(.data$disease_name, section_disease),
        entity = coalesce(.data$entity_name, .data$host_or_group, .data$vector_or_group, NA_character_),
        issue_type = "action_not_in_allowed_set",
        detail = action
      )
  }

  if ("source_access" %in% names(tables)) {
    issues[[length(issues) + 1L]] <- tables %>%
      filter(!is.na(source_access), source_access != "") %>%
      filter(!str_detect(source_access, paste(allowed_source_access, collapse = "|"))) %>%
      transmute(
        batch_id,
        table_type,
        disease_name = coalesce(.data$disease_name, section_disease),
        entity = coalesce(.data$host_or_group, .data$vector_or_group, .data$entity_name, NA_character_),
        issue_type = "source_access_not_standard",
        detail = source_access
      )
  }

  if ("role_claim" %in% names(tables)) {
    issues[[length(issues) + 1L]] <- tables %>%
      filter(table_type == "host_role_evidence", !is.na(role_claim), role_claim != "") %>%
      filter(!role_claim %in% allowed_host_roles) %>%
      transmute(
        batch_id,
        table_type,
        disease_name = coalesce(.data$disease_name, section_disease),
        entity = host_or_group,
        issue_type = "host_role_claim_not_in_current_repo_vocab",
        detail = role_claim
      )

    issues[[length(issues) + 1L]] <- tables %>%
      filter(table_type == "vector_role_evidence", !is.na(role_claim), role_claim != "") %>%
      filter(!role_claim %in% allowed_vector_roles) %>%
      transmute(
        batch_id,
        table_type,
        disease_name = coalesce(.data$disease_name, section_disease),
        entity = vector_or_group,
        issue_type = "vector_role_claim_not_in_current_repo_vocab",
        detail = role_claim
      )
  }

  if ("role_assignment" %in% names(tables)) {
    issues[[length(issues) + 1L]] <- tables %>%
      filter(table_type == "assignment_staging", !is.na(role_assignment), role_assignment != "") %>%
      filter(str_detect(role_assignment, " |confirmed_reservoir|suspected_reservoir|maintenance vector|positive vector role|vector role|cryptic|issue|no role found|negative against")) %>%
      transmute(
        batch_id,
        table_type,
        disease_name = coalesce(.data$disease_name, section_disease),
        entity = coalesce(.data$entity_name, NA_character_),
        issue_type = "assignment_value_needs_normalization_or_is_issue_not_role",
        detail = role_assignment
      )
  }

  confidence_columns <- c("confidence", "assignment_confidence")
  for (column in confidence_columns) {
    if (column %in% names(tables)) {
      issues[[length(issues) + 1L]] <- tables %>%
        filter(!is.na(.data[[column]]), .data[[column]] != "") %>%
        filter(!.data[[column]] %in% allowed_confidence) %>%
        transmute(
          batch_id,
          table_type,
          disease_name = coalesce(.data$disease_name, section_disease),
          entity = coalesce(.data$entity_name, .data$host_or_group, .data$vector_or_group, NA_character_),
          issue_type = paste0(column, "_not_in_allowed_set"),
          detail = .data[[column]]
        )
    }
  }

  bind_rows(issues)
}

write_clean_markdown <- function(lines, output_path) {
  cleaned <- clean_deep_research_markup(lines)
  cleaned <- str_replace_all(cleaned, "^## +", "## ")
  write_text(cleaned, output_path)
}

# ------------------------------------------------------------------------------|
#      Check candidate presence in local layers --------------------------------|
# ------------------------------------------------------------------------------|
present_anywhere <- function(path, value) {
  if (!file.exists(path) || is.na(value) || value == "") {
    return(NA)
  }

  data <- read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.character), clean_deep_research_markup))
  char_cols <- names(data)[vapply(data, is.character, logical(1))]
  if (length(char_cols) == 0) {
    return(FALSE)
  }

  value_clean <- str_to_lower(str_squish(value))
  any(vapply(data[char_cols], function(column) {
    any(str_to_lower(str_squish(column)) == value_clean, na.rm = TRUE)
  }, logical(1)))
}

candidate_presence_check <- function(tables, spec) {
  required <- c(
    "table_type", "action", "disease_name", "entity_type", "entity_name",
    "role_assignment", "assignment_confidence", "review_reason"
  )
  for (column in required) {
    if (!column %in% names(tables)) {
      tables[[column]] <- NA_character_
    }
  }

  candidates <- tables %>%
    filter(table_type == "assignment_staging") %>%
    filter(action %in% c(
      "candidate_add_after_source_check",
      "defer_vocabulary_or_taxonomy",
      "defer_insufficient_evidence",
      "reject_or_negative_evidence_only"
    )) %>%
    transmute(
      batch_id = spec$batch_id,
      disease_name,
      entity_type = if_else(is.na(entity_type) | entity_type == "", "host", entity_type),
      entity_name,
      role_assignment,
      assignment_confidence,
      action,
      review_reason
    ) %>%
    filter(!is.na(entity_name), entity_name != "") %>%
    distinct()

  if (nrow(candidates) == 0) {
    return(tibble())
  }

  batch_dir <- file.path(input_root, spec$batch_id)
  host_path <- file.path(batch_dir, "host_role_candidates_batch.csv")
  vector_path <- file.path(batch_dir, "vector_candidates_batch.csv")
  roster_path <- file.path(batch_dir, "species_host_vector_roster_batch.csv")

  candidates %>%
    rowwise() %>%
    mutate(
      exact_in_host_candidates = present_anywhere(host_path, entity_name),
      exact_in_vector_candidates = present_anywhere(vector_path, entity_name),
      exact_in_roster = present_anywhere(roster_path, entity_name),
      join_note = case_when(
        entity_type == "host" & isTRUE(exact_in_host_candidates) ~ "exact host candidate match",
        entity_type == "vector" & isTRUE(exact_in_vector_candidates) ~ "exact vector candidate match",
        isTRUE(exact_in_roster) ~ "exact roster match",
        str_detect(entity_name, "spp\\.|group|/|issue|arthropods|bats|rodents|livestock|primates|Aedes|Hyalomma|Ixodes|Pteropus") ~
          "group, complex, issue row, or non-exact taxon; needs manual join decision",
        TRUE ~ "no exact batch candidate/roster match found"
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------------------------|
#      Write per-batch reformatted outputs ------------------------------------|
# ------------------------------------------------------------------------------|
write_batch_outputs <- function(spec) {
  if (!file.exists(spec$source_path)) {
    stop("Missing source report: ", spec$source_path)
  }

  batch_dir <- file.path(input_root, spec$batch_id)
  output_dir <- file.path(batch_dir, "reformatted")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  lines <- readLines(spec$source_path, warn = FALSE, encoding = "UTF-8")
  tables <- extract_tables(lines, spec)
  issues <- issue_table(tables)
  presence <- candidate_presence_check(tables, spec)

  write_clean_markdown(lines, file.path(output_dir, "DEEP_RESEARCH_REFORMATTED.md"))
  write_csv(tables, file.path(output_dir, "extracted_all_tables.csv"), na = "")

  for (table_type in sort(unique(tables$table_type))) {
    out <- tables %>% filter(table_type == !!table_type)
    out <- drop_empty_columns(out)
    write_csv(out, file.path(output_dir, paste0("extracted_", table_type, ".csv")), na = "")
  }

  write_csv(issues, file.path(output_dir, "extraction_quality_issues.csv"), na = "")
  write_csv(presence, file.path(output_dir, "candidate_presence_check.csv"), na = "")

  summary <- tables %>%
    count(table_type, name = "rows") %>%
    arrange(table_type)

  action_summary <- tables %>%
    filter("action" %in% names(.), !is.na(action), action != "") %>%
    count(action, name = "rows") %>%
    arrange(desc(rows), action)

  readme_lines <- c(
    paste0("# Reformatted Deep Research Output: ", spec$batch_title),
    "",
    paste0("Source report: `", spec$source_path, "`"),
    "",
    "These files are staging artifacts only. They are not official role evidence or assignment tables.",
    "",
    "Extracted table counts:",
    if (nrow(summary) == 0) "- No tables extracted." else paste0("- `", summary$table_type, "`: ", summary$rows),
    "",
    "Action counts:",
    if (nrow(action_summary) == 0) "- No action rows found." else paste0("- `", action_summary$action, "`: ", action_summary$rows),
    "",
    "Known cleanup boundaries:",
    "- Deep Research citation and URL wrappers were stripped where possible.",
    "- Source URLs, DOI/PMID/PMCID fields, source-access labels, and exact evidence spans still require source checking before import.",
    "- Rows in `extraction_quality_issues.csv` flag vocabulary drift or values that need normalization before any official CSV update.",
    "- Rows in `candidate_presence_check.csv` show whether staged entities have exact local candidate/roster matches.",
    "",
    "Generated by:",
    "",
    "`Rscript scripts/associations/role_annotation/6_5_Reformat_Deep_Research_Reports.R`"
  )
  write_text(readme_lines, file.path(output_dir, "README.md"))

  tibble(
    batch_id = spec$batch_id,
    source_report = spec$source_path,
    output_dir = output_dir,
    table_rows = nrow(tables),
    quality_issue_rows = nrow(issues)
  )
}

# ------------------------------------------------------------------------------|
#      Run all batches and write manifest -------------------------------------|
# ------------------------------------------------------------------------------|
run_summary <- pmap_dfr(report_specs, function(source_path, batch_id, batch_title, phase) {
  write_batch_outputs(tibble(
    source_path = source_path,
    batch_id = batch_id,
    batch_title = batch_title,
    phase = phase
  ))
})

write_csv(
  run_summary,
  file.path(input_root, "deep_research_reformat_manifest.csv"),
  na = ""
)

message("Wrote reformatted Deep Research artifacts.")
print(run_summary)
