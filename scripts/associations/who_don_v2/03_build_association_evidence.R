library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_extraction.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_native_association.R"))

who_don_v2_ensure_dirs()
association_mode <- Sys.getenv("WHO_DON_V2_ASSOCIATION_MODE", unset = "native")
valid_association_modes <- c("native", "contract")
if (!association_mode %in% valid_association_modes) {
  stop(
    "WHO_DON_V2_ASSOCIATION_MODE must be one of: ",
    paste(valid_association_modes, collapse = ", "),
    call. = FALSE
  )
}

adoption_decisions <- v2_read_csv(
  who_don_v2_output_dir("review", "v2_disease_candidate_adoption_decisions.csv"),
  c("record_key", "disease_standard", "adoption_decision")
)
country_adoption_decisions <- v2_read_csv(
  who_don_v2_output_dir("review", "v2_country_candidate_adoption_decisions.csv"),
  c("record_key", "country_standard", "adoption_decision")
)

if (association_mode == "native") {
  evidence <- v2_build_native_association_evidence()
} else {
  association_contract <- v2_read_association_contract()
  evidence <- v2_association_evidence_from_reviewed_candidates(association_contract, adoption_decisions) %>%
    v2_apply_country_adoption_decisions(country_adoption_decisions)
}

v2_write_csv(evidence, who_don_v2_output_dir("evidence", "who_don_association_evidence.csv"))

records <- v2_read_records()
qa <- tibble::tibble(
  metric = c(
    "association_mode",
    "association_rows",
    "records_in_evidence",
    "records_without_association",
    "multi_country_records",
    "multi_disease_records",
    "native_reviewed_adoption_rows",
    "seeded_reviewed_adoption_rows",
    "native_country_adoption_rows",
    "legacy_country_exception_rows"
  ),
  value = c(
    association_mode,
    nrow(evidence),
    n_distinct(evidence$record_key),
    nrow(anti_join(records, evidence, by = c("record_id"))),
    evidence %>% count(record_key, country_standard) %>% count(record_key) %>% filter(n > 1) %>% nrow(),
    evidence %>% count(record_key, disease_standard) %>% count(record_key) %>% filter(n > 1) %>% nrow(),
    sum(evidence$source_method %in% c(
      "v2_native_reviewed_adoption",
      "v2_native_country_disease_reviewed_adoption"
    )),
    sum(grepl("seeded_reviewed_disease|seeded_reviewed_adoption", evidence$source_method)),
    sum(evidence$country_adoption_decision == "accept_native"),
    sum(evidence$country_adoption_decision == "accept_legacy_exception")
  )
)
v2_write_stage_diagnostic(qa, "v2_association_evidence_summary.csv")

message(
  "Wrote v2 association evidence (",
  association_mode,
  " mode): ",
  nrow(evidence)
)
