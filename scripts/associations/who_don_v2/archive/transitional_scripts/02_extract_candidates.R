library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_disease_rules.R"))

who_don_v2_ensure_dirs()
clean_final <- v2_read_clean_final()

country_candidates <- v2_country_candidates_from_clean(clean_final)
disease_candidates <- v2_disease_candidates_from_clean(clean_final)

v2_write_csv(country_candidates, who_don_v2_output_dir("candidates", "who_don_country_candidates.csv"))
v2_write_csv(disease_candidates, who_don_v2_output_dir("candidates", "who_don_disease_candidates.csv"))

country_aliases <- country_candidates %>%
  transmute(
    alias = country_raw,
    country_standard,
    alias_type = "clean_final_seed",
    is_ambiguous = FALSE,
    priority = 100L,
    notes = "Seeded from accepted clean final country evidence."
  ) %>%
  filter(!is.na(alias), alias != "", !is.na(country_standard), country_standard != "") %>%
  distinct()

disease_aliases <- disease_candidates %>%
  v2_safe_disease_aliases()

disease_resolution_seed <- v2_disease_resolution_seed_from_clean(clean_final)
disease_rule_validation <- v2_validate_disease_aliases(disease_aliases)

v2_write_csv(country_aliases, who_don_v2_rules_dir("country_aliases.csv"))
v2_write_csv(disease_aliases, who_don_v2_rules_dir("disease_aliases.csv"))
v2_write_csv(disease_resolution_seed, who_don_v2_rules_dir("disease_resolution_seed_from_clean.csv"))
v2_write_csv(disease_rule_validation, who_don_v2_output_dir("qa", "v2_disease_rule_validation.csv"))

candidate_qa <- tibble::tibble(
  metric = c("country_candidate_rows", "disease_candidate_rows", "records_with_country", "records_with_disease"),
  value = c(
    nrow(country_candidates),
    nrow(disease_candidates),
    dplyr::n_distinct(country_candidates$record_key),
    dplyr::n_distinct(disease_candidates$record_key)
  )
)
v2_write_csv(candidate_qa, who_don_v2_output_dir("qa", "v2_candidate_summary.csv"))

message("Wrote v2 candidates: ", nrow(country_candidates), " country, ", nrow(disease_candidates), " disease")
