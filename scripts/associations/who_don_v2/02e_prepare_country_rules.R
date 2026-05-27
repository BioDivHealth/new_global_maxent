source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_country_rules.R"))

who_don_v2_ensure_dirs()

aliases <- v2_prepare_country_aliases(force = TRUE)
policy <- v2_prepare_country_policy_decisions()
validation <- v2_validate_country_aliases(aliases)

v2_write_csv(validation, who_don_v2_output_dir("qa", "v2_country_rule_validation.csv"))
v2_write_stage_diagnostic(
  aliases %>%
    count(country_standard, name = "aliases"),
  "v2_country_alias_summary.csv"
)

if (nrow(validation) > 0 && any(validation$severity == "blocking")) {
  stop("Blocking country rule validation issues found.", call. = FALSE)
}

message("Prepared native country aliases: ", nrow(aliases), " rows")
message("Prepared native country policy decisions: ", nrow(policy), " rows")
