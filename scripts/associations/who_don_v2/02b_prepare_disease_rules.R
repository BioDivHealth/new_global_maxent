source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_disease_rules.R"))

prepared <- v2_prepare_disease_rules()

message(
  "Prepared safe disease rules: ",
  nrow(prepared$disease_aliases),
  " aliases; ",
  nrow(prepared$disease_rule_model),
  " unified rule-model rows; ",
  nrow(prepared$resolution_seed),
  " clean-resolution provenance rows; ",
  nrow(prepared$validation),
  " validation issues"
)
