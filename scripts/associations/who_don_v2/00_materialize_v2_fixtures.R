source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))

who_don_v2_ensure_dirs()
v2_seed_rule_tables()
v2_materialize_reference_inputs()

association_contract <- v2_read_clean_final()
records_source <- v2_read_clean_records_seed()

v2_write_csv(association_contract, who_don_v2_rules_dir("accepted_association_contract.csv"))
v2_write_csv(records_source, who_don_v2_output_dir("records", "who_don_records_source.csv"))

message(
  "Materialized v2 fixtures: ",
  nrow(association_contract),
  " association rows; ",
  nrow(records_source),
  " records"
)
