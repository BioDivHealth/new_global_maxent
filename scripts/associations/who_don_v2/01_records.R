library(dplyr)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))

who_don_v2_ensure_dirs()
v2_seed_rule_tables()

records <- v2_read_records_source()
association_contract <- v2_read_association_contract()

records <- v2_materialize_records(records, association_contract)

v2_write_csv(records, who_don_v2_output_dir("records", "who_don_records_clean.csv"))
v2_write_stage_diagnostic(
  records %>%
    count(record_key_source, name = "records"),
  "v2_record_key_repair_summary.csv"
)
v2_write_csv(
  records %>%
    select(record_key, record_key_source, source_record_id, DonId, record_id, Title, article_url),
  who_don_v2_output_dir("records", "who_don_record_key_map.csv")
)

message("Wrote v2 records: ", nrow(records))
