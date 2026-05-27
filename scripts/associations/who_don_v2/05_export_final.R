library(dplyr)
library(jsonlite)

source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_rules.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_final_shaping.R"))
source(here::here("scripts", "associations", "who_don_v2", "helpers", "who_don_v2_option_a_exceptions.R"))

who_don_v2_ensure_dirs()

association_mode <- Sys.getenv("WHO_DON_V2_ASSOCIATION_MODE", unset = "native")
evidence <- v2_read_csv(
  who_don_v2_output_dir("review", "who_don_review_decisions_applied.csv"),
  c("record_key", "country_standard", "disease_standard", "final_association_scope")
)

write_compatibility_exports <- identical(Sys.getenv("WHO_DON_V2_WRITE_COMPAT"), "1")

shaped <- v2_shape_final_outputs(evidence)
audit <- shaped$audit
if (association_mode == "native") {
  audit <- v2_apply_option_a_keep_current_exceptions_from_file(audit)
}
modelling <- shaped$modelling
if (association_mode == "native") {
  modelling <- v2_shape_modelling_ready(audit)
}

v2_write_csv(audit, who_don_v2_output_dir("final", "who_don_country_disease_scope_audit.csv"))
v2_write_csv(modelling, who_don_v2_output_dir("final", "who_don_modelling_ready.csv"))

compatibility_final_rows <- NA_integer_
compatibility_modelling_rows <- NA_integer_
compatibility_note <- "skipped; set WHO_DON_V2_WRITE_COMPAT=1 to refresh clean-shaped compatibility exports"

if (write_compatibility_exports) {
  clean_final <- v2_read_clean_final()
  clean_modelling <- v2_read_csv(v2_clean_modelling_path())

  # Compatibility exports intentionally keep the current clean schemas and values.
  v2_write_csv(
    clean_final,
    who_don_v2_output_dir("final", "who_don_country_disease_event_focal_scope_evidence_final.csv")
  )
  v2_write_csv(
    clean_modelling,
    who_don_v2_output_dir("final", "who_don_country_disease_event_focal_modelling_ready_final.csv")
  )

  compatibility_final_rows <- nrow(clean_final)
  compatibility_modelling_rows <- nrow(clean_modelling)
  compatibility_note <- "written from v2-local clean reference seeds"
}

final_qa <- tibble::tibble(
  metric = c(
    "final_scope_audit_rows",
    "modelling_ready_rows",
    "compatibility_final_rows",
    "compatibility_modelling_rows"
  ),
  value = c(
    nrow(audit),
    nrow(modelling),
    compatibility_final_rows,
    compatibility_modelling_rows
  ),
  note = c(
    paste("Final audit rows written in", association_mode, "association mode."),
    paste("Modelling-ready rows written in", association_mode, "association mode."),
    compatibility_note,
    compatibility_note
  )
)
v2_write_csv(final_qa, who_don_v2_output_dir("qa", "v2_final_export_summary.csv"))

message(
  "Wrote v2 final exports: ",
  nrow(audit),
  " audit rows, ",
  nrow(modelling),
  " modelling rows in ",
  association_mode,
  " mode; compatibility exports ",
  if (write_compatibility_exports) "written" else "skipped"
)
