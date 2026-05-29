# GenBank-Simple Evidence

This folder contains the active GenBank-simple country-evidence outputs for WHO
disease modelling readiness.

## Main Files

- `genbank_readiness_disease_country_summary_standardized.csv`: current
  disease-country evidence table used by downstream readiness scripts.

## Related Folders

- `qa/`: review and run-control tables, including manifest QA, search logs,
  target QA, and country-standardization QA.
- `../../manual/genbank_simple/`: manual query/taxid overrides used by the
  expanded readiness manifest.
- `../../staged/genbank_simple/manifests/`: generated standard and readiness
  manifests.
- `../../staged/genbank_simple/intermediate/`: derived lower-level or
  unstandardized summary tables. Record-level country CSVs are intentionally
  local/ignored because they are bulky and can be regenerated from per-target
  retrieval checkpoints.
- `../../staged/genbank_simple/maps/readiness/`: compact readiness map CSV
  outputs. Generated PNG maps are ignored.
- `../../staged/genbank_simple/local_runs/`: ignored per-target retrieval
  checkpoints and logs.

## Legacy Files

Ignored legacy standard-mode files are kept, when needed locally, under
`../../archive/genbank_simple/legacy_19_target/`. Current modelling-readiness
handoffs should use the readiness evidence output above.
