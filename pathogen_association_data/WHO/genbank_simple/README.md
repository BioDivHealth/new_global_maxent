# GenBank-Simple Outputs

This folder contains the current GenBank-simple country-evidence outputs for
WHO disease modelling readiness.

## Main Files

- `genbank_readiness_query_overrides.csv`: manual query/taxid overrides used by
  the expanded readiness manifest.
- `genbank_simple_manifest.csv`: original 19-target GenBank-simple manifest,
  retained as a reference/control surface.
- `genbank_simple_readiness_manifest.csv`: one row per approved future/current
  retrieval target.
- `genbank_readiness_disease_country_summary_standardized.csv`: current
  disease-country evidence table used by downstream readiness scripts.

## Subfolders

- `qa/`: review and run-control tables, including manifest QA, search logs,
  target QA, and country-standardization QA.
- `intermediate/`: derived lower-level or unstandardized tables. The
  record-level country CSVs are intentionally local/ignored because they are
  bulky and can be regenerated from per-target retrieval checkpoints.
- `maps_readiness/`: readiness map CSV outputs and local generated PNG maps.
  The map PNGs are ignored; the compact map CSVs are suitable to keep when
  needed for review.

## Legacy Files

Flat `genbank_*` files without the `readiness` prefix are from the older
19-target GenBank-simple run. They are retained for comparison and fallback
support, but current modelling-readiness handoffs should use the readiness
outputs above.
