# Disease Modelling Readiness

This folder contains generated modelling-readiness handoff files. These outputs
are built from the WHO disease master list, role-annotation QA summaries,
species host/vector rosters, GenBank country summaries, WHO Disease Outbreak
News country evidence, and SDM availability manifests.

When the expanded GenBank-simple readiness run has been summarized and
standardized, the readiness build uses
`pathogen_association_data/evidence/genbank_simple/genbank_readiness_disease_country_summary_standardized.csv`.
That file combines the original 19-target GenBank-simple run with the expanded
readiness run. The older
`genbank_disease_country_summary_standardized.csv` is retained as a fallback for
historical standard-mode reruns only and is treated as local/archive material.

Regenerate from the repository root with:

```sh
Rscript scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R
```

## Files

- `disease_modelling_readiness.csv` is the lean planning table for all non-held
  analysis units.
- `disease_modelling_pilot.csv` is the WHO-focused pilot handoff subset from
  the same build.
- `disease_modelling_pilot_package/` is the generated pilot package folder. The
  pilot table is the spine, and companion CSVs expose host, vector, country,
  SDM-species, and evidence-summary layers keyed by `analysis_unit_id`.
  `pilot_sdm_species.csv` treats
  `sdms/outputs/catalog/accessible_sdm_species.csv` as the SDM availability
  source of truth.
- `disease_modelling_pilot_package.rds` and
  `disease_modelling_pilot_package.xlsx` are convenience versions of the same
  pilot package tables.
- `disease_modelling_readiness_full.csv` is the wider audit companion with join
  diagnostics, upstream source fields, and provenance columns retained for
  debugging the lean planning table.
- `disease_modelling_readiness_v1.csv` is a frozen snapshot of the earlier
  wider 53-column planning table, kept so the refined table can be compared
  against the previous layout.

These files are workflow control surfaces, not final biological evidence
sources. Direct vector evidence remains limited to curated vector rows, SDM
availability is name-matched availability only, and role assignments preserve
manual-review uncertainty.
