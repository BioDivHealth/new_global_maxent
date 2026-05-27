# WHO Role Annotation

This folder separates biological role review from the core WHO disease-pathogen-host and vector evidence tables.

The files here should not replace the canonical network, vector-screening, competence, GenBank, or WHO DON outputs. Those upstream files provide evidence layers. Role annotation is an interpretation layer that should preserve uncertainty and provenance.

## Files

- `host_role_candidates.csv` is a generated table seeded from `combined_who_network_canonical_zoonotic.csv` for the current active role scope.
- `host_role_candidates_summary.csv` is a generated QA summary for the host candidate table.
- `species_host_vector_roster.csv` is a generated collaborator-facing roster with one row per disease plus species or vector taxon. It covers both vectored and non-vectored diseases by combining host rows from the canonical WHO backbone with vector rows from the curated disease-vector table, plus host-vector observation and competence flags where available.
- `species_host_vector_roster_summary.csv` is a generated per-disease QA summary for the species roster.
- `species_host_vector_roster.xlsx` is a generated two-sheet workbook for sharing: `roster` contains the same rows as the CSV, and `column_descriptions` explains each column.
- `host_role_evidence.csv` is the manually curated evidence table for host role claims.
- `host_role_assignments.csv` is the reviewed or draft host-role assignment table.
- `vector_role_evidence.csv` is the manually curated evidence table for vector role claims.
- `vector_role_assignments.csv` is the reviewed or draft vector-role assignment table.
- `qa/` contains generated role-annotation QA summaries that make the host,
  vector, GenBank-simple, WHO DON, and role-review evidence layers easier to
  inspect before assigning final roles.
- Modelling-readiness handoff tables are generated under
  `pathogen_association_data/readiness/`, not in this role-annotation QA folder.
  That folder contains the lean all-disease readiness table, the WHO-focused
  pilot subset, the pilot package CSV/RDS/XLSX bundle, and the full audit
  companion with join diagnostics. Its SDM-species package table uses the
  broader accessible SDM inventory in `sdms/outputs/catalog/` as the SDM
  availability source of truth.

## Role-Review Scope

The current role-review scope is defined as:

```r
(in_gibb_etal == TRUE | in_empres_i == TRUE) &
  Pathogens != "Genus Vesiculovirus"
```

`priority_prototype_status` is retained as descriptive metadata rather than
used as an exclusion filter for role annotation. The GenBank-simple evidence
layer now has a readiness-combined output that binds the original 19-target run
with the expanded readiness run. Evidence-readiness and modelling-readiness
scripts should prefer
`genbank_simple/genbank_readiness_disease_country_summary_standardized.csv`
when it exists, falling back to the older 19-target
`genbank_disease_country_summary_standardized.csv` only for historical reruns.
Rows outside the current role-review scope can remain in the upstream
disease/pathogen tables for provenance and later review. They should not be
deleted just because they are deferred.

The collaborator-facing species roster also excludes broad genus/subgenus source rows with low taxonomic focus:

- `Genus Vesiculovirus`
- `Subgenus Merbecovirus`
- `Subgenus Sarbecovirus`

## Interpretation Boundary

`host_role_candidates.csv` is deliberately conservative. A candidate row means that a host appears in the canonical WHO disease-pathogen-host backbone for a role-review-scope disease. It does not mean the host is a reservoir, amplifier, incidental host, or dead-end host.

Final host or vector role labels should be assigned only after source-backed evidence is captured in the corresponding `*_role_evidence.csv` table.

## Manual Curation Conventions

Disease-level role reviews live in `reviews/`. Use `reviews/_role_review_template.md` for new disease reviews, and keep the markdown source log current before adding rows to the manual CSVs.

Group-level taxa can live in the evidence tables when the source only supports a group-level claim, such as `Aves`, `Rodentia`, `Culex spp.`, or `Aedes spp.`. Keep those rows marked for manual review when a downstream species-level interpretation would be tempting. Assignment tables may contain group-level rows only when the assignment is explicitly group-level; do not propagate group evidence to every member species.

For Phase N non-vectored disease reviews, vector findings should normally be recorded in markdown as `not_applicable_non_vectored_scope`. Do not add placeholder vector evidence or vector assignment rows solely to indicate that vectors were not part of the reviewed disease system.

## Regeneration

Regenerate host role candidates from the repository root with:

```sh
Rscript scripts/associations/role_annotation/6_1_Derive_Host_Role_Candidates.R
```

Regenerate the collaborator-facing host/vector species roster with:

```sh
Rscript scripts/associations/role_annotation/6_2_Derive_Species_Host_Vector_Roster.R
```

Regenerate evidence-readiness QA tables with:

```sh
Rscript scripts/associations/role_annotation/6_3_Build_Evidence_Readiness_QA.R
```

Regenerate the disease modelling readiness table with:

```sh
Rscript scripts/associations/role_annotation/6_10_Build_Disease_Modelling_Readiness.R
```
