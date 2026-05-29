# Pathogen Association Data Layout

This folder is the data workspace for the pathogen-host-vector association
pipeline. Treat it as a mix of raw source exports, generated pipeline outputs,
manual review inputs, and local archive material.

Run scripts from the repository root so `here::here()` resolves paths
consistently. New scripts should source
`scripts/associations/working_inputs.R` and use its shared path constants before
adding new hard-coded `pathogen_association_data/...` roots.

For remaining proposed splits, see `DATA_LAYOUT_PROPOSAL.md`. That file is
planning guidance, not the source of truth for already-moved folders.

## Active Pipeline Folders

- `WHO/`: Compatibility root for remaining WHO-centred pipeline outputs that
  have not yet been split. WHO diseases, WHO Disease Outbreak News v2, and WHO
  network data now use the split `source_data/`, `manual/`, `staged/`,
  `evidence/`, and `archive/` roots below.
- `readiness/`: Generated modelling-readiness handoff files. These are planning
  and collaborator handoff surfaces, not final biological evidence claims.
- `evidence/host_vector/`: Combined VectorMap + MapVEu host-vector evidence
  outputs. Prefer this integrated evidence surface for downstream host-vector
  joins.
- `evidence/role_annotation/`: Active role-annotation evidence, assignment,
  roster, and QA outputs. Manual review materials, generated prompt staging, and
  source PDFs/text are split out under `manual/`, `staged/`, and `source_data/`.
- `evidence/genbank_simple/`: Active GenBank-simple readiness disease-country
  evidence and QA outputs. Generated manifests, intermediate summaries, and map
  controls live under `staged/genbank_simple/`; manual query overrides live
  under `manual/genbank_simple/`.
- `evidence/vector_screening/`: Active disease/pathogen-vector evidence and
  vector-competence annotation outputs. QA companions, including competence
  unmatched-review files, live under `evidence/vector_screening/qa/`.
- `evidence/who_networks/`: Active WHO host-pathogen and WHO-only host-vector
  network evidence. Host-pathogen backbones live under
  `evidence/who_networks/host_pathogen/`; current WHO-only host-vector
  integrations live under `evidence/who_networks/host_vector/who_only/`; join
  QA lives under `evidence/who_networks/qa/`. Derived WHO vector-role candidate
  surfaces live under `evidence/role_annotation/`.
- `evidence/who_diseases/`: Active WHO diseases backbone, master expansion, host
  species, and QA evidence surfaces. See `WHO_DISEASES_DATA.md` for the
  full source/manual/staged/evidence split.
- `evidence/who_don_v2/`: Active WHO Disease Outbreak News v2 evidence
  outputs. `final/who_don_modelling_ready.csv` is the tracked downstream
  country-evidence surface; larger final audit tables, intermediate evidence,
  QA summaries, and web JSON exports are generated locally and ignored.

## Raw And Staged Source Folders

- `source_data/vectormap/`, `manual/vectormap/`, and `staged/vectormap/`:
  Split VectorMap source family. Raw exports live under
  `source_data/vectormap/raw/`, reviewed crosswalks under `manual/vectormap/`,
  and VectorMap-only generated outputs under `staged/vectormap/outputs/`.
- `source_data/mapveu/`, `manual/mapveu/`, and `staged/mapveu/`: Split MapVEu
  source family. Raw exports live under
  `source_data/mapveu/raw/`, reviewed crosswalks under `manual/mapveu/`, and
  MapVEu-only generated outputs under `staged/mapveu/outputs/`.
- `source_data/clover/`: Ignored raw/vendor CLOVER source checkout. Generated
  WHO-specific CLOVER outputs live under `staged/clover/outputs/`.
- `source_data/virion/raw/`: Ignored raw VIRION download material. Generated
  WHO-specific VIRION outputs live under `staged/virion/outputs/`.
- `source_data/role_annotation/`, `manual/role_annotation/`, and
  `staged/role_annotation/`: Split role-annotation source PDFs/OCR text, manual
  reviews/source checks, and generated Deep Research prompt/report staging.
- `staged/who_networks/`: Generated WHO network staging surfaces. Source-family
  CLOVER/VIRION network components live under
  `staged/who_networks/source_components/`; canonicalization support files live
  under `staged/who_networks/canonicalization/`.
- `source_data/vector_screening/`, `manual/vector_screening/`, and
  `staged/vector_screening/`: Split Vector Screening source family. EFSA raw
  workbooks live under `source_data/vector_screening/efsa/raw/`, manual
  screening/crosswalk/taxonomy decisions under `manual/vector_screening/`, and
  generated source-specific/intermediate outputs under
  `staged/vector_screening/`. VecTraits API probe outputs remain ignored under
  `staged/vector_screening/vectraits/` until promoted.
- `manual/who_don_v2/` and `staged/who_don_v2/`: Split WHO Disease Outbreak
  News v2 review inputs and generated staging layers. Durable review decisions
  live under `manual/who_don_v2/review/`; generated records, reference seeds,
  and candidate tables live under ignored `staged/who_don_v2/` subfolders.
- `source_data/who_diseases/`, `manual/who_diseases/`, and
  `staged/who_diseases/`: Split WHO diseases source family. Raw WHO regional
  tables live under `source_data/who_diseases/regional_tables/`, manual
  name-resolution/transmission/pathogen-matching/broad-taxa decisions live under
  `manual/who_diseases/`, and generated backbone/master-expansion/host-query/
  broad-taxa staging files live under `staged/who_diseases/`.

## Archive Or Local Comparison Material

- `archive/outputs_v1/`: Legacy local comparison snapshots moved out of active
  source roots. They are not active inputs and should not be referenced by
  current scripts.
- `archive/loose_files/`: Unclassified local material moved out of active data
  roots. Do not use these files as pipeline inputs until their contents are
  reviewed and moved to a named active folder.
- `archive/genbank_simple/legacy_19_target/`: Ignored local archive of the
  older standard-mode GenBank-simple outputs.
- `archive/vector_screening/`: Inactive Vector Screening snapshots retained for
  comparison only.
- `archive/who_don_v2/`: Ignored WHO Disease Outbreak News v2 archive and QA
  comparison material. Active scripts should not read from this archive except
  for explicitly documented historical QA or migration checks.
- Loose PDFs or dragged files found at this level should be moved under
  `archive/loose_files/` unless they are explicitly documented by the relevant
  script or README.

## Shared Path Constants

`scripts/associations/working_inputs.R` defines the shared roots future scripts
should prefer:

- `pathogen_association_data_dir`
- `who_data_dir`
- `source_data_dir`, `manual_data_dir`, `staged_data_dir`,
  `evidence_data_dir`
- `vectormap_raw_dir`, `vectormap_dir`, `vectormap_outputs_dir`,
  `vectormap_manual_dir`
- `mapveu_raw_dir`, `mapveu_dir`, `mapveu_outputs_dir`, `mapveu_manual_dir`
- `vector_host_dir`, `vector_host_outputs_dir`
- `clover_source_dir`, `virion_source_dir`, `virion_source_version_dir`,
  `who_clover_dir`, `who_virion_dir`
- GenBank-simple helpers for evidence, manual overrides, staged manifests,
  staged intermediates, staged maps, ignored local runs, QA, and legacy
  compatibility locations
- Vector Screening helpers for raw EFSA source workbooks, manual screening and
  taxonomy decisions, staged EFSA/intermediate outputs, active evidence, QA, and
  legacy compatibility locations
- role-annotation helpers for evidence, manual review/source-check, staged
  Deep Research, source PDF/text, roster, and QA locations
- WHO DON v2 helpers for generated records/reference/candidates, durable manual
  review decisions, active evidence/final/web/QA outputs, and ignored archives
- WHO diseases helpers for raw regional tables, manual name resolution and
  transmission rules, staged backbone/master-expansion/pathogen-matching/
  host-query/broad-taxa outputs, and active evidence/QA surfaces
- WHO network helpers for staged source components, canonicalization support,
  active host-pathogen evidence, current WHO-only host-vector evidence, join QA,
  and WHO-scoped vector-role candidate surfaces
- `readiness_dir`

The same helper also keeps the current WHO working network/pathogen accessors,
including `who_working_network_path()` and `who_working_pathogens_path()`.
