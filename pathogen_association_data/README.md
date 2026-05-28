# Pathogen Association Data Layout

This folder is the data workspace for the pathogen-host-vector association
pipeline. Treat it as a mix of raw source exports, generated pipeline outputs,
manual review inputs, and local archive material.

Run scripts from the repository root so `here::here()` resolves paths
consistently. New scripts should source
`scripts/associations/working_inputs.R` and use its shared path constants before
adding new hard-coded `pathogen_association_data/...` roots.

## Active Pipeline Folders

- `WHO/`: Main WHO-centred pipeline output root. Active subfolders include
  `who_diseases/`, `networks/`, `vector_screening/`, `role_annotation/`,
  `genbank_simple/`, and `disease_outbreak_news_v2/`.
- `readiness/`: Generated modelling-readiness handoff files. These are planning
  and collaborator handoff surfaces, not final biological evidence claims.
- `vector_host/`: Combined VectorMap + MapVEu host-vector evidence outputs.
  Prefer `vector_host/outputs/` for downstream host-vector joins.

## Raw And Staged Source Folders

- `vectormap/`: Raw VectorMap downloads, manual crosswalks, and staged
  VectorMap-only host-vector outputs. Use these when changing VectorMap
  extraction or debugging VectorMap taxonomy/host filtering.
- `mapveu/`: Raw MapVEu exports, manual crosswalks, and staged MapVEu-only
  host-vector outputs. Use these when changing MapVEu extraction or debugging
  MapVEu taxonomy/host filtering.
- `viralemergence-clover-2604d22/`: Local CLOVER source checkout/vendor export.
  Generated WHO-specific CLOVER outputs live under `WHO/clover/`.

## Archive Or Local Comparison Material

- `*/outputs_v1/`: Legacy local comparison snapshots. They are not active
  inputs, not referenced by current scripts, and are ignored where relevant.
- `WHO/untitled folder/`: Unclassified local material. Do not use as a pipeline
  input until its contents are reviewed and moved to a named folder.
- Loose PDFs or dragged files at this level should be treated as local source
  material until explicitly documented by the relevant script or README.

## Shared Path Constants

`scripts/associations/working_inputs.R` defines the shared roots future scripts
should prefer:

- `pathogen_association_data_dir`
- `who_data_dir`
- `vectormap_dir`, `vectormap_outputs_dir`, `vectormap_manual_dir`
- `mapveu_dir`, `mapveu_outputs_dir`, `mapveu_manual_dir`
- `vector_host_dir`, `vector_host_outputs_dir`
- `clover_source_dir`, `who_clover_dir`
- `readiness_dir`

The same helper also keeps the current WHO working network/pathogen accessors,
including `who_working_network_path()` and `who_working_pathogens_path()`.
