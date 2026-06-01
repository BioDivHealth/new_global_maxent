# new_global_maxent

This repository is an R analysis workspace for zoonotic and vector-borne
disease modelling. It now contains three related but distinct work surfaces:

- legacy manuscript modelling assets from the original global MaxEnt analysis
- active pathogen-host-vector association workflows under `scripts/associations/`
- local SDM pilot and modelling-readiness work under `scripts/sdms/`,
  `sdms/`, and `pathogen_association_data/readiness/`

Run scripts from the repository root so `here::here()` resolves paths
consistently.

## Current Active Work

The actively maintained workflow is the pathogen association pipeline. It builds
WHO-linked disease/pathogen analysis units, attaches host, vector, competence,
country, role-review, and SDM-availability evidence layers, and writes
modelling-readiness handoff files.

Start with these files:

- `AGENTS.md`: repository working rules and current path conventions.
- `DATA_DECISIONS.md`: evidence-interpretation boundaries for host, vector,
  competence, country, role, and readiness layers.
- `pathogen_association_data/README.md`: current data lifecycle layout.
- `scripts/associations/README.md`: association workflow overview.
- `scripts/associations/working_inputs.R`: shared path helpers for active
  scripts.

Conceptual workflow order:

1. `scripts/associations/network_building/` builds WHO disease/pathogen
   backbones and host-pathogen networks from WHO, CLOVER, and VIRION sources.
2. `scripts/associations/vector_screening/` curates disease/pathogen-vector
   evidence and vector-competence annotations.
3. `scripts/associations/host_vector_sources/` prepares VectorMap and MapVEu
   host-vector evidence.
4. `scripts/associations/host_vector_integration/` joins disease/pathogen,
   host, vector, competence, and host-vector evidence for WHO-scoped outputs.
5. `scripts/associations/genbank_simple/` builds GenBank disease-country
   evidence for readiness workflows.
6. `scripts/associations/who_don_v2/` builds WHO Disease Outbreak News
   disease-country evidence.
7. `scripts/associations/role_annotation/` builds role-review candidates,
   source-check surfaces, QA summaries, and modelling-readiness handoffs.

## Data Layout

Versioned pathogen association data lives under `pathogen_association_data/`.
The current layout is lifecycle-based:

- `source_data/`: raw or near-raw source/vendor files.
- `manual/`: hand-edited curation, review, crosswalk, and control files.
- `staged/`: generated intermediates, prompts, manifests, and candidate tables.
- `evidence/`: active analysis-ready evidence outputs and QA surfaces.
- `readiness/`: generated modelling-readiness handoff files.
- `archive/`: inactive snapshots and local historical comparison material.

New scripts should source `scripts/associations/working_inputs.R` and use helper
functions instead of hard-coding data paths.

## Legacy Manuscript Modelling

The original manuscript work studied climate-change impacts on zoonotic and
vector-borne disease risk using host/vector ecological models, exposure, and
vulnerability layers. Some historical scripts, figures, and references remain in
the repository, but the top-level modelling order from the original manuscript
is not the current operating guide for this worktree.

The `data/` tree is sparse or ignored in normal local checkouts, and large
historical modelling inputs may not be present. Treat old modelling scripts and
point-data artifacts as legacy unless a current workflow explicitly references
them.

## Local SDM Pilot Work

The current SDM pilot surface is driven by readiness outputs, not by the legacy
manuscript script order. See:

- `pathogen_association_data/readiness/README.md`
- `scripts/sdms/README.md`

The readiness outputs are planning and collaborator handoff surfaces. They are
not final biological evidence claims.

## Validation Style

There is no `testthat` suite. Validate changes by running the smallest relevant
script or parse/smoke check and confirming outputs land in the current
`source_data/`, `manual/`, `staged/`, `evidence/`, or `readiness/` roots without
schema regressions.
