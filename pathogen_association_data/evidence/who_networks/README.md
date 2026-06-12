# WHO Network Evidence

This folder contains active WHO-linked network evidence surfaces. Source-family
network components and canonicalization support files are split out under
`pathogen_association_data/staged/who_networks/`.

Use the helper functions in `scripts/associations/working_inputs.R` rather than
hard-coding these paths.

## Layout

### `host_pathogen/`

Active WHO and master-plus host-pathogen backbones:

- `combined_who_network.csv`: combined CLOVER + VIRION WHO host-pathogen
  backbone.
- `combined_who_network_canonical.csv`: canonical pathogen-name version of the
  combined WHO backbone.
- `combined_who_network_canonical_zoonotic.csv`: legacy WHO-only canonical
  zoonotic network retained for provenance, compatibility checks, and explicit
  legacy QA scripts.
- `master_plus_who_host_network.csv`: active combined master-plus + WHO host
  network. It includes `in_legacy_canonical_zoonotic_pathogen_host`, which
  marks rows belonging to the old canonical zoonotic contract.

Related helpers:

- `who_raw_network_path()`
- `who_canonical_network_path()`
- `who_canonical_zoonotic_network_path()`
- `who_network_host_pathogen_path(filename)`

`who_working_network_path()` still resolves to the legacy canonical zoonotic
table for backwards compatibility. New active consumers should read
`master_plus_who_host_network.csv` directly or use the compatibility helper in
`scripts/associations/network_building/helpers/master_plus_compatibility_helpers.R`.

### `host_vector/who_only/`

Current WHO-only host-vector integration outputs:

- `disease_host_vector_links.csv`: conservative disease-host-vector
  intersection of WHO disease-host evidence, curated disease-vector evidence,
  and observed host-vector evidence.
- `disease_host_vector_links_competence_annotated.csv`: conservative table plus
  vector-competence annotations.
- `disease_host_vector_links_expanded.csv`: broader screened-disease table that
  keeps observed host-vector rows and marks whether curated disease-vector
  evidence exists.
- `disease_host_vector_links_expanded_competence_annotated.csv`: expanded table
  plus vector-competence annotations.
- `disease_host_vector_links_expanded_summary.csv`: per-disease summary of the
  expanded table.
- `pathogen_host_vector_links.csv`: pathogen-level host-vector integration.

These files are scoped to the current WHO vector-screened disease set, not the
full master-plus/readiness disease list.

Related helper:

- `who_network_host_vector_path(filename)`

### `qa/`

Host-vector integration QA outputs:

- `host_vector_join_qa_summary.csv`
- `host_vector_join_disease_coverage.csv`
- `host_vector_join_missing_host_tax_id.csv`
- `host_vector_join_unmatched_disease_vectors.csv`
- `host_vector_join_unmatched_pathogen_vectors.csv`
- `host_vector_join_taxonomy_caution_rows.csv`

Related helper:

- `who_network_qa_path(filename)`

## Staged Companions

Source-family network components live under:

```text
pathogen_association_data/staged/who_networks/source_components/
  clover_who_network.csv
  virion_who_network.csv
```

Canonicalization support files live under:

```text
pathogen_association_data/staged/who_networks/canonicalization/
  combined_who_pathogen_canonical_lookup.csv
```

Related helpers:

- `who_network_source_component_path(filename)`
- `who_network_canonicalization_path(filename)`

Source lookup files used during network construction live under:

```text
pathogen_association_data/source_data/who_networks/domesticated/
  domesticated_lab_farmed.csv
```

Related helper:

- `who_network_domesticated_path()`

## Role-Candidate Companions

Derived WHO vector-role candidate summaries live with the role-annotation
evidence layer, not in this network folder:

```text
pathogen_association_data/evidence/role_annotation/vector_role_candidates_who.csv
pathogen_association_data/evidence/role_annotation/vector_role_candidates_who_summary.csv
```

These are generated from the WHO-only host-vector evidence layer. They are not a
full master-plus/readiness vector-candidate surface.

Related helpers:

- `role_vector_candidate_path("who")`
- `role_vector_candidate_summary_path("who")`

## Practical Use

Use `master_plus_who_host_network.csv` for active WHO/master-plus
host-pathogen work. When an old WHO-only contract is required, filter
`in_legacy_canonical_zoonotic_pathogen_host == TRUE` or use
`read_legacy_compatible_master_plus_network()`, which returns the old
22-column association-grain view.

Use `disease_host_vector_links.csv` for conservative disease-level
host-vector summaries.

Use `disease_host_vector_links_expanded.csv` when you want the broader
screened-disease table that includes host-vector-only candidate rows.

Use `pathogen_host_vector_links.csv` when pathogen identity needs to remain
explicit.

Use `qa/` outputs to inspect join coverage, missing host taxids, unmatched
vectors, and taxonomy caution rows before interpreting host-vector results.
