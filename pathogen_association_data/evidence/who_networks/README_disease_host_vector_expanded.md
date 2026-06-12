# Data Sources For `disease_host_vector_links_expanded.csv`

This note describes the provenance chain behind:

```text
pathogen_association_data/evidence/who_networks/host_vector/who_only/disease_host_vector_links_expanded.csv
```

The expanded table is written by
`scripts/associations/host_vector_integration/5_9b_Derive_Disease_Host_Vector_Links_Expanded.R`.
The same script writes the companion summary:

```text
pathogen_association_data/evidence/who_networks/host_vector/who_only/disease_host_vector_links_expanded_summary.csv
```

## What The Expanded Table Represents

`disease_host_vector_links_expanded.csv` contains one row per:

- `disease + host + vector`

within the current WHO vector-screened disease subset.

It keeps observed host-vector combinations for WHO hosts in scope, then marks
whether curated disease-vector evidence is also present.

Important status fields include:

- `confirmed_by_both`
- `host_vector_only_candidate`
- `disease_vector_evidence`
- `disease_vector_evidence_status`
- `host_vector_evidence`
- `taxonomy_caution`

This table is broader than `disease_host_vector_links.csv`. It is not a raw
VectorMap/MapVEu extract and it is not a final vector-role assignment table.

## Immediate Inputs

The script reads three active evidence surfaces through helpers in
`scripts/associations/working_inputs.R`:

1. `who_network_host_pathogen_path("master_plus_who_host_network.csv")`
   - current target:
     `pathogen_association_data/evidence/who_networks/host_pathogen/master_plus_who_host_network.csv`
   - role: WHO/master-plus disease-host-pathogen backbone filtered by
     `in_legacy_canonical_zoonotic_pathogen_host`

2. `vector_screening_evidence_path("disease_vector_links_taxonomy_cleaned.csv")`
   - current target:
     `pathogen_association_data/evidence/vector_screening/disease_vector_links_taxonomy_cleaned.csv`
   - role: cleaned disease-vector evidence table for the screened diseases

3. `vector_host_outputs_dir/vector_host_links_join_ready.csv`
   - current target:
     `pathogen_association_data/evidence/host_vector/vector_host_links_join_ready.csv`
   - role: join-ready observed host-vector evidence table

## Source Group 1: WHO Disease-Host-Pathogen Backbone

The WHO host-pathogen backbone is built from staged source-family components:

```text
pathogen_association_data/staged/who_networks/source_components/clover_who_network.csv
pathogen_association_data/staged/who_networks/source_components/virion_who_network.csv
```

`scripts/associations/network_building/stages/legacy_who_compatibility/4_CombineNetworks.R`
combines those into:

```text
pathogen_association_data/evidence/who_networks/host_pathogen/combined_who_network.csv
```

`scripts/associations/network_building/stages/legacy_who_compatibility/4_1_Canonicalize_Combined_WHO_Network.R`
then writes the legacy canonical compatibility networks. Active host-vector
integration now reads `master_plus_who_host_network.csv` and filters the
legacy-compatible rows before collapsing to `disease + host` grain.

In `5_9b`, the working network is collapsed to `disease + host` grain while
retaining summaries such as:

- `pathogen_count_in_disease_host_network`
- `pathogen_examples`
- `detection_method_examples`
- `main_source_examples`

## Source Group 2: Curated Disease-Vector Evidence

The disease-vector branch comes from the Vector Screening workflow, not from
VectorMap or MapVEu.

Direct input:

```text
pathogen_association_data/evidence/vector_screening/disease_vector_links_taxonomy_cleaned.csv
```

This table is built from literature-review and EFSA-derived vector evidence via
the `scripts/associations/vector_screening/` pipeline.

Important consequence:

- if a host-vector combination is observed but its vector is absent from this
  curated disease-vector table, the row can still appear in the expanded table
  as `host_vector_only_candidate`

## Source Group 3: Observational Host-Vector Evidence

The host-vector branch comes from the integrated VectorMap + MapVEu workflow.

Direct input:

```text
pathogen_association_data/evidence/host_vector/vector_host_links_join_ready.csv
```

This table is collapsed to one row per:

- `host_tax_id`
- `vector_join_key`

It preserves provenance summaries such as:

- `source_platform_examples`
- `source_dataset_examples`
- `interaction_type_examples`
- `country_examples`
- `record_count`

Upstream source-specific staged branches live under:

```text
pathogen_association_data/staged/vectormap/outputs/
pathogen_association_data/staged/mapveu/outputs/
```

## How The Sources Are Combined

`5_9b_Derive_Disease_Host_Vector_Links_Expanded.R` does the following:

1. Reads the WHO working network and collapses it to `disease + host`.
2. Restricts the workflow to diseases present in
   `disease_vector_links_taxonomy_cleaned.csv`.
3. Expands each WHO host by matching observed host-vector rows from
   `vector_host_links_join_ready.csv`.
4. Left-joins curated disease-vector evidence by normalized disease and vector
   names.
5. Writes `disease_host_vector_links_expanded.csv`.
6. Summarizes the expanded table to one row per disease.

Because the disease-vector join is a left join, the expanded table contains:

- rows supported by both evidence layers: `confirmed_by_both`
- rows supported by host-vector evidence only within the screened disease
  subset: `host_vector_only_candidate`

## Vector Competence Annotation Layer

Vector competence evidence is not the row-inclusion gate for
`disease_host_vector_links.csv` or `disease_host_vector_links_expanded.csv`.
The inclusion gate remains the curated disease-vector table plus observed
host-vector evidence.

`scripts/associations/vector_screening/5_6c_Join_Vector_Competence_Evidence.R`
adds competence annotations to companion files, including:

```text
pathogen_association_data/evidence/who_networks/host_vector/who_only/disease_host_vector_links_expanded_competence_annotated.csv
```

Important added fields include:

- `vector_competence_status`
- `competence_statuses`
- `vector_competence_evidence_types`
- `transmission_demonstrated`
- `natural_infection_reported`
- `vector_role_hint`
- `uncertainty_reason`
- `competence_source_examples`

These fields are evidence annotations, not final ecological role labels.
