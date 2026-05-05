# Data Sources For `disease_host_vector_links_expanded.csv`

This note describes the full provenance chain behind `disease_host_vector_links_expanded.csv`.

That expanded table is written by `scripts/associations/host_vector_integration/5_9b_Derive_Disease_Host_Vector_Links_Expanded.R`. The same script also writes `disease_host_vector_links_expanded_summary.csv`, which is only a compact per-disease summary of this larger row-level table.

## What The Expanded Table Represents

`disease_host_vector_links_expanded.csv` contains one row per:

-   `disease + host + vector`

within the screened disease subset.

It keeps all observed host-vector combinations for WHO hosts in scope, then marks whether curated disease-vector evidence is also present.

Important status fields include:

-   `confirmed_by_both`
-   `host_vector_only_candidate`
-   `disease_vector_evidence`
-   `disease_vector_evidence_status`
-   `host_vector_evidence`
-   `taxonomy_caution`

It does not come from one raw source directly. It is derived by combining:

-   the WHO disease-host-pathogen backbone
-   curated disease-vector evidence
-   observational host-vector evidence

## Immediate Input Files Used By `5_9b`

The script reads these three files directly:

1.  `combined_who_network.csv` Path: `pathogen_association_data/WHO/networks/combined_who_network.csv` Role: canonical WHO disease-host-pathogen backbone.

2.  `disease_vector_links_taxonomy_cleaned.csv` Path: `pathogen_association_data/WHO/vector_screening/outputs/disease_vector_links_taxonomy_cleaned.csv` Role: cleaned disease-vector evidence table for the screened diseases.

3.  `vector_host_links_join_ready.csv` Path: `pathogen_association_data/vector_host/outputs/vector_host_links_join_ready.csv` Role: join-ready host-vector evidence table collapsed to one row per `host_tax_id + vector_join_key`.

## Source Group 1: WHO Disease-Host-Pathogen Backbone

### `combined_who_network.csv`

This is the main WHO-linked network used as the disease-host backbone for the expanded table.

-   Built by: `scripts/associations/network_building/4_CombineNetworks.R`
-   Path: `pathogen_association_data/WHO/networks/combined_who_network.csv`
-   Upstream components:
    -   `clover_who_network.csv`
    -   `virion_who_network.csv`

In `5_9b`, this table is collapsed to `disease + host` grain, while retaining summaries such as:

-   `pathogen_count_in_disease_host_network`
-   `pathogen_examples`
-   `detection_method_examples`
-   `main_source_examples`

This means the expanded output always starts from hosts already present in the WHO disease-host network.

## Source Group 2: Curated Disease-Vector Evidence

The disease-vector side does not come from VectorMap or MapVEu. It comes from the pathogen-association workflow that combines literature review evidence with EFSA-derived vector links.

### Direct disease-vector input used by `5_9b`

#### `disease_vector_links_taxonomy_cleaned.csv`

-   Built by: `scripts/associations/vector_screening/5_5b_Vector_Name_Cleanup.R`
-   Path: `pathogen_association_data/WHO/vector_screening/outputs/disease_vector_links_taxonomy_cleaned.csv`
-   Role in `5_9b`: provides the curated disease-vector table that is left-joined onto observed host-vector rows.

Important consequence:

-   if a host-vector combination is observed but its vector is absent from this curated disease-vector table, the row still appears in the expanded table as `host_vector_only_candidate`

### Upstream disease-vector staging

#### `disease_vector_links.csv`

-   Built by: `scripts/associations/vector_screening/5_5_Consolidate_Vector_Evidence.R`
-   Role: canonical disease-vector table before taxonomy cleanup
-   Input: `vector_table_with_efsa_standardized.csv`

#### `vector_table_with_efsa_standardized.csv`

-   Built by: `scripts/associations/vector_screening/5_4_Standardize_Filter_Vector_Table.R`
-   Role: standardized disease-vector evidence table used for canonical collapse

#### `vector_table_with_efsa.csv`

-   Built by: `scripts/associations/vector_screening/5_3_Combine_LitReview_EFSA_Vector_Table.R`
-   Role: merged disease-vector evidence table that appends EFSA rows to the curated literature-review table

### Raw or near-raw disease-vector evidence sources

#### `vector_table.xlsx`

-   Role: curated literature-review disease-vector table maintained in the repo
-   Used in: `5_3_Combine_LitReview_EFSA_Vector_Table.R`
-   Contribution: literature-review vector evidence rows

#### `pathogen_vector_links_efsa.csv`

-   Built by: `scripts/associations/vector_screening/5_2_EFSA_Vector_Crosswalk.R`
-   Path: `pathogen_association_data/WHO/vector_screening/efsa/outputs/pathogen_vector_links_efsa.csv`
-   Contribution: EFSA-derived pathogen/vector links crosswalked conservatively to WHO diseases/pathogens

So the disease-vector branch used in the expanded table is:

`vector_table.xlsx` plus `pathogen_vector_links_efsa.csv` -\> `vector_table_with_efsa.csv` -\> `vector_table_with_efsa_standardized.csv` -\> `disease_vector_links.csv` -\> `disease_vector_links_taxonomy_cleaned.csv`

## Source Group 3: Observational Host-Vector Evidence

The host-vector side used in the expanded table comes from a combined evidence workflow that merges VectorMap-derived and MapVEu-derived host-vector records.

### Direct host-vector input used by `5_9b`

#### `vector_host_links_join_ready.csv`

-   Built by: `scripts/associations/host_vector_integration/5_8_Prepare_Host_Vector_Join_Table.R`
-   Path: `pathogen_association_data/vector_host/outputs/vector_host_links_join_ready.csv`
-   Role in `5_9b`: provides the join-ready host-vector rows used to expand WHO disease-host records by observed vector use

This table is already collapsed to one row per:

-   `host_tax_id`
-   `vector_join_key`

It preserves host-vector provenance summaries such as:

-   `source_platform_examples`
-   `source_dataset_examples`
-   `interaction_type_examples`
-   `country_examples`
-   `record_count`

### Upstream host-vector staging

#### `vector_host_links_analysis_ready.csv`

-   Built by: `scripts/associations/host_vector_sources/5_13_Combined_Host_Vector_Evidence.R`
-   Path: `pathogen_association_data/vector_host/outputs/vector_host_links_analysis_ready.csv`
-   Role: combined record-level host-vector evidence table with source provenance preserved

#### `vector_host_links_analysis_summary.csv`

-   Built by: `scripts/associations/host_vector_sources/5_13_Combined_Host_Vector_Evidence.R`
-   Role: combined deduplicated host-vector summary table

### Two upstream evidence branches merged in `5_13`

#### VectorMap branch

Direct staged input:

-   `pathogen_association_data/vectormap/outputs/vectormap_vector_host_links_analysis_ready.csv`

This branch is built from VectorMap downloads and cleaning stages described in `pathogen_association_data/vectormap/README.md`.

Relevant raw VectorMap sources include:

-   `BloodMealMap_*.csv`
-   `TickMap_*.csv`
-   `FleaMap_*.csv`
-   `MiteMap_*.csv`

These provide direct host-associated vector records. In the current workflow, VectorMap/BloodMealMap is the first host-vector source branch.

#### MapVEu branch

Direct staged input:

-   `pathogen_association_data/mapveu/outputs/mapveu_vector_host_links_analysis_ready.csv`

This branch is built by:

-   `scripts/associations/host_vector_sources/5_11_MapVEu_Analysis_Ready.R`

It contributes blood-meal-based host-vector records from the MapVEu workflow and is merged with the VectorMap branch in `5_13_Combined_Host_Vector_Evidence.R`.

So the host-vector branch used in the expanded table is:

VectorMap staged analysis-ready records plus MapVEu staged analysis-ready records -\> `vector_host_links_analysis_ready.csv` -\> `vector_host_links_join_ready.csv`

## How These Sources Are Combined In `5_9b`

`5_9b_Derive_Disease_Host_Vector_Links_Expanded.R` does the following:

1.  Reads `combined_who_network.csv` and collapses it to `disease + host`.
2.  Restricts the workflow to the diseases present in `disease_vector_links_taxonomy_cleaned.csv`.
3.  Expands each WHO host by all matching observed host-vector rows from `vector_host_links_join_ready.csv`.
4.  Left-joins the curated disease-vector evidence from `disease_vector_links_taxonomy_cleaned.csv` using normalized disease names and normalized vector names.
5.  Writes `disease_host_vector_links_expanded.csv`.
6.  Summarises the expanded table to one row per disease and writes `disease_host_vector_links_expanded_summary.csv`.

Because the disease-vector join is a left join, the expanded table contains:

-   rows supported by both evidence layers: `confirmed_by_both`
-   rows supported by host-vector evidence only within the screened disease subset: `host_vector_only_candidate`

## Practical Interpretation

When reading `disease_host_vector_links_expanded.csv`, remember that each row depends on all three source groups:

-   WHO disease-host-pathogen structure determines which hosts are in scope
-   curated disease-vector evidence determines whether a row is `confirmed_by_both`
-   observational host-vector evidence determines which vectors are observed on those hosts at all

So the expanded table is not a simple vector table and not a simple VectorMap or MapVEu extract. It is a joined, screened, disease-host-vector integration output.

If you only want disease-level counts derived from this table, use `disease_host_vector_links_expanded_summary.csv`.

## Vector Competence Annotation Layer

The extracted vector competence evidence in `diseases/vector_competence.csv` is a newer mechanistic annotation layer. It comes largely from the same disease-specific literature base as the curated literature-review vector table, but it answers a different question:

-   the curated disease-vector branch asks whether a vector is implicated for a disease/pathogen
-   the competence layer asks what kind of evidence exists for infection, transmission, mixed results, or non-competence

This means vector competence evidence is not currently used as the row-inclusion gate for `disease_host_vector_links.csv` or `disease_host_vector_links_expanded.csv`. The row-inclusion gate remains the curated disease-vector table plus the host-vector evidence described above.

Instead, `scripts/associations/vector_screening/5_6c_Join_Vector_Competence_Evidence.R` collapses the competence evidence to disease-vector grain and writes competence-annotated companion files:

-   `pathogen_association_data/WHO/vector_screening/outputs/vector_competence_collapsed.csv`
-   `pathogen_association_data/WHO/vector_screening/outputs/vector_competence_join_unmatched.csv`
-   `pathogen_association_data/WHO/vector_screening/outputs/disease_vector_links_taxonomy_cleaned_competence_annotated.csv`
-   `pathogen_association_data/WHO/networks/disease_host_vector_links_competence_annotated.csv`
-   `pathogen_association_data/WHO/networks/disease_host_vector_links_expanded_competence_annotated.csv`

Important added fields include:

-   `vector_competence_status`
-   `competence_statuses`
-   `vector_competence_evidence_types`
-   `transmission_demonstrated`
-   `natural_infection_reported`
-   `vector_role_hint`
-   `uncertainty_reason`
-   `competence_source_examples`

These fields should be interpreted as evidence annotations, not final ecological role labels. They are useful for prioritising rows for bridge-vector, enzootic-vector, or role-candidate review, but they should not by themselves be treated as definitive proof of transmission-cycle role.
