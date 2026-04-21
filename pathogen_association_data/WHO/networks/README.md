# WHO Network Outputs

This folder contains the core WHO-linked host-pathogen network tables plus the new derived disease-host-vector and pathogen-host-vector outputs.

The files fall into three groups:

-   base host-pathogen network tables
-   derived host-vector integration tables
-   QA tables for the host-vector join workflow

## Base Network Tables

### `clover_who_network.csv`

Bacteria-focused WHO host-pathogen network derived from the CLOVER workflow.

-   Grain: one row per pathogen-host association retained from the CLOVER branch
-   Main role: bacterial component of the WHO network before combination
-   Key fields: `Pathogen`, `PathogenTaxID`, `Disease_name`, `Host_clean`, `HostTaxID`, host taxonomy, pathogen taxonomy, `DetectionMethod`, `MainSource`

### `virion_who_network.csv`

Virus-focused WHO host-pathogen network derived from the VIRION workflow.

-   Grain: one row per pathogen-host association retained from the VIRION branch
-   Main role: viral component of the WHO network before combination
-   Key fields: `Pathogen`, `PathogenTaxID`, `Disease_name`, `Host_clean`, `HostTaxID`, host taxonomy, pathogen taxonomy, `DetectionMethod`, `MainSource`

### `combined_who_network.csv`

Main combined WHO host-pathogen network used downstream across the repo.

-   Grain: one row per pathogen-host association after combining CLOVER and VIRION outputs
-   Main role: canonical disease-host-pathogen table for downstream joins
-   Key fields: `Pathogen`, `PathogenTaxID`, `Disease_name`, `Host`, `HostTaxID`, host taxonomy, pathogen taxonomy, `DetectionMethod`, `MainSource`, `PathogenType`

This is the main input used when connecting host-pathogen information to vector evidence.

## Derived Host-Vector Integration Tables

### `disease_host_vector_links.csv`

Disease-level integrated table linking WHO disease-host records to disease-vector evidence and observational host-vector evidence.

-   Grain: one row per `disease + host + vector`
-   Main role: disease-level view of which WHO hosts are linked to which vectors for diseases with curated vector evidence
-   Built from:
    -   `combined_who_network.csv`
    -   `disease_vector_links_taxonomy_cleaned.csv`
    -   `vector_host_links_join_ready.csv`
-   Important fields:
    -   disease-host context: `disease_name`, `host`, `host_tax_id`
    -   disease-vector context: `vector_species`, `vector_group`, `best_evidence_level`, `best_evidence_basis`
    -   host-vector provenance: `vector_host_record_count`, `source_platform_examples`, `interaction_type_examples`, `country_examples`
    -   caution flags: `vector_taxon_rank`, `vector_species_needs_review`, `taxonomy_caution`

This table is best for disease-level ecological summaries where pathogen identity does not need to stay separate.

### `disease_host_vector_links_expanded.csv`

Expanded disease-level host-vector table for the screened disease subset.

-   Grain: one row per `disease + host + vector`
-   Main role: keeps all observed host-vector links for WHO hosts in the 12 screened diseases, then marks whether curated disease-vector evidence is present
-   Built from:
    -   `combined_who_network.csv`
    -   `vector_host_links_join_ready.csv`
    -   `disease_vector_links_taxonomy_cleaned.csv`
-   Important fields:
    -   host-vector side: `host_vector_evidence`, `vector_host_record_count`, `source_platform_examples`, `interaction_type_examples`, `country_examples`
    -   disease-vector side: `disease_vector_evidence`, `disease_vector_evidence_status`, `best_evidence_level`, `best_evidence_basis`
    -   summary status: `link_status`
        -   `confirmed_by_both`
        -   `host_vector_only_candidate`
    -   caution flags: `vector_taxon_rank`, `vector_species_needs_review`, `taxonomy_caution`

This table is broader than `disease_host_vector_links.csv`. It includes observed host-vector combinations even when the vector is not supported in the curated disease-vector table.

### `disease_host_vector_links_expanded_summary.csv`

Compact per-disease summary of the expanded disease-host-vector table.

-   Grain: one row per disease
-   Main role: quick overview of how many rows are confirmed by both evidence layers versus host-vector-only candidates
-   Key fields:
    -   `expanded_rows`
    -   `confirmed_by_both_rows`
    -   `host_vector_only_candidate_rows`
    -   `distinct_hosts`
    -   `distinct_vectors`
    -   `taxonomy_caution_rows`

### `pathogen_host_vector_links.csv`

Pathogen-level integrated table linking WHO pathogen-host rows to pathogen-vector evidence and observational host-vector evidence.

-   Grain: one row per `disease + pathogen + host + vector`
-   Main role: more specific downstream table when pathogen identity must remain explicit
-   Built from:
    -   `combined_who_network.csv`
    -   `pathogen_vector_links_filled.csv`
    -   `vector_host_links_join_ready.csv`
-   Important fields:
    -   pathogen context: `pathogen`, `pathogen_tax_id`, `pathogen_type`, `pathogen_family`, `pathogen_genus`
    -   pathogen-vector context: `vector_species`, `vector_group`, `evidence_strength`, `vector_evidence_basis`, `assignment_basis`
    -   host-vector provenance: `vector_host_record_count`, `source_platform_examples`, `interaction_type_examples`, `country_examples`
    -   caution flags: `vector_taxon_rank`, `vector_species_needs_review`, `taxonomy_caution`

Use this table when pathogen-level differences matter, for example in later filtering or SDM-linked pathogen summaries.

## Host-Vector Join QA Tables

### `host_vector_join_qa_summary.csv`

Compact metric table summarizing the current host-vector integration run.

-   Grain: one row per metric
-   Main role: quick status check for counts, overlap, and unresolved join issues

### `host_vector_join_disease_coverage.csv`

Disease-by-disease coverage table for the disease-vector to host-vector join.

-   Grain: one row per disease
-   Main role: shows how much of the disease-vector evidence currently has host-vector overlap
-   Key fields:
    -   `total_disease_vector_rows`
    -   `total_distinct_disease_vectors`
    -   `disease_vectors_with_host_overlap`
    -   `final_disease_host_vector_rows`
    -   `final_distinct_hosts`
    -   `final_distinct_vectors`
    -   `taxonomy_caution_rows`

### `host_vector_join_missing_host_tax_id.csv`

Blocked host-vector evidence rows that could not be used in the join-ready host-vector table because a required host taxid was missing.

-   Grain: record-level host-vector evidence rows
-   Main role: review queue for host-matching problems in the host-vector pipeline
-   Important field: `block_reason`

If this file is empty, the current host-vector join table has no host-taxid blocked rows.

### `host_vector_join_unmatched_disease_vectors.csv`

Disease-vector rows from the pathogen-association pipeline that do not currently have any exact normalized vector-name overlap with the host-vector join table.

-   Grain: one row per unmatched disease-vector record
-   Main role: identifies vector taxa present in the disease-vector workflow but absent from current host-vector evidence

### `host_vector_join_unmatched_pathogen_vectors.csv`

Pathogen-vector rows from the pathogen-association pipeline that do not currently have any exact normalized vector-name overlap with the host-vector join table.

-   Grain: one row per unmatched pathogen-vector record
-   Main role: identifies pathogen-vector assignments that cannot yet be connected to host-vector observations

### `host_vector_join_taxonomy_caution_rows.csv`

Combined review table of integrated rows flagged with taxonomy caution.

-   Grain: mixed output table containing flagged rows from both disease-level and pathogen-level integrated outputs
-   Main role: review queue for rows involving non-species vector ranks, review-needed vector labels, or other taxonomy caveats
-   Key field: `output_level` indicates whether the row came from the disease-level or pathogen-level output

## Practical Use

If you want the main WHO host-pathogen backbone, use `combined_who_network.csv`.

If you want disease-level vector integration, use `disease_host_vector_links.csv`.

If you want the broader screened-disease table that includes host-vector-only candidates, use `disease_host_vector_links_expanded.csv`.

If you want pathogen-level vector integration, use `pathogen_host_vector_links.csv`.

If you want to understand what is still unresolved in the join workflow, start with:

-   `host_vector_join_qa_summary.csv`
-   `host_vector_join_disease_coverage.csv`
-   `host_vector_join_unmatched_disease_vectors.csv`
-   `host_vector_join_unmatched_pathogen_vectors.csv`
-   `host_vector_join_taxonomy_caution_rows.csv`
