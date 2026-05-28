# VectorMap Files

VectorMap is the pilot source family split across `source_data/`, `manual/`,
and `staged/`.

Raw exports live here under `source_data/vectormap/raw/`. Human-reviewed
crosswalks live under `manual/vectormap/`. VectorMap-only generated outputs live
under `staged/vectormap/outputs/`.

## Raw downloads

- `raw/BloodMealMap_Layer_-3496204453665016601.csv`: Raw mosquito
  blood-meal table with vector identity, host identity, pathogen context, and
  collection metadata.
- `raw/FleaMap_-6875364799851429947.csv`: Raw flea occurrence table with
  host-associated collection records and parasite-testing metadata.
- `raw/HostMap_9137669084057255600.csv`: Raw host occurrence/context layer from
  VectorMap, kept as a reference layer rather than used directly for
  host-vector links.
- `raw/MidgeMap_1481758679890687542.csv`: Raw biting-midge occurrence layer
  without direct host-link fields.
- `raw/MiteMap_-6324776740768397246.csv`: Raw mite occurrence table with
  host-associated collection records and parasite-testing metadata.
- `raw/MosquitoMap2_2627680870621077260.csv`: Raw mosquito occurrence layer
  without direct host-link fields.
- `raw/Sand_Fly_Map_-7115205178788551694.csv`: Raw sand fly occurrence layer
  without direct host-link fields.
- `raw/TickMap_4464597498443279194.csv`: Raw tick occurrence table with
  host-associated collection records and parasite-testing metadata.

## Manual inputs

- `manual/vectormap/vectormap_host_manual_crosswalk.csv`: Reviewed one-to-one
  host mapping file used to force approved scientific-name matches into the
  WHO-host-filtered output.
- `manual/vectormap/vectormap_vector_taxonomy_manual_map.csv`: Reviewed vector
  taxonomy mapping file used during VectorMap vector-name cleanup.

## Generated outputs

- `staged/vectormap/outputs/vectormap_vector_host_links_raw.csv`: Combined raw
  direct-evidence host-vector table built from BloodMealMap, TickMap, FleaMap,
  and MiteMap.
- `staged/vectormap/outputs/vectormap_vector_host_links_who_exact.csv`:
  WHO-host-matched subset containing only exact scientific-binomial matches
  before any manual crosswalk is applied.
- `staged/vectormap/outputs/vectormap_vector_host_links_who_filtered.csv`: Main
  species-level host-vector table filtered to hosts present in
  `combined_who_network.csv`.
- `staged/vectormap/outputs/vectormap_host_crosswalk_review.csv`: Full
  unresolved host-label review table with canonicalization fields, review
  buckets, and candidate match hints.
- `staged/vectormap/outputs/vectormap_host_manual_crosswalk_candidates.csv`:
  Ranked subset of unresolved scientific labels that are suitable candidates
  for manual review and inclusion in the crosswalk.
- `staged/vectormap/outputs/vectormap_host_package_candidates.csv`: Ranked
  subset of unresolved labels worth checking with a taxonomy package after
  structural filtering.
- `staged/vectormap/outputs/vectormap_host_taxize_review.csv`: Taxonomy-package
  review output for candidate host labels.
- `staged/vectormap/outputs/vectormap_host_taxize_who_hits.csv`: Taxonomy
  package matches that hit WHO host names.
- `staged/vectormap/outputs/vectormap_vector_host_links_who_vector_cleaned.csv`:
  WHO-host-filtered table after vector-name cleanup.
- `staged/vectormap/outputs/vectormap_vector_host_links_analysis_ready.csv`:
  VectorMap-only host-vector evidence surface ready for combined host-vector
  integration.
- `staged/vectormap/outputs/vectormap_vector_host_links_analysis_summary.csv`:
  Summary companion for the VectorMap-only analysis-ready output.
- `staged/vectormap/outputs/vectormap_vector_taxonomy_review_needed.csv`:
  Vector-name rows that need manual taxonomy review.
