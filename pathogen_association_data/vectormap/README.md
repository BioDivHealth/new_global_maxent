# VectorMap Files

This folder contains raw VectorMap downloads plus the staged host-filtering outputs used to build WHO-host-focused host-vector associations.

## Raw downloads

-   `BloodMealMap_Layer_-3496204453665016601.csv`: Raw mosquito blood-meal table with vector identity, host identity, pathogen context, and collection metadata.
-   `FleaMap_-6875364799851429947.csv`: Raw flea occurrence table with host-associated collection records and parasite-testing metadata.
-   `HostMap_9137669084057255600.csv`: Raw host occurrence/context layer from VectorMap, kept as a reference layer rather than used directly for host-vector links.
-   `MidgeMap_1481758679890687542.csv`: Raw biting-midge occurrence layer without direct host-link fields.
-   `MiteMap_-6324776740768397246.csv`: Raw mite occurrence table with host-associated collection records and parasite-testing metadata.
-   `MosquitoMap2_2627680870621077260.csv`: Raw mosquito occurrence layer without direct host-link fields.
-   `Sand_Fly_Map_-7115205178788551694.csv`: Raw sand fly occurrence layer without direct host-link fields.
-   `TickMap_4464597498443279194.csv`: Raw tick occurrence table with host-associated collection records and parasite-testing metadata.

## Manual inputs

-   `manual/vectormap_host_manual_crosswalk.csv`: Reviewed one-to-one host mapping file used to force approved scientific-name matches into the WHO-host-filtered output.

## Generated outputs

-   `outputs/vectormap_vector_host_links_raw.csv`: Combined raw direct-evidence host-vector table built from BloodMealMap, TickMap, FleaMap, and MiteMap.
-   `outputs/vectormap_vector_host_links_who_exact.csv`: WHO-host-matched subset containing only exact scientific-binomial matches before any manual crosswalk is applied.
-   `outputs/vectormap_vector_host_links_who_filtered.csv`: Main species-level host-vector table filtered to hosts present in `combined_who_network.csv`.
-   `outputs/vectormap_host_crosswalk_review.csv`: Full unresolved host-label review table with canonicalization fields, review buckets, and candidate match hints.
-   `outputs/vectormap_host_manual_crosswalk_candidates.csv`: Ranked subset of unresolved scientific labels that are suitable candidates for manual review and inclusion in the crosswalk.
-   `outputs/vectormap_host_package_candidates.csv`: Ranked subset of unresolved labels worth checking with a taxonomy package after structural filtering.
