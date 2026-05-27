# Rift Valley Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Rift Valley fever
- Source pathogen or analysis unit: Phlebovirus riftense
- Host candidate rows: 17
- Vector candidate rows: 78
- Competence-linked vector rows: 72
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, `diseases/rvf/rvf_vector_extractions.md`, and `diseases/rvf/rvf_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, cattle, sheep, goats, camels, buffalo, antelope, bats, and a rodent row.
- WHO source language supports group-level ruminant livestock amplification and human spillover/incidental infection.
- Species-level livestock assignments were not made because the source row used here is group-level and RVF host roles vary by animal species, age, pregnancy status, and outbreak context.

## Current Vector Candidate Highlights

- The local vector roster is broad, with many mosquito species and a few weak tick rows.
- WHO source language supports vertical maintenance in `Aedes` mosquitoes and ruminant amplification via local competent mosquitoes including `Culex`, `Mansonia`, and `Anopheles`.
- Species-heavy vector rows remain deferred because the official source used here is genus-level.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Rift Valley fever fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/rift-valley-fever | Yes | Supports livestock amplification, human spillover, `Aedes` vertical maintenance, and `Culex`/`Mansonia`/`Anopheles` mechanical-vector wording. |
| CDC About Rift Valley Fever | Official public health guidance | https://www.cdc.gov/rift-valley-fever/about/index.html | Background only | Confirms mosquito and infected animal tissue exposure pathways. |
| FAO Rift Valley fever page | Official veterinary source | https://www.fao.org/animal-health/animal-diseases/rift-valley-fever/en | Background only | Confirms livestock risk and mosquito vector genera. |
| Local RVF vector extraction markdown | Local curated extraction | `diseases/rvf/rvf_vector_extractions.md` | Background only | Used to identify candidate vector breadth and defer species promotion. |
| Local RVF competence extraction markdown | Local curated extraction | `diseases/rvf/rvf_vector_competence_extractions.md` | Background only | Used to avoid treating competence-only rows as final role assignments. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| ruminant livestock | amplifying_host | supports | high | yes | WHO supports amplification in naive ruminants. |
| Homo sapiens | incidental_host | supports | high | no | WHO supports human infection through animal tissues or mosquito bites; not maintenance. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aedes spp. | enzootic_maintenance_vector | supports | high | yes | WHO supports vertical maintenance in `Aedes` mosquitoes. |
| Culex spp. | mechanical_vector | supports | medium | yes | WHO lists `Culex` among local competent mosquitoes acting as mechanical vectors. |
| Mansonia spp. | mechanical_vector | supports | medium | yes | WHO lists `Mansonia` among local competent mosquitoes acting as mechanical vectors. |
| Anopheles spp. | mechanical_vector | supports | medium | yes | WHO lists `Anopheles` among local competent mosquitoes acting as mechanical vectors. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2 rows added for ruminant livestock and humans.
- Vector evidence rows: 4 rows added for `Aedes spp.`, `Culex spp.`, `Mansonia spp.`, and `Anopheles spp.`.

## Draft Assignments Added

- Host assignments: 2 rows added.
- Vector assignments: 4 group-level rows added.

## Deferred Candidates And Why

- Individual livestock species: deferred because the added source is group-level.
- Wildlife and bat candidates: deferred because the official-source pass did not support source-backed role rows.
- Individual mosquito species including `Aedes mcintoshi`, `Aedes vexans`, `Culex pipiens`, and `Culex poicilipes`: deferred until species-level source review.
- Tick rows: deferred because current evidence is weak and not part of the official-source vector role row.

## Open Questions For Collaborator Review

- Whether RVF should distinguish `enzootic_maintenance_vector` from `transovarial_maintenance_vector` in the vocabulary.
- Whether species-level RVF vector assignments should be promoted from the local extraction files in a separate, species-focused pass.
