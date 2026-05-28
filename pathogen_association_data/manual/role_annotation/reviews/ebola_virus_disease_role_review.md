# Ebola Virus Disease Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Ebola virus disease
- Source pathogen or analysis unit: Orthoebolavirus zairense
- Host candidate rows: 16
- Vector candidate rows: 0
- Local candidate snapshot: `host_role_candidates.csv` and `species_host_vector_roster.csv`.

## Current Host Candidate Highlights

- Local candidates include humans, fruit bats, non-human primates, rodents, and other mammals.
- CDC source language supports African fruit bats as likely involved in orthoebolavirus ecology.
- WHO source language supports non-human primates as infected source animals and humans as outbreak amplifying hosts.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no vector rows are present in the current role-review surface and no arthropod vector role was identified in this pass.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC How Ebola Disease Spreads | Official guidance | https://www.cdc.gov/ebola/causes/index.html | Yes | Supports likely African fruit bat source/reservoir group. |
| WHO Ebola disease fact sheet | Official factsheet | https://www.who.int/news-room/fact-sheets/detail/ebola-virus-disease | Yes | Supports animal spillover examples and human-to-human transmission. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| African fruit bats | reservoir_host_group | supports | medium | yes | Likely reservoir/source group. |
| non-human primates | spillover_host | supports | medium | yes | Source animals but not treated as natural reservoir. |
| Homo sapiens | amplifying_host | supports | high | no | Human-to-human outbreak transmission. |

## Source-Backed Vector Role Findings

No vector role evidence added; vector role evidence is not applicable in the current non-vectored scope.

## Rows Added To Evidence CSVs

- Host evidence rows: 3.
- Vector evidence rows: 0.

## Draft Assignments Added

- Host assignments: 2 rows for African fruit bats and humans.
- Vector assignments: 0.

## Deferred Candidates And Why

- Individual bat species: source is group-level.
- Individual non-human primate species: spillover/source role only and needs species-specific review.
- Rodents and other mammals: local host presence is not role evidence.

## Open Questions For Collaborator Review

- Whether Ebola should have a species-level bat assignment only after a dedicated reservoir literature pass.
