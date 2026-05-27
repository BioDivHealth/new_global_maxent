# Chikungunya Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Chikungunya fever
- Source pathogen or analysis unit: Alphavirus chikungunya
- Host candidate rows: 20
- Vector candidate rows: 21
- Competence-linked vector rows: 15
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, and `diseases/chikungunya/chikungunya_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, several non-human primates, bats, rodents, and one bird row.
- CDC Yellow Book supports humans and non-human primates as likely main amplifying reservoirs for mosquito infection.
- Only human assignment was added; non-human primates remain group-level evidence because local species-level candidates need separate source support.

## Current Vector Candidate Highlights

- Local vector candidates include confirmed `Aedes aegypti` and `Aedes albopictus`, plus multiple probable or candidate mosquito species.
- WHO source language supports `Aedes aegypti` and `Aedes albopictus` as the most common chikungunya vectors.
- Other local vector candidates remain deferred unless source-specific role language supports them.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Chikungunya fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/chikungunya | Yes | Supports `Aedes aegypti` and `Aedes albopictus` main vector rows. |
| CDC Transmission of Chikungunya Virus | Official public health guidance | https://www.cdc.gov/chikungunya/php/transmission/index.html | Background only | Confirms primary mosquito-borne transmission by `Aedes aegypti` and `Aedes albopictus`. |
| CDC Yellow Book: Chikungunya | Official travel medicine guidance | https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/chikungunya.html | Yes | Supports viremic humans and non-human primates as likely main amplifying reservoirs. |
| Local chikungunya competence extraction markdown | Local curated extraction | `diseases/chikungunya/chikungunya_vector_competence_extractions.md` | Background only | Used to keep competence-only and negative rows out of final role assignments. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Homo sapiens | amplifying_host | supports | high | no | CDC Yellow Book supports viremic humans as likely amplifying reservoirs. |
| non-human primates | reservoir_host_group | supports | medium | yes | CDC Yellow Book supports group-level primate amplification or reservoir context. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aedes aegypti | main_vector | supports | high | no | WHO identifies this as one of the most common vectors. |
| Aedes albopictus | main_vector | supports | high | no | WHO identifies this as one of the most common vectors. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2 rows added.
- Vector evidence rows: 2 rows added.

## Draft Assignments Added

- Host assignments: 1 human `amplifying_host` row added.
- Vector assignments: 2 rows added for `Aedes aegypti` and `Aedes albopictus`.

## Deferred Candidates And Why

- Non-human primate species: deferred because evidence is group-level.
- Bat, rodent, and bird candidates: deferred because local host presence is not role evidence.
- Other mosquito species: deferred unless source-backed role evidence distinguishes established, regional, or candidate roles.
- Negative or uncertain competence rows such as non-mosquito or weak mosquito evidence: not promoted to role assignments.

## Open Questions For Collaborator Review

- Whether group-level non-human primate evidence should become a group-level assignment after the Phase V vocabulary review.
- Whether local outbreak-region vectors such as `Aedes hensilli` should be promoted in a later regional subreview.
