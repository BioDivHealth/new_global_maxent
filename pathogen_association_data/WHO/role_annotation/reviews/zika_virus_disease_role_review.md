# Zika Virus Disease Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Zika virus disease
- Source pathogen or analysis unit: Orthoflavivirus zikaense
- Host candidate rows: 18
- Vector candidate rows: 56
- Competence-linked vector rows: 27
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, `diseases/zika/zika_vector_extractions.md`, and `diseases/zika/zika_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, non-human primates, bats, rodents, and other mammals.
- CDC source language supports humans in the urban Aedes-human-Aedes cycle.
- Non-human primate and other mammal rows remain deferred because this pass did not add source-backed sylvatic reservoir rows.

## Current Vector Candidate Highlights

- Local vector candidates include confirmed `Aedes aegypti`, `Aedes albopictus`, and multiple sylvatic or candidate `Aedes` species, plus many weaker non-Aedes rows.
- WHO source language supports `Aedes aegypti` as the main vector.
- CDC source language supports `Aedes albopictus` as part of urban transmission, but it remains marked for review because WHO frames `Aedes aegypti` as the main vector.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Zika virus fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/zika-virus | Yes | Supports `Aedes aegypti` as main vector. |
| CDC Transmission of Zika Virus | Official public health guidance | https://www.cdc.gov/zika/php/transmission/index.html | Yes | Supports urban transmission involving infected people and `Aedes aegypti`/`Aedes albopictus`. |
| Local Zika vector extraction markdown | Local curated extraction | `diseases/zika/zika_vector_extractions.md` | Background only | Used to identify breadth of candidate vectors. |
| Local Zika competence extraction markdown | Local curated extraction | `diseases/zika/zika_vector_competence_extractions.md` | Background only | Used to keep negative and competence-only rows out of final role assignments. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Homo sapiens | amplifying_host | supports | high | no | CDC supports infected people as sources for urban Aedes transmission. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aedes aegypti | primary_vector | supports | high | no | WHO identifies this as the main Zika vector. |
| Aedes albopictus | main_vector | supports | medium | yes | CDC includes it in urban transmission, but WHO frames `Ae. aegypti` as main. |

## Rows Added To Evidence CSVs

- Host evidence rows: 1 row added for human urban-cycle amplification.
- Vector evidence rows: 2 rows added for `Aedes aegypti` and `Aedes albopictus`.

## Draft Assignments Added

- Host assignments: 1 row added for human `amplifying_host`.
- Vector assignments: 2 rows added.

## Deferred Candidates And Why

- Non-human primates: deferred because sylvatic reservoir evidence was not added in this pass.
- Bat, rodent, and other mammal candidates: deferred because host presence is not role evidence.
- Sylvatic `Aedes` species and weaker non-Aedes vectors: deferred until source-specific review.
- Negative `Culex` evidence remains evidence for non-role or uncertainty, not a role assignment.

## Open Questions For Collaborator Review

- Whether Zika should get a separate sylvatic primate/reservoir group row after a close review of sylvatic ZIKV papers.
- Whether `Aedes albopictus` should stay `main_vector`, become `candidate_vector`, or wait for a future `secondary_vector` vocabulary.
