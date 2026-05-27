# Dengue Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Dengue
- Source pathogen or analysis unit: Orthoflavivirus denguei
- Host candidate rows: 37
- Vector candidate rows: 16
- Competence-linked vector rows: 15
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, `diseases/dengue/dengue_vector_extractions.md`, and `diseases/dengue/dengue_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, several non-human primates, bats, rodents, marsupials, carnivores, and tree shrews.
- Official-source evidence added in this pass supports humans as the amplifying host in the human-mosquito-human dengue transmission cycle.
- Non-human primate or other mammal candidates remain deferred because the current official row does not establish a global reservoir assignment.

## Current Vector Candidate Highlights

- Local vector candidates include confirmed `Aedes aegypti`, `Aedes albopictus`, several sylvatic `Aedes` species, and weaker `Culex quinquefasciatus` evidence.
- WHO source language supports `Aedes aegypti` as the primary vector and `Aedes albopictus` as a secondary-context vector.
- Because the current allowed vocabulary has no `secondary_vector`, `Aedes albopictus` is recorded as `candidate_vector` with manual review.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Dengue and severe dengue fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/dengue-and-severe-dengue | Yes | Supports human-to-mosquito transmission and primary/secondary Aedes vector wording. |
| CDC How Dengue Spreads | Official public health guidance | https://www.cdc.gov/dengue/transmission/index.html | Background only | Confirms spread by infected `Aedes` mosquitoes including `Ae. aegypti` and `Ae. albopictus`. |
| CDC Yellow Book: Dengue | Official travel medicine guidance | https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/dengue.html | Background only | Supports `Aedes aegypti` and `Aedes albopictus` vector context and human viremia cautions. |
| Local dengue vector extraction markdown | Local curated extraction | `diseases/dengue/dengue_vector_extractions.md` | Background only | Used to identify candidate complexity and defer secondary or regional species. |
| Local dengue competence extraction markdown | Local curated extraction | `diseases/dengue/dengue_vector_competence_extractions.md` | Background only | Used to avoid assigning competence-only or field-prevalence-only species as final roles. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Homo sapiens | amplifying_host | supports | high | no | WHO supports human-to-mosquito transmission during viremia. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aedes aegypti | primary_vector | supports | high | no | WHO describes this as the primary dengue vector. |
| Aedes albopictus | candidate_vector | supports | medium | yes | WHO says it can act as a vector but is normally secondary to `Ae. aegypti`. |

## Rows Added To Evidence CSVs

- Host evidence rows: 1 row added for human `amplifying_host`.
- Vector evidence rows: 2 rows added for `Aedes aegypti` and `Aedes albopictus`.

## Draft Assignments Added

- Host assignments: 1 row added for human `amplifying_host`.
- Vector assignments: 1 source-backed `primary_vector` row for `Aedes aegypti`, plus 1 `candidate_vector` row for `Aedes albopictus` marked `draft_needs_review`.

## Deferred Candidates And Why

- Non-human primates: deferred because no source-backed dengue reservoir assignment was added in this pass.
- Bat, rodent, marsupial, carnivore, and tree-shrew candidates: deferred because local host presence is not role evidence.
- Sylvatic `Aedes` vectors: deferred until source-specific sylvatic/regional role wording is reviewed.
- `Culex quinquefasciatus`: deferred as negative or weak field-detection evidence, not a role assignment.
- `Aedes mediovittatus`, `Aedes polynesiensis`, `Aedes malayensis`, and related regional vectors: deferred for later source-specific review rather than global assignment.

## Open Questions For Collaborator Review

- Whether the role vocabulary should add `secondary_vector`; until then, `Aedes albopictus` is represented as `candidate_vector` with a manual-review flag.
- Whether sylvatic dengue should receive separate group- or region-specific host/vector role rows after closer review of local extraction papers.
