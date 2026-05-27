# Plague Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Plague
- Source pathogen or analysis unit: Yersinia pestis
- Host candidate rows: 37
- Vector candidate rows: 47
- Competence-linked vector rows: 47
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, `diseases/plague/plague_vector_extractions.md`, and `diseases/plague/plague_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, many rodents, carnivores, lagomorphs, and some livestock-like or wild ungulate rows.
- Official-source evidence supports a group-level wild rodent-flea maintenance cycle.
- Human evidence supports an incidental/spillover role, with a separate pneumonic human-to-human transmission context that was not converted into a reservoir assignment.

## Current Vector Candidate Highlights

- Local vector candidates are flea-heavy, with 47 flea rows and many species-level competence/extraction rows.
- Official-source rows were added for fleas as a group and `Xenopsylla cheopis` as the named Oriental rat flea vector.
- The local extraction files contain many additional species-level flea rows, but most require source-specific geographic and efficiency caveats before assignment.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC How Plague Spreads | Official public health guidance | https://www.cdc.gov/plague/causes/index.html | Yes | Supports natural wild rodent-flea cycle and human infection routes. |
| CDC MMWR plague treatment guidelines | Official public health guidance | https://www.cdc.gov/mmwr/volumes/70/rr/rr7003a1.htm | Yes | Supports incidental host language and `Xenopsylla cheopis` vector evidence. |
| WHO Plague fact sheet | Official public health factsheet | https://www.who.int/news-room/fact-sheets/detail/plague | Background only | Supports broad flea, animal reservoir, and control context but less specific than CDC rows used here. |
| Local plague vector extraction markdown | Local curated extraction | `diseases/plague/plague_vector_extractions.md` | Background only | Used to identify species-level vector candidates and source caveats. |
| Local plague competence extraction markdown | Local curated extraction | `diseases/plague/plague_vector_competence_extractions.md` | Background only | Used to avoid treating all experimentally competent fleas as equivalent role assignments. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Rodentia | maintenance_host | supports | high | yes | CDC supports a natural maintenance cycle involving wild rodents and fleas. |
| Homo sapiens | incidental_host | supports | high | no | CDC supports spillover to incidental hosts including humans. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| fleas | main_vector | supports | high | yes | CDC supports fleas in the natural maintenance and transmission cycle. |
| Xenopsylla cheopis | primary_vector | supports | high | no | CDC MMWR explicitly names the Oriental rat flea as a plague vector. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2 rows added for rodent group maintenance and human incidental/spillover role.
- Vector evidence rows: 2 rows added for fleas as a group and `Xenopsylla cheopis`.

## Draft Assignments Added

- Host assignments: 2 rows added for group-level `Rodentia` maintenance and human `incidental_host`.
- Vector assignments: 2 rows added for group-level fleas and species-level `Xenopsylla cheopis`.

## Deferred Candidates And Why

- Individual rodent species: deferred because official-source evidence is group-level and plague host role varies by region and epizootic context.
- Carnivores, lagomorphs, ungulates, and other mammal candidates: deferred unless source-backed role evidence distinguishes susceptible, incidental, amplifying, or surveillance roles.
- Most flea species: deferred because competence or historical vector evidence needs geography, efficiency, host association, and source wording reviewed species by species.
- Mechanical-vector rows such as `Pulex irritans`: deferred until the vector-role assignment policy decides whether to assign `mechanical_vector` when biological-vector evidence is poor.

## Open Questions For Collaborator Review

- Whether `Rodentia` group-level maintenance assignment is acceptable now or should remain evidence-only until regional reservoir-host rules are defined.
- Whether flea species with strong local extraction evidence should be promoted in a later plague-specific subreview, rather than mixed into this official-source pass.
