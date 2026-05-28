# Mpox Monkeypox Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Mpox (Monkeypox)
- Source pathogen or analysis unit: Orthopoxvirus monkeypox
- Host candidate rows: 20
- Vector candidate rows: 0
- Competence-linked vector rows: 0
- Local candidate snapshot: host-only candidate set with humans primates and several mammal candidates.

## Current Host Candidate Highlights

- `Homo sapiens` is present as a local candidate and is supported as a human outbreak/amplification host.
- The animal reservoir remains unknown in official summaries; small mammals are treated as susceptible/candidate hosts rather than assigned reservoirs.
- Non-human primates are susceptible hosts but are not assigned reservoir status from the reviewed official sources.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no current disease-vector rows are present for this Phase N disease.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Mpox fact sheet | official public health factsheet | https://www.who.int/news-room/fact-sheets/detail/mpox | yes | Used for unknown reservoir small-mammal susceptibility and human-to-human spread. |
| CDC Monkeypox in Animals and Pets | official public health guidance | https://www.cdc.gov/monkeypox/about/mpox-in-animals-and-pets.html | yes | Used for non-human primate and small-mammal susceptibility context. |
| CDC How Monkeypox Spreads | official public health guidance | https://www.cdc.gov/monkeypox/causes/index.html | supporting context | Used to cross-check animal-to-human transmission context. |
| Local host role candidates | local candidate table | `pathogen_association_data/WHO/role_annotation/host_role_candidates.csv` | yes | Used to verify local candidate presence and tax_id for `Homo sapiens`. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| small mammals | `susceptible_host_only` | supports | medium | yes | WHO says the natural reservoir is unknown but various small mammals are susceptible. |
| `Homo sapiens` | `amplifying_host` | supports | high | no | WHO states mpox spreads from person to person mainly through close contact. |
| non-human primates | `susceptible_host_only` | supports | medium | yes | CDC states non-human primates can get sick with monkeypox and have signs of disease like humans. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| not applicable | `not_applicable_non_vectored_scope` | not applicable | not applicable | no | Phase N host-only review; no vector rows added. |

## Rows Added To Evidence CSVs

- Host evidence rows: 3
- Vector evidence rows: 0

## Draft Assignments Added

- Host assignments: 1
- Vector assignments: 0

## Deferred Candidates And Why

- Small mammals and rodents: deferred from reservoir assignment because WHO explicitly says the natural reservoir remains unknown.
- Non-human primates: deferred from reservoir assignment because the reviewed CDC source supports susceptibility rather than reservoir status.
- Species-level animal candidates: deferred unless future sources support exact species-level role claims.

## Open Questions For Collaborator Review

- Whether to keep susceptible animal groups as evidence-only rows or add low-confidence group assignments after a dedicated animal-reservoir review.
