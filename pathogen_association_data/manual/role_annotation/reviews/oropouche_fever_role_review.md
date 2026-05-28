# Oropouche Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Oropouche fever
- Source pathogen or analysis unit: Orthobunyavirus oropoucheense
- Host candidate rows: 6
- Vector candidate rows: 15
- Competence-linked vector rows: 15
- Local candidate snapshot: role roster, host candidates, vector competence annotated table, and `diseases/oropouche/*_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, primates, sloth, and genus-level primate rows.
- WHO supports sloths, non-human primates, and perhaps birds as sylvatic vertebrate hosts, but with uncertainty.
- PAHO source language supports humans as amplifying hosts in the transmission cycle; this is kept review-flagged.

## Current Vector Candidate Highlights

- Local vector candidates include `Culicoides paraensis`, `Culicoides sonorensis`, `Culex quinquefasciatus`, and other mosquito/midge candidates.
- WHO supports `Culicoides paraensis` as primary vector to humans.
- CDC HAN supports `Culex quinquefasciatus` only as a possible vector, so it is a review-flagged candidate.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Oropouche virus disease fact sheet | Official factsheet | https://www.who.int/news-room/fact-sheets/detail/oropouche-virus-disease | Yes | Supports primary `Culicoides paraensis` vector and uncertain sylvatic vertebrate host groups. |
| CDC Oropouche causes and spread | Official public health guidance | https://www.cdc.gov/oropouche/causes/index.html | Background only | Confirms primary biting midge transmission and possible sylvatic hosts. |
| CDC HAN Oropouche Virus Disease | Official health advisory | https://emergency.cdc.gov/han/2024/pdf/cdc_han_515.pdf | Yes | Supports possible `Culex quinquefasciatus` vector language. |
| PAHO Oropouche epidemiological update | Official regional update | https://www.paho.org/sites/default/files/2025-09/2024-mar-06-phe-update-oropouche-eng-final2.pdf | Yes | Supports human amplifying-host language. |
| Local Oropouche extraction markdowns | Local curated extraction | `diseases/oropouche/` | Background only | Used to keep mosquito/midge candidates conservative. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| sloths and non-human primates | reservoir_host_group | supports | medium | yes | WHO supports sylvatic vertebrate host group but uncertainty remains. |
| Homo sapiens | amplifying_host | supports | medium | yes | PAHO update supports human amplification; needs review. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Culicoides paraensis | primary_vector | supports | high | no | WHO identifies primary vector to humans. |
| Culex quinquefasciatus | candidate_vector | supports | medium | yes | Possible vector only. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2.
- Vector evidence rows: 2.

## Draft Assignments Added

- Host assignments: 2 review-flagged rows.
- Vector assignments: 2 rows.

## Deferred Candidates And Why

- Individual sloth and primate candidates: deferred because source evidence is group-level.
- `Culicoides sonorensis`: competence evidence exists locally but source-backed role assignment was not made in this official-source pass.
- Other mosquitoes and midges: deferred as candidate, possible, or weak evidence.

## Open Questions For Collaborator Review

- Whether Oropouche human amplification should be treated as a durable role or an outbreak-context finding.
- Whether `Culicoides sonorensis` should receive `competent_vector_only` after closer review.
