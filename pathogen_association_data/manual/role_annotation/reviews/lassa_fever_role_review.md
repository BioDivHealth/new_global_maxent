# Lassa Fever Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Lassa fever
- Source pathogen or analysis unit: Mammarenavirus lassaense
- Host candidate rows: 14
- Vector candidate rows: 0
- Competence-linked vector rows: 0
- Local candidate snapshot: host-only disease in the current role review scope.

## Current Host Candidate Highlights

- `Mastomys natalensis` is present as a local host candidate and matches the official reservoir source at species level.
- `Homo sapiens` is present as a local host candidate; human-to-human transmission is documented but context-dependent.
- Other rodent candidates remain candidate-only unless source-backed role evidence supports species-level assignment.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no current disease-vector rows are present for this Phase N disease.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Lassa fever fact sheet | official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/lassa-fever | yes | Used for Mastomys reservoir and context-dependent human transmission rows. |
| Local host role candidates | local candidate table | `pathogen_association_data/WHO/role_annotation/host_role_candidates.csv` | yes | Used to verify local candidate presence and tax_id. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| `Mastomys natalensis` | `reservoir_host` | supports | high | no | WHO describes Mastomys rats as the main reservoir of Lassa virus. |
| `Homo sapiens` | `amplifying_host` | supports | medium | yes | WHO describes human-to-human transmission prevention in health-care settings; assignment remains context-dependent. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| not applicable | `not_applicable_non_vectored_scope` | not applicable | not applicable | no | Phase N host-only review; no vector rows added. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2
- Vector evidence rows: 0

## Draft Assignments Added

- Host assignments: 2
- Vector assignments: 0

## Deferred Candidates And Why

- Other local rodent candidates: deferred because the reviewed source supports Mastomys reservoir evidence and does not support broad species-level assignments to every rodent candidate.

## Open Questions For Collaborator Review

- Whether human Lassa fever should remain `amplifying_host` with manual review or use a more outbreak-context-specific final label.
