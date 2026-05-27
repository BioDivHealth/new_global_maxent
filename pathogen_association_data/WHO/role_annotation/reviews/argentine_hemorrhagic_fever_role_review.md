# Argentine Hemorrhagic Fever Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Argentine hemorrhagic fever
- Source pathogen or analysis unit: Mammarenavirus juninense
- Host candidate rows: 8
- Vector candidate rows: 0
- Competence-linked vector rows: 0
- Local candidate snapshot: host-only disease in the current role review scope.

## Current Host Candidate Highlights

- `Calomys musculinus` is present as a local candidate and has source-backed reservoir support.
- `Homo sapiens` is present as a local candidate and is treated as an incidental/spillover host because human-to-human transmission is rare.
- Other local rodent candidates remain deferred unless Junin-virus-specific role evidence supports them.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no current disease-vector rows are present for this Phase N disease.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| UKHSA Argentine haemorrhagic fever guidance | official public health guidance | https://www.gov.uk/guidance/argentine-haemorrhagic-fever-origins-reservoirs-transmission-and-guidelines/ | yes | Used for Calomys reservoir and human transmission caveat rows. |
| Local host role candidates | local candidate table | `pathogen_association_data/WHO/role_annotation/host_role_candidates.csv` | yes | Used to verify local candidate presence and tax_id. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| `Calomys musculinus` | `reservoir_host` | supports | high | no | UKHSA guidance identifies Calomys rodents as the animal reservoir and the local candidate table includes `Calomys musculinus`. |
| `Homo sapiens` | `incidental_host` | supports | high | no | UKHSA describes human-to-human transmission as very rare and possible through infected body fluids. |

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

- Other Calomys or rodent candidates: deferred because the current reviewed row supports `Calomys musculinus` and does not justify broad species-level assignment to all listed rodents.

## Open Questions For Collaborator Review

- Whether to add separate evidence-only rows for other Calomys species after a dedicated Junin virus ecology review.
