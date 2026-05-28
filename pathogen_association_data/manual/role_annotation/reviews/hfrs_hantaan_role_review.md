# Hemorrhagic Fever With Renal Syndrome Hantaan Virus Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Hemorrhagic fever with renal syndrome (Hantaan virus)
- Source pathogen or analysis unit: Orthohantavirus hantanense
- Host candidate rows: 21
- Vector candidate rows: 0
- Competence-linked vector rows: 0
- Local candidate snapshot: host-only disease in the current role review scope.

## Current Host Candidate Highlights

- `Apodemus agrarius` is present as a local candidate and has species-level Hantaan virus reservoir support.
- `Homo sapiens` is present as a local candidate and is treated as an incidental/spillover host.
- Other rodents and mammals in the candidate table remain deferred pending species-specific Hantaan virus role evidence.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no current disease-vector rows are present for this Phase N disease.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC Emerging Infectious Diseases HFRS in US Soldiers South Korea | peer-reviewed article | https://wwwnc.cdc.gov/eid/article/15/11/09-0076_article | yes | Used for `Apodemus agrarius` reservoir evidence. |
| CDC HFRS Clinician Brief | official clinical guidance | https://www.cdc.gov/hantavirus/hcp/clinical-overview/hfrs.html | yes | Used for human incidental/spillover role. |
| Local host role candidates | local candidate table | `pathogen_association_data/WHO/role_annotation/host_role_candidates.csv` | yes | Used to verify local candidate presence and tax_id. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| `Apodemus agrarius` | `reservoir_host` | supports | high | no | CDC EID article identifies the striped field mouse as the Hantaan virus reservoir host. |
| `Homo sapiens` | `incidental_host` | supports | high | no | CDC HFRS guidance supports rodent-associated human infection rather than maintenance in humans. |

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

- Other local rodent candidates: deferred because HFRS includes multiple hantaviruses and this review is scoped to Hantaan virus.

## Open Questions For Collaborator Review

- Whether future review should split HFRS host roles by additional hantavirus analysis units rather than carrying them under this Hantaan-specific disease row.
