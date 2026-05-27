# Severe Fever With Thrombocytopenia Syndrome Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Severe fever with thrombocytopenia syndrome (SFTS)
- Source pathogen or analysis unit: Bandavirus dabieense
- Host candidate rows: 0 in current local host-role candidate surface
- Vector candidate rows: 15
- Competence-linked vector rows: 9
- Local candidate snapshot: role roster, vector competence annotated table, and `diseases/SFTSV/` extraction files.

## Current Host Candidate Highlights

- The current generated role surface has vector rows but no host candidate rows for SFTS.
- Evidence rows were added for domesticated animal natural infection and human spillover/rare secondary transmission.
- Domesticated animal evidence was assigned only as `host_presence_only`, not reservoir or amplifier.

## Current Vector Candidate Highlights

- Local vector candidates include `Haemaphysalis longicornis`, other `Haemaphysalis`, `Ixodes`, `Amblyomma`, `Dermacentor`, `Hyalomma`, and `Rhipicephalus` ticks.
- UKHSA supports `Haemaphysalis longicornis` as primary vector.
- CDC EID evidence supports possible reservoir-vector role for `Haemaphysalis longicornis`, kept review-flagged.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| UKHSA SFTS guidance | Official public health guidance | https://www.gov.uk/guidance/severe-fever-with-thrombocytopaenia-syndrome-sfts-epidemiology-outbreaks-and-guidance/ | Yes | Supports primary `Haemaphysalis longicornis` vector and rare human-to-human caveat. |
| CDC EID SFTSV among domesticated animals China | Peer-reviewed article | https://wwwnc.cdc.gov/eid/article/19/5/12-0245_article | Yes | Supports natural infection across domesticated animals. |
| CDC EID Haemaphysalis longicornis reservoir/vector paper | Peer-reviewed article | https://wwwnc.cdc.gov/eid/article/21/10/pdfs/15-0126.pdf | Yes | Supports review-flagged reservoir/vector evidence. |
| Local SFTS extraction markdowns | Local curated extraction | `diseases/SFTSV/` | Background only | Used to identify deferred tick candidates. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| domesticated animals | host_presence_only | supports | medium | yes | Natural infection only; no clear reservoir/amplifier role assigned. |
| Homo sapiens | incidental_host | supports | medium | yes | Tick-borne zoonosis with rare human-to-human transmission. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Haemaphysalis longicornis | primary_vector | supports | high | no | Primary vector in official guidance. |
| Haemaphysalis longicornis | enzootic_maintenance_vector | supports | medium | yes | Reservoir-vector evidence remains review-flagged. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2.
- Vector evidence rows: 2.

## Draft Assignments Added

- Host assignments: 1 `host_presence_only` deferred row for domesticated animals.
- Vector assignments: 1 primary vector row.

## Deferred Candidates And Why

- Other tick species: deferred unless source review supports role beyond detection.
- Animal host groups: deferred because the current row is natural-infection evidence only.
- Human role assignment: evidence row added but not assigned because rare secondary transmission complicates role language.

## Open Questions For Collaborator Review

- Whether SFTS needs an explicit animal amplifier/reservoir review beyond the current local host-candidate surface.
