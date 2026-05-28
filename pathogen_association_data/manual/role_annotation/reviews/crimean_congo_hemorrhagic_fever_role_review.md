# Crimean-Congo Hemorrhagic Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Crimean-Congo hemorrhagic fever
- Source pathogen or analysis unit: Orthonairovirus haemorrhagiae
- Host candidate rows: 13
- Vector candidate rows: 49
- Competence-linked vector rows: 30
- Local candidate snapshot: role roster, host candidates, vector competence annotated table, and `diseases/cchf/*_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, livestock-like ruminants, rodents, a hedgehog, hare, tortoise, and broad goat/sheep entries.
- Source-backed rows added here support livestock as group-level amplifying hosts and humans as spillover/incidental hosts.
- Species-level livestock assignments were deferred.

## Current Vector Candidate Highlights

- Local vector candidates are dominated by `Hyalomma` ticks, with additional `Rhipicephalus`, `Dermacentor`, `Haemaphysalis`, `Amblyomma`, `Ixodes`, and `Ornithodoros` candidates.
- WHO supports `Hyalomma` as the principal vector genus.
- Non-`Hyalomma` and individual `Hyalomma` species remain deferred until source-specific review.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO CCHF fact sheet | Official factsheet | https://www.who.int/en/news-room/fact-sheets/detail/crimean-congo-haemorrhagic-fever | Yes | Supports human spillover routes and `Hyalomma` principal-vector genus. |
| CDC EID CCHF Virus in Cattle and Ticks Israel | Peer-reviewed article | https://wwwnc.cdc.gov/eid/article/31/11/25-0622_article | Yes | Supports livestock amplifying-host language. |
| Local CCHF extraction markdowns | Local curated extraction | `diseases/cchf/` | Background only | Used to identify deferred vector-species candidates. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| livestock | amplifying_host | supports | high | yes | Group-level livestock amplification. |
| Homo sapiens | incidental_host | supports | high | no | Human spillover from ticks or infected animal tissues. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Hyalomma spp. | principal_vector_genus | supports | high | yes | Principal-vector genus, not automatic species assignment. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2.
- Vector evidence rows: 1.

## Draft Assignments Added

- Host assignments: 2.
- Vector assignments: 1 group-level row.

## Deferred Candidates And Why

- Individual livestock species: group-level source only.
- Individual tick species: deferred to species-level CCHF vector review.
- Non-livestock animal candidates: host presence is not role evidence.

## Open Questions For Collaborator Review

- Whether to add a distinct vocabulary value for tick-as-reservoir versus vector.
