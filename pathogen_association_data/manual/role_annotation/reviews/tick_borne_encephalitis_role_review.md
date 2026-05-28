# Tick-Borne Encephalitis Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Tick-borne encephalitis
- Source pathogen or analysis unit: Orthoflavivirus encephalitidis
- Host candidate rows: 0 in current local host-role candidate surface
- Vector candidate rows: 11
- Competence-linked vector rows: 13
- Local candidate snapshot: role roster, vector competence annotated table, and `diseases/tbe/*_extractions.md`.

## Current Host Candidate Highlights

- The current generated role surface has vector rows but no host candidate rows for TBE.
- CDC source-backed rows were still added for small rodents as primary amplifying hosts and humans as spillover/incidental hosts, with review flags where local host candidates are absent.

## Current Vector Candidate Highlights

- Local vector candidates include `Ixodes ricinus`, `Ixodes persulcatus`, and several other tick species.
- CDC Yellow Book supports `Ixodes ricinus` and `Ixodes persulcatus` as primary subtype-associated vectors.
- CDC transmission guidance supports group-level `Ixodes` environmental maintenance with small rodents.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC Yellow Book TBE | Official travel medicine guidance | https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/tick-borne-encephalitis.html | Yes | Supports small rodents and primary `Ixodes` vectors. |
| CDC Transmission of TBE | Official transmission guidance | https://www.cdc.gov/tick-borne-encephalitis/php/transmission/index.html | Yes | Supports maintenance between `Ixodes` ticks and small rodents. |
| ECDC TBE factsheet | Official public health factsheet | https://www.ecdc.europa.eu/en/tick-borne-encephalitis/facts/factsheet | Background only | Confirms small rodent reservoir/amplifier and broader indicator hosts. |
| Local TBE extraction markdowns | Local curated extraction | `diseases/tbe/` | Background only | Used to identify deferred vector candidates. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| small rodents | amplifying_host | supports | high | yes | Group-level row; no local host candidate rows currently exist. |
| Homo sapiens | incidental_host | supports | high | yes | Human spillover row added despite no local host candidate row. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Ixodes ricinus | main_vector | supports | high | no | Primary vector for European subtype. |
| Ixodes persulcatus | main_vector | supports | high | no | Primary vector for Far Eastern and Siberian subtypes. |
| Ixodes spp. | enzootic_maintenance_vector | supports | high | yes | Group-level maintenance evidence. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2.
- Vector evidence rows: 3.

## Draft Assignments Added

- Host assignments: 1 group-level small rodent row.
- Vector assignments: 3 rows.

## Deferred Candidates And Why

- Other tick species: deferred unless source-backed role evidence supports more than detection or competence.
- Other vertebrate hosts: local role candidate rows absent; broader indicator hosts require separate review.

## Open Questions For Collaborator Review

- Why TBE has no local host candidate rows despite clear source-backed small-rodent host ecology.
