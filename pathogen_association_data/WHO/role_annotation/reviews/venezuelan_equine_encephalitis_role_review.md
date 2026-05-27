# Venezuelan Equine Encephalitis Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Venezuelan equine encephalitis
- Source pathogen or analysis unit: Alphavirus venezuelan
- Host candidate rows: 31
- Vector candidate rows: 35
- Competence-linked vector rows: 31
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, and `diseases/vee/vee_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, domestic horse, rodents, birds, bats, marsupials, carnivores, and a primate row.
- Source-backed rows distinguish enzootic rodent reservoir evidence from epizootic equine amplification evidence.
- Human role was recorded as incidental and manual-review because VEE outbreak sources also note possible high human viremia in epidemic contexts.

## Current Vector Candidate Highlights

- Local vector candidates include enzootic `Culex` taxa, `Aedes taeniorhynchus`, `Psorophora confinnis`, and other mosquito candidates.
- Added vector rows distinguish `Culex Melanoconion spp.` as enzootic maintenance vectors from `Aedes taeniorhynchus` as an epizootic/epidemic vector.
- Additional species-level vectors remain deferred pending subtype and geography review.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC MMWR Update: Venezuelan Equine Encephalitis Colombia 1995 | Official outbreak report | https://www.cdc.gov/mmwr/preview/mmwrhtml/00039331.htm | Yes | Supports equines as most important vertebrate amplifying host in the epizootic cycle and local epizootic vector context. |
| CDC MMWR Venezuelan Equine Encephalitis Colombia 1995 | Official outbreak report | https://www.cdc.gov/mmwr/preview/mmwrhtml/00039070.htm | Yes | Supports human outbreak context and review caveat. |
| CDC Emerging Infectious Diseases: VEEV Infection of Cotton Rats | Peer-reviewed article | https://wwwnc.cdc.gov/eid/article/13/8/06-1157_article.htm | Yes | Supports enzootic circulation between `Culex (Melanoconion)` vectors and rodent reservoirs. |
| Venezuelan Equine Encephalitis - StatPearls | Peer-reviewed clinical review | https://www.ncbi.nlm.nih.gov/books/NBK559332/ | Yes | Supports sylvatic rodent reservoirs, equine amplification, and epizootic `Aedes taeniorhynchus` vector wording. |
| Local VEE competence extraction markdown | Local curated extraction | `diseases/vee/vee_vector_competence_extractions.md` | Background only | Used to avoid converting competence-only rows into role assignments. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| rodents | reservoir_host_group | supports | high | yes | Review source supports sylvatic rodents as primary reservoir hosts for enzootic strains. |
| Equus caballus | amplifying_host | supports | high | no | CDC MMWR supports equines as the most important vertebrate amplifying host in the epizootic cycle. |
| Homo sapiens | incidental_host | supports | medium | yes | Human cases occur in outbreaks, but possible high human viremia means this role needs review. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Culex Melanoconion spp. | enzootic_maintenance_vector | supports | high | yes | CDC EID supports enzootic circulation with rodent reservoirs. |
| Aedes taeniorhynchus | epidemic_vector | supports | high | yes | Review source supports epizootic primary vector wording; synonym and geography caveats remain. |

## Rows Added To Evidence CSVs

- Host evidence rows: 3 rows added.
- Vector evidence rows: 2 rows added.

## Draft Assignments Added

- Host assignments: 2 rows added for rodent group reservoir evidence and domestic horse amplification.
- Vector assignments: 2 rows added for `Culex Melanoconion spp.` and `Aedes taeniorhynchus`.

## Deferred Candidates And Why

- Individual rodent species: deferred because reservoir evidence is group-level and subtype/geography-specific.
- Human assignment: evidence row added but no assignment added because human epidemic viremia needs manual interpretation.
- Birds, bats, marsupials, carnivores, and primate candidate rows: deferred because local presence is not role evidence.
- Other mosquito species: deferred pending subtype, geography, and source wording review.

## Open Questions For Collaborator Review

- Whether `Aedes taeniorhynchus` should be stored under the local `Aedes` name or an `Ochlerotatus` synonym in final assignment rows.
- Whether humans should be `incidental_host`, `amplifying_host`, or evidence-only for epidemic VEE contexts.
