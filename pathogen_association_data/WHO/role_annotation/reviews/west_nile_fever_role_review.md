# West Nile Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: West Nile fever
- Source pathogen or analysis unit: Orthoflavivirus nilense
- Host candidate rows: 186
- Vector candidate rows: 75
- Competence-linked vector rows: 39
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, and QA outputs regenerated on 2026-05-08.

## Current Host Candidate Highlights

- The local host roster is bird-heavy: 140 bird candidate rows, including many passeriform, accipitriform, anseriform, strigiform, and charadriiform rows.
- The local roster also includes human, horse/domestic equid, livestock-like mammal, rodent, bat, primate, reptile, and amphibian candidate rows.
- Group-level bird reservoir evidence is valid for `Aves` but should not be propagated to every bird species without source-specific or review-backed support.
- Human and domestic horse/dead-end evidence is sufficiently explicit for draft source-backed assignments.

## Current Vector Candidate Highlights

- The local vector roster has 75 mosquito rows, including confirmed `Culex pipiens`, `Culex modestus`, `Culex tarsalis`, `Culex perexiguus`, `Culex annulirostris`, `Aedes albopictus`, and broader `Culex spp.` entries.
- Existing pilot rows support genus-level `Culex spp.` as `principal_vector_genus`, and `Culex pipiens` plus `Culex modestus` as `main_vector` with regional caveats.
- Many additional vector rows are candidate/probable or competence-layer rows only; these remain deferred unless a role source explicitly supports them.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO West Nile virus fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/west-nile-virus | Yes | Supports bird reservoir group, Culex principal-vector genus, and horse dead-end role. |
| CDC West Nile Virus Key Messages | Official public health guidance | https://www.cdc.gov/west-nile-virus/php/outbreak-communication/key-messages.html | Yes | Supports human, horse, and other mammal dead-end host language. |
| ECDC factsheet about West Nile virus infection | Official public health factsheet | https://www.ecdc.europa.eu/en/west-nile-fever/facts | Yes | Supports Europe-focused `Culex pipiens` and `Culex modestus` main vector rows. |
| ECDC Culex pipiens factsheet for experts | Official vector factsheet | https://www.ecdc.europa.eu/en/infectious-disease-topics/related-public-health-topics/disease-vectors/facts/mosquito-factsheets/culex-pipiens | Yes | Supports `Culex pipiens` as a major WNV vector; regional caveat retained. |
| CDC West Nile Virus Surveillance and Control Guidelines | Official public health guidance | https://www.cdc.gov/west-nile-virus/php/surveillance-and-control-guidelines/index.html | Not yet | Useful for later review of US enzootic/epidemic vector specificity. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aves | reservoir_host_group | supports | high | yes | WHO supports birds as reservoir hosts, but the claim is group-level. |
| Homo sapiens | dead_end_incidental_host | supports | high | no | CDC explicitly describes humans as dead-end hosts in WNV transmission logic. |
| Equus caballus | dead_end_host | supports | high | no | WHO supports horses as dead-end hosts; mapped to domestic horse candidate. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Culex spp. | principal_vector_genus | supports | high | yes | WHO supports genus-level principal vector language. |
| Culex pipiens | main_vector | supports | high | yes | ECDC supports a major/main vector role, with European geographic caveat. |
| Culex modestus | main_vector | supports | high | yes | ECDC supports a main vector role in Europe. |

## Rows Added To Evidence CSVs

- Host evidence rows: no new rows added in this pass; existing pilot rows were reviewed and one role claim was normalized from `reservoir_amplifying_host_group` to `reservoir_host_group`.
- Vector evidence rows: no new rows added in this pass; existing pilot `Culex pipiens` role claim was normalized from `major_vector_amplificatory_bridge` to `main_vector`.

## Draft Assignments Added

- Host assignments: no new rows added in this pass; existing human and domestic-horse draft assignments retained.
- Vector assignments: no new rows added in this pass; existing `Culex pipiens` assignment normalized to `main_vector`, and `Culex modestus` retained as `main_vector`.

## Deferred Candidates And Why

- Bird species-level assignments: deferred because current source-backed evidence is group-level.
- Additional `Culex` species and non-`Culex` mosquitoes: deferred unless source-backed role evidence distinguishes main, bridge, enzootic, epidemic, or candidate roles.
- Mammal, reptile, amphibian, and broad livestock candidate rows: deferred unless sources support role-specific claims beyond infection or host presence.

## Open Questions For Collaborator Review

- Whether `Aves` should receive a group-level assignment row now, or remain evidence-only until a later species/group assignment policy review.
- Whether regional vector assignments for Europe should be mirrored by separate North American source-backed rows for `Culex tarsalis`, `Culex quinquefasciatus`, and related vectors.
