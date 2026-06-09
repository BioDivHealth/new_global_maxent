# Chikungunya Fever Role Review

Phase: `Phase V`
Started: `2026-05-08`
Last updated: `2026-06-03`

## Disease Scope And Local Candidate Counts

- Disease name: Chikungunya fever
- Source pathogen or analysis unit: Alphavirus chikungunya
- Host candidate rows: 20
- Vector candidate rows: 21
- Competence-linked vector rows: 15
- Local candidate snapshot: `species_host_vector_roster.csv`, `host_role_candidates.csv`, `disease_vector_links_taxonomy_cleaned_competence_annotated.csv`, and `diseases/chikungunya/chikungunya_vector_competence_extractions.md`.

## Current Host Candidate Highlights

- Local host candidates include humans, several non-human primates, bats, rodents, and one bird row.
- CDC Yellow Book supports humans and non-human primates as likely main amplifying reservoirs for mosquito infection.
- Only human assignment was added; non-human primates remain group-level evidence because local species-level candidates need separate source support.

## Current Vector Candidate Highlights

- Local vector candidates include confirmed `Aedes aegypti` and `Aedes albopictus`, plus multiple probable or candidate mosquito species.
- WHO source language supports `Aedes aegypti` and `Aedes albopictus` as the most common chikungunya vectors.
- Other local vector candidates remain deferred unless source-specific role language supports them.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| WHO Chikungunya fact sheet | Official public health factsheet | https://www.who.int/en/news-room/fact-sheets/detail/chikungunya | Yes | Supports `Aedes aegypti` and `Aedes albopictus` main vector rows. |
| CDC Transmission of Chikungunya Virus | Official public health guidance | https://www.cdc.gov/chikungunya/php/transmission/index.html | Background only | Confirms primary mosquito-borne transmission by `Aedes aegypti` and `Aedes albopictus`. |
| CDC Yellow Book: Chikungunya | Official travel medicine guidance | https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/chikungunya.html | Yes | Supports viremic humans and non-human primates as likely main amplifying reservoirs. |
| Local chikungunya competence extraction markdown | Local curated extraction | `diseases/chikungunya/chikungunya_vector_competence_extractions.md` | Background only | Used to keep competence-only and negative rows out of final role assignments. |
| Althouse et al. 2018 | Primary field serology / transmission-dynamics study | `pathogen_association_data/source_data/role_annotation/papers/althouse2018_chikungunya_senegal_monkeys.pdf` | Pending CSV update | Strongest species-level support for African monkey amplification hosts in Senegal; argues monkeys alone do not maintain continuous circulation. |
| Eastwood et al. 2017 | Primary NHP serology study | `pathogen_association_data/source_data/role_annotation/papers/eastwood2017_chikungunya_kenya_primates.pdf` | Background or evidence-only | Supports Kenyan NHP CHIKV exposure and likely enzootic circulation; species do not directly match current Chikungunya roster assignments. |
| Evans et al. 2022 | Primary NHP serology study | `pathogen_association_data/source_data/role_annotation/papers/evans2022_chikungunya_myanmar_primates.pdf` | Background only | Supports Myanmar NHP exposure / possible sylvatic circulation; sampled macaque species do not directly support a current roster species assignment. |
| Vourc'h et al. 2014 | Primary vertebrate serology study | `pathogen_association_data/source_data/role_annotation/papers/vourch2014_chikungunya_indian_ocean_primates_rats.pdf` | Evidence-only candidate | Supports `Macaca fascicularis` antibody exposure only; no CHIKV RNA detected, so do not infer amplification or reservoir role. |
| Patouillat et al. 2024 | Systematic review | `pathogen_association_data/source_data/role_annotation/papers/patouillat2024_asian_primate_zoonotic_pathogens_review.pdf` | Background only | Useful for Asian primate surveillance gaps and context; not a direct host-role assignment source. |

## New NHP Source Check

Local PDF text has been extracted under `pathogen_association_data/source_data/role_annotation/pdf_text/` for the five primate-focused Chikungunya papers above.

| Source | Directly supported host evidence | Role interpretation for this review |
|---|---|---|
| Althouse et al. 2018 | `Chlorocebus sabaeus`, `Erythrocebus patas`, and `Papio papio` in Kedougou, Senegal had high CHIKV seropositivity and force-of-infection / reproductive-number support. | Use as the main source-backed basis for `amplifying_host` evidence for `Erythrocebus patas` and `Papio papio`. Do not assign reservoir status from this paper. `Chlorocebus sabaeus` is direct in the paper, but the local roster has `Chlorocebus aethiops`; any `Chlorocebus` row needs taxonomy review. |
| Eastwood et al. 2017 | Kenyan NHP sera show CHIKV-neutralizing antibodies, especially in western Kenya; 2014 positives included `Papio anubis`, `Cercopithecus mitis`, and `Cercopithecus ascanius`. | Supports group-level NHP exposure / enzootic-circulation context. Current positives do not map directly to the local Chikungunya host roster, so do not add species-level assignments from this source alone. |
| Evans et al. 2022 | Myanmar NHPs had CHIKV antibodies; sampled positives were `Macaca mulatta` and `Macaca nemestrina`, with no PCR-confirmed active infection. | Use as Asian NHP exposure context only. Do not assign `Macaca fascicularis` from this paper. |
| Vourc'h et al. 2014 | `Macaca fascicularis` and `Rattus rattus` showed antibody exposure after the Indian Ocean outbreak; CHIKV RNA was not detected in tested sera or rat organs. | Supports at most `host_presence_only` / exposure evidence for `Macaca fascicularis`. Do not infer amplification, reservoir, or maintenance role. Do not transfer rat evidence to `Mus musculus` or `Xerus erythropus`. |
| Patouillat et al. 2024 | Review notes that chikungunya is among viruses reported across multiple wild Asian primate species and emphasizes surveillance gaps. | Background context only; do not use as a direct evidence row unless paired with the underlying primary source. |

Current conservative implication: African primates can be treated as source-backed amplification hosts only where the species and geography are directly supported, with `Erythrocebus patas` and `Papio papio` the clearest additions. Reservoir or maintenance-host status remains unresolved in these papers.

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Homo sapiens | amplifying_host | supports | high | no | CDC Yellow Book supports viremic humans as likely amplifying reservoirs. |
| non-human primates | reservoir_host_group | supports | medium | yes | CDC Yellow Book supports group-level primate amplification or reservoir context. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| Aedes aegypti | main_vector | supports | high | no | WHO identifies this as one of the most common vectors. |
| Aedes albopictus | main_vector | supports | high | no | WHO identifies this as one of the most common vectors. |

## Rows Added To Evidence CSVs

- Host evidence rows: 2 rows added.
- Vector evidence rows: 2 rows added.

## Draft Assignments Added

- Host assignments: 1 human `amplifying_host` row added.
- Vector assignments: 2 rows added for `Aedes aegypti` and `Aedes albopictus`.

## Deferred Candidates And Why

- Non-human primate species: deferred because evidence is group-level.
- Bat, rodent, and bird candidates: deferred because local host presence is not role evidence.
- Other mosquito species: deferred unless source-backed role evidence distinguishes established, regional, or candidate roles.
- Negative or uncertain competence rows such as non-mosquito or weak mosquito evidence: not promoted to role assignments.

## Open Questions For Collaborator Review

- Whether group-level non-human primate evidence should become a group-level assignment after the Phase V vocabulary review.
- Whether local outbreak-region vectors such as `Aedes hensilli` should be promoted in a later regional subreview.
