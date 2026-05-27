# Role Evidence Curation Progress

Last updated: 2026-05-08

## Scope

Current plan scope from `ROLE_EVIDENCE_FULL_CURATION_PLAN.md`: 21 diseases total.

- Phase V vectored diseases: 12
- Phase N non-vectored diseases: 9
- Broader readiness table diseases outside this scope: deferred expansion/reconciliation problem.

## Batch 0: Setup And Guardrails

Status: complete for this pass.

- Created `pathogen_association_data/WHO/role_annotation/reviews/`.
- Added `reviews/_role_review_template.md`.
- Extended `scripts/associations/role_annotation/6_3_Build_Evidence_Readiness_QA.R` to write:
  - `qa/role_evidence_claim_summary.csv`
  - `qa/role_assignment_claim_summary.csv`
- Added role curation conventions to `pathogen_association_data/WHO/role_annotation/README.md`.
- Reviewed West Nile pilot rows against the plan vocabulary.
- Normalized pilot role values:
  - `reservoir_amplifying_host_group` -> `reservoir_host_group`
  - `major_vector_amplificatory_bridge` -> `main_vector`

## V Batch 1: Highest-Value Rich-Source Vectored Diseases

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| West Nile fever | `west_nile_fever_role_review.md` | 3 | 3 | 2 | 2 | Existing pilot reviewed and vocabulary normalized. |
| Yellow fever | `yellow_fever_role_review.md` | 2 | 4 | 1 | 4 | Official-source cycle rows added; species-heavy rows deferred. |
| Dengue | `dengue_role_review.md` | 1 | 2 | 1 | 2 | `Aedes albopictus` kept as `candidate_vector` pending secondary-vector vocabulary review. |
| Plague | `plague_role_review.md` | 2 | 2 | 2 | 2 | Rodent/flea group evidence kept group-level; `Xenopsylla cheopis` added as species-level vector. |

## V Batch 2: Arbovirus Systems With Strong Vector Complexity

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| Rift Valley fever | `rift_valley_fever_role_review.md` | 2 | 4 | 2 | 4 | Ruminant livestock amplification, human spillover, and genus-level vector roles added. |
| Chikungunya fever | `chikungunya_fever_role_review.md` | 2 | 2 | 1 | 2 | Human amplification and main `Aedes` vectors added; primate group evidence kept review-flagged. |
| Zika virus disease | `zika_virus_disease_role_review.md` | 1 | 2 | 1 | 2 | Human urban-cycle row and `Aedes aegypti`/`Aedes albopictus` vector rows added. |
| Venezuelan equine encephalitis | `venezuelan_equine_encephalitis_role_review.md` | 3 | 2 | 2 | 2 | Rodent reservoir, equine amplification, and enzootic/epizootic vector rows added. |

## V Batch 3: Tick And Midge Systems

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| Crimean-Congo hemorrhagic fever | `crimean_congo_hemorrhagic_fever_role_review.md` | 2 | 1 | 2 | 1 | Livestock amplification, human spillover, and `Hyalomma` principal-vector genus added. |
| Tick-borne encephalitis | `tick_borne_encephalitis_role_review.md` | 2 | 3 | 1 | 3 | Host rows added despite no local host candidate rows; `Ixodes` subtype vectors added. |
| Severe fever with thrombocytopenia syndrome (SFTS) | `sfts_role_review.md` | 2 | 2 | 1 | 1 | Host evidence is conservative host-presence/spillover; `Haemaphysalis longicornis` primary vector added. |
| Oropouche fever | `oropouche_fever_role_review.md` | 2 | 2 | 2 | 2 | Sylvatic vertebrate hosts and human amplification review rows added; `Culicoides paraensis` primary vector added. |

## N Batch 1: High-Priority Viral Zoonoses With Strong Reservoir Literature

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| Ebola virus disease | `ebola_virus_disease_role_review.md` | 3 | 0 | 2 | 0 | African fruit bat reservoir group and human amplification rows added; primate spillover evidence kept review-flagged. |
| Sudan virus disease (Ebola virus disease) | `sudan_virus_disease_role_review.md` | 2 | 0 | 1 | 0 | Human amplification row added; bat reservoir evidence low-confidence and evidence-only. |
| Marburg virus disease | `marburg_virus_disease_role_review.md` | 2 | 0 | 2 | 0 | `Rousettus aegyptiacus` reservoir and human amplification rows added. |
| Nipah virus disease | `nipah_virus_disease_role_review.md` | 3 | 0 | 3 | 0 | `Pteropus` reservoir group, pig amplification, and context-dependent human amplification rows added. |

## N Batch 2: Rodent-Borne And Arenavirus Systems

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| Lassa fever | `lassa_fever_role_review.md` | 2 | 0 | 2 | 0 | `Mastomys natalensis` reservoir and context-dependent human amplification rows added. |
| Hemorrhagic fever with renal syndrome (Hantaan virus) | `hfrs_hantaan_role_review.md` | 2 | 0 | 2 | 0 | `Apodemus agrarius` reservoir and human incidental/spillover rows added. |
| Argentine hemorrhagic fever | `argentine_hemorrhagic_fever_role_review.md` | 2 | 0 | 2 | 0 | `Calomys musculinus` reservoir and human incidental/spillover rows added. |

## N Batch 3: High-Contact / Multi-Host Systems

Status: complete for this pass.

| Disease | Review markdown | Host evidence rows | Vector evidence rows | Host assignment rows | Vector assignment rows | Notes |
|---|---|---:|---:|---:|---:|---|
| Influenza (H5N1 avian influenza) | `influenza_h5n1_avian_influenza_role_review.md` | 4 | 0 | 4 | 0 | Wild aquatic bird group reservoir row plus poultry cattle and human rows added; species-heavy avian candidates deferred. |
| Mpox (Monkeypox) | `mpox_monkeypox_role_review.md` | 3 | 0 | 1 | 0 | Human amplification assignment added; animal reservoir kept unknown and susceptible animal groups evidence-only. |

## Validation

Latest validation command:

```sh
Rscript scripts/associations/role_annotation/6_3_Build_Evidence_Readiness_QA.R
```

Latest manual checks:

- All four manual CSVs parse with `readr::read_csv()`.
- All host role evidence and assignment values are in the plan vocabulary.
- All vector role evidence and assignment values are in the plan vocabulary.
- QA summaries show substantive evidence for all 21 reviewed diseases.
- Phase N diseases have zero vector evidence and assignment rows in the role tables.

Latest role table totals:

- `host_role_evidence.csv`: 47 rows
- `vector_role_evidence.csv`: 29 rows
- `host_role_assignments.csv`: 37 rows
- `vector_role_assignments.csv`: 27 rows

## Remaining Plan Work

Completed in this pass:

- Batch 0
- V Batch 1
- V Batch 2
- V Batch 3
- N Batch 1
- N Batch 2
- N Batch 3

Remaining:

- Final full-plan QA and completion audit
