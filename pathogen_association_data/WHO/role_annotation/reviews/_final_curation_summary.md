# Role Evidence Full Curation Summary

Last updated: 2026-05-08

## Objective

Complete `ROLE_EVIDENCE_FULL_CURATION_PLAN.md` for the 21 scoped diseases:

- Batch 0 setup and guardrails.
- Phase V: 12 vectored diseases with host and vector role evidence review.
- Phase N: 9 non-vectored diseases with host role evidence review and explicit vector non-applicability.
- Batch-level progress recorded in markdown review files.
- No commit or push.

## Deliverables

| Requirement | Artifact or validation evidence | Status |
|---|---|---|
| Reviews directory exists | `pathogen_association_data/WHO/role_annotation/reviews/` | complete |
| Reusable disease template exists | `reviews/_role_review_template.md` | complete |
| Progress tracked in markdown | `reviews/_curation_progress.md` | complete |
| Final summary recorded | `reviews/_final_curation_summary.md` | complete |
| Batch 0 schema and guardrails documented | `pathogen_association_data/WHO/role_annotation/README.md` | complete |
| QA script summarizes evidence and assignment claims | `scripts/associations/role_annotation/6_3_Build_Evidence_Readiness_QA.R` writes `role_evidence_claim_summary.csv` and `role_assignment_claim_summary.csv` | complete |
| Phase V review markdowns | 12 `*_role_review.md` files for vectored diseases | complete |
| Phase N review markdowns | 9 `*_role_review.md` files for non-vectored diseases | complete |
| Source search log per disease | Each disease review has `## Sources Searched` | complete |
| Deferred candidates documented | Each disease review has `## Deferred Candidates And Why` | complete |
| Phase tags per disease | Each disease review has `Phase: Phase V` or `Phase: Phase N` | complete |
| Evidence rows parse | all four manual role CSVs parse with `readr::read_csv()` | complete |
| Role vocabulary valid | host and vector evidence/assignment values are in plan vocabulary | complete |
| Assignment rows are source-backed | zero host assignments and zero vector assignments lack matching evidence rows | complete |
| Phase V vector evidence | all 12 Phase V diseases have nonzero vector role evidence rows | complete |
| Phase N vector scope | all 9 Phase N diseases have zero vector evidence rows and markdown records `not_applicable_non_vectored_scope` | complete |

## Final QA Snapshot

Validation command:

```sh
Rscript scripts/associations/role_annotation/6_3_Build_Evidence_Readiness_QA.R
```

Final role table totals:

- `host_role_evidence.csv`: 47 rows
- `vector_role_evidence.csv`: 29 rows
- `host_role_assignments.csv`: 37 rows
- `vector_role_assignments.csv`: 27 rows

Target disease coverage:

- 21 of 21 scoped diseases have review markdowns.
- 21 of 21 scoped diseases have nonzero host role evidence.
- 12 of 12 Phase V diseases have nonzero vector role evidence.
- 9 of 9 Phase N diseases have zero vector role evidence and no placeholder vector assignment rows.

## Curation Boundaries

- Group-level official evidence was kept group-level and not propagated to every local species candidate.
- Species-heavy candidate sets were deferred unless source-backed role evidence supported direct assignment.
- Natural infection or host presence alone was not treated as reservoir or vector role evidence.
- Phase N vector findings were recorded in markdown as `not_applicable_non_vectored_scope` rather than as placeholder CSV rows.
- Several rows remain `needs_manual_review` where role wording is context-dependent or group-level.
