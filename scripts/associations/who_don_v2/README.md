# WHO DON V2 Pipeline

This folder contains the active v2 WHO Disease Outbreak News pipeline. The
pre-v2 legacy and clean-migration script/data folders were archived under
`archive/who_don_pre_v2/` after v2 was accepted as the production surface.

V2 is organized around native candidate layers, explicit claim types, a
canonical country-disease-scope evidence table, and reviewed adoption/policy
layers:

Shared implementation modules live under `helpers/`; top-level numbered scripts
are the runnable production and optional QA entrypoints.

The default production command is Option A native association mode:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R
```

Option A builds routine association evidence from native country and native
disease adoption layers. The old accepted association contract is retained as
audit/reference material and as a rollback association mode:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --association-mode contract
```

Default production writes a minimal top-level `qa/` surface. To refresh detailed
stage diagnostics for debugging, opt in explicitly:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --verbose-qa
WHO_DON_V2_VERBOSE_QA=1 Rscript scripts/associations/who_don_v2/02d_compare_disease_candidates.R
```

Verbose diagnostics are written to `qa/archive/stage_diagnostics/`, not to the
top-level `qa/` folder.

The production runner executes:

0. `00_materialize_v2_fixtures.R` is not part of routine production; run it only
   when intentionally refreshing v2-owned fixtures from accepted clean outputs.
1. `01_records.R`
2. `02e_prepare_country_rules.R`
3. `02f_extract_country_candidates_native.R`
4. `02g_compare_country_candidates.R`
5. `02b_prepare_disease_rules.R`
6. `02c_extract_disease_candidates_native.R`
7. `02d_compare_disease_candidates.R`
8. `03_build_association_evidence.R`
9. `04_classify_scope.R`
10. `04b_write_policy_review_manifest.R`
11. `05_export_final.R`
12. `06_export_web.R`

It then runs hard production checks and writes
`final/who_don_v2_output_manifest.csv`.

Run the optional clean-vs-v2 audit only when you want to inspect differences
from the accepted clean reference:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --audit-clean
```

To audit existing v2 outputs without regenerating production outputs:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --skip-production --audit-clean
```

To refresh optional post-v2 review/quality-tightening surfaces:

```sh
Rscript scripts/associations/who_don_v2/07_quality_tightening_review_surfaces.R
```

To refresh the deterministic medium native-new country review sample after the
quality-tightening surface exists:

```sh
Rscript scripts/associations/who_don_v2/08_sample_medium_native_new_country_candidates.R
```

To refresh the optional scope QA closure workpack and durable decision file
after the quality-tightening surface exists:

```sh
Rscript scripts/associations/who_don_v2/09_scope_qa_closure.R
```

The former Option B contract-dependency audit scripts have been archived under
`archive/option_b_shadow/`. They are historical migration diagnostics, not
routine production commands. To rerun the initial dependency audit from the
archive:

```sh
Rscript scripts/associations/who_don_v2/archive/option_b_shadow/10_contract_dependency_audit.R
```

This diagnostic is not part of routine production. It compares the current
contract-based final audit to a native-only association build, then writes
summary counts, row-level differences, native-only associations, and candidate
slim-contract rows under `qa/archive/contract_dependency_audit/`.

Stages `07`, `08`, and `09` are not part of routine production. They summarize
country recovery gaps, rank native-new country candidates, sample/close optional
scope adjudication candidates, and write targeted review surfaces for future
manual QA. They are not broad LLM inputs.

Option A native association is now the default production path. Historical
Option A migration/review scripts have been archived under
`archive/option_a_migration/`. For example, the residual scope-consensus
handoff for the earlier Option A scope disagreements is prepared by:

```sh
Rscript scripts/associations/who_don_v2/archive/option_a_migration/22_prepare_option_a_residual_scope_consensus_workpacks.R
```

It writes review workpacks under
`qa/archive/option_a_scope_consensus/` and blank durable reviewer templates under
`review/`. After those templates are completed by reviewers, merge consensus
decisions with:

```sh
Rscript scripts/associations/who_don_v2/archive/option_a_migration/23_merge_option_a_residual_scope_consensus_reviews.R
```

This is optional QA for the Option A transition history and is not part of
routine production.

The completed Option A full-article keep-current decisions have been
materialized into a durable explicit exception file:

```text
review/option_a_full_article_keep_current_exception_rows.csv
```

The historical shadow scripts that applied those decisions and froze the
pre-switch remaining gap are:

```sh
Rscript scripts/associations/who_don_v2/archive/option_a_migration/27_apply_option_a_full_article_scope_decisions_shadow.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/28_option_a_remaining_gap_triage.R
```

Those outputs are written under
`qa/archive/option_a_native_association/` and
`qa/archive/option_a_remaining_gap_triage/`. The production-readiness ledger is
under `qa/archive/option_a_production_readiness/` and records `0` unexplained
blockers. The final acceptance report is:

```text
scripts/associations/who_don_v2/OPTION_A_PRODUCTION_ACCEPTANCE_REPORT.md
```

The current implementation materializes accepted clean outputs into
v2-owned fixtures with `00_materialize_v2_fixtures.R`. Routine v2 production
reads those v2-owned fixture/rule files rather than the old clean folders or
`reference/` seed files. Final v2 association evidence is built through
`review/v2_disease_candidate_adoption_decisions.csv`. Native disease candidates
are adopted only when an explicit deterministic policy accepts them. Native
country candidates are extracted, compared to the accepted country layer, and
applied through `review/v2_country_candidate_adoption_decisions.csv`. Exact
native record-country matches carry native country evidence into final
association evidence; accepted countries not recovered natively are retained as
explicit `legacy_country_exception` rows. Native-only country candidates are not
adopted unless reviewed or covered by deterministic policy; currently this
includes high-confidence title-country candidates after explicit
false-positive filters and the narrow
`medium_native_new_reported_cases_policy` for explicit reported-case wording.
Seeded clean candidates are retained or removed through explicit adoption
decisions, not by an implicit exact-match contract.

Main data outputs are written to split data-layout roots:

```text
pathogen_association_data/staged/who_don_v2/
pathogen_association_data/manual/who_don_v2/
pathogen_association_data/evidence/who_don_v2/
pathogen_association_data/archive/who_don_v2/
```

The canonical v2 outputs are:

- `staged/who_don_v2/reference/who_don_clean_final_seed.csv`
- `staged/who_don_v2/reference/who_don_clean_modelling_seed.csv`
- `staged/who_don_v2/reference/who_don_clean_records_seed.csv`
- `rules/accepted_association_contract.csv`
- `staged/who_don_v2/records/who_don_records_source.csv`
- `staged/who_don_v2/records/who_don_records_clean.csv`
- `staged/who_don_v2/candidates/who_don_country_candidates_native.csv`
- `staged/who_don_v2/candidates/who_don_disease_candidates_native.csv`
- `manual/who_don_v2/review/v2_country_candidate_review_queue.csv`
- `manual/who_don_v2/review/v2_country_candidate_adoption_decisions.csv`
- `evidence/who_don_v2/qa/v2_native_country_vs_accepted_summary.csv`
- `manual/who_don_v2/review/v2_disease_candidate_adoption_decisions.csv`
- `rules/disease_rule_model.csv`
- `evidence/who_don_v2/evidence/who_don_association_evidence.csv`
- `evidence/who_don_v2/evidence/who_don_claims.csv`
- `manual/who_don_v2/review/who_don_review_queue.csv`
- `manual/who_don_v2/review/who_don_scope_adjudication_candidates.csv`
- `manual/who_don_v2/review/who_don_review_decisions_seeded_from_clean.csv`
- `manual/who_don_v2/review/who_don_review_decisions_applied.csv`
- `evidence/who_don_v2/qa/v2_policy_review_decision_manifest.csv`
- `evidence/who_don_v2/final/who_don_country_disease_scope_audit.csv`
- `evidence/who_don_v2/final/who_don_modelling_ready.csv`
- `evidence/who_don_v2/web/who_don_web.json`
- `evidence/who_don_v2/web/who_don_meta.json`
- `evidence/who_don_v2/qa/v2_final_export_summary.csv`
- `evidence/who_don_v2/qa/v2_review_queue_summary.csv`
- `evidence/who_don_v2/qa/v2_production_checks.csv`
- `evidence/who_don_v2/final/who_don_v2_output_manifest.csv`

Only selected durable surfaces are tracked in git:

- `manual/who_don_v2/review/*_decisions.csv` files that encode reviewed or
  policy-controlled decisions
- `evidence/who_don_v2/final/who_don_modelling_ready.csv`, the downstream
  country-evidence surface consumed by role-annotation/readiness scripts

Generated staging, intermediate evidence, QA summaries, web JSON, large final
audit tables, and archived comparison material are ignored unless explicitly
promoted.

The default top-level `qa/` files are:

- `evidence/who_don_v2/qa/v2_country_rule_validation.csv`
- `evidence/who_don_v2/qa/v2_disease_rule_validation.csv`
- `evidence/who_don_v2/qa/v2_final_export_summary.csv`
- `evidence/who_don_v2/qa/v2_native_country_vs_accepted_summary.csv`
- `evidence/who_don_v2/qa/v2_policy_review_decision_manifest.csv`
- `evidence/who_don_v2/qa/v2_production_checks.csv`
- `evidence/who_don_v2/qa/v2_review_queue_summary.csv`

Optional post-v2 quality tightening outputs are generated by stages `07`, `08`,
and `09`. Current completed snapshots have been moved under
`qa/archive/completed_review_surfaces/`; rerun the optional scripts to refresh
archive snapshots if more review is needed:

- `qa/archive/completed_review_surfaces/v2_country_recovery_gap_review.csv`
- `qa/archive/completed_review_surfaces/v2_country_recovery_gap_summary.csv`
- `qa/archive/completed_review_surfaces/v2_scope_adjudication_candidates_enriched.csv`
- `qa/archive/completed_review_surfaces/v2_scope_adjudication_summary.csv`
- `qa/archive/completed_review_surfaces/v2_scope_adjudication_review_sample.csv`
- `qa/archive/completed_review_surfaces/v2_native_new_country_priority_review.csv`
- `qa/archive/completed_review_surfaces/v2_native_new_country_priority_summary.csv`
- `qa/archive/completed_review_surfaces/v2_medium_native_new_country_sample.csv`
- `qa/archive/completed_review_surfaces/v2_medium_native_new_country_sample_manifest.csv`
- `qa/archive/completed_review_surfaces/v2_scope_qa_closure_workpack.csv`
- `qa/archive/completed_review_surfaces/v2_scope_qa_closure_summary.csv`
- `qa/archive/completed_review_surfaces/v2_scope_qa_closure_manifest.csv`
- `qa/archive/completed_review_surfaces/v2_quality_tightening_manifest.csv`
- `review/v2_targeted_adjudication_subset_candidates.csv`
- `review/v2_medium_native_new_country_review_decisions.csv`
- `review/v2_scope_adjudication_review_decisions.csv`

The medium native-new country sample is an optional QA/review artifact. Its
durable decision file records accepted, rejected, and closed-insufficient rows,
but it is not read directly by production. The accepted medium-country pattern
has already been promoted through the normal v2 policy layer as
`medium_native_new_reported_cases_policy`; remaining closed rows stay optional
QA unless a later full-article review finds another repeated,
evidence-supported pattern.

Audit-only clean comparison outputs are generated by `--audit-clean`. Current
snapshots have been moved under `qa/archive/clean_audit/`; rerun the audit
command to refresh archive snapshots:

- `qa/archive/clean_audit/v2_vs_clean_summary.csv`
- `qa/archive/clean_audit/v2_native_adoption_gate.csv`

`qa/archive/clean_audit/v2_vs_clean_summary.csv` is now an
intentional-difference report. Expected
non-exact categories are policy explained, for example rows added by accepted
native disease candidates, rows added by reviewed native-country policy, scope
changes from claim policy, or clean rows removed by seeded-weak policies.
After the pre-v2 archive move, audit/fixture-refresh helpers resolve accepted
clean outputs from `archive/who_don_pre_v2/data/disease_outbreak_news_clean/`
when the old active clean folder is not present.

Current production snapshot, refreshed with the default native Option A runner
on 2026-05-15:

- Association evidence rows: `9957`
- Final audit rows: `9963`
- Modelling rows: `4288`
- Native country candidates: `14459`
- Native disease candidates: `8842`
- Claim rows: `9957`
- Web export rows: `9963`
- Web strict/default rows: `4288`
- Production checks: `19` pass, `0` fail
- Production-readiness unexplained blockers: `0`
- Full-article keep-current exceptions materialized: `222` output rows from
  `216` decision rows
- Web-facing influenza `disease` and `disease_display` labels have `0`
  malformed `Influenza influenza(...)` labels and `0` legacy
  `Influenza (H...)` display labels. The final audit keeps raw
  `disease_standard` values for provenance, so use the web/display fields for
  app labels.
- Clean audit differences are optional audit outputs against the older accepted
  clean reference; use the Option A production-readiness ledger for
  pre-switch current-vs-Option-A blocker status.
- Medium native-new country sample: `1171` sampled rows from `2391` medium
  native-new candidates; `326` accept-pattern, `159` reject-pattern, and `686`
  defer-insufficient-evidence-closed decisions.
- Accepted medium reported-case policy adoption: `374` country-adoption rows
  through `country_candidate_medium_reported_cases_policy`.
- Native-country accepted misses: `157`; no actionable tail or
  `high_rule_review` rows remain.
- Scope adjudication candidates: `661`. Remaining readiness gaps are documented
  as non-blocking optional QA debt or accepted Option A policy differences in
  `qa/archive/option_a_production_readiness/option_a_production_readiness_gap_ledger.csv`.

Compatibility exports are skipped by default. To refresh the old clean-shaped
final filenames for a temporary downstream compatibility check, opt in
explicitly:

```sh
WHO_DON_V2_WRITE_COMPAT=1 Rscript scripts/associations/who_don_v2/run_who_don_v2.R
```

The canonical v2 outputs are the split `evidence/who_don_v2/`,
`manual/who_don_v2/review/`, and `staged/who_don_v2/` files listed above. The
clean-shaped compatibility files are not authoritative v2 outputs.

Web app JSON is now exported natively by `06_export_web.R` from the v2 final
audit and modelling outputs:

- `web/who_don_web.json`
- `web/who_don_meta.json`

The exporter does not copy into the separate `who_don_app` repository by
default. To copy deliberately, set both environment variables:

```sh
WHO_DON_V2_COPY_WEB_TO_APP=1 WHO_DON_APP_DATA_DIR=/path/to/who_don_app/public/data \
  Rscript scripts/associations/who_don_v2/06_export_web.R
```

`06_compare_to_clean.R` is audit-only; it is not required for routine v2
production.

Archived pre-v2 folders:

- `archive/who_don_pre_v2/scripts/who_don/`
- `archive/who_don_pre_v2/scripts/who_don_clean/`
- `archive/who_don_pre_v2/data/disease_outbreak_news/`
- `archive/who_don_pre_v2/data/disease_outbreak_news_clean/`

Archived v2 cleanup leftovers:

- `scripts/associations/who_don_v2/archive/transitional_scripts/`
- `scripts/associations/who_don_v2/archive/historical_docs/`
- `pathogen_association_data/archive/who_don_v2/candidates/`
- `pathogen_association_data/archive/who_don_v2/final_compatibility_exports/`
- `pathogen_association_data/archive/who_don_v2/qa_seeded_baseline/`
- `pathogen_association_data/archive/who_don_v2/qa/orphaned/`
- `pathogen_association_data/archive/who_don_v2/qa/clean_audit/`
- `pathogen_association_data/archive/who_don_v2/qa/completed_review_surfaces/`
- `pathogen_association_data/archive/who_don_v2/qa/stage_diagnostics/`

## LLM Policy

OpenAI Batch submission is not part of the main v2 pipeline. V2 may carry
accepted LLM-derived provenance from older clean/reference outputs, but routine
v2 runs do not call an LLM.

LLM use, if needed later, should be a targeted adjudication fallback after
deterministic rules and review tables have narrowed the problem. Do not send the
broad `review/who_don_review_queue.csv` as an LLM input. That file is a review
surface and includes audit-only rows carried forward from accepted clean
evidence. Rows that are plausible future adjudication candidates are separated
into `review/who_don_scope_adjudication_candidates.csv`; even those should be
manually subsetted before any paid/manual LLM run.

### Option B slim-layer shadow lane

Status: closed as optional QA/audit, not an active production blocker. The
current status note is
`archive/option_b_shadow/OPTION_B_STATUS.md`.

The slim-layer migration path is implemented as a non-production shadow lane.
Run it only if you want to inspect accepted-contract dependency after the
contract dependency audit exists:

```sh
Rscript scripts/associations/who_don_v2/archive/option_b_shadow/10b_contract_dependency_audit_summaries.R
Rscript scripts/associations/who_don_v2/archive/option_b_shadow/11_build_association_evidence_slim_shadow.R
Rscript scripts/associations/who_don_v2/archive/option_b_shadow/12_compare_slim_shadow.R
```

Slim decision schemas live in:

```text
scripts/associations/who_don_v2/rules/slim_country_decisions.csv
scripts/associations/who_don_v2/rules/slim_disease_decisions.csv
scripts/associations/who_don_v2/rules/slim_scope_decisions.csv
scripts/associations/who_don_v2/rules/slim_association_decisions.csv
```

Shadow outputs are written to
`pathogen_association_data/archive/who_don_v2/qa/slim_layer_shadow/`. This lane
does not write to production `evidence/who_don_v2/`,
`manual/who_don_v2/review/`, `staged/who_don_v2/`, or `web/` outputs. Active
slim decision tables are intentionally empty for now; do not copy the large
`accepted_association_contract.csv` into slim tables.
