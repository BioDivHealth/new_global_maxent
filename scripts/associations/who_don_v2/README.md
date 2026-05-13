# WHO DON V2 Pipeline

This folder contains the active v2 WHO Disease Outbreak News pipeline. The
pre-v2 legacy and clean-migration script/data folders were archived under
`archive/who_don_pre_v2/` after v2 was accepted as the production surface.

V2 is organized around native candidate layers, explicit claim types, a
canonical country-disease-scope evidence table, and reviewed adoption/policy
layers:

The default production command is:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R
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

Stages `07`, `08`, and `09` are not part of routine production. They summarize
country recovery gaps, rank native-new country candidates, sample/close optional
scope adjudication candidates, and write targeted review surfaces for future
manual QA. They are not broad LLM inputs.

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

Main outputs are written to:

```text
pathogen_association_data/WHO/disease_outbreak_news_v2/
```

The canonical v2 outputs are:

- `reference/who_don_clean_final_seed.csv`
- `reference/who_don_clean_modelling_seed.csv`
- `reference/who_don_clean_records_seed.csv`
- `rules/accepted_association_contract.csv`
- `records/who_don_records_source.csv`
- `records/who_don_records_clean.csv`
- `candidates/who_don_country_candidates_native.csv`
- `candidates/who_don_disease_candidates_native.csv`
- `review/v2_country_candidate_review_queue.csv`
- `review/v2_country_candidate_adoption_decisions.csv`
- `qa/v2_native_country_vs_accepted_summary.csv`
- `review/v2_disease_candidate_adoption_decisions.csv`
- `rules/disease_rule_model.csv`
- `evidence/who_don_association_evidence.csv`
- `evidence/who_don_claims.csv`
- `review/who_don_review_queue.csv`
- `review/who_don_scope_adjudication_candidates.csv`
- `review/who_don_review_decisions_seeded_from_clean.csv`
- `review/who_don_review_decisions_applied.csv`
- `qa/v2_policy_review_decision_manifest.csv`
- `final/who_don_country_disease_scope_audit.csv`
- `final/who_don_modelling_ready.csv`
- `web/who_don_web.json`
- `web/who_don_meta.json`
- `qa/v2_final_export_summary.csv`
- `qa/v2_review_queue_summary.csv`
- `qa/v2_production_checks.csv`
- `final/who_don_v2_output_manifest.csv`

The default top-level `qa/` files are:

- `qa/v2_country_rule_validation.csv`
- `qa/v2_disease_rule_validation.csv`
- `qa/v2_final_export_summary.csv`
- `qa/v2_native_country_vs_accepted_summary.csv`
- `qa/v2_policy_review_decision_manifest.csv`
- `qa/v2_production_checks.csv`
- `qa/v2_review_queue_summary.csv`

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

Current production snapshot:

- Final audit rows: `10738`
- Modelling rows: `6486`
- Native country candidates: `14459`
- Native disease candidates: `8246`
- Claim rows: `10738`
- Clean audit differences: `7406` exact matches, `171` v2 additions by disease
  policy, `521` v2 additions by country policy, `2639` scope changes explained
  by claim policy, and `76` clean rows removed by policy.
- Medium native-new country sample: `1171` sampled rows from `2391` medium
  native-new candidates; `326` accept-pattern, `159` reject-pattern, and `686`
  defer-insufficient-evidence-closed decisions.
- Accepted medium reported-case policy adoption: `374` country-adoption rows
  through `country_candidate_medium_reported_cases_policy`.
- Native-country accepted misses: `157`; no actionable tail or
  `high_rule_review` rows remain.
- Scope adjudication candidates: `1245`, with `0` remaining high-priority
  possible focal-event rows. The scope QA closure pass reduced the previous
  `1627` optional rows through conservative event/context claim rules; the
  current remaining closure surface has `216` context-pattern rows and `1029`
  closed-insufficient-evidence rows.

Compatibility exports are skipped by default. To refresh the old clean-shaped
final filenames for a temporary downstream compatibility check, opt in
explicitly:

```sh
WHO_DON_V2_WRITE_COMPAT=1 Rscript scripts/associations/who_don_v2/run_who_don_v2.R
```

The canonical v2 outputs are the v2 `evidence/`, `review/`, `final/`, and `qa/`
files listed above. The clean-shaped compatibility files are not authoritative
v2 outputs.

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
- `pathogen_association_data/WHO/disease_outbreak_news_v2/archive/candidates/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/archive/final_compatibility_exports/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/archive/qa_seeded_baseline/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/qa/archive/orphaned/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/qa/archive/clean_audit/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/qa/archive/completed_review_surfaces/`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/qa/archive/stage_diagnostics/`

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
