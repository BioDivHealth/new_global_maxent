# WHO DON v2 Option A Production Acceptance Report

Date: 2026-05-15

## Status

Option A is now the default WHO DON v2 production association mode.

Default command:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R
```

Rollback/reference command for the old accepted-contract association base:

```sh
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --association-mode contract
```

## Final Native Production Snapshot

Validated default native run on 2026-05-15:

| Output | Rows |
| --- | ---: |
| Association evidence | 9,957 |
| Final audit | 9,963 |
| Modelling-ready | 4,288 |
| Web rows | 9,963 |
| Web strict/default rows | 4,288 |
| Native country candidates | 14,459 |
| Native disease candidates | 8,842 |
| Scope adjudication candidates | 661 |
| Production checks | 19 pass / 0 fail |
| Unexplained production-readiness blockers | 0 |

Web-facing influenza labels are clean:

- malformed `Influenza influenza(...)` labels: `0`
- legacy `Influenza (H...)` display labels: `0`

## What Changed

Routine association evidence now comes from native country and native disease
adoption layers through `who_don_v2_native_association.R`. The old
`rules/accepted_association_contract.csv` is retained as an audit/reference
artifact and for rollback/reference comparisons, but it is not the default
association-evidence base.

Resolved full-article keep-current decisions are materialized as a durable
explicit exception file:

```text
pathogen_association_data/WHO/disease_outbreak_news_v2/review/option_a_full_article_keep_current_exception_rows.csv
```

This keeps the accepted exceptions reproducible after default production no
longer writes the old contract-based output.

Active slim decision tables remain empty:

- `rules/slim_association_decisions.csv`: `0`
- `rules/slim_country_decisions.csv`: `0`
- `rules/slim_disease_decisions.csv`: `0`
- `rules/slim_scope_decisions.csv`: `0`

No broad LLM queue is part of routine production.

## Current-vs-Previous Production Differences

The pre-switch current production snapshot had `10,834` final audit rows and
`6,381` modelling rows. Reviewed Option A production now has `9,963` final
audit rows and `4,288` modelling rows.

Pre-switch current-vs-reviewed-Option-A gap:

| Metric | Count |
| --- | ---: |
| Net final audit row gap | 871 |
| Net modelling row gap | 2,093 |
| Current semantic rows missing from reviewed Option A | 846 |
| Reviewed Option A semantic additions versus current | 121 |
| Current modelling row-instances missing from reviewed Option A | 2,599 |
| Reviewed Option A modelling row-instances added versus current | 506 |
| Full-article keep-current exceptions applied | 216 decision rows / 222 output rows |

All differences are categorized in:

```text
pathogen_association_data/WHO/disease_outbreak_news_v2/qa/archive/option_a_production_readiness/option_a_production_readiness_gap_ledger.csv
```

Weighted readiness categories:

| Category | Weighted instances | Status |
| --- | ---: | --- |
| accepted Option A non-focal scope policy | 1,861 | non-blocking |
| accepted old-contract removals | 842 | non-blocking |
| accepted Option A scope promotions | 492 | non-blocking |
| optional native disease recovery QA | 327 | non-blocking |
| optional conservative uncertain-scope QA | 176 | non-blocking |
| optional full-article exception QA | 118 | non-blocking |
| optional native country/pairing recovery QA | 109 | non-blocking |
| accepted native additions | 75 | non-blocking |
| optional native addition QA | 60 | non-blocking |
| optional duplicate/multiplicity QA | 12 | non-blocking |

Unexplained blocker rows: `0`.

## Main Artifacts

- `scripts/associations/who_don_v2/run_who_don_v2.R`
- `scripts/associations/who_don_v2/03_build_association_evidence.R`
- `scripts/associations/who_don_v2/05_export_final.R`
- `scripts/associations/who_don_v2/helpers/who_don_v2_option_a_exceptions.R`
- `scripts/associations/who_don_v2/archive/option_a_migration/29_option_a_production_readiness_gap_review.R`
- `scripts/associations/who_don_v2/archive/option_a_migration/30_materialize_option_a_keep_current_exceptions.R`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/final/who_don_country_disease_scope_audit.csv`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/final/who_don_modelling_ready.csv`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/web/who_don_web.json`
- `pathogen_association_data/WHO/disease_outbreak_news_v2/final/who_don_v2_output_manifest.csv`

## Validation Commands Run

```sh
Rscript scripts/associations/who_don_v2/archive/option_a_migration/15_build_association_evidence_option_a_shadow.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/16_run_option_a_shadow_pipeline.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/17_compare_option_a_to_current.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/18_option_a_transition_adjudication.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/19_option_a_scope_policy_review.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/20_option_a_scope_rule_experiment.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/21_option_a_residual_scope_pattern_review.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/23_merge_option_a_residual_scope_consensus_reviews.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/26_prepare_option_a_full_article_scope_review.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/27_apply_option_a_full_article_scope_decisions_shadow.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/28_option_a_remaining_gap_triage.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/29_option_a_production_readiness_gap_review.R
Rscript scripts/associations/who_don_v2/archive/option_a_migration/30_materialize_option_a_keep_current_exceptions.R
Rscript scripts/associations/who_don_v2/run_who_don_v2.R
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --audit-clean
Rscript scripts/associations/who_don_v2/run_who_don_v2.R --skip-production --audit-clean
```

Additional checks:

- production checks: `19` pass, `0` fail
- `rules/slim_*_decisions.csv`: all `0` rows
- web rows: `9,963`, strict rows: `4,288`
- web malformed influenza labels: `0`
- scripts added/changed in this pass parse successfully

## Remaining Optional QA

The remaining non-blocking QA debt is intentionally narrow and documented in
the readiness ledger. It covers native disease recovery, native country/pairing
recovery, conservative uncertain-scope rows, duplicate/multiplicity checks, and
native additions that are useful to review later but are not required for the
native production switch.

The clean-vs-v2 audit remains optional. It compares against the older accepted
clean reference, not the pre-switch current-vs-Option-A readiness ledger, so its
`unexplained` labels should not be treated as production blockers without
cross-checking the readiness artifacts above.

## Non-Goals

- No broad LLM adjudication queue was introduced.
- No large slim exception layer was populated.
- Old H1N1/global pandemic country-list rows were not broadly preserved.
- Compatibility exports remain opt-in through `WHO_DON_V2_WRITE_COMPAT=1`.
