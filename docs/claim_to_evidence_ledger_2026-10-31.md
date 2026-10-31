# Claim-to-evidence ledger for the dark-siren synthetic manuscript

**Version:** 2026-10-31  
**Purpose:** submission-readiness control. This ledger maps manuscript-level statements to the strongest currently documented evidence and defines the permitted interpretation. It is not a new analysis or independent scientific review.

## Executive disposition

The project supports a substantial, reproducible **synthetic-methodology study** with several audited simulation campaigns. It does not support an observational Hubble-constant result, a validated real-survey completeness correction, a claim of superiority over modern galaxy-catalogue reconstruction, or a new general theorem. The manuscript remains an internal draft until a domain expert reviews the likelihood and the release artifacts are independently rechecked.

## Claim ledger

| Manuscript claim | Evidence currently documented | Audit/reproduction level | Permitted wording | Do not claim |
|---|---|---|---|---|
| In the fixed synthetic host-density Stage B model, correct weighting reduced absolute posterior-mean error relative to catalogue-only analysis. | `docs/selection_misspecification_stage_b_report_2026-10-29.md`; 500 paired replicates, six arms, source commit `12ed103e83bc398df37db8eba8584a9b2b9bc9ad`; artifact ID 12864010259. | Artifact hashes, raw-row aggregates, posterior-grid reconstruction and primary-estimand audit reported PASS. | A prespecified paired contrast favored the known-law weighting arm in this specific toy model. | Weighting is universally better, or the result validates an observational pipeline. |
| Host-density results depend on density-estimation and catalogue-selection settings. | Canonical report: `docs/host_density_sensitivity_screening_results_2026-10-29.md`; supersedes `docs/host_density_validation_report_2026-10-29.md`. | Report documents a 13-setting, 40-paired-replicate-per-setting campaign and its CI replay/audit. | Describe as a finite exploratory sensitivity screen under the declared simulator. | Use the superseded local snapshot's values as current results or call the screen a high-precision confirmatory study. |
| Selection-curve misspecification changes the error penalty in the tested conditions. | `docs/catalogue_selection_confirmatory_replication_2026-10-29.md`; three conditions × 250 paired replicates; source commit `9aadda0ffde113ebb2584aa5425424dbe154dd2b`. | All three condition auditors reconstructed raw rows, arm aggregates and posterior grids; an additional independent reanalysis recomputed estimates from downloaded CSVs. | Report the three condition-specific paired contrasts and their stated Monte Carlo intervals. | Call the conditions independent replications: their seed ranges overlap. Do not pool them as independent conditions. |
| The selection-study artifacts are internally consistent with saved posterior grids and exported replicate records. | `docs/independent_reproduction_confirmatory_2026-10-29.md`. | Hash checks, posterior-grid normalization and aggregate recomputation passed. The report explicitly distinguishes this from seed-level regeneration. | “Independent artifact and numerical reanalysis of saved records.” | “Independent reproduction of the full simulation from seeds” for Study C. |
| The H0-dependent-noise/selection-conditioned model has been exercised in simulation and reveals calibration limitations. | Study B section in `docs/draft_manuscript_2026-10-29.md`, associated fixed-truth/SBC audit reports and workflow artifacts. | The manuscript reports a 1,200-dataset audit/reconstruction; individual claims must be mapped to their exact report and run before release. | State the precise simulated design, audit operation and observed calibration limits. | General calibration claims outside the tested model, or infer that a passing numerical audit proves model correctness. |
| The linear finite-population bootstrap software can generate and independently audit a small paired multi-frame smoke artifact. | `docs/multiframe_v3_smoke_audit_report_2026-10-31.md`; successful run 41176574199 at source head `f869524c3372c1285e6daa212cbfbfa6af4498ea`. | In-job smoke audit reconstructed six frames, 36 outer sample keys and 108 method rows. The first successful run did not archive the smoke directory; a later workflow correction is awaiting exact-head CI/artifact verification. | Software-validation milestone, with the archival caveat stated. | A confirmatory coverage result, method ranking, or verified durable archive from the corrected workflow until the artifact is downloaded and checked. |
| The v1 linear bootstrap pilot gives descriptive coverage/width observations in three fixed finite frames. | `docs/audited_linear_finite_population_bootstrap_pilot_results_2026-10-31.md`; source commit `74e6dacd6690875fdeb6c6d6fc0a1fcaa40fde65`. | Independent audit PASS for its documented outer artifact scope; v1 does not serialize all inner-bootstrap attempt diagnostics. | Historical, descriptive pilot results with denominators and limitations. | Treat as v2/v3 results, claim a winner, or imply complete inner-attempt audit. |

## Mandatory language controls

1. **Synthetic provenance:** use “synthetic,” “toy model,” or “controlled simulation” wherever results are introduced. No real LVK/GWTC observation is analyzed in the manuscript's reported experiments.
2. **Known selection law:** distinguish the oracle/known-inclusion arm from estimated-inclusion and misspecified arms. A known law is a diagnostic control, not a practical solution to survey completeness.
3. **Audit vocabulary:** say exactly which layer was checked: file integrity, reconstruction of posterior grids, aggregate recomputation, seed-level regeneration, or independent expert review. These are not interchangeable.
4. **Coverage:** report numerator, denominator, interval level and Wilson/other justified Monte Carlo interval. Do not infer calibration from lower RMSE or a smaller absolute error.
5. **Overlapping seeds:** the three Study C conditions are not independent replications. Preserve this qualification in the abstract, figures, captions and conclusions.
6. **Model scope:** the low-redshift one-dimensional model does not include realistic sky localization, angular clustering, a three-dimensional galaxy field, magnitude distributions, detector-calibrated selection, or a complete LVK event likelihood.
7. **External comparison:** do not claim superiority over Gaussian-process, clustering-based, LIGHT, or other modern catalogue-completion methods. The models and estimands are not matched.
8. **Publication status:** internal draft only; no submission, peer review, or external faculty sign-off is established by these records.

## Release checklist

- [ ] Freeze the intended manuscript source commit and record the exact file hash.
- [ ] Download every artifact used in a table or figure; verify byte counts and SHA-256 against its manifest.
- [ ] Recompute every reported table cell from raw replicate rows using a separately maintained analysis script.
- [ ] For every figure/table, record source run, source commit, artifact ID, hash, script/config, seed range, replicate denominator and audit scope.
- [ ] Verify the current manuscript uses the canonical host-density screen report, not the superseded snapshot.
- [ ] Check every Study C mention for seed-overlap and no-seed-replay qualifications.
- [ ] Reconcile all Study B calibration statements to exact fixed-truth/SBC artifacts.
- [ ] Obtain domain-expert review of likelihood conditioning, selection normalization, host-population assumptions and interpretation.
- [ ] Keep the CHIMERA v2.2 external benchmark marked NO-GO until the exact compatible inputs and a complete independently auditable execution exist.
- [ ] Run exact-head CI and independently reproduce the final artifact set after all manuscript/code changes.

## Recommended next highest-value scientific step

Do not add another exploratory arm until the manuscript's primary claims are frozen and its raw artifact set is assembled. The strongest next scientific upgrade is an independently seeded, preregistered replication of the **nonlinear dark-siren** result with a fixed estimand, fully replayed data-generation seeds, realistic selection normalization checks, failure accounting, and a separate validation of coverage. The linear SRSWOR bootstrap benchmark is useful statistical software work but cannot validate the nonlinear dark-siren posterior or H0 inference.
