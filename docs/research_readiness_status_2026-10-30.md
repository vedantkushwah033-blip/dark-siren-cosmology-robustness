# Research readiness and completion audit — 30 October 2026

## Executive status

This repository is a substantial, multi-stage synthetic statistical-methods project for simplified dark-siren inference. It is not a real LVK/GWTC analysis, does not produce an observational Hubble-constant measurement, and does not resolve the Hubble tension.

The repository contains several completed software and simulation studies, including a 500-replicate host-density benchmark, a 13-setting sensitivity screen, a 1,200-dataset H0-dependent-noise study, and paired selection-function misspecification experiments. Those results are conditional on their declared toy assumptions. They do not establish performance on real detector data.

The current focused manuscript candidate is:

**Finite-Validation Uncertainty and Selection-Model Misspecification in Inverse-Probability Host Weighting: A Controlled Dark-Siren Simulation Benchmark**

Its core question is whether an ordinary paired-record bootstrap adequately estimates uncertainty when external validation records are actually sampled without replacement from a finite parent galaxy frame, and how this uncertainty propagates into host-density and toy H0 inference.

## What is already implemented in the wider repository

1. **Core low-redshift inference primitives**
   - Synthetic parent galaxy populations and separate catalogue/GW selection.
   - A selection-aware toy likelihood with a fixed-detected-sample correction.
   - Bounded, reflected weighted KDE host-density estimation with explicit bandwidth semantics.
   - Posterior-grid normalization and numerical-stability safeguards.

2. **Earlier controlled simulation studies**
   - Host-density-to-H0 benchmark with 500 paired replicates and multiple inference arms; raw artifacts and posterior-grid reconstruction audits are documented in the corresponding research PR.
   - Exploratory 13-setting host-density sensitivity screen with 40 paired replicates per setting and independent reconstruction of the recorded posterior fits and aggregate fields.
   - H0-dependent-noise simulation with 1,200 synthetic datasets and a recorded independent replay/audit.
   - Selection-function misspecification experiments with paired conditions and independent raw-artifact reconstruction.
   - Additional distance-uncertainty, posterior-grid, calibration, hidden-host, and host-population sensitivity diagnostics.

3. **Estimated-inclusion validation workflow**
   - Logistic inclusion-probability estimation from external validation records.
   - Oracle, catalogue-only, and estimated-weighting comparison arms.
   - Explicit fit failures, boundary flags, effective sample size, and extreme-weight diagnostics.
   - A persisted paired-record bootstrap workflow with saved resample indices, success masks, and an independent refit/KDE/posterior auditor.
   - Important limitation: the ordinary bootstrap resamples records with replacement and currently does not apply a finite-population correction. It is not yet demonstrated to be calibrated for the finite-parent SRSWOR design.

4. **Current finite-parent SRSWOR work (draft PR #16)**
   - A scenario fixes the parent population, incomplete catalogue, and GW distance sample.
   - Validation draws are exact samples without replacement.
   - Nested sample sizes use prefixes of one random permutation per draw, allowing paired comparisons across validation sizes.
   - Seeds, validation indices, hashes, estimator diagnostics, and failed attempts are persisted.
   - The independent auditor regenerates parent/GW inputs and thinning, checks provenance hashes, and replays nested SRSWOR indices.
   - A second-stage independent reconstruction has now been implemented for successful records: bounded logistic fit, inverse-weight ESS/extremes, reflected weighted KDE, density ISE, selected toy likelihood, posterior mean/median/68% endpoints, interval width, absolute error, and coverage indicator.
   - The auditor also independently re-derives pooled and scenario-balanced summary fields and has adversarial tests for tampered posterior metrics and missing grid metadata.
   - These latest auditor changes are committed but not yet confirmed by a green exact-head CI run.

## Current verification state

- The last completed workflow run inspected for the pre-fix branch head was not green. It failed in the persisted-bootstrap artifact test because the summary manifest omitted the n_validation field, which the adversarial resample-index test required.
- The bootstrap runner now writes n_validation to its summary manifest.
- The finite-parent test assertions and fixture handling were corrected.
- An earlier exact-head run, #41036271780, completed successfully on Python 3.10, 3.11, and 3.12 for the initial finite-parent structural-audit implementation.
- The expanded numerical reconstruction auditor was added after that run. The subsequent workflow attempt was cancelled while additional commits were pushed; a green run for the expanded auditor has not yet been observed. Do not call the latest branch fully validated until all configured Python versions and the smoke artifact audit pass.

## What remains before the focused paper can be considered complete

1. Obtain a green exact-head CI run for all configured Python versions.
2. Extend the finite-parent auditor to independently reconstruct:
   - validation logistic coefficients and fitted probabilities;
   - inverse-probability weights and weight diagnostics;
   - reflected weighted KDE and density error;
   - selected toy likelihood and posterior summaries;
   - scenario-balanced and pooled summaries.
3. Implement the prespecified paired comparison of ordinary bootstrap uncertainty against the empirical repeated-sampling distribution under SRSWOR. Clearly distinguish:
   - posterior credible intervals for H0;
   - repeated-sampling uncertainty in posterior summaries;
   - bootstrap estimates of that repeated-sampling uncertainty.
4. Run deterministic smoke tests, then a pilot with enough independent parent scenarios and nested validation draws to estimate runtime, failure rates, and variance components.
5. Freeze a confirmatory design before inspecting confirmatory outcomes. Set precision targets for coverage and uncertainty-calibration metrics; count parent scenarios as top-level independent units and validation draws as nested observations.
6. Evaluate nonlinear selection misspecification, positivity/overlap stress, and extreme-weight behavior under predeclared regimes.
7. Complete sensitivity to grid resolution, KDE bandwidth, validation sample size, parent-frame size, and GW event count.
8. Reconcile licensing metadata and produce a clean environment/reproduction specification.
9. Finish the manuscript only after the results are independently reconstructed; obtain external domain/statistical review before making publication-readiness claims.

## Completion estimate (engineering judgement, not a measured metric)

- **Repository-wide software and exploratory simulation foundation:** approximately 70–80% implemented.
- **Focused finite-parent benchmark implementation:** approximately 55–65% implemented, pending exact-head CI validation of the expanded numerical auditor.
- **Bootstrap-versus-SRSWOR experiment:** approximately 10–15% complete (a detailed prespecified protocol exists; implementation and results remain outstanding).
- **Focused manuscript as a scientifically defensible, independently audited paper:** approximately 40–50% complete.
- **Real-data/CHIMERA v2.2 validation:** 0% complete for this focused paper; the required exact benchmark data and end-to-end reproduction have not been established.

These percentages describe task completion, not scientific quality, publication probability, or admissions value. The remaining work is disproportionately important: independent estimator reconstruction, bootstrap-vs-SRSWOR calibration, confirmatory simulation, and external review can change the scientific conclusion.

## Claim discipline

Do not claim a new H0 measurement, real LVK validation, superiority over established dark-siren pipelines, or broad methodological novelty. Catalogue incompleteness, host weighting, selection normalization, and population misspecification are established topics. The candidate contribution is narrower: a reproducible finite-parent validation-sampling and uncertainty-calibration benchmark under explicit synthetic assumptions.


## Update after expanded auditor work

The current branch adds an independent numerical reconstruction of successful finite-parent records and an independent re-derivation of summary metrics. It also adds the protocol document `docs/bootstrap_vs_srswor_comparison_protocol_2026-10-30.md`, defining the finite-frame target, ordinary paired bootstrap, FPC-adjusted comparison, primary calibration ratios, hierarchical design, staged pilot/confirmatory plan, and audit requirements. This protocol is not an implemented experiment and contains no numerical finding. Latest recorded code/document commits include `ce4532a`, `9cb2667`, and `d621490`; the latest exact-head CI outcome is not yet confirmed.


## Update: related-work and novelty-boundary audit

A targeted literature-positioning document was added at `docs/related_work_and_novelty_boundary_audit_2026-10-30.md`. It records established work on catalogue incompleteness, host weighting, GP reconstruction, and population inference, and narrows the candidate contribution to design-based calibration of validation-sample uncertainty under finite-parent SRSWOR. This is not a systematic review and does not establish priority or a first-of-its-kind claim. The literature audit does not change the completion estimate: the bootstrap-vs-SRSWOR experiment and confirmatory evidence remain outstanding.

The latest CI run for the expanded auditor and literature audit is run **41040303109**, head commit **7e6e166140600b932a7402fd946b9a2419de1a15**. At the time of this update it is still running; no green exact-head result is claimed.
