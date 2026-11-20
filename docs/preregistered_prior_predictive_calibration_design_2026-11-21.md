# Preregistered prior-predictive calibration study — design freeze v1

**Status:** design proposal; not executed; not a calibration result.  
**Prepared:** 2026-11-21  
**Base code head:** `4c512e03aa8d01a57b8afacf8e5ac20bc6f21577` (PR #140).  
**Scope:** synthetic low-redshift Euclidean dark-siren benchmark only. No LVK/GWTC data, observational H0 estimate, real-catalogue validation, or universal estimator claim.

## 1. Confirmatory question

Under a declared prior-predictive distribution for H0 and the existing synthetic population/event model, how do catalogue-only inference, oracle inclusion weighting, known-probability IPW, degree-one estimated IPW, and degree-two estimated IPW compare in posterior calibration and weight stability? Does coefficient integration improve calibration relative to MAP plug-in, and does any apparent gain persist across catalogue-inclusion curvature and validation-sample size?

The primary target is calibration under the *declared synthetic generative model*. Density recovery, fixed-truth error, and calibration are separate outcomes.

## 2. Frozen generative design

Unless a pre-run engineering audit finds a code-level inconsistency, freeze:

- H0 prior: Uniform(30, 110) km s^-1 Mpc^-1.
- Draw one independent H0 truth per replicate from this prior. Use the same truth and same underlying population/event random streams for all paired curvature and model-degree comparisons within that replicate.
- Parent size: 1,200; detected events: 30.
- Curvature: gamma in {-1, 0, +1}; slope beta=-3.9.
- Validation frame: 500; nested validation sizes: 100, 300, 500.
- Five arms: A catalogue-only; B full-parent oracle; C known-probability IPW; D MAP estimated IPW; E coefficient-integrated estimated IPW.
- Compare inclusion model degrees one and two using the existing centered basis [1, x, x^2-c_g] and the prior-matched alpha/beta marginal priors. Keep gamma's prior and all prior scales frozen from the matched-degree audit.
- Keep the catalogue inclusion process distinct from GW detection selection. Do not apply an additional missing-host mixture correction on top of pure IPW.
- Keep oracle arm B and known-probability arm C labelled simulation diagnostics, not deployable procedures.

**Replication unit:** an independent H0 truth, parent population, detected-event realization, and validation-frame realization. Curvature and model-degree arms are paired within this unit. Do not count validation sizes, arms, curvature settings, or degree settings as independent replicates.

## 3. Sample size and execution stages

### Stage 0 — implementation and numerical gate (not inferential)

Before generating confirmatory outcomes, extend the replicate orchestrator so that it:
1. draws H0 truth from the declared prior using a dedicated, reproducible seed stream;
2. records the truth draw and all stream seeds in the top-level config and per-replicate ledger;
3. passes the identical H0 truth to all paired curvature/degree children for that replicate;
4. rejects a rerun into a non-empty output directory;
5. marks every planned cell SUCCESS, FAILED, or NON_EVALUABLE, with no silent dropping;
6. independently reconstructs seed schedule, truth draws, child summaries, paired contrasts, and aggregates.

Run a 10-replicate non-confirmatory engineering smoke. Verify all five arms, both degrees, all curvatures, all validation sizes, and audit outputs. Do not inspect comparative outcomes to change the design.

### Stage 1 — confirmatory prior-predictive run

Target **300 independent replicates**, each with a distinct prior draw of H0 truth. This gives 300 independent prior-predictive draws per curvature/validation-size combination through paired reuse of each replicate's underlying streams. At nominal 68% coverage, 300 replicates give a binomial Monte Carlo standard error of about 2.7 percentage points (roughly ±5.3 points at 95% confidence); 100 replicates would give about 4.7 points (±9.1 points), which is too imprecise for a serious coverage comparison. Treat 100 replicates only as an interim engineering/pilot checkpoint, not the final calibration study.

Planned workload:
- 300 independent replication units.
- 3 curvature settings.
- 3 validation sizes.
- 5 arms.
- 2 inclusion-model degrees.
- 27,000 arm-level lifecycle cells (300 x 3 x 3 x 5 x 2).
- 1,800 degree/replicate/curvature child executions (300 x 3 x 2).

If the engineering gate shows that this workload exceeds the available CI budget, first run a 100-replicate interim checkpoint and then continue with disjoint replicate IDs to reach 300 total. Split execution into deterministic batches of 10 or 20 replication IDs. Batches must use disjoint replicate IDs and the same frozen master seed schedule; aggregate only after verifying uniqueness and full planned-cell coverage. Do not reduce grid or quadrature tolerances to make the job fit.

## 4. Numerical accuracy gate

- Main H0 grid: at least 2,241 equally spaced points over [30,110].
- Quadrature: order 9 for the integrated inclusion posterior, with order 7 versus 9 sensitivity on a prespecified subset.
- Before confirmatory launch, verify posterior-density convergence on a deterministic subset of at least 20 replication units, covering all three curvature values, both degrees, and all validation sizes.
- Compare 2,241 versus 4,481 H0 grids on that subset. Report L1 posterior-density distance, median shift, 68% and 95% interval endpoint shifts, and boundary posterior mass.
- Compare quadrature orders 7 and 9 on the same subset; report changes in integrated posterior density and posterior summaries.
- Freeze acceptance thresholds from the numerical-method tolerances already established in the repository before viewing calibration outcomes. A failed numerical check makes the affected result NON_EVALUABLE; it must not be silently retained as a successful posterior.
- Archive H0 and redshift grids, likelihood arrays, quadrature nodes/weights, and configuration needed for independent reconstruction.

The earlier 561-to-1,121 grid comparison had L1 distance 0.00191816, while 1,121-to-2,241 was 0.00012453. This motivates a finer main grid; it does not itself certify the proposed 2,241-point grid.

## 5. Outcomes and estimands

### Primary calibration outcomes

For each arm, degree, curvature and validation size, report:
- prior-predictive rank histogram for the true H0, using the posterior CDF at the truth;
- empirical 68% and 95% credible-interval coverage;
- rank-uniformity diagnostics with prespecified tests and multiplicity control;
- calibration curves or rank-bin counts with Monte Carlo uncertainty;
- posterior interval width and boundary-mass frequency.

A uniform mean rank alone is not sufficient evidence of calibration. Show the complete rank histogram and interval coverage.

### Secondary estimation outcomes

- Bias and RMSE of posterior median minus true H0.
- Mean absolute error.
- Paired differences D2-D1, E-D, and each estimated method versus catalogue-only and known-probability IPW.
- Paired confidence intervals and Monte Carlo standard errors.
- Failure-inclusive bounds in addition to success-only summaries.

### Weight-stability outcomes

For D and E:
- Kish ESS and entropy ESS;
- maximum normalized weight;
- top 1% and top 5% normalized weight mass;
- support-gap diagnostics;
- for E, posterior-mass-weighted quantiles of node-wise ESS and node-wise top-5% mass, plus worst-node extrema and quadrature ESS.

Do not use a single global ESS threshold to discard results. Treat diagnostics as explanatory measurements, not predictive proof.

## 6. Prespecified comparisons and multiplicity

Primary comparisons:
1. E versus D within each degree/curvature/validation-size cell (does coefficient integration help?).
2. D1 versus D2 and E1 versus E2 in paired cells (does the flexible inclusion model help?).
3. Estimated IPW versus catalogue-only (does correction help?).

Report effect sizes and 95% paired intervals. Apply Holm correction within each explicitly named family of confirmatory comparisons. Label all other contrasts exploratory. Do not select a preferred model based on the smallest unadjusted p-value.

## 7. Failure and non-evaluable accounting

Every planned replicate/curvature/degree/validation-size/arm cell remains in the ledger.

- SUCCESS: runner and independent child audit both pass, with finite and valid posterior output.
- FAILED: software, optimizer, or audit failure.
- NON_EVALUABLE: numerical convergence, support, or posterior validity criterion fails under the frozen protocol.
- NOT_RUN: only permitted when a batch has not yet been executed; final analysis must disclose all such cells.

Publish counts and reasons by arm, degree, curvature and validation size. Success-only calibration is potentially selected; always show failure-inclusive coverage bounds and RMSE bounds under the declared prior support. Do not code failures as covered or as zero error.

## 8. Reproducibility and audit requirements

- Deterministic, separately named seed streams for H0 truth, parent population, event generation, validation frame, and any within-frame sample.
- Pairing checks must verify the intended shared inputs across degrees and curvatures.
- Independent audit must reconstruct the truth draw from its seed and declared prior, all child configs, manifest hashes/byte counts, posterior summaries, rank values, paired contrasts, and aggregate metrics.
- Preserve immutable raw artifacts and SHA-256 manifests.
- Repeated runs with the same seed schedule are reproducibility reruns, not additional independent replicates.
- Keep PRs draft/unmerged until code, numerical gates, independent audits, and interpretation are reviewed.

## 9. Decision rules and claim limits

A method is not called calibrated merely because its average bias is small, RMSE is low, mean rank is near one-half, or a coverage estimate is near its nominal value. Calibration requires coherent rank histograms and coverage across the declared prior-predictive design, with Monte Carlo uncertainty and failures disclosed.

The study can support conclusions only about this synthetic generative model and frozen parameter ranges. It cannot establish observational H0 performance, real survey completeness correction, LVK/GWTC applicability, or universal superiority of an estimator.

## 10. Known implementation gap — do not launch Stage 1 yet

The current `run_quadratic_five_arm_replicates.py` accepts one fixed `--h0-true` for the entire run. It does not yet draw a distinct prior-predictive H0 truth per replication unit or audit that truth draw. Therefore, the existing fixed-truth runner cannot be relabelled as the prior-predictive confirmatory experiment. Stage 0 requires a code change and tests before Stage 1 is authorized.
