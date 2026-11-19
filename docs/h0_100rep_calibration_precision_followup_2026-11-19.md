# 100-replicate fixed-truth H0 calibration precision follow-up

**Status:** completed fixed-truth precision pilot; workflow 44411667480 succeeded at runner SHA `652021943b31490302d586b5b9c77c7ccdfcc8f7`. The 600-cell artifact passed the independent orchestration audit. The independent known-π/full-parent oracle decomposition completed successfully in workflow 44433623231 (PR #97).

## Motivation

The preceding 90-cell pilot used 10 replicates per curvature-by-validation-size cell. It exposed a substantial positive median bias even at the correctly specified reference curvature, but its coverage estimates have very large Monte Carlo uncertainty. The next experiment increases the number of paired replicates to 100 per cell while retaining the same low-redshift Euclidean toy model and fixed-truth target. It is a precision follow-up, not a change in the scientific scope.

## Fixed design

- 100 replicate seeds, crossed with catalogue-inclusion curvature \(\gamma\in\{-1,0,+1\}\) and validation sample size \(n_{\rm val}\in\{100,500\}\): 600 planned cells.
- Fixed generating \(H_0=70\;\mathrm{km\,s^{-1}\,Mpc^{-1}}\), 1,200 parent galaxies, and exactly four detected events per replicate.
- Validation frame size 500. The 100-unit sample is a prefix of the 500-unit simple random sample without replacement.
- The analysis seed initializes the parent-population, temporary catalogue-thinning, and GW-event generation sequence; this combined analysis sequence is reused across curvature conditions within replicate. The validation frame uses a separate `validation_seed`. This runner does not currently expose a third, independent event seed. Catalogue retention outcomes vary with curvature as intended.
- The DGP uses the quadratic-logit inclusion model; inference fits a degree-one logistic inclusion model. Thus \(\gamma=0\) is the correctly specified reference within the fitted family, while \(\gamma=\pm1\) are misspecified conditions.
- The fitted inclusion-coefficient posterior uses the declared independent Gaussian prior, adaptive Laplace-proposal Gauss-Hermite quadrature at orders 3, 5, 7, and 9, and a shared coefficient node across all four events.
- \(H_0\) prior support is uniform on [20, 300], represented on a 281-point grid (spacing 1 km s\(^{-1}\) Mpc\(^{-1}\)). A paired grid-refinement audit has now completed on 27 fixed-truth pilot cells: 281→561 passed both predeclared tolerances in 27/27 comparisons, while 141→281 passed only 1/27. Thus 281 points is supported relative to 561 on this subset; this remains a numerical sensitivity audit, not a calibration result.
- The completed run contains all 600 planned cells in its ledger. A separate code review found that the runner persists the ledger incrementally after attempted children, so an interrupted run could omit later unattempted cells; a pre-populated NOT_RUN/RUNNING ledger is still a follow-up hardening item. In this completed run, all planned cells are accounted for.

## Endpoints and analysis

### Primary descriptive endpoints

For each \((\gamma,n_{\rm val})\) cell, report:

1. Mean posterior-median bias relative to 70.
2. RMSE of posterior medians.
3. Empirical 68%, 90%, and 95% equal-tailed interval coverage.
4. Mean 68% interval width.
5. Median minimum catalogue Kish ESS and median maximum normalized catalogue weight.

Coverage is accompanied by a binomial Monte Carlo standard error and a 95% Wilson interval. These quantify simulation precision; they do not correct model misspecification or establish nominal frequentist coverage.

### Paired contrasts

- At each validation size, compare the posterior-median bias under \(\gamma=+1\) with \(\gamma=-1\), paired by replicate.
- Within each curvature condition, compare \(n_{\rm val}=500\) with \(n_{\rm val}=100\), paired by replicate.
- Report paired means, standard deviations, and Monte Carlo standard errors. Do not interpret the paired contrasts as universal monotonic effects.

### Integrity gates

Before interpreting numerical results:

- Run the focused auditor regression tests.
- Independently verify every child manifest and child audit report.
- Reconcile all 600 planned cell IDs and all status counts against the replicate ledger.
- Verify nested validation-ID prefixes, shared parent/event/frame identities, and cross-curvature validation ID/redshift pairing.
- Recompute all aggregate summary metrics from successful ledger rows.
- Preserve the complete artifact, including failed/non-evaluable attempts, precision summaries, and the independent orchestration audit.

## Interpretation limits

This is **fixed-truth repeated-sampling coverage**, not joint simulation-based calibration. The true inclusion coefficients are fixed by the DGP rather than jointly drawn from the fitted coefficient prior. It conditions on exactly four detected events and therefore does not validate event-rate or detection-count inference. The model remains a synthetic low-redshift Euclidean benchmark; it is not LVK/GWTC analysis, a real galaxy-survey validation, or an observational \(H_0\) measurement.

A 100-replicate cell still has non-negligible Monte Carlo uncertainty (for example, the standard error of a coverage proportion near 0.68 is about 0.047). Treat the results as a materially improved precision diagnostic, not definitive calibration evidence. Any apparent baseline bias at \(\gamma=0\) must be investigated separately rather than attributed automatically to curvature misspecification.

## Reproducibility

The branch-scoped workflow is `.github/workflows/h0-100rep-calibration-precision.yml`. It uses seed 2026111927, uploads the raw 600-cell artifact, and runs the independent orchestration auditor. No production likelihood code or `main` branch is changed by this follow-up.


## Updated numerical-grid evidence

The second-stage paired grid experiment has now also completed on 27 fixed-truth pilot cells (run 44417364702). All 54 child runs at 561 and 1,121 points passed their independent raw-data audits; the seven observable/truth input tables were byte-identical within each pair. The 561→1,121 comparison passed both predeclared tolerances in all 27 cells (maximum posterior-density difference 5.93×10^-5; maximum reported-quantile change 0.01072 km s^-1 Mpc^-1). The 141-point baseline failed the quantile tolerance in all 27 cells. This supports the 281-point grid used in the planned 600-cell follow-up relative to finer grids on this 27-cell subset; it is not a calibration result and does not establish global convergence.


## Fixed-truth results (600-cell run)

The 600-cell run completed with 599 SUCCESS and one NON_EVALUABLE cell; there were no FAILED or NOT_RUN cells. The non-evaluable cell was retained in the ledger and excluded only from fitted-arm posterior summaries. The independent orchestration audit passed all 600 child audits, 300 within-curvature nested-validation blocks, 200 paired-validation blocks, and 14,106 aggregate manifest entries; all summary metrics were recomputed from the ledger.

| Curvature γ | Validation n | Successful / planned | 68% coverage | 90% coverage | 95% coverage | Mean median bias | RMSE | Mean 68% width | Median min Kish ESS | Median max normalized weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 99 / 100 | 0.545 | 0.727 | 0.798 | +2.020 | 17.455 | 23.389 | 19.70 | 0.1552 |
| -1 | 500 | 100 / 100 | 0.510 | 0.740 | 0.780 | +2.255 | 17.215 | 22.725 | 60.89 | 0.0776 |
| 0 | 100 | 100 / 100 | 0.510 | 0.790 | 0.860 | +6.428 | 19.930 | 24.376 | 13.91 | 0.1976 |
| 0 | 500 | 100 / 100 | 0.510 | 0.770 | 0.850 | +6.768 | 19.843 | 23.582 | 42.41 | 0.0954 |
| +1 | 100 | 100 / 100 | 0.560 | 0.820 | 0.880 | +10.396 | 21.587 | 24.517 | 13.30 | 0.1623 |
| +1 | 500 | 100 / 100 | 0.540 | 0.720 | 0.880 | +11.114 | 22.071 | 23.493 | 34.53 | 0.0912 |

Bias and RMSE are in km s\(^{-1}\) Mpc\(^{-1}\). Coverage Monte Carlo standard errors are approximately 0.05 in every cell, so differences of a few percentage points are not precisely resolved. The 68% coverage values (0.51–0.56) are materially below nominal in this fixed-truth pilot; 90% coverage is 0.72–0.82 and 95% coverage is 0.78–0.88. These are descriptive repeated-sampling results, not joint SBC.

### Paired contrasts

- Mean posterior-median bias difference, γ=+1 minus γ=−1: +8.273 (MCSE 0.512; 99 paired replicates) at validation n=100, and +8.859 (MCSE 0.466; 100 paired replicates) at n=500.
- Mean bias difference, validation n=500 minus n=100: +0.152 (MCSE 0.229) for γ=−1, +0.340 (MCSE 0.266) for γ=0, and +0.718 (MCSE 0.258) for γ=+1.

### Interpretation

Larger validation samples clearly improve weight stability: median minimum catalogue Kish ESS rises in every curvature condition, while the median maximum normalized weight falls. However, this did not produce a clear reduction in posterior-median bias or an improvement in coverage. In particular, the γ=0 reference condition still has mean bias around +6.4 to +6.8 and 68% coverage of 0.51. Therefore, the undercoverage and baseline bias cannot be attributed solely to the omitted quadratic curvature in the fitted inclusion model.

The γ=+1 versus γ=−1 paired contrast is large relative to its Monte Carlo SE and supports a curvature-associated difference in this synthetic setting. It does not establish a universal causal or estimator-superiority claim. The full-parent and known-π oracle decomposition on these exact 100 event realizations is the next required analysis before assigning the baseline bias to catalogue estimation, finite-catalogue noise, or the selected-event likelihood.

The single NON_EVALUABLE cell is γ=−1, validation n=100; its fitted-arm result is excluded from fitted-arm means and paired fitted-arm contrasts. Oracle posteriors can still be reconstructed from its archived parent/catalogue/event inputs and should retain the full 100-replicate oracle denominator.

Artifact: workflow 44411667480, artifact 14485703115, ZIP SHA-256 `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`. The full raw artifact contains 14,708 files and is retained for 90 days.


## Same-realization oracle decomposition

The independent oracle workflow (run 44433623231, commit `871d120bf3c68ed5e8d2d893de1117edcff00685`) downloaded the immutable 600-cell artifact, verified its SHA-256, re-hashed all 14,106 aggregate manifest entries through the orchestration audit, and reconstructed known-inclusion-probability IPW and full-parent posteriors independently from the raw parent, catalogue, event, and validation records. The oracle reconstruction used a 561-point H₀ grid. All 600 oracle cells were evaluated, including the archived raw inputs for the one NON_EVALUABLE fitted cell.

| γ | Validation n | Estimated-IPW mean bias | Known-π IPW mean bias | Full-parent mean bias | Estimated 68% coverage | Known-π 68% coverage | Full-parent 68% coverage | Known-π 90% | Full-parent 90% | Known-π 95% | Full-parent 95% |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +2.020 (99 fitted) | +5.770 | +7.048 | 0.545 | 0.49 | 0.54 | 0.74 | 0.83 | 0.79 | 0.91 |
| -1 | 500 | +2.255 | +5.770 | +7.048 | 0.510 | 0.49 | 0.54 | 0.74 | 0.83 | 0.79 | 0.91 |
| 0 | 100 | +6.428 | +6.474 | +7.048 | 0.510 | 0.50 | 0.54 | 0.78 | 0.83 | 0.86 | 0.91 |
| 0 | 500 | +6.768 | +6.474 | +7.048 | 0.510 | 0.50 | 0.54 | 0.78 | 0.83 | 0.86 | 0.91 |
| +1 | 100 | +10.396 | +6.851 | +7.048 | 0.560 | 0.51 | 0.54 | 0.82 | 0.83 | 0.89 | 0.91 |
| +1 | 500 | +11.114 | +6.851 | +7.048 | 0.540 | 0.51 | 0.54 | 0.82 | 0.83 | 0.89 | 0.91 |

Bias is in km s\(^{-1}\) Mpc\(^{-1}\). Monte Carlo standard errors for mean bias are approximately 1.7–1.9 in the estimated and known-π arms and 1.81 for the full-parent oracle. Each oracle coverage estimate uses 100 paired event realizations and has binomial MCSE of roughly 0.05.

### Paired oracle contrasts

| γ | Validation n | Estimated − known-π median | MCSE | Known-π − full-parent median | MCSE |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | -3.594 | 0.312 | -1.278 | 0.502 |
| -1 | 500 | -3.516 | 0.271 | -1.278 | 0.502 |
| 0 | 100 | -0.046 | 0.295 | -0.574 | 0.354 |
| 0 | 500 | +0.294 | 0.163 | -0.574 | 0.354 |
| +1 | 100 | +3.545 | 0.382 | -0.197 | 0.252 |
| +1 | 500 | +4.263 | 0.260 | -0.197 | 0.252 |

The estimated-minus-known-π contrast is strongly curvature-dependent: negative under γ=−1, close to zero under γ=0, and positive under γ=+1. The positive-curvature misspecified fitted arm therefore has a several-unit median shift relative to the known-π oracle, while the negative-curvature fitted arm shifts in the opposite direction. This is evidence of a sign-dependent degree-one inclusion-model misspecification effect in this synthetic design, not a universal ranking of estimated weights.

Known-π versus full-parent median differences are much smaller: −1.278 (MCSE 0.502) at γ=−1, −0.574 (0.354) at γ=0, and −0.197 (0.252) at γ=+1. In this 100-replicate artifact, finite-catalogue realization under the true inclusion probabilities is not the dominant source of the curvature-dependent fitted-arm contrast.

Crucially, the full-parent oracle itself has mean posterior-median bias +7.048 (MCSE 1.809) and coverage 0.54/0.83/0.91 at 68/90/95% on these exact 100 event realizations. Known-π coverage is 0.49–0.51 at 68%, 0.74–0.82 at 90%, and 0.79–0.89 at 95%; estimated-arm coverage is 0.51–0.56, 0.72–0.82, and 0.78–0.88 respectively. All three full-parent coverage estimates are below nominal, but each uses only 100 parent/event realizations: the 68% deficit (0.54 versus 0.68) is the clearest signal, while the 90% and 95% deficits (0.83 and 0.91) are less precisely resolved at this Monte Carlo size. The 68% undercoverage is also present without catalogue incompleteness or estimated inclusion. The baseline undercoverage cannot be assigned solely to the inclusion model; it points to a separate finite-event/selected-likelihood calibration question that requires the larger repeated-sampling baseline and, ultimately, joint SBC. The oracle result is conditional on this fixed-truth, four-event-per-replicate synthetic experiment and does not establish joint SBC.

The oracle artifact is workflow 44433623231, artifact 14488270129, SHA-256 `3bd8aeaf57e604b200d6da646f9ba0430f1e8995f9bdcecb3c1869c8c3b139fa`. Its audit status is PASS; it contains the six-cell oracle summary, paired contrasts, 600-row cell readout, grid configuration, and audit status.


## 20,000-draw full-parent repeated-sampling baseline

A separate workflow reran the full-parent selected-event likelihood on the **same 100 archived parent populations**, generating 200 fresh four-event realizations per parent (20,000 total). This holds the finite parent population fixed within each cluster and isolates event-sample variation conditional on those parent populations.

**Provenance**
- Workflow: `44435303396`
- Artifact: `14489935029`
- Artifact ZIP SHA-256: `c60f5cd64f489a19fc0595a7a394394fed2dfede6fd0a359d1036700d9d38748`
- 100 parent clusters × 200 fresh event realizations; 20,000 total.
- H0 truth 70; 561-point H0 grid; four detected events per realization.
- The workflow independently reconstructed all 100 parent summaries from event-level rows, verified all source manifest entries, and passed.

### Pooled event-level result

| Metric | Estimate | Event-level MCSE / 95% interval |
|---|---:|---:|
| Mean posterior-median bias (km s^-1 Mpc^-1) | +6.4388 | MCSE 0.1301; t 95% [6.1838, 6.6937] |
| 68% coverage | 0.6535 | MCSE 0.00336; Wilson 95% [0.6469, 0.6601] |
| 90% coverage | 0.8794 | MCSE 0.00230; Wilson 95% [0.8748, 0.8838] |
| 95% coverage | 0.9335 | MCSE 0.00176; Wilson 95% [0.9300, 0.9369] |
| Mean 68% interval width (km s^-1 Mpc^-1) | 23.9233 | — |

### Parent-cluster uncertainty

The independent sampling units for unconditional generalization are the 100 parent populations, not the 20,000 event draws.

| Metric | Mean across parent clusters | Parent-cluster MCSE | Parent-cluster t 95% interval |
|---|---:|---:|---:|
| Mean posterior-median bias | +6.4388 | 0.1408 | [6.1595, 6.7181] |
| 68% coverage | 0.6535 | 0.00293 | [0.6477, 0.6593] |
| 90% coverage | 0.8794 | 0.00255 | [0.8743, 0.8845] |
| 95% coverage | 0.9335 | 0.00187 | [0.9298, 0.9372] |

The event-level MCSE is conditional on the archived parent populations; the parent-cluster MCSE is the appropriate uncertainty for generalizing over parent-population draws.

### Reconciliation and interpretation

The original 100 parent/event realizations had full-parent coverage 0.54/0.83/0.91. Repeatedly resampling detected events on those **same** 100 parents raises it to 0.6535/0.8794/0.9335. Thus the original four-event realization set was not representative of the event-sampling distribution, especially at the 68% level. The +7.048 median-bias estimate from the original 100 event samples is also compatible with the repeated-draw estimate +6.439 once its original MCSE 1.809 is considered.

Two subsequent complete reruns reproduced these 20,000-draw summary values exactly:
- Workflow 44435303348, artifact 14490785294, ZIP SHA-256 `f44ecb66cd541f4cda376a26534eb31255feeec4493ff4eb133d52ee63caac42`.
- Workflow 44436218047, artifact 14490011968, ZIP SHA-256 `204f32af9860743e9964745fa5e6f97fcb097bf92f10d6c8fa579343a32e7d1c`.

Both had 100 parent clusters × 200 event draws, and both reproduced mean median bias +6.4388, coverage 0.6535/0.8794/0.9335, and the parent-cluster intervals above. Their archive hashes differ because they are separate workflow artifacts; the reported numerical summaries agree. These reruns are deterministic computational replications of the same script and seed schedule, not an independent reimplementation of the 20,000-event sampler. The workflow's independent checks verify source manifests and recompute parent-cluster aggregates from the saved event-level rows. A separate raw-seed event-level auditor would be needed to independently reconstruct all 20,000 event samples and posterior curves.

The paired 500-new-parent event-count experiment in PR #101 compared nested 4-event and 30-event samples on the same 500 parent/event realizations. Four-event coverage was 0.696/0.910/0.954 with mean median bias +5.060 (MCSE 0.717); the 30-event arm had coverage 0.662/0.860/0.920 and mean median bias +0.866 (MCSE 0.182). The paired 30-minus-4 median-bias contrast was −4.194 (MCSE 0.699; 95% interval [−5.567,−2.820]), and mean 68% interval width fell from 23.44 to 7.35 km s\(^{-1}\) Mpc\(^{-1}\). This supports a strong event-count effect on median bias and precision, while modest undercoverage remains in the 30-event arm. The 500-parent sample is different from the 100 archived parent populations used in the 20,000-draw conditional baseline, so these results should not be pooled or treated as a paired contrast with that baseline.

The main conclusion is now clearer: the full-parent fixed-count likelihood has positive median bias around +6.4 km s^-1 Mpc^-1 in this four-event setup, and conditional repeated-sampling coverage is below nominal (about 65.4%, 87.9%, and 93.4%). The original 54% 68%-coverage estimate overstated the deficit because it used only one four-event realization per parent. This remains fixed-truth coverage, not joint SBC. The paired 30-event experiment separately shows that more detected events greatly reduce median bias and interval width, but does not by itself resolve the remaining coverage behavior.


## Cross-experiment synthesis and current scientific position

The completed synthetic evidence now separates three effects that should not be conflated:

1. **Baseline selected-likelihood / finite-event behavior.** The full-parent oracle has mean posterior-median bias about +6.44 km s^-1 Mpc^-1 in 20,000 repeated four-event draws on 100 parent populations, with 68/90/95% coverage 0.6535/0.8794/0.9335. This bias and undercoverage persist without catalogue thinning or estimated inclusion probabilities.
2. **Finite detected-event count.** In 500 new parent populations, moving from four to thirty nested detected events reduced mean median bias from +5.06 to +0.87, RMSE from 16.79 to 4.16, and mean 68% interval width from 23.44 to 7.35 km s^-1 Mpc^-1. The paired median-bias change was −4.19 (MCSE 0.70; 95% interval [−5.57,−2.82]). This is strong evidence of a finite-event effect in the synthetic full-parent model; 30-event coverage remains somewhat below nominal.
3. **Catalogue-inclusion model misspecification.** On the first 100-replicate H0 artifact, the estimated-minus-known-π posterior-median contrast is sign-dependent: about −3.5 under γ=−1, near zero under γ=0, and +3.5 to +4.3 under γ=+1. The known-π oracle and full-parent oracle are much closer than the estimated-versus-known-π contrast. This supports a curvature-associated misspecification effect in this synthetic setup, not a general estimator ranking.

The validation-size increase from 100 to 500 improves weight stability (higher Kish ESS, lower maximum normalized weight) but does not clearly improve H0 bias or coverage in the first pilot. Weight diagnostics alone therefore do not guarantee posterior calibration.

### Active follow-ups

- PR #96: independent second 100-replicate/600-cell calibration replication, seed 2026111930, workflow #44437373259; currently running.
- PR #104: independent raw-seed reconstruction audit of the 500-replicate event-count experiment; workflow #44437893958 passed, regenerating all 500 parent populations, 15,000 detected events, and 1,000 posterior curves.
- PR #107: H0-grid refinement on the first 100-replicate artifact, 17 or 18 successful cells depending on whether the non-evaluable cell falls in the selected first-three-replicate subset; 281/561/1,121 grids, workflow #44438462754.
- PR #110: pooled summary across two independent 100-replicate calibration replications; wait for the second complete artifact and both orchestration audits before interpreting pooled estimates.
- PR #111: 1,000-replicate full-parent prior-predictive SBC extension; queued, with a 1,121-point H0 grid.
- PR #108: quadratic-logit five-arm H0 smoke with 30 detected events and gamma = −1, 0, +1; queued, one paired realization per curvature.

### Claim boundary

All results are from a synthetic low-redshift Euclidean fixed-truth model with exactly four or thirty detected events and finite simulated parent catalogues. They are not joint SBC, not real galaxy-survey validation, not LVK/GWTC analysis, and not observational H0 measurement. The defensible contribution remains a controlled benchmark of how finite validation, inclusion-model curvature misspecification, IPW weight concentration, and finite detected-event count interact in this specific synthetic likelihood.
