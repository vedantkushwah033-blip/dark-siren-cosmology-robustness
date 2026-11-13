# CI artifact review: five-arm engineering pilot and numerical diagnostics

**Reviewed:** 2026-11-13  
**Workflow:** GitHub Actions run #4896 (`43456923945`)  
**Reviewed commit:** `e2cfd85ae3ee103804f627735db199881396a09e`  
**Artifact:** `numerical-diagnostics-py3.12`  
**Artifact SHA-256:** `16fd3c3e1067efd26c8d5786e224a67f7435993d1e89bd3a0d30d8e24010a639`

## Scope

This note records an inspection of the Python 3.12 CI artifact from the matched five-arm branch. The workflow passed on Python 3.10, 3.11 and 3.12 for this exact commit. The outputs are engineering smoke tests and numerical-integrity diagnostics, not a confirmatory simulation, calibration study, estimator-ranking result, or observational measurement.

## 1. Multi-replicate five-arm smoke

The saved configuration used 2 independent analysis replicates, 5 arms, nested validation sizes 80 and 160, parent population 800, 3 detected events, validation frame 250, H0 truth 70, 31 H0 grid points, quadrature order 3, distance uncertainty 140 Mpc and detection threshold 250 Mpc. The Scott multiplier 0.7 is explicitly an engineering fixture, not a training/holdout-selected bandwidth.

The independent orchestration report returned **PASS** and reconstructed all 20 planned replicate × validation-size × arm lifecycle cells: 2 successful replicates, 0 failed replicates, 20 posterior-summary rows, 10 aggregate-metric rows and 32 paired-contrast rows. The per-replicate audits independently reconstructed input joins, inclusion probabilities, A–D/E weights, selected likelihoods, posterior summaries and file hashes. The reported maximum absolute likelihood difference was 0.0 in both replicate audits.

**Audit boundary:** the replicate auditors reused stored coefficient-posterior quadrature nodes. They did not independently reconstruct the logistic posterior or establish the correctness of its Laplace approximation and quadrature rule. That remains separate validation work.

The two-replicate smoke confirms lifecycle orchestration and artifact reconstruction only. Two replicates are far too few to support bias, RMSE, coverage, interval-width or estimator-ranking conclusions. A/B/C summaries being unchanged across validation sizes is consistent with those arms not using the validation-fitted inclusion model; D/E can change because they use estimated inclusion information. This is an implementation consistency observation, not evidence of improvement.

## 2. H0 posterior grid stability

| Case | Grid points | Prior support | Posterior median | 68% interval |
|---|---:|---:|---:|---:|
| Baseline | 1,601 | [50, 90] | 70.415478 | [59.723130, 82.573676] |
| Finer baseline | 3,201 | [50, 90] | 70.415477 | [59.723136, 82.573676] |
| Wider numerical grid, same prior | 6,401 | [50, 90] | 70.415477 | [59.723136, 82.573676] |
| Wider grid and wider prior | 6,401 | [30, 110] | 68.200082 | [40.935988, 89.985520] |

For 1,601 versus 3,201 points, the largest reported absolute change among the median and interval endpoints is approximately (5.0\times10^{-6}) km/s. Widening the numerical grid while holding the prior fixed leaves the summaries unchanged to displayed precision. Widening the prior changes the posterior materially; this is prior sensitivity, not a grid-convergence failure.

A targeted deterministic stress case evaluated 1,001 through 64,001 grid points. Its upper 95% endpoint change fell from 0.1061 between 1,001 and 2,001 points to 0.000093 between 32,001 and 64,001 points. This is a numerical stress case, not repeated-sampling evidence.

## 3. Inclusion-coefficient quadrature convergence

The independent validation-inclusion H0 smoke evaluated tensor quadrature orders 5, 7 and 9 (25, 49 and 81 nodes). The saved convergence report declares convergence at order 9 after two consecutive passing comparisons.

From order 7 to order 9:
- maximum absolute H0 posterior-density difference: (1.7442\times10^{-7});
- maximum absolute reported quantile difference: (4.8821\times10^{-5}) km/s.

This demonstrates numerical order stability for one fixture only. It cannot be generalized to other validation sizes, near-separation cases, misspecified inclusion models or tail-sensitive regimes. At order 9, quadrature ESS was about 12.1 of 81 nodes and the largest node mass was about 0.164; these diagnostics should be monitored in harder settings.

## 4. Prior sensitivity

Holding observable data and other analysis settings fixed, the inclusion-logistic slope-prior scale was multiplied by 0.5, 1 and 2. Relative to the reference scale, the H0 median shifted by about -1.251 km/s at multiplier 0.5 and +0.391 km/s at multiplier 2. The reported minimum catalogue Kish ESS values also changed (about 85.5, 52.0 and 44.5). This indicates that prior sensitivity and weight concentration matter in this fixture; it does not establish calibration or an optimal prior.

## 5. Limited conclusions

The inspected artifacts support only that:
1. The five-arm smoke persisted all planned lifecycle cells for its two-replicate fixture.
2. The independent orchestration auditor reconstructed lifecycle and aggregate arithmetic.
3. The replicate audits reproduced likelihood arrays exactly conditional on stored coefficient-posterior nodes.
4. H0 grid summaries were stable under refinement with fixed prior support in the inspected examples.
5. Inclusion quadrature was numerically stable at orders 5, 7 and 9 for one validation-inclusion fixture.
6. The prior-sensitivity check showed nonzero H0 shifts and changes in weight-concentration diagnostics.

They do **not** establish joint SBC, frequentist coverage, robustness under nonlinear misspecification or severe positivity violations, reliable ranking of the five arms, performance at production event counts, real galaxy-catalogue validity, LVK/GWTC data validity, or an observational H0 result.

## 6. Required next gates

1. Obtain exact-head CI for the lifecycle-auditor and truth-dependent failure-bound PRs after the latest fixes.
2. Validate the inclusion posterior independently against high-accuracy low-dimensional integration, including near-separation and tail-sensitive fixtures.
3. Expand quadrature and H0-grid checks across validation sizes, catalogue floors, curvature/misspecification settings and weight-concentration regimes.
4. Complete bandwidth training/holdout selection before fixing the confirmatory KDE bandwidth.
5. Freeze the matched five-arm design, failure rules, Monte Carlo precision targets and replication count before confirmatory execution.
6. Run and independently audit the confirmatory repeated-sampling study, retaining failures and non-evaluable attempts in the denominator.

**Bottom line:** the inspected CI artifact is useful evidence of engineering and numerical reproducibility for selected synthetic fixtures. It is not scientific evidence that estimated catalogue-inclusion correction improves dark-siren inference.
