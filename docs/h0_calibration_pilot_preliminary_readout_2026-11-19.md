# Preliminary fixed-truth H0 calibration pilot readout (2026-11-19)

## Scope and provenance

This is a preliminary readout of the archived 90-cell fixed-truth pilot from GitHub Actions run `44402834830`, artifact `14473894295`. The workflow's original aggregate audit failed, so this readout distinguishes raw simulation completion from independent artifact verification.

The pilot design is:
- 10 replicates × 3 curvature values (γ = -1, 0, +1) × 3 validation sizes (100, 300, 500).
- Fixed synthetic H0 = 70 km s^-1 Mpc^-1.
- 1,200 parent galaxies, 500-object validation frame, four detected events per replicate.
- Degree-one fitted inclusion model under quadratic-logit catalogue inclusion.
- Uniform H0 prior [20, 300], 141-point grid (spacing 2 km s^-1 Mpc^-1).
- Inclusion-coefficient quadrature orders 3, 5, 7, 9; two consecutive posterior-stability passes required.

This is a synthetic low-redshift Euclidean benchmark, not an observational H0 measurement, not an LVK/GWTC analysis, and not a real galaxy-survey validation.

## Artifact and pairing checks independently reproduced

The fresh ZIP contains 2,206 files. I independently verified:
- All 2,114 entries in the fresh aggregate manifest match their SHA-256 digests.
- All 90 child audit reports say PASS.
- All 90 planned cells are present and marked SUCCESS; no failed, non-evaluable, or not-run cells.
- All 90 inclusion-quadrature runs report CONVERGED: 72 converged at order 7 and 18 at order 9.
- All 30 replicate × curvature nested-validation blocks were checked.
- All 60 adjacent validation-size prefix comparisons were verified.
- Parent populations, observed GW events, and validation-frame redshifts are invariant across curvature within replicate.
- Catalogue, parent, GW-event, and invariant validation-frame records are identical across validation sizes within replicate/curvature.
- The validation ID/redshift samples are paired across curvature for all 30 replicate × validation-size combinations.
- Every reported aggregate summary field (bias, RMSE, 68/90/95% coverage, interval width, Kish ESS, maximum normalized weight, and status counts present in the archive) was independently recomputed from the replicate ledger and matched.

The historical run `44402834830` failed because its first aggregate auditor compared the size-specific `selected_for_validation` flag across validation sizes. Those flags are supposed to change as the nested sample grows. The corrected replay passed in run `44409048613`; the fresh run `44410603590` then passed the full 90-cell simulation and independent orchestration audit. PR #83 adds stronger invariant-frame checks, including validation-frame inclusion probabilities and retention outcomes; its dedicated workflow remains separate.

## Fixed-truth pilot results

Mean median-bias uncertainty across the 10 independent replicates (95% t intervals, km s^-1 Mpc^-1):

| γ | n=100 | n=300 | n=500 |
|---:|---:|---:|---:|
| -1 | 3.95 [-3.21, 11.11] | 3.95 [-3.15, 11.05] | 3.46 [-3.80, 10.72] |
| 0 | 8.21 [2.13, 14.29] | 8.22 [1.92, 14.52] | 7.50 [1.05, 13.95] |
| +1 | 13.29 [6.04, 20.54] | 13.62 [6.30, 20.94] | 12.86 [5.60, 20.12] |

These are small-pilot descriptive intervals, not confirmatory uncertainty bounds.

All nine cells have 10/10 successful replicates.

| Curvature γ | Validation n | Mean median-bias | RMSE | 68% coverage | 90% coverage | 95% coverage | Mean 68% width | Median minimum Kish ESS | Median maximum normalized weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 3.95 | 10.29 | 0.70 | 0.90 | 1.00 | 24.69 | 14.49 | 0.221 |
| -1 | 300 | 3.95 | 10.21 | 0.70 | 0.90 | 1.00 | 24.23 | 37.57 | 0.127 |
| -1 | 500 | 3.46 | 10.23 | 0.70 | 0.90 | 1.00 | 24.12 | 63.83 | 0.074 |
| 0 | 100 | 8.21 | 11.51 | 0.70 | 1.00 | 1.00 | 25.15 | 10.04 | 0.255 |
| 0 | 300 | 8.22 | 11.72 | 0.80 | 0.90 | 1.00 | 24.60 | 20.43 | 0.165 |
| 0 | 500 | 7.50 | 11.38 | 0.90 | 0.90 | 1.00 | 24.66 | 30.49 | 0.127 |
| +1 | 100 | 13.30 | 16.40 | 0.50 | 0.90 | 0.90 | 25.39 | 9.67 | 0.229 |
| +1 | 300 | 13.63 | 16.73 | 0.50 | 0.80 | 1.00 | 24.66 | 16.28 | 0.152 |
| +1 | 500 | 12.87 | 16.07 | 0.50 | 0.80 | 1.00 | 24.88 | 24.61 | 0.117 |

Bias is in km s^-1 Mpc^-1; interval width is in the same units.

## Paired diagnostics across curvature and validation size

Because the design reuses the same parent/event seeds across curvature and nested validation prefixes within replicate, paired differences are more informative than comparing unrelated cell means. Averaged over the 30 replicate × validation-size pairs:

- γ=0 minus γ=-1 median-bias difference: +4.19 km s^-1 Mpc^-1 (replicate-averaged paired SD 3.50; 95% t interval [1.69, 6.69]).
- γ=+1 minus γ=0: +5.29 km s^-1 Mpc^-1 (SD 4.04; 95% t interval [2.40, 8.17]).
- γ=+1 minus γ=-1: +9.47 km s^-1 Mpc^-1 (SD 5.35; 95% t interval [5.65, 13.30]).

These intervals use the 10 replicate-level averages (each average pools the three validation sizes), not 30 falsely independent replicate×size rows.

By contrast, changing validation size while holding replicate and curvature fixed has much smaller median-bias movement: n=300 minus 100 averages +0.11 (SD 1.32; 95% t interval [-0.83, 1.06]), n=500 minus 300 averages -0.66 (SD 0.99; [-1.37, 0.06]), and n=500 minus 100 averages -0.54 (SD 1.21; [-1.41, 0.32]) km s^-1 Mpc^-1. These are descriptive paired contrasts, not causal decomposition. They suggest that this four-event realization and curvature condition matter more than increasing validation size from 100 to 500 in this particular pilot.

Across all 90 prior-sensitivity records, changing the slope-prior SD multiplier from 1 to 0.5 shifts the H0 posterior median downward by a median 1.21 km s^-1 Mpc^-1 (range -6.61 to -0.24); changing it to 2 shifts the median upward by a median 0.34 (range 0.07 to 2.48) km s^-1 Mpc^-1. This confirms non-negligible prior sensitivity in some realizations, but does not by itself explain the full positive-bias pattern.

## Oracle decomposition of the positive H0 bias

A separate independent post-hoc calculation on the corrected 90-cell artifact (documented in PR #89) compared estimated-inclusion IPW with two simulation-only references:

- **Full-parent oracle:** use every parent galaxy with the declared merger weight. Its mean posterior-median bias is +7.58 km s^-1 Mpc^-1 across the 10 unique paired event realizations (Monte Carlo SE 3.15; paired-t 95% interval [+0.46, +14.70]); its 68% interval covered 6/10 truths.
- **Known-π IPW:** use the observed catalogue but weight each galaxy by the exact generating inclusion probability. Its mean median biases are +6.45, +7.43, and +8.84 at γ=-1, 0, and +1, respectively.
- **Estimated-π IPW:** the results in the main table above. Relative to known-π IPW, its mean bias contrast is -2.50 to -2.99 at γ=-1, +0.79 to +0.07 at γ=0, and +4.45 to +4.02 at γ=+1 (n=100 to 500).

Known-π IPW minus full-parent median contrasts are -1.13 (95% paired-t interval [-5.30, +3.03]) at γ=-1, -0.16 ([-2.97, +2.66]) at γ=0, and +1.25 ([-0.16, +2.67]) at γ=+1. These intervals are wide and do not establish that finite-catalogue noise is negligible.

This is a useful diagnostic decomposition, not an additive causal decomposition: the full-parent oracle shows that the four-event finite realization already produces positive median error before catalogue incompleteness is introduced; known-π versus full-parent isolates catalogue-realization effects; estimated-π versus known-π isolates the additional effect of learning the inclusion function. The 10-replicate evidence is too small to identify a systematic likelihood defect. The 100-replicate follow-up in PR #89 is intended to improve precision.

## Interpretation

1. Posterior-median bias is positive in all nine cells. Mean bias ranges from +3.46 to +13.62 km s^-1 Mpc^-1.
2. The positive-curvature (γ=+1) cells have the largest bias and RMSE, and the lowest observed 68% coverage (5/10 in every validation-size cell).
3. Increasing validation size improves the weight-concentration diagnostics: median minimum Kish ESS rises and maximum normalized weight falls in every curvature condition.
4. Those improved weight diagnostics do not translate into a clear H0 calibration improvement in this small pilot. For γ=+1, mean bias remains about +13 km s^-1 Mpc^-1 across validation sizes.
5. The correctly specified γ=0 reference also has notable positive bias (+7.5 to +8.2). Therefore, the observed bias cannot be attributed solely to quadratic misspecification.
6. The pilot's 10 replicates per cell are far too few for confirmatory coverage statements. Wilson 95% intervals for observed coverage counts are: 5/10, [0.24, 0.76]; 7/10, [0.40, 0.89]; 8/10, [0.49, 0.94]; 9/10, [0.60, 0.98]; and 10/10, [0.72, 1.00]. Thus even 10/10 observed coverage is compatible with substantially lower repeated-sampling coverage.
7. The archive does not contain PIT values, so PIT calibration cannot be assessed from this pilot. It should be added to a future preregistered run rather than inferred from interval coverage.
8. The prior-slope sensitivity records show non-negligible posterior movement when the slope-prior scale is changed. This deserves a paired sensitivity audit before attributing bias to the IPW correction itself.

## Mandatory limitations and next gates

- The H0 grid has spacing 2 km s^-1 Mpc^-1. Coefficient-quadrature convergence is not H0-grid convergence.
- The design conditions on exactly four detected events; it omits count-rate information.
- This is fixed-truth repeated-sampling coverage, not joint SBC.
- The validation sample sizes and event count are small; no confirmatory or universal method claim is justified.
- The fitted coefficient prior is independent Normal with SD (5, 2.5), while the generating coefficients are fixed rather than prior-drawn.
- The historical archive passed corrected replay in workflow `44409048613` (PR #87). The latest post-hoc 90-cell run `44410603590` passed with 90/90 SUCCESS cells, 90/90 child audits, 30 nested blocks, 30 paired validation blocks, and 2,114 manifest entries; I downloaded artifact `14477511669`, verified its ZIP digest and every manifest entry, and separately verified all 90 child manifests. The post-hoc report adds Wilson coverage intervals, t intervals for bias, bootstrap intervals for RMSE and paired contrasts, and explicit failure accounting. PR #83's additional invariant-frame inclusion-probability and retention checks remain a separate, stronger audit.

Next: (i) finish the corrected replay and regression workflows; (ii) run paired H0-grid refinement on the same seeds at finer grids; (iii) diagnose the positive bias in the γ=0 reference using prior, event-count, and likelihood-component sensitivity; and only then (iv) enlarge repeated-sampling calibration and plan joint SBC.


## Paired curvature contrasts in posterior-median bias

Because the parent population, GW event realizations, and validation samples are paired across curvature, I also computed replicate-paired differences in posterior-median error relative to the correctly specified gamma=0 reference. The table reports mean paired differences and ordinary 95% paired-t intervals over 10 replicate seeds.

| Contrast | n_validation=100 | n_validation=300 | n_validation=500 |
|---|---:|---:|---:|
| gamma=-1 minus gamma=0 | -4.26 [-7.20, -1.32] | -4.27 [-6.67, -1.87] | -4.04 [-6.33, -1.75] |
| gamma=+1 minus gamma=0 | +5.08 [+2.32, +7.84] | +5.40 [+2.38, +8.42] | +5.36 [+2.38, +8.34] |

These are exploratory paired contrasts with no multiplicity adjustment. They suggest a curvature-associated shift that is asymmetric in sign: negative curvature shifts median error downward relative to gamma=0, while positive curvature shifts it upward. However, gamma=0 itself has a substantial positive baseline bias, and the pilot has only 10 replicates per cell. The contrasts therefore do not establish that the entire bias is caused by catalogue-inclusion misspecification.


## Exploratory finite-parent oracle decomposition

I independently re-evaluated the latest audited 90-cell artifact (run `44410603590`, artifact `14477511669`) using a direct implementation of the same one-sided selected-event likelihood, but with the complete simulated parent population and its known merger weights. This is a diagnostic oracle, not an available-data estimator and not a superpopulation truth calculation.

For parent galaxy redshifts `z_j`, define
`a_j = 1 + 1.5 (z_j-z_min)/(z_max-z_min)` and normalize `w_j ∝ a_j` over the complete finite parent. The oracle likelihood used was

`L(H0) = product_i [sum_j w_j Normal(d_i; c z_j/H0, sigma^2)] / beta(H0)^N`,

where

`beta(H0) = sum_j w_j Phi((c z_j/H0 - d_threshold)/sigma)`.

This matches the generator's host-weight and one-sided distance-threshold convention. The oracle uses truth-only parent redshifts solely for this diagnostic. Parent IDs/redshifts and observed GW events were checked invariant across the three curvature settings for each replicate; validation-frame IDs/redshifts were also invariant across curvature. The 90 child audits in artifact `14477511669` all report PASS, and I independently rechecked all 2,114 aggregate-manifest SHA-256 entries.

Across the 10 paired event realizations, the finite-parent oracle posterior-median bias has mean **+7.58 km s^-1 Mpc^-1**, replicate SD **9.96**, and RMSE **12.11**. Thus a substantial positive bias is already present when catalogue-inclusion estimation is bypassed and the exact finite parent host mixture is supplied. This is consistent with the very small four-event sample being highly variable; it does not establish the cause of the bias.

The table compares the production validation-posterior result with that same-replicate oracle. Intervals are ordinary paired-t 95% intervals over the 10 replicate seeds, without multiplicity adjustment.

| Curvature γ | Validation n | Mean production median bias | Mean oracle median bias | Production − oracle median (95% paired CI) |
|---:|---:|---:|---:|---:|
| -1 | 100 | +3.95 | +7.58 | -3.63 [-7.06, -0.20] |
| -1 | 300 | +3.95 | +7.58 | -3.63 [-7.18, -0.09] |
| -1 | 500 | +3.46 | +7.58 | -4.12 [-7.62, -0.62] |
| 0 | 100 | +8.21 | +7.58 | +0.63 [-1.99, +3.25] |
| 0 | 300 | +8.22 | +7.58 | +0.64 [-2.60, +3.88] |
| 0 | 500 | +7.50 | +7.58 | -0.08 [-3.33, +3.17] |
| +1 | 100 | +13.30 | +7.58 | +5.72 [+3.45, +7.97] |
| +1 | 300 | +13.63 | +7.58 | +6.05 [+3.95, +8.13] |
| +1 | 500 | +12.87 | +7.58 | +5.29 [+3.34, +7.21] |

Interpretation:
1. The correctly specified γ=0 production result is close to the finite-parent oracle in this paired diagnostic; its positive bias is not primarily explained by estimating inclusion probabilities in this pilot.
2. Negative and positive curvature shift the production posterior median below and above the oracle, respectively. The shift changes little between validation n=100 and n=500 in this 10-replicate pilot.
3. This is an exploratory post-hoc diagnostic, not a preregistered endpoint. The direct oracle implementation has not yet been packaged as a repository audit script or CI test. It should be independently reviewed and made reproducible before being used as a manuscript result.
4. The oracle uses a finite parent population, so it is a matched simulation benchmark rather than an exact superpopulation posterior. The 141-point H0 grid is still unverified for convergence, and no coverage or SBC conclusion follows from this decomposition.
