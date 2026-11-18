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

## Separate four-arm feasibility smoke (not pooled with the 90-cell pilot)

A separate two-replicate run in PR #82 compared catalogue-only, degree-one, degree-two, and oracle arms at γ=0 and +1, n=100/300/1000. It used a narrower uniform H0 prior [55,85], a 61-point grid, quadrature orders 5/7, and a different event-generation setup; therefore its metrics must not be pooled with the 10-replicate [20,300] pilot above. Run 44410778481 passed the full 947-test suite and completed all 48 planned arm rows. Artifact 14477209501 (ZIP SHA-256 41dcb04d9c8ce66441b4d404a7d9fb9984017e28fb75f521e56a5f0d3d1f8c6f) has 102 aggregate manifest entries and 12 child manifests; I verified all hashes, paired parent/event/frame inputs, nested validation prefixes, and independently reconstructed all 48 posterior summaries from the saved curves.

At n=1000, the two-replicate RMSEs were:

| Curvature | Catalogue-only | Degree one | Degree two | Oracle |
|---:|---:|---:|---:|---:|
| γ=0 | 7.59 | 4.20 | 3.19 | 3.63 |
| γ=+1 | 7.21 | 3.93 | 4.52 | 4.31 |

Degree two is not uniformly better in this tiny smoke: it has lower RMSE at γ=0 but higher RMSE than degree one at γ=+1. Its quadrature ESS is higher (about 29 versus 9.4–9.6 for degree one), but this is a numerical integration diagnostic, not evidence of better calibration. With only two replicates, observed coverage is either 1/2 or 2/2 and cannot establish nominal coverage. This smoke is a feasibility check for the multi-arm implementation, not an estimator ranking.

## Oracle decomposition of the positive H0 bias

A separate independent post-hoc calculation on the corrected 90-cell artifact (documented in PR #89) compared estimated-inclusion IPW with two simulation-only references:

- **Full-parent oracle:** use every parent galaxy with the declared merger weight. Its mean posterior-median bias is +7.58 km s^-1 Mpc^-1 across the 10 unique paired event realizations (Monte Carlo SE 3.15; paired-t 95% interval [+0.46, +14.70]); its 68% interval covered 6/10 truths.
- **Known-π IPW:** use the observed catalogue but weight each galaxy by the exact generating inclusion probability. Its mean median biases are +6.45, +7.43, and +8.84 at γ=-1, 0, and +1, respectively.
- **Estimated-π IPW:** the results in the main table above. Relative to known-π IPW, its mean bias contrast is -2.50 to -2.99 at γ=-1, +0.79 to +0.07 at γ=0, and +4.45 to +4.02 at γ=+1 (n=100 to 500).

Known-π IPW minus full-parent median contrasts are -1.13 (95% paired-t interval [-5.30, +3.03]) at γ=-1, -0.16 ([-2.97, +2.66]) at γ=0, and +1.25 ([-0.16, +2.67]) at γ=+1. These intervals are wide and do not establish that finite-catalogue noise is negligible.

This is a useful diagnostic decomposition, not an additive causal decomposition: the full-parent oracle shows that the four-event finite realization already produces positive median error before catalogue incompleteness is introduced; known-π versus full-parent isolates catalogue-realization effects; estimated-π versus known-π isolates the additional effect of learning the inclusion function. The 10-replicate evidence is too small to identify a systematic likelihood defect. The 100-replicate follow-up in PR #89 is intended to improve precision.

## H0-grid refinement — corrected workflow completed

The corrected paired grid-refinement workflow completed SUCCESS in run 44411660552. It used the audited 90-cell baseline artifact, preserved the archived observable/truth inputs byte-for-byte, and evaluated 27 fixed-truth cells (three replicates per γ∈{-1,0,+1} and n∈{100,300,500}) on 281- and 561-point H0 grids. All 54 refined runs succeeded, all 54 independent child audits passed, and all 998 child-manifest entries verified. The baseline's 2,110 aggregate manifest entries were also verified. Artifact 14478766488 (ZIP SHA-256 e1f926d666fcf9b9bc6769f282e23aa1728005264f99d6f283a1a1829c710d12) was downloaded and independently checked.

The final 81 posterior-grid comparisons are:

| Grid comparison | Within both predeclared tolerances |
|---|---:|
| 141 → 281 points | 1/27 |
| 141 → 561 points | 0/27 |
| 281 → 561 points | 27/27 |

All maximum absolute density differences were below 9.35e-4 (density tolerance 1e-3). The failures were due to posterior-quantile shifts: maximum 0.189 km s^-1 Mpc^-1 against the 0.05 tolerance. Therefore the original 141-point grid is too coarse for the declared quantile tolerance, while 281 points is stable relative to 561 points on these 27 cells.

Important status distinction: the workflow report says PASS because all 54 reruns, audits and artifact checks completed. It separately reports all_comparisons_within_tolerance=false (28/81 comparisons pass both tolerances, 53/81 fail). This is a completed numerical-sensitivity study, not evidence that every grid pair converged, and not a calibration study or SBC. The 100-replicate follow-up uses the 281-point grid, which is supported by the 281→561 comparison on this subset.

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


## H0 grid-refinement stage 1: 141 vs 281 vs 561 points

Workflow run `44411660552`, artifact `14478766488`, re-ran the first three paired fixed-truth replicates across all nine curvature × validation-size cells at 281 and 561 H0 grid points. The archived 90-cell baseline manifest contained 2,110 entries; all were verified. The resulting ZIP manifest contained 1,271 entries, and an independent SHA-256 pass found **zero mismatches**. All 54 refined child runs succeeded, all 54 child raw-data audits passed, and the seven observable/truth input tables were byte-identical between the 281- and 561-point runs in every cell.

The numerical comparisons use the predeclared thresholds: maximum posterior-density difference ≤ 0.001 and maximum absolute change in the reported posterior quantiles ≤ 0.05 km s^-1 Mpc^-1.

| Comparison | Cells | Passed | Maximum density difference | Maximum quantile difference |
|---|---:|---:|---:|---:|
| 141 → 281 | 27 | 1/27 | 0.0009346 | 0.15857 |
| 141 → 561 | 27 | 0/27 | 0.0009345 | 0.18872 |
| 281 → 561 | 27 | **27/27** | **0.0002352** | **0.03813** |

Thus the original 141-point grid is too coarse under the predeclared quantile tolerance for these cells. The 281-point grid agrees with 561 points in all 27 tested cells, while the maximum 281→561 quantile change is 0.0382 km s^-1 Mpc^-1, below the 0.05 threshold. The maximum posterior-density difference is also below 0.001 in every adjacent comparison.

This is evidence for local numerical stability between 281 and 561 points on these 27 selected cells, not proof of global convergence or calibration. It does not validate all 90 pilot cells, and the selection of only three replicates limits the scope. A second stage now runs the same archived inputs at 561 and 1,121 points to check whether the apparent 281→561 stability persists at a finer resolution. No final H0 grid is declared until that stage is independently audited.


## H0 grid-refinement stage 2: 561 vs 1,121 points

Workflow run `44417364702`, artifact `14480782008`, repeated the same 27 fixed-truth cells on 561- and 1,121-point grids. The archived 90-cell baseline manifest (2,110 entries) was verified. The stage-2 ZIP manifest contained 1,273 entries; an independent SHA-256 check found zero mismatches. All 54 refined runs succeeded, all 54 child raw-data audits passed, and the seven observable/truth input tables were byte-identical between 561- and 1,121-point runs in all 27 cells.

| Comparison | Cells | Passed | Maximum posterior-density difference | Maximum quantile difference |
|---|---:|---:|---:|---:|
| 141 → 561 | 27 | 0/27 | 0.0009345 | 0.18872 |
| 141 → 1,121 | 27 | 0/27 | 0.0009344 | 0.19647 |
| 561 → 1,121 | 27 | **27/27** | **0.00005931** | **0.01072** |

Combined with stage 1 (281→561 passed 27/27; maximum density difference 0.0002352 and maximum quantile difference 0.03813), this supports adjacent-grid stability at 281→561 and 561→1,121 on the selected 27 cells. The original 141-point grid fails the predeclared quantile tolerance in every cell. The 1,121-point grid adds little relative to 561 in this subset.

For this pilot subset, **561 points is a conservative numerically stable reference**, and 281 points also passed against 561. This is not a proof for every one of the 90 calibration cells or for all possible event realizations. The 100-replicate follow-up uses 281 points; a subset-specific grid check on its own event realizations remains desirable before treating the entire follow-up as grid-validated. No coverage or SBC claim follows from this numerical experiment.


## Second-stage H0 grid refinement: 561 vs 1,121 points

Workflow run `44417364702`, artifact `14480782008`, completed 54 refined inference runs across the same 27 paired fixed-truth cells at 561 and 1,121 H0 grid points.

Independent post-download checks:
- All 1,273 entries in the stage-2 artifact manifest passed SHA-256 verification.
- All 2,110 entries in the archived 90-cell baseline manifest were verified by the workflow.
- All 54 refined runs are marked SUCCESS and all 54 child raw-data audits report PASS.
- The seven observable/truth CSV inputs are byte-identical between 561- and 1,121-point runs in each of the 27 cells.
- Coefficient quadrature converged at order 7 in 40 runs and order 9 in 14 runs.

| Comparison | Cells | Within predeclared tolerances | Maximum posterior-density difference | Maximum quantile difference |
|---|---:|---:|---:|---:|
| 141 → 561 | 27 | 0/27 | 0.0009345 | 0.18872 |
| 141 → 1,121 | 27 | 0/27 | 0.0009344 | 0.19647 |
| 561 → 1,121 | 27 | **27/27** | **0.00005931** | **0.010715** |

The predeclared limits remain 0.001 for posterior-density differences and 0.05 km s^-1 Mpc^-1 for the maximum absolute change in reported quantiles. The maximum 561→1,121 quantile difference is 0.01072, and the maximum posterior-density difference is 5.93×10^-5. Thus the 561-point grid is stable against this finer 1,121-point check on all 27 selected cells. The 141-point grid fails the quantile tolerance in every cell.

This strengthens the numerical-grid result but does not establish calibration, validate all 90 original pilot cells, or replace an SBC/PIT study. The 27 cells use only the first three paired replicates. A final numerical-grid recommendation should be based on a wider paired subset or all pilot cells if computationally feasible.
