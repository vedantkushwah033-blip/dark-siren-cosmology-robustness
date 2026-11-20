# Paired degree-one / degree-two H0 engineering pilot readout

**Run:** GitHub Actions workflow `44741808673`  
**Artifact:** `paired-degree-multireplicate-pilot` (artifact ID `14635611763`)  
**Code head used by the run:** `ba08cddb281feff2370aa23f4552dad5b3bff18f`  
**Current PR #140 head at readout preparation:** `8d435a7ff2f1ac0f7e5887aa2e95228806d08a46`  
**Scope:** finite paired synthetic engineering pilot only. This is not confirmatory calibration, estimator-superiority evidence, or observational inference.

## Design

- Parent galaxies: 1,200; detected events: 30.
- Fixed synthetic truth: (H_0=70); prior support: [30, 110].
- H0 grid: 281 points; quadrature order: 5.
- Curvature values: (gammain{-1,0,+1}).
- Validation frame: 500; validation sizes: 100 and 500.
- Replicates: 20, with matched analysis and validation seeds between degree-one and degree-two fits.
- Degree one uses ([1,x]); degree two uses centered ([1,x,x^2-c_g]). Marginal priors on ((\alpha,\beta)) are matched across degrees.

Each degree has 600 planned lifecycle cells (20 replicates × 3 curvatures × 2 validation sizes × 5 arms). All 600/600 cells were SUCCESS for each degree. There were 60 successful child runs total, and all child audits passed.

## Independent artifact audit

The root independent audit reports **PASS**:
- 2,190 outer-manifest entries and 3,736 nested-manifest entries verified.
- 60 paired child blocks verified.
- 600 paired posterior-summary rows reconstructed.
- 30 paired aggregate rows reconstructed.
- Common-arm likelihood invariance verified.
- The two degree-specific aggregate audits each report 600/600 SUCCESS and reconstruct 30 aggregate rows and 960 paired-contrast rows.

The auditor explicitly limits this to artifact integrity, lifecycle, pairing, and arithmetic. It does not establish calibration.

## Selected H0 results

Values below are median-error bias, RMSE, and nominal 68% interval coverage. All estimates use 20 replicates per cell. Arms A/B/C are common to both degrees; their values are identical by construction and independent audit.

| Curvature | n_validation | Arm | D1 bias | D2 bias | D1 RMSE | D2 RMSE | D1 coverage | D2 coverage |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| -1 | 100 | D: MAP IPW | -2.644 | -0.325 | 4.771 | 4.191 | 0.60 | 0.60 |
| -1 | 100 | E: integrated IPW | -1.949 | +0.029 | 3.989 | 4.202 | 0.65 | 0.75 |
| -1 | 500 | D: MAP IPW | -1.365 | +1.650 | 3.523 | 5.104 | 0.80 | 0.60 |
| -1 | 500 | E: integrated IPW | -1.219 | +1.505 | 3.376 | 4.639 | 0.80 | 0.60 |
| 0 | 100 | D: MAP IPW | -0.007 | +0.879 | 4.196 | 4.139 | 0.60 | 0.55 |
| 0 | 100 | E: integrated IPW | +0.373 | +0.936 | 3.852 | 4.355 | 0.65 | 0.75 |
| 0 | 500 | D: MAP IPW | +1.343 | +1.957 | 3.786 | 4.588 | 0.70 | 0.55 |
| 0 | 500 | E: integrated IPW | +1.395 | +1.928 | 3.809 | 4.434 | 0.65 | 0.55 |
| +1 | 100 | D: MAP IPW | +3.054 | +2.238 | 5.156 | 5.167 | 0.40 | 0.35 |
| +1 | 100 | E: integrated IPW | +3.027 | +1.916 | 4.950 | 4.831 | 0.45 | 0.70 |
| +1 | 500 | D: MAP IPW | +4.212 | +2.224 | 5.605 | 4.976 | 0.30 | 0.45 |
| +1 | 500 | E: integrated IPW | +4.103 | +2.202 | 5.530 | 4.810 | 0.40 | 0.45 |

For context, the common parent-oracle arm has bias +1.816, RMSE 4.318 and coverage 0.65 in all six cells. The common known-probability IPW arm has:
- (gamma=-1): bias +1.377, RMSE 3.920, coverage 0.60.
- (gamma=0): bias +1.259, RMSE 3.548, coverage 0.65.
- (gamma=+1): bias +1.242, RMSE 3.371, coverage 0.70.
These common-arm values do not depend on validation size.

## Paired uncertainty: D2 minus D1 mean absolute error

The aggregate table reports paired t-based 95% intervals across the 20 matched replicates. Positive values mean degree two has larger absolute median error.

| Curvature | n_validation | D: MAP IPW mean difference [95% CI] | E: integrated IPW mean difference [95% CI] |
|---:|---:|---:|---:|
| -1 | 100 | -0.469 [-1.564, +0.626] | +0.074 [-1.087, +1.235] |
| -1 | 500 | +1.296 [-0.229, +2.820] | +1.287 [+0.048, +2.526] |
| 0 | 100 | -0.014 [-1.160, +1.132] | +0.383 [-0.740, +1.505] |
| 0 | 500 | +0.440 [-0.690, +1.570] | +0.462 [-0.453, +1.377] |
| +1 | 100 | +0.018 [-1.569, +1.605] | -0.153 [-1.290, +0.984] |
| +1 | 500 | -0.789 [-2.032, +0.454] | -0.717 [-1.823, +0.390] |

These are 12 exploratory comparisons without a familywise multiplicity adjustment. The one interval excluding zero (integrated IPW, (gamma=-1,n=500)) points toward worse absolute error under degree two. Most comparisons are inconclusive. At positive curvature and n=500, point estimates favor degree two, but paired intervals include zero.

## Weight-stability diagnostics

For D (MAP plug-in IPW), the table below gives median Kish ESS and median top-5%-weight mass across the 20 replicates. For E (coefficient-posterior-integrated IPW), it gives the median across replicates of the *minimum over coefficient-quadrature nodes* Kish ESS and the *maximum over nodes* top-5%-weight mass. These are different summaries and should not be compared as if they were the same estimator-level ESS.

| Curvature | n_validation | D1 D Kish ESS | D2 D Kish ESS | D1 D top-5% mass | D2 D top-5% mass | D1 E minimum-node Kish ESS | D2 E minimum-node Kish ESS | D1 E maximum-node top-5% mass | D2 E maximum-node top-5% mass |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 223.7 | 170.1 | 27.2% | 32.2% | 46.9 | 2.9 | 56.9% | 96.4% |
| -1 | 500 | 169.7 | 110.3 | 31.2% | 36.2% | 73.5 | 8.6 | 45.3% | 78.4% |
| 0 | 100 | 177.7 | 136.0 | 31.6% | 35.1% | 38.2 | 2.9 | 64.1% | 97.6% |
| 0 | 500 | 146.8 | 149.6 | 35.8% | 35.6% | 65.7 | 12.2 | 52.0% | 77.2% |
| +1 | 100 | 136.4 | 164.9 | 37.0% | 35.5% | 27.9 | 4.7 | 69.8% | 96.0% |
| +1 | 500 | 106.1 | 172.1 | 40.0% | 33.8% | 44.0 | 20.7 | 56.6% | 68.3% |

The integrated degree-two arm has particularly severe node-wise IPW concentration in several cells: median minimum-node Kish ESS is about 2.9 at (gamma=-1,n=100) and (gamma=0,n=100), while the median maximum-node top-5% mass is 96–98%. This is a plausible instability mechanism worth investigating; it is not a universal ESS cutoff or proof that degree two is intrinsically inferior. It also explains why coefficient-posterior integration alone should not be assumed to regularize inverse-probability weights.

## Interpretation and limits

1. **No general degree-two advantage is demonstrated.** Its effect changes with curvature, validation size, and IPW arm.
2. At (gamma=-1,n=500), degree two worsens both D and E RMSE point estimates; for E, the paired absolute-error interval is above zero.
3. At (gamma=+1,n=500), degree two reduces D/E RMSE point estimates, but paired absolute-error intervals include zero.
4. At (gamma=+1), both degrees show substantial positive bias and poor coverage in several cells. More flexible inclusion fitting does not resolve all inference failures.
5. The 20-replicate pilot is too small for definitive rankings; 281 H0 grid points and quadrature order 5 are engineering-pilot settings, not a final convergence guarantee.

**Next validation required:** a pre-registered larger paired study (at least 100 replicates per curvature/validation-size cell), higher-resolution H0-grid and quadrature convergence checks, and explicit calibration/coverage analysis. Do not infer an astrophysical H0 result from this synthetic benchmark.
