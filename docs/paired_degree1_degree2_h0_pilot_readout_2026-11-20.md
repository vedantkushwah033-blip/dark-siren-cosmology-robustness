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

Paired 68% coverage differences (D2 minus D1) were also reconstructed from the 20 matched replicate indicators. The only unadjusted 95% t-interval excluding zero was E at \\(\gamma=+1,n=100\\): +0.25 coverage (14/20 versus 9/20), 95% CI [+0.042,+0.458], unadjusted p=0.021. This does not survive a 12-comparison Bonferroni threshold of 0.0042; the other 11 intervals include zero. Treat this as a hypothesis for follow-up, not evidence of improved calibration.

## Weight-stability diagnostics

For E (coefficient-posterior-integrated IPW), the original stored diagnostics report the minimum node-wise Kish ESS and maximum node-wise top-5%-weight mass across all quadrature nodes. Those are worst-case summaries. I independently reconstructed the node-specific inverse-inclusion weights from the archived posterior nodes, posterior quadrature weights, catalogue redshifts, and common host-rate base weights. The reconstructed extrema match the stored diagnostics to floating-point precision across all 240 degree × replicate × curvature × validation-size records.

To avoid treating negligible-posterior-mass tail nodes as representative, the table also reports posterior-mass-weighted quantiles across coefficient nodes. “Weighted q05/q50 ESS” are the 5th percentile and median of node-wise Kish ESS under the validation-posterior quadrature weights. “Weighted q50 top-5%” is the posterior-weighted median node-wise top-5%-weight mass.

| Curvature | n_validation | D1 min ESS | D2 min ESS | D1 weighted q05 / q50 ESS | D2 weighted q05 / q50 ESS | D1 max top-5% | D2 max top-5% | D1 weighted median top-5% | D2 weighted median top-5% |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 46.9 | 2.9 | 103.2 / 223.7 | 21.1 / 170.1 | 56.9% | 96.4% | 27.2% | 32.2% |
| -1 | 500 | 73.5 | 8.6 | 114.3 / 169.7 | 39.4 / 110.3 | 45.3% | 78.4% | 31.2% | 36.2% |
| 0 | 100 | 38.2 | 2.9 | 78.3 / 177.7 | 19.1 / 133.3 | 64.1% | 97.6% | 31.6% | 35.1% |
| 0 | 500 | 65.7 | 12.2 | 99.5 / 146.8 | 59.6 / 149.6 | 52.0% | 77.2% | 35.8% | 35.6% |
| +1 | 100 | 27.9 | 4.7 | 60.1 / 136.4 | 26.5 / 152.5 | 69.8% | 96.0% | 37.0% | 35.5% |
| +1 | 500 | 44.0 | 20.7 | 68.2 / 106.1 | 93.3 / 172.1 | 56.6% | 68.3% | 40.0% | 33.8% |

The distinction matters. Degree two has much more extreme worst-case nodes, especially at n=100, and a lower posterior-weighted 5th-percentile ESS in every curvature/size cell. But its posterior-weighted median ESS is not uniformly worse, and its weighted-median top-5% mass is broadly comparable to degree one. Therefore the earlier worst-node warning should not be read as saying that most posterior mass has catastrophic ESS. It identifies a low-posterior-mass tail that still deserves scrutiny; it does not establish a universal ESS cutoff or estimator failure.

## Exploratory weight-diagnostic / error association

Within each of the 12 degree × curvature × validation-size cells (20 replicates each), I computed Spearman rank correlations between E-arm absolute H0 median error and (i) the posterior-mass-weighted 5th percentile of node-wise Kish ESS and (ii) the posterior-mass-weighted 90th percentile of node-wise top-5%-weight mass. These are post-hoc exploratory checks, not a prespecified diagnostic-prediction test.

The largest unadjusted associations were:
- Degree one, \(\gamma=+1,n=100\): \(\rho(q_{05}(ESS),|e|)=-0.430\), p=0.058; \(\rho(q_{90}(top5),|e|)=+0.483\), p=0.031.
- Degree one, \(\gamma=+1,n=500\): \(\rho(q_{05}(ESS),|e|)=-0.480\), p=0.032; top-5 mass correlation +0.454, p=0.044.
- Degree two, \(\gamma=+1,n=500\): \(\rho(q_{05}(ESS),|e|)=-0.531\), p=0.016; top-5 mass correlation +0.514, p=0.020.

None of the 24 correlation tests survives a Bonferroni threshold of \(0.05/24\\approx0.0021\). Other cells show mixed directions. This pilot therefore does **not** establish that pre-posterior weight diagnostics predict posterior error. A larger preregistered study with held-out replicates and multiplicity-controlled prediction is needed.

## Interpretation and limits

1. **No general degree-two advantage is demonstrated.** Its effect changes with curvature, validation size, and IPW arm.
2. At (gamma=-1,n=500), degree two worsens both D and E RMSE point estimates; for E, the paired absolute-error interval is above zero.
3. At (gamma=+1,n=500), degree two reduces D/E RMSE point estimates, but paired absolute-error intervals include zero.
4. At (gamma=+1), both degrees show substantial positive bias and poor coverage in several cells. More flexible inclusion fitting does not resolve all inference failures.
5. The 20-replicate pilot is too small for definitive rankings; 281 H0 grid points and quadrature order 5 are engineering-pilot settings, not a final convergence guarantee.

**Next validation required:** a pre-registered larger paired study (at least 100 replicates per curvature/validation-size cell), higher-resolution H0-grid and quadrature convergence checks, and explicit calibration/coverage analysis. Do not infer an astrophysical H0 result from this synthetic benchmark.

## Reproduction after the posterior-mode convergence repair

A fresh paired run was completed on the convergence-repair branch after the earlier pilot:
- Workflow: [44752343250](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44752343250)
- Code commit: `8d435a7ff2f1ac0f7e5887aa2e95228806d08a46`
- Artifact: `14646356614`, SHA-256 `0485475ea62f41ef29c25d495eec63a25a5af882efd208b08f182cbec06ac91a`
- Workflow job completed SUCCESS; all run steps, including paired orchestration, independent audit, and artifact upload, completed successfully.

The downloaded ZIP was independently checked against its outer, degree-level, and 120 child manifests: 2,190 + 1,088 + 1,088 + 1,560 = **5,926 file entries verified**, with zero missing files, byte-count mismatches, or SHA-256 mismatches. The root audit reports PASS, 60 paired child blocks, 600 paired posterior-summary rows, 30 paired aggregate rows, and common-arm likelihood invariance. Both degree ledgers contain 600/600 SUCCESS cells and all child audit statuses are PASS.

The paired aggregate CSV reproduces the numerical tables above (20 matched replicates per cell). This is a successful rerun of the same frozen seed/design, not an independent Monte Carlo replication and not an increase in effective replicate count. It supports reproducibility of the archived run and successful completion after the numerical repair; it does not strengthen calibration or estimator-ranking claims.

The full repository test-matrix workflow triggered on the subsequent documentation-only branch update was cancelled before jobs were created. Therefore, this paired run's embedded experiment/audit checks passed, but a fresh full Python 3.10–3.12 repository test-matrix result on the current branch head is still outstanding.

## Clarification of nested manifest counts

The post-repair artifact contains 2,314 ZIP members. Its outer manifest covers 2,190 payload paths; the degree-level manifests each cover 1,088 paths; and 120 child manifests cover 1,560 child-level entries. These scopes are nested and overlap, so their counts must **not** be added and described as a unique-file total. Independent verification found all 2,190 outer-manifest payloads present with matching byte counts and SHA-256 hashes; all 1,088 entries in each degree manifest also matched; all 120 child manifests reported PASS and all 1,560 child-manifest entries matched. The archive's additional members are the nested manifest files and the independent root audit report.
