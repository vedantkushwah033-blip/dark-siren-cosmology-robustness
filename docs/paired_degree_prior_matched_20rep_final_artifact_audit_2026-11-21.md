# Final prior-matched D1/D2 paired pilot audit — 2026-11-21

## Provenance and scope

- PR: #140, branch `research/strict-inclusion-posterior-mode-20261120`.
- Workflow: [44773150331](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44773150331), completed SUCCESS.
- Code head: `489234ec26ef9f505047dccba5f7ab35f456bab4`.
- Artifact: `paired-degree-multireplicate-pilot`, ID `14657642461`.
- Artifact SHA-256: `e92cdbf589548a2113c6db8f9a764d577bdcbc410554722a36e1cffec43b76f1`.
- Scope: fixed-truth, paired synthetic low-redshift Euclidean engineering pilot. Not prior-predictive calibration, an estimator-superiority result, observational H0 inference, or LVK/GWTC analysis.

Design: 20 matched replicates; gamma in {-1, 0, +1}; N_parent=1,200; 30 detected events; validation frame=500; nested validation sizes 100 and 500; H0 truth=70 km/s/Mpc; prior [30,110]; 281-point H0 grid; quadrature order 5; master seed 2026112020. Degree one uses [1,x]; degree two uses centered [1,x,x^2-c_g] with c_g=0.2845005579191241 and matched marginal priors on (alpha,beta).

## Independent artifact verification

The artifact's independent audit reports PASS. I independently rechecked the ZIP payloads against the recorded manifests:

- Root manifest: 2,190/2,190 entries present; zero byte-count or SHA-256 mismatches.
- Degree-one manifest: 1,088/1,088 entries verified; zero mismatches.
- Degree-two manifest: 1,088/1,088 entries verified; zero mismatches.
- Both degree ledgers: 600/600 planned lifecycle cells SUCCESS; 60 successful child runs per degree.
- Paired summary: 600 rows; paired aggregate table: 30 rows; all aggregate metrics independently reconstructed from replicate-level paired contrasts with maximum absolute discrepancy 1.78e-15.
- All 60 paired input/common-likelihood checks passed; maximum common-arm log-likelihood discrepancy 1.71e-13.
- The root auditor independently confirms 60 paired child blocks, 3,736 nested-manifest entries and common-arm invariance.

A separate run on PR #133 (workflow 44754584784, artifact 14660473745) used the same frozen configuration and master seed. Its aggregate bias/RMSE values differed from this PR #140 artifact by at most 1.1e-14 and its coverage values were identical. This is a reproducibility rerun of the same fixed-seed experiment, not an additional independent Monte Carlo sample.

## Selected paired results

RMSE in km/s/Mpc; coverage is empirical inclusion of the fixed truth in the nominal 68% interval. Each cell contains 20 paired replicates.

| gamma | n_validation | Arm | D1 RMSE | D2 RMSE | D1 coverage | D2 coverage |
|---:|---:|---|---:|---:|---:|---:|
| -1 | 100 | D: MAP IPW | 4.771 | 4.191 | 0.60 | 0.60 |
| -1 | 100 | E: integrated IPW | 3.989 | 4.202 | 0.65 | 0.75 |
| -1 | 500 | D: MAP IPW | 3.523 | 5.104 | 0.80 | 0.60 |
| -1 | 500 | E: integrated IPW | 3.376 | 4.639 | 0.80 | 0.60 |
| 0 | 100 | D: MAP IPW | 4.196 | 4.139 | 0.60 | 0.55 |
| 0 | 100 | E: integrated IPW | 3.852 | 4.355 | 0.65 | 0.75 |
| 0 | 500 | D: MAP IPW | 3.786 | 4.588 | 0.70 | 0.55 |
| 0 | 500 | E: integrated IPW | 3.809 | 4.434 | 0.65 | 0.55 |
| +1 | 100 | D: MAP IPW | 5.156 | 5.167 | 0.40 | 0.35 |
| +1 | 100 | E: integrated IPW | 4.950 | 4.831 | 0.45 | 0.70 |
| +1 | 500 | D: MAP IPW | 5.605 | 4.976 | 0.30 | 0.45 |
| +1 | 500 | E: integrated IPW | 5.530 | 4.810 | 0.40 | 0.45 |

The paired D2-minus-D1 mean absolute-error differences (95% paired t intervals) are:

| gamma | n_validation | D: MAP IPW | E: integrated IPW |
|---:|---:|---:|---:|
| -1 | 100 | -0.469 [-1.564, +0.626] | +0.074 [-1.087, +1.235] |
| -1 | 500 | +1.296 [-0.229, +2.820] | +1.287 [+0.048, +2.526] |
| 0 | 100 | -0.014 [-1.160, +1.132] | +0.383 [-0.740, +1.505] |
| 0 | 500 | +0.440 [-0.690, +1.570] | +0.462 [-0.453, +1.377] |
| +1 | 100 | +0.018 [-1.569, +1.605] | -0.153 [-1.290, +0.984] |
| +1 | 500 | -0.789 [-2.032, +0.454] | -0.717 [-1.823, +0.390] |

Positive differences mean larger absolute median error for D2. Only the E, gamma=-1, n=500 interval excludes zero, in the direction of worse D2 absolute error. These are exploratory unadjusted comparisons; they do not support a universal ranking. In particular, the degree effect changes with curvature and validation size.

Common arms A/B/C are identical across degree by construction and independently verified. The full-parent oracle B has bias +1.816, RMSE 4.318 and 68% coverage 0.65 in all six cells. Known-probability IPW C has RMSE 3.920, 3.548 and 3.371 for gamma -1, 0 and +1 respectively (independent of validation size).

## Weight stability

For the integrated-IPW E arm, the median across 20 replicates of the worst quadrature-node Kish ESS and maximum node top-5% weight mass is:

| gamma | n_validation | D1 min-node ESS | D2 min-node ESS | D1 max node top-5% | D2 max node top-5% |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 46.86 | 2.86 | 56.9% | 96.4% |
| -1 | 500 | 73.53 | 8.64 | 45.3% | 78.4% |
| 0 | 100 | 38.24 | 2.93 | 64.1% | 97.6% |
| 0 | 500 | 65.73 | 12.16 | 52.0% | 77.2% |
| +1 | 100 | 27.88 | 4.67 | 69.8% | 96.0% |
| +1 | 500 | 44.04 | 20.72 | 56.6% | 68.3% |

These are worst-node diagnostics, not posterior-mass-weighted summaries. They flag severe tail-node concentration for D2, especially at n=100. They do not imply that every quadrature node or most posterior mass has catastrophic ESS. They do show that a more flexible inclusion model does not automatically stabilize inverse-probability weights.

## Conclusions and next gate

1. The paired pipeline and artifact pass independent integrity, lifecycle, pairing, and arithmetic audits.
2. Degree two does not consistently improve H0 RMSE, coverage, or paired absolute error over degree one.
3. The strongest warning is D2's extreme worst-node IPW concentration at small validation size.
4. This experiment is fixed-truth, uses only 20 replicates per cell, a 281-point H0 grid, and quadrature order 5. It is not a calibration study and cannot support a degree ranking.
5. The next scientific gate is a predeclared prior-predictive study with at least 100 independent population/event realizations per curvature/validation-size cell, finer H0 grids, quadrature convergence checks, explicit lifecycle accounting, and calibration/rank diagnostics. Weight stability must be reported alongside posterior error and coverage.

Do not pool repeated runs of this frozen seed as independent replicates. Do not promote oracle weights as deployable. Keep catalogue inclusion separate from GW detection selection, and do not add a second missing-host mixture correction on top of pure IPW.
