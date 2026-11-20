# 30-event degree-matched inclusion quadrature: post-run readout (2026-11-21)

## Status and provenance

- Pull request: #141, draft; branch `research/30event-degree-quadrature-convergence-20261121`.
- Workflow: [run 44806809662](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44806809662), completed successfully.
- The Python 3.12 job executed the 30-event experiment, independently audited the child bundles, and uploaded the experiment artifact. Python 3.10 and 3.11 ran the general test matrix successfully; the 30-event experiment is intentionally gated to Python 3.12 and was skipped on those two jobs.
- Artifact: `degree-matched-30event-quadrature-sensitivity` (GitHub artifact ID `14669447416`).
- Design: 1,200 parent galaxies; 30 detected events; fixed truth H0=70; prior support [30,110]; 281 H0 grid points; gamma in {-1,0,+1}; validation frame 500 with nested validation sizes 100 and 500; inclusion degrees 1 and 2; Gauss-Hermite orders 5, 7, 9.
- Seeds: master 2026112101; validation 2026112102.

The downloaded artifact was inspected independently after the workflow. All 18 child audit reports say PASS. Each of the three gamma-level manifests lists 114 payload files; all 342 listed payload entries were present and matched both declared byte counts and SHA-256 hashes. The three gamma-level manifests themselves report PASS. This verifies artifact integrity and the declared reconstruction checks; it is not an independent scientific validation of the model.

## Numerical quadrature findings

For the uncertainty-integrated IPW arm (E), the largest differences between quadrature orders 7 and 9 across all 12 gamma × degree × validation-size cells were:

| Quantity | Maximum absolute difference |
|---|---:|
| Posterior density L1 distance | 0.0006544 |
| Posterior median | 0.002071 km/s/Mpc |
| 2.5th percentile | 0.003339 km/s/Mpc |
| 16th percentile | 0.001244 km/s/Mpc |
| 84th percentile | 0.000352 km/s/Mpc |
| 97.5th percentile | 0.000831 km/s/Mpc |

The largest order-5 versus order-9 differences were materially larger: maximum density L1 0.004818, median shift 0.01674 km/s/Mpc, and 2.5th-percentile shift 0.01736 km/s/Mpc. The order-7 to order-9 comparison is therefore the more defensible numerical stability comparison for this pilot. It does not establish convergence beyond order 9 or across other realizations.

## Paired degree comparison at order 9

Posterior medians in the uncertainty-integrated IPW arm E (km/s/Mpc):

| Curvature gamma | Validation n | Degree 1 | Degree 2 | Degree 2 − Degree 1 |
|---:|---:|---:|---:|---:|
| -1 | 100 | 73.983 | 76.001 | +2.018 |
| -1 | 500 | 72.450 | 76.563 | +4.113 |
| 0 | 100 | 77.707 | 77.267 | -0.440 |
| 0 | 500 | 76.064 | 78.298 | +2.234 |
| +1 | 100 | 77.552 | 76.376 | -1.177 |
| +1 | 500 | 76.381 | 76.750 | +0.369 |

The sign of the degree-2 minus degree-1 contrast changes across curvature and validation size. This one fixed-truth realization does not support a general degree-2 advantage. Nor should these six contrasts be interpreted as calibration or bias estimates.

## Weight-instability warning

Degree-2 uncertainty-integrated IPW has severe node-level weight concentration in this design. Across its saved quadrature nodes and the three curvature settings:

- Minimum node-level Kish ESS was 1.238 at validation size 100.
- Maximum node-level top-5%-weight mass was 0.999923 at validation size 100.
- At validation size 500, minimum node-level Kish ESS was 2.175 and maximum node-level top-5%-weight mass was 0.988894.

These are diagnostic extremes, not universal thresholds. They show that a large nominal validation sample does not by itself guarantee stable inverse-probability weights in this fitted model. The coefficient-posterior quadrature ESS is a different quantity and must not be confused with host-weight Kish ESS.

## Interpretation and next gate

This run is a numerical quadrature and paired-model sensitivity pilot using one fixed synthetic realization. It is not prior-predictive calibration, simulation-based calibration, estimator-superiority evidence, observational H0 inference, or validation on a real galaxy catalogue or public gravitational-wave data.

Next, retain the current run as the numerical convergence check and proceed to an independently seeded, multi-truth prior-predictive design. That design should use at least 100 independent truths per curvature × validation-size cell, preserve paired arms, predeclare posterior-grid and quadrature tolerances, and report lifecycle failures rather than dropping them. Continue to report host-weight concentration alongside posterior outcomes. Do not promote degree 2 or uncertainty integration as a default based on this pilot.
