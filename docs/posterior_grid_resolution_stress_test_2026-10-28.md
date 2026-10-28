# Posterior-grid resolution stress test — 2026-10-28

## Status and scope

This is a **targeted numerical stress test of the posterior-grid calculation**, not a new astrophysical result. The initial table below came from an independent reimplementation; a follow-up audit, `scripts/audit_extreme_grid_resolution.py`, now performs the same comparison using the actual `posterior_h0` and `posterior_summary` package functions. Regression tests are in `tests/test_extreme_grid_resolution_audit.py`, and the test workflow has been updated to run the test and archive its JSON report. The latest commit's CI must complete before this repository-level audit is called validated.

## Question

Can a fixed (H_0) grid with spacing (0.025\ \mathrm{km\,s^{-1}\,Mpc^{-1}}) accurately summarize a posterior that is substantially narrower than that spacing?

## Setup

Use the existing simplified single-host likelihood, with host redshift (z=0.05), the low-redshift relation (d_L = cz/H_0), and observed distance chosen to put the likelihood maximum at (H_0=70):

`d_obs = c z / 70`.

The distance likelihood is Gaussian with standard deviation `sigma_d`. Compare:
- **coarse grid:** 1,601 equally spaced points from (H_0=50) to 90 (spacing 0.025);
- **fine reference grid:** 400,001 equally spaced points over the same interval (spacing 0.0001).

Both grids use the same prior support and trapezoidal normalization/quantile calculation. This isolates numerical grid resolution within the toy model; it does not represent a realistic GW distance-error distribution.

## Results

| Distance uncertainty `sigma_d` (Mpc) | Coarse posterior SD | Fine-grid posterior SD | Absolute SD discrepancy | Absolute error in 16th percentile |
|---:|---:|---:|---:|---:|
| 2.00 | 0.65430 | 0.65430 | < 0.00001 | 0.00018 |
| 0.20 | 0.06538 | 0.06538 | < 0.00001 | 0.00184 |
| 0.05 | 0.01628 | 0.01634 | 0.00006 | 0.00502 |
| 0.02 | 0.00091 | 0.00654 | 0.00562 | 0.01051 |
| 0.01 | 0.00000002 | 0.00327 | 0.00327 | 0.01375 |

Units for posterior SD and quantile errors are `km s⁻¹ Mpc⁻¹`. Values are rounded; the underlying grid calculations use full floating-point values.

## Interpretation

When the posterior width is comfortably resolved, the coarse and fine summaries agree closely. Once the likelihood becomes narrower than the grid spacing, the coarse grid can materially misrepresent uncertainty: for `sigma_d=0.02` Mpc, its reported posterior SD is about seven times smaller than the fine-grid value; for `sigma_d=0.01` Mpc, it collapses almost to zero. Quantile errors also become material.

This is a general numerical limitation of summarizing a sharply peaked density from a grid that does not resolve it. It does **not** establish that existing project headline experiments use an inadequate grid; their relevant posterior widths and convergence evidence must be checked separately. Nor does the extreme-distance-precision stress setup establish a realistic GW regime.

## Recommended engineering response

1. The direct-package deterministic stress audit and regression test have now been added. Check their GitHub Actions result for the exact head SHA before treating them as validated.
2. Define an explicit resolution criterion and behavior. A possible criterion is the number of grid intervals spanning the central 68% interval, but the threshold must be tested against reference-grid errors across peak locations, prior boundaries, and mixture/multimodal cases before being adopted.
3. Prefer a warning or explicit convergence diagnostic over silently returning a precise-looking but unresolved summary. Avoid changing default behavior until backward compatibility and downstream workflow effects are assessed.
4. Re-run grid-convergence checks on the exact posterior summaries used in headline experiments, and record the tested commit, grid, and maximum discrepancies.

## Limitations

- The stress test used one host and a deliberately narrow likelihood to probe numerical failure, not a population of realistic gravitational-wave events.
- The fine grid is a numerical reference, not an analytic posterior.
- The original table is from an independent reimplementation; the newly added direct-package audit should reproduce it closely, but that has not yet been confirmed by the pending CI run.
- No production inference code, historical artifact, or headline experiment output was changed by this report.
