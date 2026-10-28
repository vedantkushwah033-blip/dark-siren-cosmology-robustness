# Posterior-grid resolution stress test — 2026-10-28

## Status and scope

This is a **targeted numerical stress test of the posterior-grid calculation**, not a new astrophysical result. The initial table below came from an independent reimplementation. A follow-up audit, `scripts/audit_extreme_grid_resolution.py`, evaluates the same simplified single-host setup using the package's actual `posterior_h0` and `posterior_summary` functions. Regression tests are in `tests/test_extreme_grid_resolution_audit.py`.

The direct-package audit and tests passed on PR head `847fb6ffd46c628e1243bc279e9296f78da48c5d` in the full Python 3.10, 3.11, and 3.12 test workflow: [run 40638997215](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40638997215). The workflow archived `extreme_grid_resolution_audit.json`; the Python 3.12 artifact was inspected after the run.

## Question

Can a fixed (H_0) grid with spacing (0.025\ \mathrm{km\,s^{-1}\,Mpc^{-1}}) accurately summarize a posterior that is substantially narrower than that spacing?

## Setup

Use the existing simplified single-host likelihood, with host redshift (z=0.05), the low-redshift relation (d_L = cz/H_0), and observed distance chosen to put the likelihood maximum at (H_0=70):

`d_obs = c z / 70`.

The distance likelihood is Gaussian with standard deviation `sigma_d`. Compare:
- **coarse grid:** 1,601 equally spaced points from (H_0=50) to 90 (spacing 0.025);
- **fine reference grid:** 400,001 equally spaced points over the same interval (spacing 0.0001).

Both grids use the same prior support and trapezoidal normalization/quantile calculation. This isolates numerical grid resolution within the toy model; it does not represent a realistic GW distance-error distribution.

## Initial independent-reimplementation results

| Distance uncertainty `sigma_d` (Mpc) | Coarse posterior SD | Fine-grid posterior SD | Absolute SD discrepancy | Absolute error in 16th percentile |
|---:|---:|---:|---:|---:|
| 2.00 | 0.65430 | 0.65430 | < 0.00001 | 0.00018 |
| 0.20 | 0.06538 | 0.06538 | < 0.00001 | 0.00184 |
| 0.05 | 0.01628 | 0.01634 | 0.00006 | 0.00502 |
| 0.02 | 0.00091 | 0.00654 | 0.00562 | 0.01051 |
| 0.01 | 0.00000002 | 0.00327 | 0.00327 | 0.01375 |

Units for posterior SD and quantile errors are `km s⁻¹ Mpc⁻¹`. Values are rounded. This table reports the 16th-percentile discrepancy specifically; it is not the maximum error over all reported quantiles.

## Direct-package audit results

The archived JSON compares the production package's coarse-grid summaries with the fine-grid reference. It reports the maximum absolute discrepancy over the 16th, 84th, 2.5th, and 97.5th percentiles plus the median, alongside posterior-SD discrepancy. The predeclared tolerances are 0.001 in H0 for both maximum quantile discrepancy and SD discrepancy.

| `sigma_d` (Mpc) | Coarse posterior SD | Fine-grid posterior SD | Absolute SD discrepancy | Maximum absolute discrepancy across five quantiles | Status |
|---:|---:|---:|---:|---:|---|
| 2.00 | 0.65429946 | 0.65429946 | 2.22e-16 | 0.000330 | Within tolerance |
| 0.20 | 0.06537908 | 0.06537908 | 2.64e-16 | 0.003175 | Underresolved by quantile criterion |
| 0.05 | 0.01628481 | 0.01634465 | 0.0000598 | 0.012297 | Underresolved |
| 0.02 | 0.00091326 | 0.00653786 | 0.00562460 | 0.010955 | Underresolved |
| 0.01 | 0.00000002 | 0.00326893 | 0.00326891 | 0.017343 | Underresolved |

The direct-package audit flags four of five cases. In the `sigma_d=0.20` case, the SD is stable but the maximum discrepancy across all five quantiles exceeds the declared 0.001 tolerance. This illustrates why a check of posterior SD alone is insufficient.

The two tables are **not numerically identical metrics**: the first is from an independent reimplementation and reports the 16th-percentile error; the second uses the repository package and reports the maximum discrepancy over five quantiles. Do not describe the first table as exactly reproduced by the second.

## Interpretation

When the posterior width is comfortably resolved, the coarse and fine summaries agree closely. Once the likelihood becomes narrow relative to the grid spacing, the coarse grid can materially misrepresent uncertainty: for `sigma_d=0.02` Mpc, its reported posterior SD is about seven times smaller than the fine-grid value; for `sigma_d=0.01` Mpc, it collapses almost to zero. Quantile discrepancies also become material. Even when posterior SD agrees, some quantiles can exceed the chosen numerical tolerance.

This is a general numerical limitation of summarizing a sharply peaked density from a grid that does not resolve it. It does **not** establish that existing project headline experiments use an inadequate grid; their relevant posterior widths and convergence evidence must be checked separately. Nor does the extreme-distance-precision stress setup establish a realistic GW regime.

## Recommended engineering response

1. The direct-package deterministic stress audit and regression test passed on the exact head SHA cited above.
2. Define an explicit resolution criterion and behavior. A possible criterion is the number of grid intervals spanning the central 68% interval, but the threshold must be tested against reference-grid errors across peak locations, prior boundaries, and mixture/multimodal cases before being adopted.
3. Prefer a warning or explicit convergence diagnostic over silently returning a precise-looking but unresolved summary. Avoid changing default behavior until backward compatibility and downstream workflow effects are assessed.
4. Re-run grid-convergence checks on the exact posterior summaries used in headline experiments, and record the tested commit, grid, and maximum discrepancies.

## Limitations

- The stress test used one host and a deliberately narrow likelihood to probe numerical failure, not a population of realistic gravitational-wave events.
- The fine grid is a numerical reference, not an analytic posterior.
- The original independent-reimplementation table and direct-package audit use different quantile-error summaries, so they must remain separately identified.
- No production inference code, historical artifact, or headline experiment output was changed by this report.
