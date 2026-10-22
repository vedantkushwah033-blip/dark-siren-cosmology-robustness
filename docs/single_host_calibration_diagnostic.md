# Single-host calibration diagnostic — initial result

**Status:** initial diagnostic reproduced in an execution environment on 2026-10-22 using the same parameter settings as `scripts/run_single_host_calibration.py`. This is not yet a GitHub Actions artifact. The code and result should be rerun through the committed workflow and compared before treating this as a release result.

## Configuration

- Replicates: 500
- Random seed: 20261022
- True H0: 70 km s^-1 Mpc^-1
- Host redshift: 0.02
- Distance uncertainty: 10 Mpc
- H0 grid: 50–90 km s^-1 Mpc^-1, 4,001 points
- Interval: equal-tail 68%
- Distance model: d_L = c z / H0

## Initial metrics

| Metric | Result |
|---|---:|
| Covered intervals | 349 / 500 |
| Empirical coverage | 0.698 (69.8%) |
| Approximate coverage Monte Carlo SE | 0.0205 (2.05 percentage points) |
| Bias of posterior-grid mode | +0.959 km s^-1 Mpc^-1 |
| RMSE of posterior-grid mode | 8.092 km s^-1 Mpc^-1 |
| Mean interval width | 14.797 km s^-1 Mpc^-1 |

The approximate Monte Carlo standard error is sqrt(p_hat * (1-p_hat)/n). With 500 replicates, observed coverage can fluctuate by a few percentage points. A rough two-standard-error range around the observed coverage is about 65.7%–73.9%; this is only a quick approximation, not an exact confidence interval.

## Interpretation

Under this particular fixed-parameter, single-host setup, observed 68% equal-tail interval coverage was 69.8%, which is close to nominal given Monte Carlo uncertainty. This does **not** explain or reproduce the original project's reported approximately 94.35% coverage. The original result comes from a different, multi-host simulation design and must be audited separately.

This is a repeated-sampling coverage check at a fixed true parameter, not simulation-based calibration (SBC). It does not validate the realism of the dark-siren model, catalogue-incompleteness assumptions, host weights, or GW detection selection.

## Reproducibility checklist

1. Run the manually triggered `single-host-calibration` GitHub Actions workflow with the settings above.
2. Download the `single-host-calibration-results` artifact.
3. Confirm the artifact reports the same settings and metrics (minor platform-level floating-point differences are acceptable).
4. Save the run URL and commit SHA before citing these values in a manuscript.

## Next diagnostic

Reproduce the original project's coverage calculation using its own generative process, likelihood, interval definition, prior range, and random-seed protocol. Do not infer a cause for the original high coverage from this single-host control.
