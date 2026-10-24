# Validation milestone: joint cosmology experiment and grid audit

**Date:** 24 October 2026

## CI status

The following GitHub Actions runs completed successfully:

- [Full test suite](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117473575): pytest jobs succeeded on Python 3.10, 3.11, and 3.12.
- [Joint distance-model × host-mixture experiment](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117477689): completed successfully; focused tests passed (11 passed), experiment and independent artifact audit passed.
- [Posterior-grid resolution audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117473534): completed successfully; focused tests passed (12 passed).

## Joint experiment audit

The workflow independently recalculated the outputs from raw CSV files and reported:

- 32,400 per-run inference rows
- 108 summary rows
- 135 paired-contrast rows
- 27 experimental cells with 300 replicates per cell
- Maximum absolute discrepancy in independently recomputed summary values: 0.0
- Maximum absolute discrepancy in paired-contrast means: 0.0
- Maximum absolute discrepancy in bootstrap interval endpoints: 0.0

Artifact: [joint-distance-host-mixture-results](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117477689/artifacts/12569556869).

The zero discrepancies mean the audit reproduced the stored summaries and contrasts exactly under the same data and computation definitions. They do not independently validate the scientific assumptions or establish astrophysical realism.

## Grid-resolution audit

The resolution audit tested 810 observations across 27 configurations and four model arms at grid sizes 2,001, 4,001, and 8,001. The maximum absolute difference between the 4,001-point and 8,001-point posterior median was 0.00894 km s⁻¹ Mpc⁻¹; the maximum interval-endpoint difference was below 0.009 km s⁻¹ Mpc⁻¹. This supports using 4,001 points for the tested synthetic setup, while leaving other posterior shapes and likelihoods untested.

Artifact: [joint-grid-resolution-audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117473534/artifacts/12570341647).

## Scope and next scientific step

These are controlled synthetic experiments, not measurements of H0 from real gravitational-wave events. The next scientific step is to evaluate catalogue incompleteness separately from detection selection and host-model misspecification, keeping the simulation-only oracle missing-host correction clearly separated from methods that could be used with real data. Any new conclusions should be based on completed runs, preserved raw outputs, and independent checks.
