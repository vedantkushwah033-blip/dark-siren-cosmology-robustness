# Independent raw-artifact audit: broad-transition selection experiment

**Audit date:** 2026-10-29  
**Scope:** independent reconstruction from the saved broad-transition workflow artifact. This is an audit of exported posterior grids and aggregate summaries, not a fresh rerun of the simulation generator or a physical validation of gravitational-wave inference.

## Artifact provenance

- Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`
- Source experiment commit recorded in the artifact manifest: `9aadda0ffde113ebb2584aa5425424dbe154dd2b`
- Workflow run: [#40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231)
- Artifact: `catalogue-selection-confirmatory-broad-transition-9aadda0ffde113ebb2584aa5425424dbe154dd2b`
- Artifact SHA-256 reported by GitHub: `3ce2d085cc6669e694a953190fe995b2fd75363efaad3164a55907ef5cf2a79b`
- Manifest SHA-256 checks verified against the downloaded bytes:
  - `replicate_metrics.csv`: `d7e2c8028780ed69ec58f1b4f785b727bb5f0f59753dda756c293d0976a94d65` (3,046,144 bytes)
  - `summary.json`: `eb5f52758159b9b72830e47bdd7a74275283b25bde57af65fc2db61c7cc4bfa5` (9,239 bytes)

The manifest declares 250 replicates, 4,000 parent galaxies, 30 detected events per replicate, seed start 20271200, 81 H0 grid points and 401 redshift grid points. The true catalogue selection parameters are floor 0.10, midpoint 0.22 and width 0.08.

## Checks performed

1. **File integrity:** both files' byte counts and SHA-256 digests exactly match the artifact manifest.
2. **Record structure:** 1,500 records are present: 250 per each of six arms (`analytic_truth`, `oracle`, `catalogue_only`, `correct_weighting`, `shifted_transition`, `shallow_slope`). There are no duplicate (arm, replicate) keys.
3. **Pairing:** each replicate index maps to one seed and one true H0 across all six arms. This supports within-replicate paired contrasts; it does not make separate selection-condition experiments independent.
4. **Detection count:** every record reports exactly 30 detected events.
5. **Posterior-grid reconstruction:** for every row, the saved 81-point density grid was integrated using the declared H0 grid; posterior mean and 16th, 50th and 84th percentiles were independently reconstructed. Maximum absolute discrepancies from saved values were:
   - integral versus 1: **4.45e-16**
   - posterior mean: **2.85e-14 km s^-1 Mpc^-1**
   - posterior median: **1.43e-14 km s^-1 Mpc^-1**
   - lower 68% endpoint: **1.43e-14 km s^-1 Mpc^-1**
   - upper 68% endpoint: **1.43e-14 km s^-1 Mpc^-1**
6. **Aggregate reconstruction:** recomputed bias, RMSE, empirical 68% coverage, mean interval width, mean density ISE, mean effective-sample-size fraction, and mean maximum normalized weight for all six arms. Every field agrees with `summary.json` to the precision stored there (maximum displayed discrepancy zero at 12 decimal places).
7. **Paired primary contrast:** for each replicate, define
   `D = abs(error_catalogue_only) - abs(error_correct_weighting)`.
   Positive D means correct weighting had lower absolute posterior-mean error for that paired replicate. The raw table gives mean D = **2.3537 km s^-1 Mpc^-1**; weighting has lower absolute error in **181/250 (72.4%)** pairs. A fresh percentile bootstrap with 20,000 resamples and audit seed 20261029 gives a 95% interval **[1.7767, 2.9286] km s^-1 Mpc^-1**. This is close to, but not bit-for-bit identical with, the workflow summary's reported interval [1.779, 2.923], as expected when bootstrap seeds/resampling conventions differ.

## Reconstructed arm-level results

All values are conditional on this synthetic configuration and fixed truth H0 = 70 km s^-1 Mpc^-1.

| Arm | Mean bias | RMSE | Empirical 68% coverage | Mean density ISE | Mean ESS fraction |
|---|---:|---:|---:|---:|---:|
| Analytic truth | +1.077 | 4.043 | 69.2% | 0.0000 | 1.000 |
| Oracle | +1.298 | 4.110 | 69.2% | 0.0197 | 1.000 |
| Catalogue only | -4.994 | 6.153 | 33.6% | 0.4724 | 1.000 |
| Correct weighting | +1.304 | 4.154 | 70.8% | 0.0329 | 0.808 |
| Shifted transition | +0.549 | 3.931 | 71.6% | 0.0550 | 0.878 |
| Shallow slope | -0.545 | 3.864 | 66.0% | 0.0692 | 0.907 |

Coverage is a Monte Carlo proportion, not proof of nominal frequentist calibration beyond this setup. Correct weighting substantially reduces catalogue-only bias in this configuration, but it is not best on every metric: the shallow-slope arm has lower RMSE here, while the correctly weighted arm has lower density ISE. This is a concrete example of bias-variance trade-offs and model/configuration dependence, not evidence that one weighting rule is universally optimal.

## Limits of this audit

- The posterior summaries were independently reconstructed from the saved posterior-density grids, and reported aggregate statistics were independently recomputed from the exported replicate table.
- The original simulation was **not** rerun from source during this audit; latent data generation and likelihood evaluation were not independently regenerated from the raw random seeds.
- The experiment is synthetic and simplified. It does not validate CHIMERA, a real LVK/GWTC analysis, detector selection modelling, or a real-world H0 estimate.
- This artifact is one of three selection-condition experiments. Its paired replicate structure must not be confused with independent replication across conditions.

**Audit outcome: PASS for file integrity, posterior-summary reconstruction, within-artifact pairing, and reported aggregate consistency. No conclusion is made here about physical validity or external generalization.**
