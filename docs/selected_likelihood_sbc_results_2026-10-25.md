# Selected-likelihood PIT calibration results — 2026-10-25

## Status and provenance

A successful GitHub Actions run completed for commit `9600f60eee88f7e8efa0809be6ccd641669e6f9f`:

- Workflow run: [40216862319](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40216862319)
- Conclusion: success across Python 3.10, 3.11, and 3.12.
- Each job generated 1,000 replicates, ran the independent PIT audit, completed the remaining numerical diagnostics, and uploaded its artifact.
- Each uploaded `audit_report.json` reported `PASS`, with 1,000 replicates audited and both PIT reconstruction and summary reconstruction passed.

The successful run predates a subsequent test-only change that added extra tampering cases for host index, distance, and rejection-attempt count. That later commit exposed one incorrect expected error message in the new host-index tampering test. The test expectation was corrected in commit `3c3823f9e5756df58a7cdb4f7f2e5264c28d344e`; the full workflow for that latest commit is being rerun. The numerical results below are from the completed, successful run above, not from the pending rerun.

## Experiment specification

- Replicates: 1,000
- Seed: `20261025`
- Declared prior: uniform (H_0 in [40,100])
- Host redshifts: 0.02 and 0.04
- Host weights: 0.7 and 0.3
- Distance noise: 15 Mpc Gaussian standard deviation
- Detection interval: (50 < d_{\\rm obs} < 260) Mpc
- Posterior grid: 2,401 equally spaced points on [40,100], trapezoidal normalization
- Generation: draw (H_0) once per replicate; rejection-sample host and distance conditional on that fixed truth until the observation passes the detection window.
- Audit: replay the seeded generator protocol and independently reconstruct each PIT from the Gaussian host mixture and the selection-efficiency denominator.

## Results

The archived summary values were identical across all three Python versions.

| Quantity | Result |
|---|---:|
| Mean PIT | 0.4965598 |
| Monte Carlo standard error of mean PIT | 0.0088742 |
| PIT sample variance | 0.0787507 |
| Uniform KS statistic | 0.0248506 |
| KS p-value (descriptive only) | 0.5587675 |
| Maximum absolute PIT-decile frequency deviation from 0.1 | 0.023 |
| Mean simulated truth (H_0) | 69.6394553 |
| Expected prior mean | 70.0 |
| Truth-prior uniform KS statistic | 0.0255953 |
| Truth-prior KS p-value (descriptive only) | 0.5205726 |

PIT decile counts, from [0.0, 0.1) through [0.9, 1.0], were:

`[95, 96, 99, 108, 102, 110, 111, 88, 114, 77]`

The independent audit reported:

```json
{
  "status": "PASS",
  "n_replicates_audited": 1000,
  "independent_pit_reconstruction": "passed",
  "summary_reconstruction": "passed",
  "interpretation": "Synthetic consistency only; not astrophysical validation."
}
```

## Interpretation

These results are **consistent with the expected PIT behavior in this finite synthetic matched-model experiment**. The mean PIT is about 0.39 Monte Carlo standard errors below 0.5, and the descriptive KS test does not indicate a conspicuous departure from uniformity in this run. The same values and audit outcome across Python 3.10, 3.11, and 3.12 support computational reproducibility for the declared environment and seed.

This is not proof of exact calibration: 1,000 replicates have finite Monte Carlo power, decile checks involve multiple comparisons, and a non-significant KS result is not a model validation certificate. The audit establishes internal consistency of the recorded synthetic experiment, not astrophysical validity.

## Remaining validation work

1. Confirm the expanded tampering-test suite and complete numerical-diagnostics workflow pass on the latest commit.
2. Repeat calibration at multiple posterior-grid resolutions (for example, 1,201, 2,401, and 4,801 points) and quantify changes in PIT summaries.
3. Repeat with independent seeds and report Monte Carlo variability across runs.
4. Add stress cases near the detection-window boundaries and across distance-noise levels.
5. Only then treat this toy result as a calibrated numerical component. Real-event selection functions, catalogue completeness, and observational-systematics validation remain separate research tasks.
