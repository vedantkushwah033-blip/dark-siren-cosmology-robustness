# Prior-support × detection-threshold factorial audit: CI artifact review

## Status and scope

The compact GitHub Actions experiment completed successfully on Python 3.10, 3.11, and 3.12 in run [40082401431](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40082401431). The full repository test suite also passed on the subsequent commit [faf58bdd](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/commit/faf58bddaed7a6390eb43cec73df16c0c61629e3), with runs [40082435746](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40082435746) and [40082435744](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40082435744).

This report reviews the **compact CI artifacts**: 100 replicates per truth and a 1,001-point posterior grid. It is a structural/numerical smoke test, not the full 300-replicate, 4,001-point experiment specified in the protocol. The reported estimates have appreciable Monte Carlo uncertainty and must not be presented as publication-grade final results.

## Artifact integrity

The Python 3.10 artifact was downloaded and inspected. It contained:

- `per_run.csv`: 7,200 rows;
- `summary.csv`: 72 rows;
- `paired_contrasts.csv`: 396 rows.

The per-run row count agrees with (3\text{ truths}\times100\text{ replicates}\times3\text{ priors}\times4\text{ thresholds}\times2\text{ methods}=7,200). The summary and contrast counts agree with the documented design. All numeric cells checked in the three CSVs were finite, and no expected truth, prior, threshold, or method level was missing.

The uncorrected posterior is invariant to the assumed detection threshold, as expected: that method does not use the threshold in its likelihood. Its rows are repeated across thresholds by design; they are not independent datasets. Do not count these repeated rows as additional independent evidence.

## Selected compact-run results

All values below use the compact CI sample of (n=100) repetitions per truth. Coverage is the fraction of simulated intervals containing the fixed truth.

### Correctly specified threshold (120 Mpc), baseline prior [50,90]

| True (H_0) | Method | Bias | RMSE | 68% coverage | 95% coverage |
|---:|---|---:|---:|---:|---:|
| 52 | Uncorrected | 12.950 | 13.602 | 0.00 | 0.85 |
| 52 | Corrected | 10.868 | 11.761 | 0.05 | 0.93 |
| 70 | Uncorrected | 1.969 | 5.964 | 0.95 | 1.00 |
| 70 | Corrected | 0.508 | 6.206 | 0.96 | 1.00 |
| 88 | Uncorrected | -12.757 | 14.415 | 0.05 | 0.98 |
| 88 | Corrected | -13.865 | 15.691 | 0.04 | 0.85 |

The compact run shows why multiple metrics must be retained. At (H_0=70), correction reduces the magnitude of bias but slightly increases RMSE. Near the prior boundaries, 68% coverage is very poor for both methods. These small-sample numbers should be treated as diagnostics, not precise performance estimates.

### Prior-support sensitivity under the 120 Mpc threshold

For corrected inference, the compact-run 68% coverage was:

| True (H_0) | Prior [50,90] | Prior [40,100] | Prior [30,110] |
|---:|---:|---:|---:|
| 52 | 0.05 | 0.79 | 0.80 |
| 70 | 0.96 | 0.88 | 0.79 |
| 88 | 0.04 | 0.74 | 0.83 |

Widening the prior strongly changes interval coverage in this compact sample, especially for truths near the original boundaries. It also increases interval width and therefore trades apparent precision for better boundary coverage in some configurations. A wider prior is not automatically preferable: the full decision must consider calibration, interval width, and the scientific prior information available before seeing the data.

### Threshold sensitivity is method-dependent

At fixed prior [50,90], the uncorrected posterior results are identical across assumed thresholds because its likelihood ignores the threshold. The corrected posterior changes with the assumed threshold. In the compact sample at true (H_0=52), corrected RMSE was 11.761 at the correctly specified 120 Mpc threshold and 12.601 at the 210 Mpc assumption. At true (H_0=88), the corresponding RMSEs were 15.691 and 15.019. The direction is therefore not uniform across truth values.

These comparisons are descriptive for this simulation. The full paired-bootstrap artifact provides uncertainty intervals for configuration differences; threshold assumptions must not be selected retrospectively to optimize outcomes.

## Interpretation guardrails

1. The observations are synthetic draws from a simplified low-redshift Gaussian distance model, truncated to (0<d_{obs}<120\) Mpc.
2. This is not an LVK detection pipeline. It omits SNR-dependent efficiency, detector duty cycles, waveform dependence, multiple-event hierarchical inference, redshift uncertainty, and realistic galaxy-catalogue selection.
3. The prior-boundary truths are intentionally close to the baseline bounds. Severe undercoverage is not by itself proof of a bug; it can arise from boundary truncation and the model's repeated-sampling behavior.
4. The compact CI sample is designed to catch coding and output-shape failures. Do not quote it as a final scientific result.
5. The next research run should use the protocol's full 300 repetitions and 4,001-point grid, retain the exact artifacts, and independently validate paired contrast definitions and Monte Carlo uncertainty before conclusions are finalized.

## Reproduction

Run the compact validation:

```bash
python scripts/run_prior_threshold_factorial.py --replicates 100 --seed 20261028 --grid-size 1001 --output-dir results/prior-threshold-factorial
python -m pytest tests/test_prior_threshold_factorial.py -q
```

Run the protocol's full experiment:

```bash
python scripts/run_prior_threshold_factorial.py --replicates 300 --seed 20261028 --grid-size 4001 --output-dir results/prior-threshold-factorial-full
```
