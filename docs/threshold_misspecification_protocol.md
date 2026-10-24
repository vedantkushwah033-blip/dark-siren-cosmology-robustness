# Detection-threshold misspecification audit

## Question

How sensitive is selection-corrected inference to using the wrong detection threshold in the selection normalization?

## Design

- Synthetic host redshifts: 0.02 and 0.04.
- Host mixture weights: 0.7 and 0.3.
- Distance noise: Gaussian, sigma = 20 Mpc.
- Data-generation detection interval: 0 < observed distance < 120 Mpc.
- True values: H0 = 52, 60, 70, 80, 88.
- Inference upper thresholds: 120 (correct), 150, 180, 210 Mpc (deliberately misspecified).
- Prior: [50,90].
- Repetitions per truth: 300.
- Random seed: 20261027.
- Posterior grid: 4,001 points.

For each truth, one detected observation set is generated and reused across every assumed threshold. The output has 6,000 per-run rows and 20 summary rows.

## Metrics

For corrected inference under each assumed threshold, the summary reports bias, mean absolute error, RMSE, empirical 68% and 95% coverage, and mean interval widths. This permits comparison of the correctly specified threshold with deliberately incorrect alternatives.

## Reproduction

```bash
python scripts/run_threshold_misspecification.py --replicates 300 --seed 20261027 --grid-size 4001 --output-dir results/threshold-misspecification
python -m pytest tests/test_threshold_misspecification.py -q
```

GitHub Actions runs the audit on Python 3.10, 3.11, and 3.12 and uploads the CSV outputs.

## Limitations

This is a deliberately simplified Gaussian distance-threshold model, not a physical gravitational-wave detection pipeline. The test isolates sensitivity to a single assumed threshold; it does not cover SNR selection, detector duty cycles, waveform efficiency, multiple events, redshift uncertainty, or realistic galaxy-catalogue incompleteness. A threshold that performs best in this finite simulation must not be selected after seeing outcomes and then presented as the correct real-world threshold.
