# Matched fixed-H0 coverage audit: first results

**Successful run:** [GitHub Actions run 39995536840](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39995536840)  
**Artifact:** `matched-fixed-h0-coverage` (7,000 posterior fits; 500 replicates at each of seven fixed truths, for each of two scenarios).  
**Design:** same independently drawn distance-noise scale as the SBC diagnostic; catalogue and host generated with the pinned simulator revision.  
**Important status:** these are diagnostic results, not evidence that the model is physically realistic or ready for observational inference.

## Results

All errors, RMSE values and interval widths are in km/s/Mpc. Coverage values are fractions. Monte Carlo standard errors (MCSEs) are in `summary.csv` in the workflow artifact.

| True H0 | Scenario | Mean signed median error | RMSE | 68% coverage | 95% coverage | Mean 68% width |
|---:|---|---:|---:|---:|---:|---:|
| 55 | Known host | +2.632 | 5.452 | 71.6% | 95.6% | 8.978 |
| 55 | Unknown host | +10.110 | 13.483 | 43.0% | 99.8% | 15.772 |
| 60 | Known host | +2.072 | 5.944 | 68.8% | 95.8% | 11.097 |
| 60 | Unknown host | +7.723 | 11.197 | 55.4% | 99.6% | 17.926 |
| 65 | Known host | +1.386 | 5.952 | 69.0% | 95.0% | 12.186 |
| 65 | Unknown host | +5.129 | 8.637 | 92.0% | 99.0% | 19.715 |
| 70 | Known host | +1.212 | 6.216 | 67.8% | 96.6% | 12.852 |
| 70 | Unknown host | +1.947 | 6.300 | 90.0% | 98.6% | 20.908 |
| 75 | Known host | +0.102 | 5.807 | 73.6% | 98.0% | 13.456 |
| 75 | Unknown host | −1.087 | 4.828 | 88.8% | 98.8% | 22.037 |
| 80 | Known host | −1.908 | 6.040 | 78.6% | 98.2% | 13.021 |
| 80 | Unknown host | −4.604 | 5.583 | 91.0% | 99.6% | 22.831 |
| 85 | Known host | −4.157 | 6.277 | 77.8% | 98.2% | 11.994 |
| 85 | Unknown host | −8.886 | 9.223 | 78.8% | 98.6% | 23.126 |

## First interpretation

1. **The unknown-host model has substantial conditional bias at the lower end of the tested H0 range.** Mean signed median error is +10.11 at true H0=55 and +7.72 at H0=60, compared with +2.63 and +2.07 for known-host inference.
2. **The unknown-host 68% coverage is highly non-uniform across fixed truths.** It is 43.0% at H0=55 and 55.4% at H0=60, but 90.0–92.0% at H0=65–70. Those departures are much larger than the approximate Monte Carlo uncertainty from 500 replicates.
3. **95% coverage is often higher than nominal**, particularly in the unknown-host case. This is consistent with its wider posterior intervals, but does not resolve the poor lower-H0 68% coverage.
4. Known-host coverage also worsens toward the upper prior boundary: 68% coverage reaches 77.8% at H0=85. This makes prior-boundary and conditional-coverage diagnostics important even before interpreting host marginalization.
5. The prior-predictive SBC results were broadly nominal, while these fixed-truth results vary substantially with H0. This is not logically impossible: aggregate prior-predictive calibration can conceal conditional failures. It means the SBC pass is not sufficient to certify the estimator's conditional performance.

## Critical validation caveat

This experiment is a first diagnostic, not yet a final scientific conclusion. Before attributing the discrepancy to host marginalization, inspect the exact implementation of `posterior_h0`, the catalogue's redshift-generating distribution, the host-selection rule, and whether the inference likelihood conditions on the observed catalogue in the same way as the simulator. Verify explicitly that the true host is present in the candidate catalogue and that the likelihood for a fixed catalogue is the correctly normalized host mixture.

The fixed-H0 script currently evaluates known-host and unknown-host scenarios from the same simulated catalogue and observation. Its output should be preserved, but no causal claim should be made until the generative and inference models are audited line by line.

## Next action

The workflow now also runs a paired bootstrap comparison of unknown-host minus known-host performance by replicate. Once that updated workflow passes, add its paired intervals to the report. Then implement small deterministic tests where a catalogue contains exactly one candidate (the unknown-host posterior must reduce to the known-host posterior) and where a multi-host catalogue has a hand-computable mixture likelihood. Those tests will help isolate implementation errors before changing the scientific model.
