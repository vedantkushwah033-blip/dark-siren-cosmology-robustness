# Lower-threshold sensitivity audit

**Status:** Completed successfully by GitHub Actions on 2026-10-24.  
**Workflow:** [Toy lower-threshold sensitivity audit, run 40072322309](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40072322309)  
**Artifact:** `lower-threshold-sensitivity-results`.  
**Scope:** deterministic synthetic model only.

## Question and design

We varied the lower boundary of the toy detection interval while holding the upper boundary fixed. The tested intervals were (0, 180), (10, 180), (20, 180), and (40, 180) Mpc. Other settings were fixed: host redshifts 0.02 and 0.04 with weights 0.70 and 0.30, Gaussian distance uncertainty 20 Mpc, prior H0 in [50, 90] km s^-1 Mpc^-1, and a 4,001-point H0 grid from 25 to 125. Three fixed observations (50, 90, 140 Mpc) were each analyzed with the uncorrected and selection-normalized likelihoods.

## Detection efficiency over the prior support

| Lower threshold (Mpc) | beta minimum on H0 in [50,90] | beta maximum on H0 in [50,90] |
|---:|---:|---:|
| 0 | 0.699484 | 0.996789 |
| 10 | 0.699484 | 0.995468 |
| 20 | 0.699484 | 0.990214 |
| 40 | 0.699461 | 0.952365 |

The minimum is almost unchanged because the low-efficiency edge within the prior support is dominated by the same parameter/host configuration. The maximum falls as the lower threshold rises, since more low-distance observations are excluded by the detection rule. These are efficiencies of the simplified Gaussian-mixture toy model, not detector efficiencies.

## Posterior medians

| Lower threshold (Mpc) | 50 Mpc observation: uncorrected | 50 Mpc: corrected | 90 Mpc observation: uncorrected | 90 Mpc: corrected | 140 Mpc observation: uncorrected | 140 Mpc: corrected |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 81.368 | 80.955 | 69.608 | 67.858 | 75.771 | 72.815 |
| 10 | 81.368 | 80.959 | 69.608 | 67.861 | 75.771 | 72.823 |
| 20 | 81.368 | 80.974 | 69.608 | 67.873 | 75.771 | 72.858 |
| 40 | 81.368 | 81.139 | 69.608 | 68.050 | 75.771 | 73.306 |

All H0 values are in km s^-1 Mpc^-1. The uncorrected posterior is unchanged by design because its likelihood does not use the detection interval. The corrected posterior changes with the interval because beta(H0) changes. The effects are modest for the 50 Mpc observation and larger for the 140 Mpc observation in this selected set of examples.

## Interpretation

1. **The selection boundary matters.** The corrected median is not invariant to the lower threshold, especially for the 140 Mpc example.
2. **The tested efficiencies are not near zero within the analysis prior.** The earlier value beta_min about 0.00097 came from evaluating the broader computational grid [25,125], including values outside the prior. For this audit the efficiency range is explicitly restricted to [50,90], where its minimum is about 0.699.
3. **These fixed observations do not assess estimator bias or coverage.** No true H0 was assigned to these three observations; posterior shifts alone cannot establish which method is more accurate.
4. **Implementation check only.** The audit shows sensitivity to the selection interval in a controlled model. It does not validate a realistic gravitational-wave detection model.

## Validation and limitations

The audit produced 24 finite posterior summaries and its focused test passed. This is a deterministic sensitivity analysis; repeated-sampling coverage is addressed separately by the fixed-truth coverage-map experiment.

The model lacks a real search statistic, injection campaign, detector network sensitivity, duty cycle, source-population effects, cosmological volume terms, and a realistic galaxy catalogue. Do not interpret the results as an LVK result or as a measurement of the real Hubble constant.
