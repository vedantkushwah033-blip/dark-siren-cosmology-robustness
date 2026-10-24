# Cosmological distance-model robustness: results and independent audit (2026-10-24)

## Status

**Experiment completed; artifact independently audited.** The GitHub Actions run succeeded and produced the raw results, summary table, paired contrasts, and manifest. An independent recomputation from the downloaded artifact reproduced every reported summary statistic and every paired-bootstrap interval to floating-point precision.

- Source commit: `66430038b3d77102c9cf51714296c00b64b6dd14`
- Workflow run: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40115709336
- Artifact: `cosmological-model-robustness-results`
- Artifact SHA-256: `0cb8c5fbfe234ce75d6752460e7708b73f8ae8708d5a29baccc5aa66299b8fde`
- Seed: 20261024; paired bootstrap seed: 20273369
- Runtime environment: Python 3.11.17, NumPy 2.4.6, SciPy 1.17.1

## Question and design

This experiment asks how inference changes when synthetic luminosity-distance estimates are generated under flat ΛCDM but analyzed either with the matched flat-ΛCDM distance relation or with the deliberately mismatched low-redshift approximation (d_L\approx cz/H_0).

The design has three true values of (H_0) (60, 70, 80 km s⁻¹ Mpc⁻¹), six redshifts (0.02, 0.05, 0.10, 0.20, 0.40, 0.60), and 200 replicates per truth/redshift cell. Each generated observation is analyzed under both inference models, giving 7,200 posterior summaries (36 model-by-truth-by-redshift groups, 200 each). Matter density is fixed at Ωm = 0.3. The Gaussian distance-noise scale is max(20 Mpc, 5% of the true distance). The uniform (H_0) prior is restricted to [40, 100], evaluated on a 2,001-point grid spanning [25, 125].

This is a single-host Gaussian-distance diagnostic. It does **not** include a detector-level gravitational-wave likelihood, a realistic galaxy catalogue, a calibrated detection-selection function, uncertainty in Ωm, or the full hierarchical structure of an observational dark-siren analysis.

## Main findings

### 1. The deliberate low-redshift approximation becomes increasingly misspecified as redshift increases

Averaged across the three truths and 200 replicates per truth/redshift, the mismatched model's absolute error and coverage were:

| Redshift | Mean absolute error in H₀ | 68% interval coverage | 95% interval coverage |
|---:|---:|---:|---:|
| 0.02 | 8.43 | 78.8% | 98.3% |
| 0.05 | 5.09 | 67.8% | 94.5% |
| 0.10 | 4.78 | 36.7% | 74.8% |
| 0.20 | 8.37 | 7.0% | 33.3% |
| 0.40 | 14.21 | 0.17% | 2.17% |
| 0.60 | 18.64 | 0.0% | 0.17% |

Units for absolute error are km s⁻¹ Mpc⁻¹. These empirical coverage rates are Monte Carlo estimates from this finite design, not universal probabilities.

### 2. Matched flat-ΛCDM inference remains much better calibrated across the tested redshifts

For the matched model, averaged across truths and replicates at each redshift:

| Redshift | Mean signed bias in H₀ | Mean absolute error in H₀ | 68% interval coverage | 95% interval coverage |
|---:|---:|---:|---:|---:|
| 0.02 | +2.17 | 8.46 | 78.0% | 98.2% |
| 0.05 | +1.27 | 5.05 | 66.5% | 95.0% |
| 0.10 | +0.56 | 2.91 | 66.0% | 93.7% |
| 0.20 | +0.52 | 2.80 | 69.3% | 95.2% |
| 0.40 | +0.72 | 2.89 | 68.8% | 94.5% |
| 0.60 | +0.46 | 2.99 | 64.3% | 94.0% |

The relatively poor precision at (z=0.02) reflects the intentionally large minimum noise scale of 20 Mpc and the weak distance information at very low redshift in this setup. It should not be read as evidence that matched inference is generally biased upward.

### 3. Paired absolute-error contrasts show a redshift-dependent deterioration

The contrast is defined as mean absolute error under mismatched low-(z) inference minus mean absolute error under matched flat-ΛCDM inference. Positive values therefore favor matched inference. Each cell uses 200 paired observations and a percentile interval from 5,000 paired bootstrap resamples.

| Redshift | True H₀ | Difference in mean absolute error | 95% paired-bootstrap interval |
|---:|---:|---:|---:|
| 0.02 | 60 | -0.252 | [-0.349, -0.154] |
| 0.02 | 70 | -0.027 | [-0.126, 0.073] |
| 0.02 | 80 | +0.180 | [0.094, 0.267] |
| 0.05 | 60 | -0.079 | [-0.371, 0.222] |
| 0.05 | 70 | +0.052 | [-0.283, 0.386] |
| 0.05 | 80 | +0.122 | [-0.209, 0.451] |
| 0.10 | 60 | +1.363 | [0.886, 1.799] |
| 0.10 | 70 | +1.956 | [1.421, 2.497] |
| 0.10 | 80 | +2.277 | [1.713, 2.835] |
| 0.20 | 60 | +5.098 | [4.622, 5.549] |
| 0.20 | 70 | +5.137 | [4.547, 5.717] |
| 0.20 | 80 | +6.478 | [5.790, 7.138] |
| 0.40 | 60 | +9.559 | [9.024, 10.049] |
| 0.40 | 70 | +11.857 | [11.318, 12.358] |
| 0.40 | 80 | +12.537 | [11.818, 13.232] |
| 0.60 | 60 | +12.711 | [12.199, 13.210] |
| 0.60 | 70 | +16.029 | [15.461, 16.550] |
| 0.60 | 80 | +18.211 | [17.523, 18.842] |

The bootstrap intervals are descriptive, per-cell intervals; no multiplicity adjustment was applied across the 18 contrasts. The very small differences at low redshift should not be overinterpreted.

## Independent artifact audit

The audit was run separately after downloading the artifact. It checked:

- 7,200 raw rows, 36 summary rows, and 18 paired contrasts;
- no missing values or non-finite numeric raw fields;
- all posterior medians inside the stated prior support;
- exactly two model rows for every truth/redshift/replicate key;
- identical observed distances for the two models within each pair;
- all summary values recomputed from raw rows;
- all paired mean contrasts recomputed from raw rows;
- all 5,000-draw percentile bootstrap intervals reproduced using the documented seed.

All checks passed. Maximum absolute discrepancy between stored and independently recomputed summary metrics was below (8\times10^{-15}); the maximum discrepancy for paired contrasts and interval endpoints was below (4\times10^{-15}).

## Interpretation and limitations

The experiment demonstrates, within this simulator, that using (cz/H_0) outside its low-redshift regime can induce large systematic error and severe credible-interval undercoverage. It does not establish that every low-redshift approximation is harmful in every setting, and it does not validate a real dark-siren analysis. Fixed Ωm, one host, Gaussian distance noise, a finite prior, and a finite grid are important simplifications. The paired intervals are not multiplicity-adjusted.

The next scientifically useful extension is not merely to add more repetitions. It is to integrate this distance-model comparison with the already-developed catalogue-completeness and host-weighting mechanisms, while keeping matched seeds, documenting shared random draws, and testing the joint model with a pre-specified set of primary contrasts. A later phase should also assess grid convergence and consider marginalizing over cosmological parameters before making claims relevant to real observations.

## Reproducibility

Run the workflow above, download the artifact, and use the manifest to reproduce the design. The source implementation is in `scripts/run_cosmological_model_robustness.py` and `src/dark_siren_robustness/cosmological_inference.py`; the prespecified design and limitations are in `docs/cosmological_model_robustness_protocol.md`.
