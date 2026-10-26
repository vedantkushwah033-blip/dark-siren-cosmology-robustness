# Full distance-uncertainty sensitivity run: results and numerical checks

## Run provenance

The default experiment was run in GitHub Actions on Python 3.12.15 with NumPy 2.5.4 and SciPy 1.18.1, seed `20261026`, 300 replicates per truth/fractional-error cell, a 1,601-point grid over (H_0\in[50,90]), and 2,000 paired-bootstrap resamples. The run produced 32,400 raw posterior rows, 108 summary cells, 10,800 paired A/B replicate comparisons, and 180 paired-bootstrap contrasts.

The independent auditor returned `PASS`: it reconstructed all 32,400 posterior rows, all 108 cell summaries, the A/B pair structure, and all 180 paired-bootstrap intervals from the raw CSV. It also replayed the recorded random seed and loop order. The audit identifies the tested workflow revision as `3b7efd516d7c9ca1a096c529369416f250747e2b`.

- Raw CSV SHA-256: `576d0cf61f9d89362133f9c607f3476d775b2477525a3f02e8d039d9c816c722`
- Summary JSON SHA-256: `5578b4b15101fde0f4d9f27457d0addb0434760efc6cee6af67a7a5a4b52b077`

The files are in the `numerical-diagnostics-py3.12` Actions artifact under `distance_uncertainty_sensitivity_full/`.

## Descriptive results

The following averages are computed by first summarizing each cell, then equally averaging the 7 *interior* truth cells (H0 = 55, 60, 65, 70, 75, 80, 85) separately at each fractional error. Boundary truths H0 = 50 and 90 are excluded from this table. These are descriptive averages over simulation cells, not a pooled estimator or a universal ranking of the arms.

| Fractional error | Arm | Mean 68% interval coverage | Mean 95% interval coverage | Mean absolute median bias (km/s/Mpc) | Mean cell RMSE of posterior median (km/s/Mpc) |
|---:|---|---:|---:|---:|---:|
| 5% | A: fixed realized sigma | 0.7057 | 0.9581 | 0.588 | 3.257 |
| 5% | B: fractional sigma(H0) | 0.7090 | 0.9595 | 0.478 | 3.231 |
| 5% | C: fixed reference sigma | 0.7048 | 0.9586 | 0.582 | 3.209 |
| 10% | A: fixed realized sigma | 0.7067 | 0.9681 | 2.020 | 6.089 |
| 10% | B: fractional sigma(H0) | 0.7105 | 0.9681 | 1.864 | 6.058 |
| 10% | C: fixed reference sigma | 0.7081 | 0.9671 | 1.914 | 5.929 |
| 15% | A: fixed realized sigma | 0.7533 | 0.9781 | 3.256 | 7.448 |
| 15% | B: fractional sigma(H0) | 0.7567 | 0.9810 | 3.236 | 7.465 |
| 15% | C: fixed reference sigma | 0.7562 | 0.9767 | 3.154 | 7.339 |
| 20% | A: fixed realized sigma | 0.7657 | 0.9805 | 4.502 | 8.344 |
| 20% | B: fractional sigma(H0) | 0.7686 | 0.9824 | 4.549 | 8.354 |
| 20% | C: fixed reference sigma | 0.7657 | 0.9776 | 4.435 | 8.334 |

In this setup, interior-cell mean coverage is above nominal at the larger noise fractions, and average absolute bias generally grows with noise. The differences between A and B are modest and vary by truth/error cell; they do not support a universal claim that either likelihood is superior. The 300-replicate cell coverage estimates have Monte Carlo uncertainty, which is why the cell-level Wilson intervals and paired-bootstrap contrasts should be consulted rather than interpreting tiny differences as decisive.

## Independent grid-resolution check

As a separate post-run check, all 32,400 raw records were recomputed on a 3,201-point H0 grid and compared with the stored 1,601-point summaries:

| Posterior summary | Median absolute change | 95th percentile absolute change | Maximum absolute change |
|---|---:|---:|---:|
| Median | 0.0000028 | 0.0000141 | 0.0001238 |
| Lower 68% bound | 0.0000080 | 0.0000294 | 0.0001121 |
| Upper 68% bound | 0.0000069 | 0.0000333 | 0.0001562 |
| Lower 95% bound | 0.0000157 | 0.0000580 | 0.0001244 |
| Upper 95% bound | 0.0000076 | 0.0000570 | 0.0001764 |

The 68% interval coverage proportion across all truths, including the boundaries, was unchanged at the reported precision for each arm when increasing the grid from 1,601 to 3,201 points. This supports numerical stability of this grid calculation for the simulated records; it does **not** validate the scientific realism of the model.

## Interpretation limits

- **Arm A is oracle-like:** its realized sigma is generated from the simulation's true H0, then held fixed over candidate H0. It is a conditional benchmark, not automatically a usable uncertainty prescription when the truth is unknown.
- **Arm C changes the generator away from H0 = 70:** comparisons with A/B combine a change in noise generation with a change in likelihood treatment.
- The model has one known host, (d=cz/H_0), Gaussian noise, and a uniform prior. It omits realistic gravitational-wave distance likelihoods, selection effects, catalogue incompleteness, host-population uncertainty, and cosmological-model complexity.
- This is a synthetic sensitivity study, not LVK-data analysis and not a new H0 measurement. The results only describe the specified simulation and should not be generalized to real dark-siren inference without a separate, realistic validation program.
