# Selection-function misspecification: audited pilot and three-condition stress screen (2026-10-29)

## Status and provenance

The 200-replicate Stage A pilot and all three predeclared stress-screen conditions completed successfully on commit `4a95ff40b68849e0ff33c9cae8916a3e96eeaa7b`. The workflows each ran the independent artifact auditor and uploaded raw replicate metrics, a summary, and a manifest.

- Stage A pilot workflow: [run #40811157211](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811157211); artifact ID `12869859795`.
- Three-condition stress screen: [run #40811159599](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811159599).
- Low-inclusion-floor artifact ID `12870297736`; summary SHA-256 `3b60f6151ab2d34a3cc8fbb18c5fcc4cd72664bb49bed944962f349c5c297b0f`.
- High-inclusion-floor artifact ID `12870541126`; summary SHA-256 `[recorded in artifact manifest]`.
- Broad-transition artifact ID `12870855101`; summary SHA-256 `ec2856409994a36a33a2698e851c511ae530d36f9cb4de68e54fd6cd6841638c`.

These are controlled synthetic, low-redshift Euclidean toy experiments. They are not real LVK observations, not a real Hubble-constant measurement, and not evidence that the method is validated for detector data.

## Design

Each stress condition used 50 paired replicates, 4,000 parent galaxies per replicate, 30 detected events, true H0 = 70 km s^-1 Mpc^-1, an 81-point H0 grid over [50, 90], and 401 redshift-grid points over [0.01, 0.4]. The detection threshold was 250 Mpc and the Gaussian distance-noise scale was 140 Mpc. The six arms were analytic-truth density (diagnostic control), oracle KDE using the simulated parent catalogue, incomplete catalogue without correction, inverse-inclusion weighting under the declared true selection law, a shifted transition midpoint, and a shallower assumed transition.

The three stress settings were:
- Low floor: true inclusion floor 0.05, midpoint 0.22, width 0.05; seed start 20270000.
- High floor: true inclusion floor 0.25, midpoint 0.22, width 0.05; seed start 20270100.
- Broad transition: true inclusion floor 0.10, midpoint 0.22, width 0.08; seed start 20270200.

The Stage A pilot used 200 replicates, seed start 20261030, and the baseline true selection parameters (floor 0.10, midpoint 0.22, width 0.05). All conditions use the corrected weighted Scott bandwidth by default. Explicit bandwidth factors have a separate documented meaning and were not used here.

## Main results

Bias is mean posterior mean minus true H0. Coverage is empirical coverage of nominal 68% posterior intervals. RMSE is in km s^-1 Mpc^-1. Coverage estimates from 50-replicate stress conditions are imprecise; the 200-replicate pilot is more informative but remains a toy-model experiment.

### Stage A pilot: 200 replicates

| Arm | Bias | RMSE | 68% coverage | Mean 68% interval width | Mean weight ESS fraction |
|---|---:|---:|---:|---:|---:|
| Analytic truth | +0.069 | 3.754 | 136/200 = 68.0% | 7.223 | 1.000 |
| Oracle KDE | +0.290 | 3.771 | 137/200 = 68.5% | 7.245 | 1.000 |
| Catalogue only | -9.468 | 9.984 | 15/200 = 7.5% | 8.654 | 1.000 |
| Correct weighting | +0.351 | 3.785 | 142/200 = 71.0% | 7.265 | 0.658 |
| Shifted transition | -0.409 | 3.756 | 139/200 = 69.5% | 7.653 | 0.785 |
| Shallow slope | -1.824 | 4.091 | 127/200 = 63.5% | 7.737 | 0.809 |

The preregistered paired absolute-error reduction, defined as absolute error for catalogue-only minus absolute error for correct weighting, was **6.536 km s^-1 Mpc^-1** (10,000-resample bootstrap 95% interval **[5.920, 7.119]**). This is a within-model comparison under the known synthetic selection law.

### Three-condition stress screen: 50 replicates per condition

| True selection setting | Arm | Bias | RMSE | 68% coverage | Mean 68% interval width | Mean weight ESS fraction |
|---|---|---:|---:|---:|---:|---:|
| Low floor (0.05) | Catalogue only | -11.096 | 11.470 | 1/50 = 2% | 8.380 | 1.000 |
| Low floor (0.05) | Correct weighting | +0.830 | 4.301 | 33/50 = 66% | 7.440 | 0.557 |
| Low floor (0.05) | Shifted transition | -0.884 | 4.199 | 34/50 = 68% | 8.005 | 0.742 |
| Low floor (0.05) | Shallow slope | -2.347 | 4.655 | 27/50 = 54% | 7.988 | 0.770 |
| High floor (0.25) | Catalogue only | -4.910 | 5.900 | 19/50 = 38% | 8.520 | 1.000 |
| High floor (0.25) | Correct weighting | +0.760 | 3.680 | 36/50 = 72% | 7.270 | 0.830 |
| High floor (0.25) | Shifted transition | +0.840 | 3.700 | 37/50 = 74% | 7.420 | 0.880 |
| High floor (0.25) | Shallow slope | -0.250 | 3.540 | 34/50 = 68% | 7.540 | 0.900 |
| Broad transition (width 0.08) | Catalogue only | -5.928 | 6.740 | 12/50 = 24% | 8.208 | 1.000 |
| Broad transition (width 0.08) | Correct weighting | +0.155 | 3.641 | 37/50 = 74% | 7.242 | 0.807 |
| Broad transition (width 0.08) | Shifted transition | -0.568 | 3.612 | 37/50 = 74% | 7.422 | 0.877 |
| Broad transition (width 0.08) | Shallow slope | -1.640 | 3.841 | 33/50 = 66% | 7.511 | 0.907 |

Analytic-truth and oracle controls were also recorded for every condition. Their purpose is diagnostic: analytic truth bypasses density estimation; oracle KDE uses the full simulated parent population, which would not be observable in practice.

The paired absolute-error reduction from catalogue-only to correct weighting was:
- Low floor: **7.944 km s^-1 Mpc^-1** (bootstrap 95% interval **[6.490, 9.265]**).
- High floor: **2.280 km s^-1 Mpc^-1** (bootstrap 95% interval **[1.122, 3.391]**).
- Broad transition: **3.486 km s^-1 Mpc^-1** (bootstrap 95% interval **[2.341, 4.479]**).

Positive values favor correct weighting under this estimand. These bootstrap intervals describe variation across the paired synthetic replicates; they do not include model-form uncertainty or make the synthetic assumptions realistic.

## Interpretation

1. In all three stress settings, ignoring the known catalogue-inclusion mechanism was associated with negative H0 bias and larger RMSE than the corrected-weighting arm.
2. The apparent benefit varied materially with the inclusion floor: the paired absolute-error reduction was largest in the low-floor condition and smallest in the high-floor condition. This is consistent with the severity of catalogue incompleteness affecting the value of correction in this toy setup.
3. Misspecifying the selection transition was not harmless. The shallow-slope arm had worse RMSE than correct weighting in the low-floor and broad-transition conditions; the direction and size of the effect varied with the true setting.
4. A finite 50-replicate stress screen is for robustness triage, not a precise coverage study. The coverage proportions should not be overinterpreted; Wilson intervals are wide at n=50.
5. The pilot suggests that the analytic-truth and oracle arms behave similarly to correct weighting under the baseline toy law, but they are controls rather than practical alternatives. Correct weighting assumes the true selection law is known; real applications must estimate or otherwise justify it.
6. These results do not address realistic gravitational-wave detector selection, event-rate marginalization, peculiar velocities, redshift uncertainty, non-Euclidean distance-redshift relations, or real galaxy-catalogue systematics.

## Reproduction

The exact commands and configurations are preserved in the workflow manifests. The pilot command was:

```bash
python scripts/run_selection_misspecification.py --replicates 200 --n-parent 4000 --n-detected 30 --seed 20261030 --h0-points 81 --z-points 401 --catalogue-floor 0.10 --catalogue-midpoint 0.22 --catalogue-width 0.05
```

The stress-screen workflow used the same grid and sample sizes, with 50 replicates per condition and the seeds/settings listed above. Download the workflow artifacts to inspect raw replicate records and manifests. Each workflow's independent audit step succeeded.

## Current verification boundary

At report preparation time, the selection pilot and three stress conditions had succeeded with audit steps. The full multi-version test workflow for the same commit was still running; do not describe the full suite as passed until its final result is successful. The external CHIMERA benchmark remains blocked pending verification of a version-compatible input bundle. No real-data inference has been performed.
