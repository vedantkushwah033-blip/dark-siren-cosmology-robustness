# Selection-function misspecification: matched grid and bandwidth sensitivity (2026-10-29)

## Status and scope

**The matched sensitivity workflow completed successfully, and all seven configuration artifacts passed the independent selection-misspecification auditor.** The three-version numerical-diagnostics CI matrix also passed on Python 3.10, 3.11, and 3.12 at follow-up head \`f1173c81a37eefed21092d46d2c177b80d597515\`.

- Sensitivity workflow: [run 40795743899](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40795743899).
- Artifact: \`selection-misspecification-sensitivity-243e1bdc3210951dadcf5b2d05d2b97a6ebba4b9\`, artifact ID \`12863708515\`, 90-day retention.
- Sensitivity-run code commit: \`243e1bdc3210951dadcf5b2d05d2b97a6ebba4b9\`. The later \`f1173c8...\` commit fixes a missing parenthesis in a separate supplementary Stage B paired-analysis script; it does not change the sensitivity runner, auditor, or sensitivity results.
- Full test matrix on the follow-up head: [run 40795868301](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40795868301).

This is a controlled synthetic-methodology study. It is not a real LVK analysis, a realistic detector-selection validation, or a measurement of the real Hubble constant.

## Matched design

Seven configurations each used 30 paired replicates, 4,000 simulated parent galaxies per replicate, 30 detected GW-like events, and the same seed range **20262000–20262029**. All six arms were retained in each replicate, giving 180 raw rows per configuration. Seeds do not overlap Stage A (20261030–20261229) or the planned Stage B confirmation (20262030–20262529).

| Configuration | H₀ grid points | Redshift grid points | KDE bandwidth |
|---|---:|---:|---|
| Scott baseline | 161 | 801 | Weighted Scott rule |
| H₀ refinement | 321 | 801 | Weighted Scott rule |
| H₀ refinement | 641 | 801 | Weighted Scott rule |
| Redshift refinement | 161 | 1,601 | Weighted Scott rule |
| Bandwidth sensitivity | 161 | 801 | \(h=0.15\,SD_w\) |
| Bandwidth sensitivity | 161 | 801 | \(h=0.25\,SD_w\) |
| Bandwidth sensitivity | 161 | 801 | \(h=0.35\,SD_w\) |

Explicit factors multiply the weighted standard deviation of the **original observed sample**. They are not SciPy covariance multipliers and are not multiples of Scott's bandwidth.

Each configuration's raw CSV, summary, and manifest passed the repository's artifact auditor. I independently checked the SHA-256 hashes in all seven manifests, verified 180 rows per configuration, confirmed 30 replicates and all six arms per replicate, and recomputed the paired primary contrast from the raw rows. The comparison CSV hashes also match the sensitivity summary.

## Grid-resolution result: GO under the prespecified tolerances

The preregistered numerical acceptance thresholds are a maximum absolute posterior-mean change of 0.1 km s⁻¹ Mpc⁻¹ and a maximum absolute 68% interval-endpoint change of 0.2 km s⁻¹ Mpc⁻¹. These checks compare each refined grid with the 161-point H₀ / 801-point redshift baseline on matched seeds.

| Configuration | Maximum absolute posterior-mean change | Maximum absolute 68% endpoint change | Gate |
|---|---:|---:|---|
| H₀ grid: 321 points | 0.000232 km s⁻¹ Mpc⁻¹ | 0.003786 km s⁻¹ Mpc⁻¹ | PASS |
| H₀ grid: 641 points | 0.000290 km s⁻¹ Mpc⁻¹ | 0.004133 km s⁻¹ Mpc⁻¹ | PASS |
| Redshift grid: 1,601 points | 0.000123 km s⁻¹ Mpc⁻¹ | 0.000142 km s⁻¹ Mpc⁻¹ | PASS |

All observed grid changes are far below the prespecified tolerances. For this 30-seed sensitivity subset and the declared toy model, the reported posterior summaries are not materially affected by these H₀/redshift grid refinements. This does not establish convergence for every model, every event configuration, or a realistic detector likelihood.

## KDE bandwidth result: material sensitivity remains

The primary contrast is
\[
D_r=|\hat H_{0,\mathrm{catalogue\ only},r}-H_{0,\mathrm{true}}|
-|\hat H_{0,\mathrm{correct\ weighting},r}-H_{0,\mathrm{true}}|.
\]
Positive values indicate lower absolute posterior-mean error under correct weighting than under the unweighted incomplete-catalogue arm.

The Scott baseline has mean \(D=7.1263\) km s⁻¹ Mpc⁻¹, with an exploratory 10,000-resample paired-bootstrap 95% interval of [5.9766, 8.2126]. This is a small sensitivity subset, not the 500-replicate Stage B confirmation.

| Configuration | Mean primary contrast | Change vs Scott baseline | Paired-bootstrap 95% interval for change | Maximum absolute H₀ posterior-mean change across arms | Maximum absolute 68% endpoint change |
|---|---:|---:|---:|---:|---:|
| Scott baseline | 7.1263 | Reference | — | — | — |
| \(h=0.15\,SD_w\) | 7.2343 | +0.1080 | [+0.0151, +0.2077] | 0.5314 | 0.6567 |
| \(h=0.25\,SD_w\) | 7.0248 | −0.1015 | [−0.1474, −0.0596] | 0.4821 | 0.5343 |
| \(h=0.35\,SD_w\) | 6.6754 | −0.4510 | [−0.6821, −0.2374] | 1.3997 | 1.7955 |

These results expose a real methodological distinction:

- **Grid refinement is stable under the chosen tolerances.**
- **Bandwidth changes are not negligible.** The broadest tested factor shifts some posterior means by as much as 1.40 km s⁻¹ Mpc⁻¹ and some interval endpoints by as much as 1.80 km s⁻¹ Mpc⁻¹ relative to the Scott baseline.
- In the correct-weighting arm, mean RMSE changed by −0.0091, +0.0054, and +0.0737 km s⁻¹ Mpc⁻¹ for factors 0.15, 0.25, and 0.35, respectively. Mean density ISE changed by −0.0171, +0.0073, and +0.0566, respectively.
- The paired bootstrap intervals for the primary contrast's change exclude zero for these three factor comparisons in this small subset. This is a warning that downstream results can depend on bandwidth, not a universal ranking of bandwidth choices. These are exploratory comparisons across several settings, and no confirmatory claim should rest on them alone.

The bandwidth sensitivity therefore **does not pass a “results are bandwidth-invariant” interpretation**. The corrected Scott default remains the declared primary convention, but it must not be described as uniquely optimal. The sensitivity results should be shown alongside it, and the broad-factor changes should be investigated before any publication-level performance claim.

## Decision and next gates

1. **Numerical grid gate: PASS** for the tested matched 30-seed subset.
2. **Artifact-integrity gate: PASS** for all seven configurations and the paired comparison outputs.
3. **Bandwidth-robustness gate: NOT PASSED as invariance.** Report the observed dependence; do not select a bandwidth after looking for the most favorable H₀ result.
4. Continue with a predeclared calibration across multiple true H₀ values, event counts, catalogue-completeness strengths, and host-population shapes.
5. Keep the external realistic-pipeline benchmark separate. The synthetic results do not validate realistic detector selection or replace the missing-input CHIMERA benchmark.

The results support a narrow claim: the grid choices tested here are numerically stable under the stated thresholds, while plausible KDE bandwidth alternatives can materially change the toy inference. They do not support a real-data H₀ claim or a general guarantee for inverse-probability weighting.
