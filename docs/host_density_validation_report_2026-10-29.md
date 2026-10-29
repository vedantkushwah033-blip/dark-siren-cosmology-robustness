# Host-density validation campaign — 29 October 2026

**Status:** completed exploratory CI campaign; artifact integrity and aggregate reconstruction checked. This is a synthetic screening study, not a real-data result or a confirmatory universal claim.

## Provenance and design

- Workflow: [host-density validation campaign run #40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794)
- Source commit: `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`
- Artifact ID: `12865051001`
- Design: 13 configurations × 40 paired replicates × three arms.
- Toy settings: fixed H0 = 70 km s⁻¹ Mpc⁻¹; 6,000 parent galaxies; 40 detected events; distance-noise scale 140 Mpc; 250 Mpc detection threshold.
- Integrity: all 14 manifest-listed campaign files matched their recorded SHA-256 hashes and byte counts. The saved CSV aggregates were independently recomputed from raw rows to numerical precision (maximum discrepancy about 1.8×10⁻¹⁵). A separate independent posterior-reconstruction smoke audit passed on two fresh seeds. It did **not** independently regenerate every posterior in all 520 campaign replicates.

## Baseline Scott-bandwidth results

| Arm | Bias (km s⁻¹ Mpc⁻¹) | RMSE | Nominal 68% coverage | Density L1 |
|---|---:|---:|---:|---:|
| Complete-parent oracle | — | 3.936 | — | — |
| Catalogue only, no correction | −8.634 | 9.165 | 10.0% | 0.537 |
| Known inclusion probabilities | +1.180 | 3.896 | 55.0% | 0.082 |

These screening figures are separate from the locked 500-replicate Stage B experiment and must not be pooled with it.

The paired known-inclusion versus oracle absolute-error difference was −0.057 km s⁻¹ Mpc⁻¹, bootstrap interval [−0.199, +0.083], an inconclusive contrast. The catalogue-only arm was much worse than the oracle on this constructed population. The weighted arm's coverage was only 55% in this 40-replicate screen, so the RMSE similarity must not be presented as proof of calibrated intervals.

## Sensitivity findings

- Grid-resolution changes were tiny over the tested settings: posterior-mean changes were at most 8.51×10⁻⁴ km s⁻¹ Mpc⁻¹ for the coarse grid and 1.03×10⁻⁴ for the fine grid.
- Bandwidth changes were material. Factors 0.5–1.2 changed density and posterior summaries; at factor 1.2, catalogue-only RMSE was about 4.238 while oracle RMSE worsened to about 5.339 in the campaign summary. This is consistent with accidental bias compensation under oversmoothing, not reliable correction of catalogue incompleteness.
- Selection-curve width assumptions also changed bias and coverage. The campaign does not establish that midpoint misspecification is generally harmless or that inverse-inclusion weighting is robust when inclusion probabilities are small.
- A separate coarse-grid oracle calibration workflow reported 300 prior/fixed-truth replicates, but did not preserve row-level seeds and raw records. Treat it as weaker aggregate evidence, not a fully replayable audit.

## Interpretation and release decision

**Supported:** in this particular synthetic setup, ignoring the declared catalogue-inclusion mechanism badly distorts the estimated host density and degrades H0 error; known-inclusion weighting substantially improves the density estimate.

**Not supported:** a universal completeness correction, calibrated real-world H0 inference, superiority over Gaussian-process methods or external GW pipelines, or a new measurement of H0.

The experiment conditions on a fixed detected-event count and uses a simplified low-redshift Gaussian-distance model. It omits realistic detector/network selection, uncertain redshifts, peculiar velocities, event-rate marginalization and survey-level hierarchical population inference. Correct weighting assumes the inclusion probabilities are known exactly.

## Remaining gates

1. The exact-head Python 3.10/3.11/3.12 numerical-diagnostics workflow passed on commit `e9a512d46aeeb7070649ac8f851c1d86103f4379`: [run #40797969025](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797969025).
2. CHIMERA v2.2 remains **NO-GO** until its exact required mock and completeness inputs and compatible execution are available and verified. No external benchmark run is claimed.
3. The literature positioning is focused, not a systematic review; no external subject-matter review or peer review is claimed.

**Bottom line:** suitable as a transparent exploratory synthetic-methodology report, with these limitations explicit. Not ready to be represented as a real-data astrophysics result.
