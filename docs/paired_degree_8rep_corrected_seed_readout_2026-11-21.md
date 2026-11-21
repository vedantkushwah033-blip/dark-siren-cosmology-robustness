# Corrected-seed eight-replicate paired-degree follow-up — independent readout

**Date:** 2026-11-21  
**Workflow:** 44860088934  
**Artifact:** 14691734869  
**Artifact ZIP SHA-256:** `30c629433a6e815487eebc3fcf55640314019915fad930542e8ba2204e69c07f`  
**PR:** #151, draft/open; no merge.

## Why this is the controlling artifact

An earlier artifact, 14691457237, used master seed `2026112121`. That seed duplicated the batch-2 seed used by PR #150, so it was not an independent follow-up and must not be used as independent replication evidence. The controlling rerun uses master seed `2026112122`, as specified on the corrected PR head. The earlier artifact and its performance summary are superseded.

The corrected artifact is fixed-truth (`H0=70)), not prior-predictive. Its design is 8 replicates, gamma {-1, 0, +1}, validation sizes {100, 300, 500}, N_parent=1200, 30 detected events, validation frame=500, H0 prior support [50, 90], 281-point H0 grid, quadrature order 5, and Scott multiplier 0.7.

## Independent artifact integrity and pairing audit

I downloaded the corrected ZIP and independently checked its contents:

- ZIP SHA-256 matches the artifact digest.
- Each degree root manifest lists 464 files; all 464 byte counts and SHA-256 hashes verified for each degree (928 root-manifest entries total).
- Each degree contains 24 child runs. All 48 child audit reports say PASS; all 48 child manifests independently verified, with no hash or byte-count mismatch.
- Each degree has 360/360 planned cells SUCCESS, 360 posterior-summary rows, and 45 aggregate metric rows.
- Each degree auditor reports 576 paired-contrast rows reconstructed and 24 successful child runs.
- Across degrees, 216 shared observable/truth files (9 per paired child block) have identical hashes.
- Child auditors reconstructed 360 posterior-summary rows per degree. Maximum reported absolute likelihood reconstruction error was 1.71e-13.
- Common arms A (catalogue-only), B (full-parent oracle), and C (known-probability IPW) are exactly identical across degrees, as expected.

These are artifact, lifecycle, pairing, and numerical-reconstruction checks. They are not calibration evidence.

## Corrected-seed paired performance comparison

I independently joined degree-one and degree-two posterior summaries by replicate, curvature, validation size, and arm. For each estimated-IPW arm I calculated the paired difference in squared median error, D2 minus D1, and a two-sided 95% paired t interval across the 8 replicates.

| gamma | n | Arm | RMSE D1 | RMSE D2 | Mean paired Δ squared error | 95% t interval |
|---:|---:|---|---:|---:|---:|---:|
| -1 | 100 | MAP plug-in IPW | 5.395 | 5.194 | -2.131 | [-19.672, 15.411] |
| -1 | 100 | Integrated IPW | 4.545 | 4.882 | 3.174 | [-10.257, 16.605] |
| -1 | 300 | MAP plug-in IPW | 4.332 | 4.780 | 4.082 | [-23.619, 31.782] |
| -1 | 300 | Integrated IPW | 4.084 | 4.304 | 1.842 | [-18.976, 22.661] |
| -1 | 500 | MAP plug-in IPW | 4.045 | 4.558 | 4.415 | [-15.477, 24.306] |
| -1 | 500 | Integrated IPW | 3.917 | 4.425 | 4.235 | [-14.187, 22.656] |
| 0 | 100 | MAP plug-in IPW | 4.618 | 5.316 | 6.939 | [-14.347, 28.225] |
| 0 | 100 | Integrated IPW | 4.397 | 5.256 | 8.292 | [-8.832, 25.415] |
| 0 | 300 | MAP plug-in IPW | 3.603 | 4.137 | 4.132 | [-10.539, 18.803] |
| 0 | 300 | Integrated IPW | 3.727 | 3.991 | 2.038 | [-6.540, 10.616] |
| 0 | 500 | MAP plug-in IPW | 3.426 | 4.423 | 7.826 | [-7.857, 23.509] |
| 0 | 500 | Integrated IPW | 3.522 | 4.213 | 5.350 | [-5.891, 16.591] |
| +1 | 100 | MAP plug-in IPW | 4.972 | 5.763 | 8.497 | [-27.557, 44.552] |
| +1 | 100 | Integrated IPW | 5.258 | 6.076 | 9.279 | [-25.476, 44.035] |
| +1 | 300 | MAP plug-in IPW | 4.946 | 4.839 | -1.050 | [-25.742, 23.642] |
| +1 | 300 | Integrated IPW | 5.067 | 4.984 | -0.832 | [-25.227, 23.564] |
| +1 | 500 | MAP plug-in IPW | 4.420 | 4.054 | -3.103 | [-21.154, 14.948] |
| +1 | 500 | Integrated IPW | 4.568 | 4.283 | -2.526 | [-20.044, 14.993] |

All 18 paired intervals include zero. D2 RMSE is lower in 3/9 MAP cells and 2/9 integrated-IPW cells, and higher in the remaining cells. This is not a degree-two performance advantage; the corrected seed changes the exploratory pattern relative to the superseded duplicated-seed artifact.

## Weight concentration remains the principal warning

For integrated IPW, median minimum node-wise Kish ESS by degree two was:

| gamma | n=100 | n=300 | n=500 |
|---:|---:|---:|---:|
| -1 | 2.258 | 3.944 | 6.092 |
| 0 | 2.968 | 5.439 | 7.079 |
| +1 | 6.677 | 10.972 | 18.627 |

The corresponding degree-two maximum top-5%-weight-mass medians ranged from 0.737 to 0.956. Degree-one integrated-IPW minimum node-wise Kish ESS medians ranged from 40.313 to 75.775, with maximum top-5% mass from 0.467 to 0.619. Degree-two integrated IPW therefore remains substantially more concentrated in this pilot.

## Seed-dependence and pooling restrictions

The separate broader-curvature workflow #44866644584 also uses master seed `2026112122`. Its corresponding replicate stream seeds are therefore reused under the common SeedSequence schedule. Treat the eight-replicate follow-up and broader-curvature run as related/shared-random-number sensitivity analyses, not as independent replications of one another.

Do not pool this fixed-truth [50, 90] prior-support run with PR #150's 40-replicate readout: PR #150 used a different H0 prior support [30, 110] and a different frozen design. Report the runs separately unless a protocol explicitly harmonizes the support and re-runs a matched analysis.

## Interpretation and claim limit

This corrected eight-replicate result does not establish calibration, estimator superiority, or observational validity. It does reinforce the need to investigate degree-two IPW concentration and finite-validation instability. No confirmatory prior-predictive run should be launched until the four-state cell-level lifecycle ledger, interruption-safe recovery, frozen evaluability rules, and independent reconstruction gates described in the PR's design note are implemented.
