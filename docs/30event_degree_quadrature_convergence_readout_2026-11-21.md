# 30-event inclusion quadrature-convergence readout (2026-11-21)

## Status and provenance

- Workflow: [44806809662](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44806809662) — completed successfully.
- Source commit: `6ea5216516f84777e21a156fb6362bde18790c79` on `research/30event-degree-quadrature-convergence-20261121`.
- Artifact: `degree-matched-30event-quadrature-sensitivity`, ID `14669447416`; GitHub-reported SHA-256 `d0a8817b0716b2e6168ad4d69ad27b939d1359033fe49e6567b75c75077b77d0`.
- All three Python matrix jobs (3.10, 3.11, 3.12) completed successfully. The 30-event numerical experiment is intentionally gated to Python 3.12; the 3.10/3.11 jobs skip that heavy step.
- Independent post-download check: all three curvature bundle manifests passed; 114 files per bundle (342 total) matched both declared byte counts and SHA-256 hashes. All 18 degree × quadrature-order child audit reports were PASS.

## Frozen design

One fixed synthetic realization: 1,200 parent galaxies, 30 detected events, H0 truth 70 km s^-1 Mpc^-1, prior support [30, 110], 281-point H0 grid, validation frame 500, validation sizes 100 and 500, curvature gamma in {-1, 0, +1}. Degree-one and centered degree-two inclusion models were paired at Gauss-Hermite orders 5, 7, and 9. The same master and validation seeds were used across the curvature conditions. Each curvature bundle independently checked 35 paired input files and audited six child runs.

This is a single-realization numerical sensitivity pilot. It is not repeated-sampling calibration, SBC, estimator-ranking evidence, or observational H0 inference.

## Quadrature convergence of integrated-IPW posterior E

For the more informative adjacent comparison, order 7 versus order 9:

| Inclusion degree | Maximum posterior-density L1 difference | Maximum absolute median shift (km s^-1 Mpc^-1) | Maximum absolute shift among reported 2.5/16/84/97.5% quantiles |
|---:|---:|---:|---:|
| 1 | 0.0006544 | 0.002071 | 0.003339 |
| 2 | 0.0003966 | 0.001670 | 0.000831 |

Across the six degree × curvature × validation-size strata per degree, the 7→9 differences are small in this realization. The 5→7 comparison is less stable in some n=100 strata: for example, degree one at gamma=0 has density L1 0.004187 and median shift 0.01467; degree two at gamma=0 has density L1 0.003105 and median shift 0.008866. Thus order 5 should not be assumed adequate from this pilot; orders 7 and 9 agree much more closely here. This does not establish a universal quadrature order or a formal convergence guarantee.

## Degree-one versus degree-two integrated-IPW posterior medians

Order-9 results (km s^-1 Mpc^-1):

| gamma | Validation n | Degree 1 median | Degree 2 median | Degree 2 − degree 1 |
|---:|---:|---:|---:|---:|
| -1 | 100 | 73.9830 | 76.0009 | +2.0178 |
| -1 | 500 | 72.4497 | 76.5627 | +4.1130 |
| 0 | 100 | 77.7068 | 77.2672 | -0.4396 |
| 0 | 500 | 76.0637 | 78.2980 | +2.2343 |
| +1 | 100 | 77.5525 | 76.3757 | -1.1768 |
| +1 | 500 | 76.3814 | 76.7503 | +0.3690 |

The degree-two minus degree-one shift changes sign across curvature and validation-size settings and ranges from about -1.18 to +4.11. This single realization therefore gives no consistent direction of change and cannot support a claim that degree two improves H0 inference. The true H0 is 70, but these fixed-realization posterior medians are not coverage or bias estimates.

## Weight-concentration warning

At degree two, order 9, integrated-IPW arm E, the minimum node-wise Kish ESS across coefficient-posterior nodes is:

| gamma | n=100 | n=500 |
|---:|---:|---:|
| -1 | 2.84 | 4.50 |
| 0 | 1.24 | 2.18 |
| +1 | 1.57 | 3.99 |

The maximum node-level top-5%-weight mass ranges from 0.961 to 0.9999. Coefficient-posterior quadrature ESS is approximately 41.6–42.6 of 729 nodes. The very low minimum node-wise Kish ESS and near-total top-5% mass concentration are substantive instability warnings: numerical convergence of the quadrature integral does not imply stable IPW weights or reliable H0 inference.

## Interpretation and next gate

1. The workflow and artifact-integrity checks pass, and the degree-one/degree-two quadrature comparison is now extended from the four-event engineering smoke to the 30-event synthetic design.
2. Orders 7 and 9 give close E-arm posterior curves in this one realization; order 5 can differ more noticeably at n=100.
3. Degree-two posterior shifts are curvature- and validation-size-dependent, not uniformly favorable.
4. Severe node-wise IPW concentration persists even after quadrature-order stabilization. This is a more important inferential warning than the small 7→9 numerical differences.
5. The 281-point H0 grid has not been independently refined in this experiment. The Scott multiplier is fixed at 0.7. There is only one truth and one paired realization per curvature condition.

**Next scientific gate:** first run an independent H0-grid refinement (at least 281→561→1,121, with a stricter adjacent-grid comparison) conditional on the archived KDE inputs; then run multiple independent prior-predictive truths with the refined H0 grid and a predeclared lifecycle/failure ledger. Preserve the A–E arms, independent child audits, and IPW tail diagnostics. Do not promote degree two, claim calibration, or set universal ESS cutoffs based on this pilot.
