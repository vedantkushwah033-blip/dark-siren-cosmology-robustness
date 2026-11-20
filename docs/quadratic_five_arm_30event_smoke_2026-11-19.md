# Quadratic-logit five-arm 30-event H0 smoke

## Design

This is the first end-to-end H0 propagation smoke that joins the quadratic-logit inclusion DGP to the matched A-E KDE likelihood.

- Curvature: gamma = -1, 0, +1; beta = -3.9.
- Parent size: 1,200; 30 detected events.
- Fixed H0 truth: 70 km s^-1 Mpc^-1; posterior prior support [50, 90].
- Validation frame: 500; nested SRS sizes 100, 300, 500.
- Five arms: catalogue-only, full-parent oracle, known-inclusion IPW, degree-one MAP plug-in IPW, and coefficient-posterior-integrated IPW.
- All arms share event realization, host-rate weighting convention, redshift grid, H0 grid, reflected KDE rule, and selected-event likelihood.
- Exact true quadratic inclusion probabilities enter only the diagnostic known-π arm. The fitted MAP and integrated arms use validation observations only.
- Every gamma cell is independently audited from saved raw inputs and persisted coefficient posterior nodes.

## Limitations

One paired parent/event/validation realization per curvature is an engineering smoke, not a repeated-sampling experiment. The 281-point H0 grid and Scott multiplier 0.7 are fixed engineering settings, not a training/holdout-selected bandwidth protocol. The independent auditor reconstructs the A-D/E KDEs, selected likelihoods, and posterior summaries but uses stored coefficient posterior nodes; it does not independently validate the inclusion-coefficient posterior quadrature itself.

This is synthetic low-redshift Euclidean methodology only. It is not an observational H0 measurement, LVK/GWTC analysis, real galaxy survey validation, joint SBC, or evidence of universal estimator superiority.


## First completed 30-event smoke readout

Workflow `44440251039` completed successfully; artifact `14489939342` contains three gamma cells and 45 A-E posterior summaries. The aggregate manifest and all three child manifests were checked in the artifact, and each cell's independent likelihood reconstruction passed (maximum absolute log-likelihood difference below `2.0e-13`). This is still one realization per curvature, not a repeated-sampling result.

Posterior-median errors (km s^-1 Mpc^-1), with truth 70:

| Curvature | Validation n | A catalogue only | B parent oracle | C known-π IPW | D MAP IPW | E integrated IPW |
|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | -15.01 | +1.25 | -1.60 | -4.79 | -3.97 |
| -1 | 300 | -15.01 | +1.25 | -1.60 | -5.56 | -5.25 |
| -1 | 500 | -15.01 | +1.25 | -1.60 | -5.45 | -5.26 |
| 0 | 100 | -14.72 | +1.25 | +3.13 | +2.69 | +3.09 |
| 0 | 300 | -14.72 | +1.25 | +3.13 | +1.94 | +2.29 |
| 0 | 500 | -14.72 | +1.25 | +3.13 | +1.98 | +2.21 |
| +1 | 100 | -13.88 | +1.25 | +2.37 | +4.59 | +4.32 |
| +1 | 300 | -13.88 | +1.25 | +2.37 | +4.24 | +4.16 |
| +1 | 500 | -13.88 | +1.25 | +2.37 | +4.09 | +4.07 |

The catalogue-only arm is strongly low in all three curvature cells. The known-π arm varies with curvature, and the degree-one fitted arms move toward or away from the truth depending on curvature. In this one realization, D/E are closer to C under gamma 0, but have larger absolute errors than C under gamma +1. Under gamma -1, they are worse than C but remain substantially closer than A. This is a useful stress-test illustration, not an estimator ranking.

Weight concentration is severe in several cells. For example, at gamma 0 the known-π Kish ESS is about 101 out of 621 catalogue objects; MAP ESS ranges about 117–140, while the integrated arm's minimum node-wise Kish ESS ranges about 17–54. Its coefficient quadrature ESS is only about 6.8–7.0 of 25 nodes, and the maximum node-level top-5% weight mass ranges about 0.50–0.73. These diagnostics make it essential to report weight concentration alongside posterior summaries; coefficient integration does not by itself remove IPW instability.

The A/B/C results are identical across validation sizes because those arms do not use the fitted validation inclusion model. D/E change with validation size. The 281-point H0 grid, 501-point redshift grid, and Scott multiplier 0.7 remain smoke settings, not confirmed numerical or bandwidth choices. A multi-replicate 30-event study should wait for a dedicated grid and bandwidth protocol.
