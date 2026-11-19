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
