# Prospective held-out weight-diagnostic prediction readout — 2026-11-20

## Design and integrity

The prospective evaluator used two frozen independent training streams (master seeds 2026111929 and 2026111940; 200 replicates per cell) and one untouched holdout stream (seed 2026111950; 100 replicates per cell). It evaluated 18 prespecified combinations of curvature, validation size, and arm. The primary outcome was holdout MAE difference for a diagnostic-only OLS predictor versus a training-mean baseline; negative values favor the diagnostic.

Workflow run: `44676975586`; artifact: `14609519527`; artifact ZIP SHA-256: `0048b1f1cccf8e9b243a6e3f62a5c778c187bd1ed0474ef612804d4dd43d2cb6`.

The workflow completed successfully. It generated the untouched third stream, independently audited all 4,500 holdout lifecycle cells, ran the frozen evaluator and uploaded the source and prediction output. The evaluator's own result JSON reports PASS and confirms 18 prespecified tests, 3,000 bootstrap resamples per cell, and Holm adjustment across the 18-cell family.

## Result

**Zero of 18 cells met the predeclared predictive-value criterion after the familywise procedure.**

Some unadjusted cells show apparent MAE improvement. For example, degree-one MAP-IPW at gamma=-1 and n=300 had delta MAE -0.266 (bootstrap interval [-0.483, -0.059]); gamma=+1 and n=100 had delta MAE -0.301 (interval [-0.552, -0.049]). However, neither survives the 18-test Holm adjustment at the predeclared level. Other cells show negligible or adverse changes, especially near zero curvature.

This means the diagnostic has not demonstrated reliable prospective predictive value in this held-out synthetic setting. It does not prove that all weight diagnostics are useless; the test is limited to one DGP family, one fixed H0 truth, the selected diagnostic predictors, and this train/holdout split.

## Interpretation

Do not turn this result into a universal ESS or top-weight threshold. The earlier exploratory binary screen's moderate discrimination was not enough to justify an operational rule, and this stricter continuous-outcome holdout does not validate a general-purpose predictor.

The scientifically appropriate conclusion is a negative validation result for the prespecified diagnostic-only prediction procedure in this benchmark. Future diagnostic work should be motivated by a new mechanistic hypothesis and evaluated on a new untouched DGP/seed family, not by tuning thresholds on this holdout.
