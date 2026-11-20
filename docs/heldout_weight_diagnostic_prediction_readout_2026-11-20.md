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


## Independent source-artifact re-audit (2026-11-20)

The completed prospective workflow was run on source commit `c4867f9546e71d17555d3a4f94f425a3e20b7cf1`:
- Workflow: `44676975586`.
- Combined source/readout artifact: `14609519527`.
- Artifact ZIP SHA-256: `0048b1f1cccf8e9b243a6e3f62a5c778c187bd1ed0474ef612804d4dd43d2cb6`.

I independently checked the two frozen training artifacts (`14520551336`, seed 2026111929; `14527343221`, seed 2026111940) and the holdout artifact (`14609519527`, seed 2026111950):
- Each root manifest declared 5,708 files; every file's byte count and SHA-256 matched, with zero mismatches.
- Each stream had 4,500/4,500 SUCCESS lifecycle cells, 4,500/4,500 audit-PASS cells, and 300 successful child runs.
- Ledger, posterior-summary, and weight-diagnostic tables each contained 4,500 rows. Their exact `(replicate, gamma, validation_size, arm)` key sets matched, with no duplicate keys.
- For every child and validation size, I independently recomputed E-arm quadrature ESS as `1 / sum(q_i**2)` from the archived normalized coefficient weights. There were 900 checks per stream; the maximum absolute difference from the stored diagnostic was `8.88e-16`, with zero discrepancies above `1e-10`. All 300 child audit reports per stream were PASS.

This closes the row-key integrity and archived-weight quadrature-ESS concerns for the exact frozen artifacts. The check was independently reproduced from the archived inputs; the full current-head workflow replay remains pending separately.

## Detailed held-out result

Across the 18 predeclared cells:
- 12/18 point estimates had lower holdout MAE than the training-mean baseline.
- 5/18 unadjusted 95% bootstrap intervals excluded zero.
- 0/18 passed the predeclared criterion after Holm adjustment.

The most favorable unadjusted cell was degree-one MAP IPW at gamma=+1 and validation size 100: holdout MAE difference `-0.301`, 95% bootstrap interval `[-0.552, -0.049]`, Holm-adjusted p-value `0.198`. This does not survive the prespecified familywise procedure.

A same-seed completed replay (workflow `44675034340`, source artifact `14633605475`) produced the same qualitative readout: zero of 18 cells passed the predeclared criterion. It is a reproducibility replay, not an independent holdout.
