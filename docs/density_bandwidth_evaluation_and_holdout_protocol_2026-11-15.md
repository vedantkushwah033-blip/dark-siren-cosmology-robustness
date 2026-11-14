# Density bandwidth evaluation and holdout protocol

Status: design note for review; no bandwidth selected and no confirmatory run authorized by this document alone.

## Why this gate exists

The five-width and expanded eleven-width pilots show that the best observed physical KDE width can differ across estimator arms, catalogue-selection conditions, and validation sizes. Several minima occur at the edge of the tested grid. A boundary minimum is evidence that the candidate range has not bracketed the empirical risk curve; it is not a valid selected bandwidth.

The experiments so far are exploratory synthetic density-recovery studies. They do not justify selecting a width from the same replicates used to report performance, and they do not justify advancing to H0 propagation.

## Estimand and comparison rules

1. The target remains the normalized parent latent-intensity times the separately specified host/merger weight. Catalogue inclusion and host/merger weighting remain distinct.
2. The primary Arm E estimand remains loss of the posterior-mean density. Posterior expected node-specific loss is a secondary estimand and must not be substituted for the primary one.
3. Compare A–E at the same physical bandwidth values. Do not give an arm a more favorable evaluation grid.
4. Retain the complete risk curve for each arm, gamma/floor condition, validation size, and replicate. Report IAE, ISE, and boundary IAE.
5. Keep oracle Arm B and known-probability Arm C explicitly labelled as simulation diagnostics. They are not deployable estimators.
6. Report weight concentration and inclusion-model discrepancy separately. ESS alone is not a model-adequacy test.
7. Do not average likelihoods or introduce an H0 result in this density stage.

## Candidate grid expansion

The next grid should be broad enough to diagnose both undersmoothing and oversmoothing. A candidate engineering grid is 25 log-spaced physical redshift widths from 0.003 to 0.12, with the existing pilot grid points retained exactly where they overlap. This is a candidate for a bracketing experiment, not a selected grid or a claim that the endpoints are sufficient.

Grid adequacy must be checked separately for every primary arm and DGP condition. If the empirical minimum is at either endpoint, expand that side and rerun the bracketing experiment. Do not select an endpoint as an optimum.

The numerical redshift integration grid must also be checked for resolution. Repeat a subset at a finer redshift grid and verify that IAE/ISE rankings and values are stable to a predeclared numerical tolerance. KDE width-grid convergence and redshift-grid convergence are separate checks.

## Replicate separation

The existing 20-replicate pilots are development evidence only. They must not be recycled as confirmatory holdout replicates.

A future selection study should use disjoint, reproducible seed blocks:

- **Training block:** used only to compare candidate bandwidth rules and freeze any selected width or selector.
- **Holdout block:** untouched until the selector, estimand, candidate grid, failure handling, and reporting table are frozen.
- Within each replicate, use common random numbers across gamma conditions so gamma contrasts remain paired.
- Across replicate IDs, use disjoint population and independent-validation seed pairs.
- Validation sizes remain nested prefixes of the same independent validation frame within a replicate.
- Record every planned replicate × gamma × validation-size cell in a lifecycle ledger. Never silently omit a failed or non-evaluable cell.

A reasonable planning starting point is 40 training replicates and 100 holdout replicates per gamma condition. This is a planning proposal, not a power guarantee. Before holdout is opened, use training-only MCSE and between-replicate variability to decide whether the planned holdout count is adequate. If it is not, increase the holdout count before unblinding results.

## Separate two questions

### Question A: estimator risk curves

This is the primary scientific comparison. It does not require choosing one bandwidth. Report paired replicate-level risk curves over the complete candidate grid, with mean differences and uncertainty intervals for:

- D − A: estimated correction versus catalogue-only baseline;
- D − C: estimated correction versus known-probability oracle;
- C − B: known-probability IPW versus parent oracle;
- E − D: posterior-mean density versus MAP plug-in density.

Also report the corresponding contrasts for ISE and boundary IAE. Preserve per-replicate values so intervals can be reconstructed independently.

### Question B: operational bandwidth selection

This is a separate secondary study. If a data-driven bandwidth selector is evaluated, it must be specified before holdout evaluation and must use only observable training information. It must not minimize the known truth-density IAE on the training set; that would be an oracle selector and must be labelled separately.

A candidate selector may use a weighted leave-one-out predictive score, but it requires a separate implementation review before adoption. In particular, the scoring rule for the posterior-mean Arm E density, normalization of importance weights, boundary treatment, and comparability across A–E must be specified and tested. Do not silently assume that one cross-validation score is valid for all five arms.

If an oracle training selector is included as an upper-bound diagnostic, report it separately from every observable-data selector. It must never be described as an implementable bandwidth-selection procedure.

## Frozen analysis plan required before holdout

Before any holdout run, add a versioned machine-readable analysis plan recording:

- exact gamma/floor conditions and replicate seed schedule;
- target density and normalization convention;
- primary and secondary estimands;
- exact candidate bandwidth list and physical units;
- redshift integration grid and numerical convergence tolerance;
- whether the analysis reports risk curves only or also selects a bandwidth;
- selector formula and all tie-breaking rules, if any;
- primary contrasts and uncertainty method;
- treatment of failed, non-evaluable, and incomplete cells;
- required artifact and independent-audit fields;
- code commit and environment versions.

The analysis plan and source commit must be frozen before holdout artifacts are generated or inspected.

## Independent audit requirements

The holdout artifact must be independently reconstructed from raw observable/truth-separated child inputs. The audit must verify:

- all planned replicate × condition keys are present exactly once;
- seed pairing and separation are correct;
- observable validation records contain no truth inclusion probabilities;
- every successful child can be independently regenerated and refitted;
- aggregate metric and diagnostic tables equal the reconstructed child rows;
- paired contrasts and uncertainty summaries equal independently reconstructed values;
- manifests and file hashes match;
- all failures and non-evaluable cells remain visible in the ledger;
- no training or holdout data were reused across the split.

A refreshed manifest hash must not be sufficient to make a corrupted aggregate pass.

## Stop/go gates

Do not start H0 propagation until all of the following are true:

1. The current seven-gamma pilot is complete and its independent audit passes.
2. The expanded physical-width grid brackets the relevant risk curves, or the endpoint behavior is explicitly reported as unresolved.
3. Redshift-grid and quadrature numerical convergence checks pass.
4. The density estimands and any operational selector are frozen before holdout.
5. The holdout evaluation is independently audited.
6. The inclusion-model discrepancy and weight-instability diagnostics are both reported.
7. The conclusions remain within the synthetic low-redshift Euclidean scope.

Passing these gates would authorize a separate H0-stage design review, not an automatic H0 run. The later H0 study must use a shared H0 parameter across events and the complete joint selection-conditioned likelihood, with catalogue inclusion and GW detection selection treated as distinct mechanisms.

## Current decision

**NO BANDWIDTH SELECTED. NO HOLDOUT OPENED. NO H0 PROPAGATION.**

The immediate task remains to complete and audit the current seven-gamma paired density pilot. The 3-replicate pilot is an engineering and design diagnostic, not confirmatory evidence.
