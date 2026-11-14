# Quadratic-logit curvature pilot: paired contrast protocol

**Date:** 2026-11-15  
**Stage:** synthetic density recovery only  
**Status:** implementation and focused tests in progress; no completed seven-gamma result is claimed here.

## Scientific question

How does a degree-one logistic inclusion model behave when the true catalogue-inclusion mechanism contains a quadratic term, and can the effect be separated from changes in the average catalogue-retention rate?

The DGP uses
[
x=(z-0.205)/0.195,qquad
pi_\gamma(z)=\operatorname{expit}\{\alpha_\gamma-3.9x+\gamma(x^2-c)\},
]
where (c=E_g[x^2]) is defined under the parent latent-galaxy intensity measure. The intercept is recalibrated for every (gamma) so that the parent-average inclusion rate is held fixed. The planned values are
[
\gamma\in\{-2,-1,-0.5,0,0.5,1,2\}.
]
At (gamma=0), the degree-one logistic analysis model is correctly specified. Nonzero values create structural curvature misspecification for that analysis model.

## Pairing and independent replication

Within a replicate, the same population and validation seeds are used for every gamma value. This common-random-number design pairs the parent redshift frame, inclusion uniforms, validation redshift frame, and validation retention uniforms across gamma. Replicate seed pairs are disjoint. The gamma-specific inclusion probabilities and retained indicators still differ, as required by the DGP.

The outer ledger has one row per replicate × gamma cell. A child is counted as successful only after its independent raw-input auditor passes. Failed cells remain in the lifecycle ledger and are not silently omitted from the planned design. The outer auditor independently checks the ledger, seed schedule, child artifacts, aggregate metrics, and paired contrasts.

## Primary estimands

For each replicate, validation size, physical bandwidth, and metric:

1. **Gamma-versus-zero paired contrast**
   [
   \Delta_{\gamma,0}=L(\widehat f_\gamma)-L(\widehat f_0).
   ]
   This is computed separately for the catalogue-only, parent-oracle, known-probability IPW, MAP estimated-IPW, and posterior-mean-density arms. Negative values mean lower loss in that gamma cell than in its paired gamma-zero cell; they do not establish that misspecification is beneficial in general.

2. **Within-gamma estimator contrasts**
   - D − A: MAP estimated IPW minus catalogue-only;
   - D − C: MAP estimated IPW minus known-probability IPW;
   - C − B: known-probability IPW minus parent oracle;
   - E − D: loss of posterior-mean density minus MAP plug-in density.

The primary Arm E estimand is loss of the posterior-mean density (\bar f_E=\sum_k q_k f_{E,k}). The posterior expected node-specific loss (\sum_k q_k L(f_{E,k})) is a separate secondary estimand and is not substituted for the primary Arm E result. Node-specific rows are retained for diagnostics, but are not treated as extra independent replicates.

Metrics are IAE, ISE, and boundary IAE. The paired Monte Carlo standard error is (s_\Delta/\sqrt{n_{\rm pairs}}), with sample standard deviation across replicate-level paired differences. The summary reports the actual number of available pairs. With only three pilot replicates per gamma, uncertainty summaries are preliminary and should not be interpreted as high-precision inference.

## Bandwidth and analysis restrictions

The physical bandwidth grid is shared across arms and gamma values. No bandwidth is selected from the same replicate set used for evaluation. The pilot does not include a disjoint training/holdout selection stage, so all bandwidth-specific contrasts are descriptive sensitivity results.

Arm C uses truth-only quadratic inclusion probabilities and is an oracle diagnostic, not an implementable estimator. Arms D/E fit the degree-one logistic model from independent validation outcomes. Catalogue inclusion, host/merger weighting, and gravitational-wave detection selection remain distinct.

## Interpretation boundaries and next gates

- Results apply only to this synthetic low-redshift Euclidean setup.
- No real galaxy survey, LVK/GWTC data, or observational (H_0) inference is involved.
- No H0 posterior, PIT, coverage, or calibration result is produced by this density pilot.
- A successful artifact audit establishes internal reproducibility of the declared simulation and aggregation, not model validity or astrophysical validity.
- After the pilot, inspect gamma-by-validation-size-by-bandwidth surfaces and paired contrasts. Freeze any bandwidth protocol before running independent held-out seeds. Only after density recovery and its audit gates pass should density-to-H0 propagation be considered.
