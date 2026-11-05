# Independent mathematical adjudication of the Consensus dossier

**Date:** 2026-11-05  
**Scope:** Adjudication of the supplied Consensus synthesis against the inspected implementation and the project's declared synthetic design. This is a methodological review, not a new literature search and not a claim that every cited source has been independently verified.

## Executive verdict

The Consensus synthesis is useful as a checklist, but it mixes correct mathematical statements with several overstatements and a few implementation-inaccurate specifications. Its most important useful conclusion is the narrowed novelty framing: the project is an underexplored application/combination of established two-phase sampling, logistic estimation, IPW, and dark-siren posterior calibration—not a new IPW estimator or a new Bayesian method.

The project should **not** adopt the dossier's P0 instruction to add an unspecified "Method C" likelihood marginalization. Nor should it adopt a universal hard probability floor, assert that a constrained boundary solution is automatically an invalid numerical artifact, or describe V2 as ridge penalized.

## 1. Logistic likelihood and optimization

For binary validation outcomes and design matrix X, the unpenalized logistic negative log likelihood is

\[
L(\beta)=\sum_i\{\log(1+e^{x_i^T\beta})-y_i x_i^T\beta\},
\quad g=X^T(p-y),\quad H=X^T W X,
\]

where W has diagonal entries p_i(1-p_i). For finite coefficients and finite predictors, every diagonal entry is positive. Therefore H is positive semidefinite, and it is positive definite exactly when X has full column rank. Convexity does not imply existence of a finite unconstrained minimizer: under complete or quasi-complete separation the infimum can be approached along a diverging coefficient sequence.

Newton's method and IRLS are equivalent representations of the same local quadratic update for ordinary unpenalized logistic regression when the relevant system is nonsingular. This equivalence should not be generalized to arbitrary box-constrained projected algorithms or to penalized estimators.

For box constraints l <= beta <= u, a KKT-consistent coordinate sign convention is:
- interior: g_j = 0;
- at lower bound: g_j >= 0;
- at upper bound: g_j <= 0.
A projected-gradient mapping is a useful convergence diagnostic, but active-set sign checks provide a valuable independent test.

### Important correction: V2 is not ridge logistic regression

The inspected V2 implementation minimizes the **ordinary unpenalized logistic negative log likelihood subject to coefficient box constraints**. Its diagonal damping is a numerical linear-system safeguard, not an L2 penalty. The V2 objective therefore has no \(\lambda\|\beta\|_2^2/2\) term, and its gradient/Hessian have no \(\lambda\beta\) or \(\lambda I\) additions.

A box-constrained optimum can be a mathematically valid optimum of the constrained problem. If separation means the unconstrained MLE does not exist, a boundary solution must be labelled as a constrained estimate and must not be presented as a finite unconstrained MLE. Calling every boundary optimum an "artificial artifact" is too broad: the box is part of the declared estimator, although it may be a scientifically undesirable constraint.

Firth bias reduction is a distinct adjusted-score/penalized-likelihood estimator. It must not be conflated with V2, ridge, or an automatically "proper MAP" interpretation. Any Firth arm would need its own specification, implementation, validation, and preregistered role.

## 2. Positivity and inverse weights

The population positivity condition is a property of the data-generating/validation design: the inclusion probability must be positive on the target population's relevant support. A stronger uniform-positivity condition requires a positive lower bound over that support. Neither condition can be proved merely by taking the minimum of finitely many fitted scores.

For the current V2 coefficient boxes and standardized x in [-1,1], the linear predictor is bounded in absolute value by 75 for the linear model and 125 for the quadratic model. Thus the permitted parameter region allows probabilities as small as expit(-75) and expit(-125), respectively. Those are theoretical extrema over the coefficient box, not observed fitted minima and not evidence that the benchmark truth has such probabilities.

The right analysis should report separately:
1. true DGP inclusion-probability distribution (available only to the simulation scorer);
2. fitted probability distribution over catalogue objects and over the declared redshift grid;
3. validation support/density diagnostics;
4. inverse-weight concentration, ESS, and tail shares;
5. reconstructed-density error and downstream H0 calibration.

A finite validation sample cannot establish positivity outside its support. Min/max range overlap is only a crude proxy; local validation density or prespecified redshift-bin support diagnostics are more informative.

No universal ESS, top-1% mass, or probability threshold proves an estimator usable. If clipping, truncation, or stabilization is studied, it must be a separately named estimator arm, with the uncapped estimator retained as a reference and the rule fixed before looking at H0 outcomes.

## 3. IPW normalization and the missing-host mixture

The supplied first-order Hájek expansion is sound under its small, bounded relative-error approximation:
\[
\hat f(z)/f(z)-1 \approx -\{\delta(z)-E_f[\delta(Z)]\}.
\]
It correctly highlights that constant multiplicative inclusion-probability error can cancel after normalization, whereas redshift-varying error changes the reconstructed shape. It does not describe all finite-sample, KDE, bandwidth, or positivity effects.

The statement "IPW requires f=0" is not a universal mathematical rule. The correct condition is model-specific: if the density supplied to the event numerator and selection normalizer is already the reconstructed full parent-host density, an additional hidden-host mixture correction for the same catalogue incompleteness would double-correct that incompleteness. The inspected estimated-inclusion pilot path has no separate hidden-host mixture fraction. This is a scoped source finding, not proof that every likelihood in the repository is free of such a mixture.

The numerator and selection normalization must use the same population density under the declared generative model. This consistency is more fundamental than a universal instruction to set a parameter named f to zero.

## 4. Bootstrap and uncertainty scope

The ordinary validation bootstrap in the inspected code resamples paired (redshift, retention) validation records with replacement, refits the inclusion model, and rebuilds the host density while holding the catalogue and GW distances fixed. It is an approximate IID validation-sample sensitivity analysis. It does not apply a finite-population correction and is not a full repeated-sampling uncertainty interval for H0.

The repository also has an explicitly exploratory pseudo-population SRSWOR bootstrap. Its documentation correctly says that finite-population mean guarantees do not automatically transfer to the nonlinear logistic/KDE estimator. It should not be described as design-valid until that specific estimator's properties are established.

"Method C likelihood marginalization" is not an actionable statistical specification. Before implementing any integration over inclusion-model uncertainty, the target joint model, conditioning structure, nuisance prior/design likelihood, and whether feedback from GW data to the inclusion stage is allowed must be written down. A cut posterior, plug-in bootstrap, likelihood averaging, and full joint Bayesian marginalization are different procedures.

## 5. Downstream H0 likelihood and design limitations

The supplied dossier correctly emphasizes that the selection normalization must match the population density used in the event integral. It should also preserve the current design's exact limitations:
- fixed detected-event count; no event-rate likelihood;
- simplified one-dimensional redshift/distance model;
- the current toy generator/likelihood retain events when observed distance satisfies `d_obs >= tau`, i.e. they retain the far-distance side of the threshold. This is the reverse of the common simplified intuition that detects events below a maximum distance, and it does not match the Consensus dossier's displayed `I(d <= 250 Mpc)` equation. The generator and likelihood are internally aligned to their toy rule, but that rule is not detector-realistic;
- synthetic-only validation, not LVK data or a measured H0 result.

A fixed truth H0=70 does not create a generating-truth stratum of 55<=H0<=85. If that interval is used, call it a posterior-grid interior region or boundary-sensitivity region, not a truth stratum. A separate truth-stratified design would require varying H0_truth across replicates and preregistering the strata.

The proposed R=2,000 and Ndet=300 are candidate design values, not justified sample sizes. Estimate paired outcome variance, coverage MCSE, and failure-rate precision from pilots, then freeze replication counts and stopping rules.

## 6. Novelty and manuscript language

The defensible position is an underexplored intersection of established methods. Two-phase sampling already supplies estimated inclusion probabilities, weighting, and replication machinery; dark-siren studies already study catalogue incompleteness, galaxy weighting, and H0 calibration. The claimed contribution is the controlled synthetic bridge between these domains, including finite validation, nonlinear misspecification, weight-instability diagnostics, downstream calibration, and failure-aware software auditing.

Use a qualified gap statement such as:

> In the literature reviewed for this study, we did not identify a single dark-siren study that jointly evaluates finite validation samples, estimated catalogue-inclusion probabilities, nonlinear inclusion-model misspecification, inverse-probability host reconstruction, weight-instability diagnostics, downstream H0 posterior calibration, and explicit optimizer-failure accounting. The contribution is a reproducible synthetic benchmark connecting established two-phase-sampling and dark-siren methodologies, not a new weighting estimator or a new Bayesian method.

Do not use an absolute "first" claim unless a documented systematic review, inclusion criteria, search dates, and primary-source verification support it.

## 7. Implementation-grounded status and next work

Source inspection of the current branch shows:
- the production wrapper `estimate_host_density_from_validation` still dispatches only to V1 linear/quadratic fitters;
- the pilot calls that wrapper without a solver-version option;
- the ordinary validation bootstrap refits through that same wrapper;
- the preposterior diagnostics currently accept V1 model identifiers only;
- V2 remains standalone, despite successful focused adversarial CI.

Therefore, V2 integration is a real outstanding engineering task. It should be added as an explicit opt-in solver path, preserving V1 as the default reference. The solver identity, optimizer status, boundary flag, projected-gradient/KKT diagnostics, and failure reason must survive into pilot rows and bootstrap artifacts. Tests must prove that the V2 path was actually invoked end-to-end; a unit test of the standalone fitter is not enough.

### Priority order
1. Complete exact-head Python 3.10/3.11/3.12 CI for the current V2 branch.
2. Integrate opt-in V2 through the validation-to-density wrapper and IID bootstrap, with diagnostics support.
3. Extend pilot configuration and artifact schema to retain solver identity and V2-specific convergence fields.
4. Add end-to-end artifact auditor checks that reject solver metadata mismatches.
5. Calibrate any proposed positivity/weight warning thresholds using simulations; do not hard-code universal thresholds.
6. Freeze estimands, prior/grid policy, truth strata, operational failure denominators, and MC precision before confirmatory simulations.
7. Recover historical D-01 inputs if possible; otherwise retain "historical root cause unavailable" and distinguish prospective stress tests from a historical replay.

**Release status remains NO-GO for confirmatory claims.** No main-branch edits, merges, or release declaration are made by this audit.
