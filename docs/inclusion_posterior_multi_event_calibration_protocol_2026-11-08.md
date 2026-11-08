# Multi-event inclusion-posterior calibration protocol

**Status:** pre-execution protocol; no calibration result is claimed.  
**Scope:** synthetic low-redshift Euclidean distance model only.  
**Inference target:** modular/cut posterior in which the independent validation sample determines the inclusion-coefficient posterior and GW events do not update those coefficients.

## 1. Purpose and separation of claims

This experiment asks whether finite validation information and inclusion-model misspecification propagate into calibration and repeated-sampling coverage of a shared-coefficient, selection-conditioned multi-event H0 posterior.

Keep the following outcomes separate:

1. **Inclusion-model recovery:** held-out probability calibration and prediction error.
2. **IPW stability:** Kish/entropy ESS, maximum normalized host weight, and top-tail weight concentration.
3. **Quadrature convergence:** stability of coefficient moments and validation posterior predictions.
4. **Downstream numerical convergence:** stability of the shared-node selected H0 log-likelihood shape and normalized posterior across quadrature orders and H0 grids.
5. **Conditional H0 PIT:** prior-predictive PIT over H0 conditional on fixed inclusion-truth parameters. This is not joint SBC.
6. **Fixed-truth coverage:** repeated-sampling coverage at predeclared H0 values.
7. **Failure/evaluability:** proportion of attempted replicates that produce a valid, numerically converged posterior.

A pass in one category does not establish any of the others.

## 2. Data-generating design

Use the same low-redshift Euclidean model already implemented:

- Parent redshift support: [0.01, 0.40].
- Parent redshifts sampled from the declared smooth latent intensity.
- Merger/host base weight: (q(z)=1+1.5(z-z_{min})/(z_{max}-z_{min})).
- GW distance model: (d_{obs}sim N(cz/H_0,sigma^2)).
- Common distance noise: (sigma=140) Mpc.
- GW detection: (d_{obs}geq250) Mpc.
- Fixed detected-event count per replicate; the generator records every candidate, rejection, and selected event.
- Analysis catalogue inclusion and validation retention are separate Bernoulli draws from the same declared inclusion-probability function.

### Inclusion truth regimes

Use the implemented floor-plus-sigmoid truth:

[
pi(z)=f+rac{1-f}{1+exp((z-z_m)/w)},
quad z_m=0.22,quad w=0.05.
]

Predeclare at least these floor regimes:

| Floor (f) | Interpretation |
|---:|---|
| 0.00 | Exactly logistic in redshift; degree-one logistic fit is correctly specified. |
| 0.10 | Moderate floor-induced nonlinear-logit misspecification. |
| 0.25 | Stronger floor-induced nonlinear-logit misspecification. |

The floor-zero case is a necessary reference: it distinguishes finite-validation uncertainty under a correctly specified logistic model from error caused by the positive-floor misspecification. Do not describe the positive-floor cases as generic or universal misspecification.

### Validation design

For each replicate:

1. Generate one analysis parent population, its thinned observed catalogue, and the multi-event GW sample.
2. Generate a separate validation parent frame with a distinct seed.
3. Draw validation observations by SRS without replacement.
4. Reuse the same frame permutation to create nested validation sizes (n_Vin{100,300,1000}).
5. Hold the analysis catalogue and GW events fixed across those three validation sizes.

This is the fixed-analysis-catalogue design. Do not remove validation units from the analysis catalogue: the validation frame is independently generated. Keep any disjoint finite-parent design as a separate experiment and do not pool its results with this one.

## 3. Inference arms

### Primary inclusion model

- Logistic polynomial degree 1.
- Proper independent Gaussian coefficient prior with intercept SD 5 and slope SD 2.5 on a fixed domain-standardized covariate, (x=(z-0.205)/0.195) for the declared redshift domain [0.01, 0.40]. The coordinate origin and scale are fixed before validation sampling and are identical across replicates and validation sizes.
- Slope-prior sensitivity multipliers: 0.5, 1.0, 2.0; hold the intercept prior fixed.
- No weight clipping in the primary analysis.
- Same node-specific normalized IPW host weights in every event mixture and the one-sided beta selection factor.
- One shared inclusion-coefficient node across all events, followed by one joint log-sum-exp marginalization.

Degree 2 may be included as a predeclared secondary model, not substituted after inspecting the degree-one results.

### Numerical settings

Pilot:
- 20 replicates per design cell.
- H0 grid: 41 or 81 points.
- Quadrature order ladder: 5, 7, 9, 11.

Confirmatory:
- Target 1,000 replicates per design cell where computationally feasible.
- H0 grid refinement: at least 81 versus 161 points.
- Quadrature order ladder: 7, 11, 15, 19, with additional orders only under a predeclared convergence rule.
- A replicate is not numerically evaluable for the primary result if downstream H0 log-likelihood shape, posterior density, or quantiles fail the declared convergence criteria.

These are proposed starting settings. Freeze the final configuration and tolerances before launching confirmatory runs.

## 4. Two distinct calibration experiments

### A. Conditional H0 prior-predictive PIT (fixed inclusion truth)

For each replicate, draw (H_0^{true}sim Uniform(50,90)), then generate the analysis population, validation sample, and detected events conditional on that same (H_0^{true}). Hold the inclusion-truth parameters fixed within each catalogue-floor cell, as in the current runner.

Compute the posterior CDF at the generating value:

[
u_r=F_r(H_0^{true}).
]

Report:
- PIT histogram with fixed equal-width bins.
- Empirical PIT mean and variance.
- Uniformity diagnostic with effect size and uncertainty, not only a p-value.
- Quantile-based calibration summaries.
- Failure/non-evaluable fraction.

Because positive-floor inclusion truth is intentionally misspecified relative to the fitted logistic model, uniform PITs are not guaranteed. The floor-zero scenario is the correctly specified inclusion reference, subject to finite catalogue approximation and all other model assumptions.

### B. Fixed-truth repeated-sampling coverage

Use (H_0^{true}=70) Mpc-based km/s/Mpc units (70 km s(^{-1}) Mpc(^{-1})) for every replicate. Independently regenerate the parent population, catalogue, validation frame, and detected events.

For 68%, 90%, and 95% equal-tail credible intervals, report:
- Empirical coverage.
- Binomial Monte Carlo standard error and confidence interval.
- Mean and median interval width.
- Posterior median bias.
- Evaluability and failure rates.

Do not call conditional H0 PIT a joint SBC experiment, do not call PIT a coverage experiment, and do not call coverage a calibration proof.

## 5. Replication and uncertainty

Use a pilot of 20 replicates only to test orchestration, failure handling, runtime and auditability. Do not interpret its coverage/PIT as confirmatory.

For a coverage probability near 0.68, 1,000 evaluable replicates give approximate Monte Carlo standard error

[
sqrt{0.68(1-0.68)/1000}approx0.0147,
]

or roughly ±2.9 percentage points for a normal-approximation 95% interval. If the final number of evaluable replicates is lower, report the wider uncertainty and do not imply the original target was achieved.

For paired validation-size contrasts, retain replicate identity and use the same analysis catalogue/events and nested validation frame across (n_V). Bootstrap at the replicate level, not at the event or host level.

## 6. Failure ledger

Every attempted replicate and design cell must produce one explicit status:

- `SUCCESS`: fit, likelihood, posterior and required numerical checks completed.
- `NON_EVALUABLE`: data generated, but predeclared numerical convergence/evaluability criteria failed.
- `FAILED`: execution or validation error prevented a valid result.
- `NOT_RUN`: planned but not attempted; retain the planned row and reason.

Record stage, exception class/message, seeds, configuration hash, event-generation attempts, quadrature order, posterior-grid size, and audit status. Never drop failed replicates from the ledger. Report coverage among evaluable replicates together with the fraction evaluable; optionally provide conservative all-attempt bounds.

## 7. Required raw artifacts and independent audit

For every replicate retain:

- Observable analysis catalogue: object IDs and redshifts only.
- Observable GW events: event IDs and measured distances only.
- Observable validation sample: validation IDs, redshifts, retention labels only.
- Truth-only parent inclusion probabilities, host IDs/redshifts, H0 truth, and complete GW candidate-generation trace in separate files.
- Coefficient quadrature nodes/weights, node-specific host weights or enough raw inputs to reconstruct them, nodewise event likelihoods, nodewise beta, marginalized log likelihood, H0 posterior curve and summaries.
- Configuration, RNG stream identifiers, file hashes, and failure ledger.

The independent auditor must:
1. Verify the inclusion DGP parameters match between analysis and validation frames.
2. Reconstruct the validation sample from the validation truth frame and SRS IDs.
3. Reconstruct every GW selection decision and selected event from the generation trace.
4. Reconstruct coefficient quadrature from validation-only observations and declared prior.
5. Reconstruct node-specific IPW weights, event-mixture likelihoods, beta and shared-node marginalization.
6. Reconstruct H0 posterior curves and credible intervals.
7. Verify checksums and preserve any discrepancy as an audit failure.

Truth-only columns must never enter non-oracle fitting or likelihood calculation.

## 8. Predeclared interpretation boundaries

- A higher Kish ESS does not guarantee low host-density error or calibrated H0 inference.
- A converged quadrature ladder does not imply posterior calibration.
- A numerically normalized H0 posterior does not imply correct coverage.
- A floor-zero success does not validate positive-floor misspecification cases.
- A positive-floor failure does not imply all logistic inclusion models fail.
- No universal IPW superiority claim follows.
- No real galaxy survey, detector, LVK/GWTC event, or observational H0 claim follows from this synthetic experiment.

**Go/no-go:** run a small smoke and independent audit first. Start any confirmatory conditional-PIT, joint-SBC, or coverage experiment only after the generator/inference selection match, raw-artifact reconstruction, prior sensitivity, downstream quadrature convergence, and failure ledger are all verified.


## 9. Prior-coordinate correction (2026-11-08)

The earlier sample-standardized coefficient prior made its implied prior over
physical-redshift inclusion curves depend on the realized validation sample
mean and standard deviation. That is unsuitable for clean comparisons across
replicates and validation sizes.

The H0 propagation runner now uses fixed domain-based coordinates:
`x=(z-z_mid)/z_half_range`, with center 0.205 and scale 0.195 for the
declared [0.01, 0.40] redshift domain. The coordinate specification is
recorded in the artifact configuration and checked by the independent auditor.
The posterior fitter and quadrature ladder retain sample-standardization as
their backward-compatible default, but accept an explicit fixed coordinate
center and scale; the H0 propagation and prior-sensitivity path passes these
fixed values explicitly.

This correction is a prerequisite for interpretable repeated-sampling
Bayesian comparisons. It does not itself establish calibration or coverage.


## 10. Distinguish the implemented conditional PIT from future joint SBC

The current orchestration runner now labels its randomized-H0 arm
`H0_PIT_FIXED_INCLUSION_TRUTH`. It draws H0 from the inference H0 prior but
holds the inclusion DGP parameters fixed. This is a conditional H0
prior-predictive PIT, not joint SBC for (H0, inclusion coefficients).

In the floor-zero reference, the fixed inclusion truth in the declared
fixed-coordinate basis is (theta_0, theta_1)=(0.30, -3.90), while inference
uses independent Gaussian coefficient priors with standard deviations
(5.0, 2.5). The coefficient truth is therefore not drawn from the inference
prior. Positive catalogue floors additionally introduce functional-form
misspecification.

A genuine joint or modular-SBC design requires an explicit prior draw for
the generating inclusion coefficients, a matched logistic DGP, and a precise
definition of the target posterior. Because the current method is modular/cut
and deliberately prevents GW data from updating inclusion coefficients,
ordinary full-joint SBC uniformity must not be presumed for that target.
See `conditional_h0_pit_vs_joint_sbc_audit_2026-11-08.md` for the detailed
derivation and proposed requirements.
