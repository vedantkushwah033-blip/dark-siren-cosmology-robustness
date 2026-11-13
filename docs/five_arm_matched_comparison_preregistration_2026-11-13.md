# Five-arm matched comparison: preregistration and implementation contract

**Date:** 2026-11-13. **Status:** design contract only; no five-arm experiment has been executed. **Scope:** synthetic low-redshift Euclidean dark-siren methodology, not LVK/GWTC analysis or observational H0 measurement.

## 1. Research question

Under what conditions does estimating galaxy-catalogue inclusion probabilities improve dark-siren inference, and under what conditions do model misspecification, finite validation samples, and inverse-probability-weight instability eliminate or reverse that benefit?

## 2. Five distinct comparison arms

The informal labels “oracle inclusion weighting” and “known-probability IPW” can refer to the same calculation. They must not be counted as two methods. This protocol makes the five arms mathematically distinct:

| Arm | Host-population input | Interpretation |
|---|---|---|
| A. Catalogue-only | Observed thinned catalogue, equal galaxy weights | Uncorrected practical baseline |
| B. Parent-population oracle | Full simulated parent population before thinning | Infeasible finite-parent oracle benchmark |
| C. Known-probability IPW | Observed catalogue weighted by true 1/pi(z) | Oracle inclusion-probability benchmark |
| D. Plug-in estimated IPW | Observed catalogue weighted by 1/pi(z; theta_MAP) | Point-estimated inclusion correction |
| E. Uncertainty-integrated IPW | Shared-node selected likelihood integrated over p(theta|validation) | Bayesian modular uncertainty propagation |

B and C are not duplicates: B removes thinning by supplying the full parent population; C retains the thinned catalogue but supplies true inclusion probabilities. B is a simulation-only reference, not an implementable observational estimator. C is also an oracle diagnostic, not a deployable method.

## 3. Matched data-generating design

For each replicate, generate once and reuse across every arm and validation size:
- One parent galaxy population on z in [0.01, 0.40].
- One thinned analysis catalogue using the inclusion truth.
- One detected GW-event set from the existing Euclidean distance model and separate GW detection rule.
- One independent validation frame and a single random permutation defining nested samples n=100, 300, 1000.

The analysis catalogue, GW events, and parent population must be identical across all five arms. Do not regenerate events by arm. Catalogue inclusion and GW detection remain separate processes. Keep the independent validation frame separate from the analysis catalogue, consistent with the existing Design A implementation.

Use the existing floor-plus-sigmoid truth:
pi_true(z) = f + (1-f) / (1 + exp((z-z_m)/w)), z_m=0.22, w=0.05.
Primary floors: f=0 (degree-one logistic correctly specified), 0.10 (moderate floor-induced misspecification), and 0.25 (stronger floor-induced misspecification). Fit the declared degree-one logistic model using validation data only, with fixed domain-standardized coordinates and the existing proper Gaussian coefficient prior.

## 4. Exact arm definitions

**A — Catalogue-only.** Estimate host-redshift density from the observed analysis catalogue with equal weights.

**B — Parent-population oracle.** Estimate host-redshift density from all galaxies in the simulated parent population before thinning. This is a finite-parent reference, not the exact analytic population density unless the latter is explicitly used as a separate benchmark.

**C — Known-probability IPW.** For each observed catalogue galaxy use raw weight w_i=1/pi_true(z_i). Do not clip weights in the primary analysis. Report effective sample size and concentration diagnostics.

**D — Plug-in estimated IPW.** Fit the same degree-one logistic model to the validation data and use the posterior mode theta_MAP: w_i=1/pi(z_i;theta_MAP). This defines the point-estimate arm. Do not silently replace it with posterior-mean probabilities; that would be a distinct sensitivity estimator.

**E — Uncertainty-integrated estimated IPW.** Use the existing proper validation posterior quadrature. At each coefficient node theta_b, calculate node-specific IPW host weights, evaluate every event numerator and the GW detection-efficiency normalization using that same theta_b, and marginalize the complete joint likelihood once:
L_E(H0) = sum_b omega_b exp{ell(H0,theta_b)}.
Normalize the H0 posterior once. Do not average separately normalized posterior curves or average event likelihood and selection efficiency independently. This is the existing modular/cut target: validation data determine p(theta|V), and GW data do not feed back into theta.

## 5. Common numerical choices and fairness

- Same H0 prior support/density, H0 grid, redshift grid, merger-rate/base-population function, event likelihood, and GW selection convention in all arms.
- Same fixed KDE bandwidth in A-E, chosen in a separate parent-density training/holdout study without looking at H0 outcomes. Use the corrected reflected-KDE Scott implementation. Bandwidth sensitivity is a separate factorial analysis.
- B uses the same KDE and bandwidth on the full parent sample. C and D use the same weighted KDE with true and MAP inclusion probabilities. E uses that weighted KDE at each posterior node.
- No second hidden-host mixture correction on top of IPW.
- No clipping/trimming/stabilization in the primary arms. Such variants must be separately named and preregistered.
- Any comparison using different event sets, KDE bandwidths, H0 priors, or selection conventions across arms is invalid.

## 6. Outcomes and estimands

For fixed-truth repeated sampling, report by arm and design cell:
- Posterior-median bias and RMSE.
- 68%, 90%, 95% equal-tail interval coverage.
- Mean/median interval width.
- Posterior evaluability and numerical-failure rate.
- Paired within-replicate contrasts against A and C.

For conditional prior-predictive experiments, report PIT distributions and interval-inclusion rates separately. Do not call conditional PIT joint SBC.

For C, D and each E coefficient node, record Kish and entropy ESS, maximum normalized weight, top 1% and 5% mass, minimum/lower-quantile inclusion probability, support gaps, and redshift concentration. Separately assess fitted inclusion-model prediction and calibration on an independent holdout not used to fit the model.

The primary comparison is paired replicate-level differences, not differences of unpaired aggregate means. Report complete-pair counts and preserve FAILED/NON_EVALUABLE attempts in lifecycle tables and failure-inclusive bounds. Report Monte Carlo uncertainty for coverage and paired contrasts.

## 7. Replication and precision

The existing 20-replicate pilot is engineering-scale only. For nominal 68% coverage, approximate MCSE is sqrt(0.68*0.32/1000)=0.0148 at 1,000 replicates, versus about 0.104 at 20. Thus 20 replicates cannot support fine coverage comparisons. Before confirmatory execution, freeze replicate count using a target MCSE or confidence-interval half-width, including precision for paired contrasts. About 1,000 replicates per cell is a planning benchmark, not a substitute for precision calculation.

## 8. Required gates before confirmatory run

1. Execute and independently audit corrected-KDE training/holdout bandwidth selection.
2. Compare inclusion posterior quadrature with an independent reference, including weakly identified and near-separated validation samples.
3. Demonstrate H0-grid convergence at representative difficult cells, at least 81/161/321 points.
4. Demonstrate coefficient-quadrature convergence after propagation through the full H0 likelihood.
5. Independently verify identical parent/catalogue/validation/event inputs across arms.
6. Independently reconstruct arm-specific weights, posterior summaries, paired contrasts, and failure-inclusive bounds from raw artifacts.
7. Freeze configuration, seed derivation, software commit, and artifact hashes.

Software tests are not calibration results; smoke experiments are not confirmatory experiments; synthetic results are not observational validation.

## 9. Interpretation rules

- C better than A: known-probability IPW reduces catalogue-inclusion bias under this DGP; it does not show estimated probabilities work.
- D worse than C: quantifies the cost of estimating inclusion probabilities under the tested validation design/model.
- E better than D: uncertainty integration may mitigate finite-validation uncertainty; it does not automatically correct structural misspecification.
- D and E both poor under positive-floor misspecification: larger validation samples alone may not solve model error.
- C itself unstable: limitation may be intrinsic to IPW and catalogue support, not merely probability estimation.
- B versus C differences require inspection of finite-parent variation, KDE smoothing and support before attributing effects to inclusion correction.
- No universal method ranking is permitted; conclusions are conditional on this DGP, validation design, KDE and GW selection model.

## 10. Implementation sequence

1. Freeze this contract.
2. Add a deterministic fixture where all five arms can be independently reconstructed.
3. Implement one matched-replicate runner with common inputs and explicit lifecycle statuses.
4. Implement an independent raw-artifact auditor that does not import runner aggregation.
5. Run a tiny orchestration smoke test.
6. Complete bandwidth, coefficient-quadrature and H0-grid convergence gates.
7. Preregister confirmatory replication/precision.
8. Execute and independently audit the confirmatory synthetic experiment.

**No five-arm results are claimed by this document.**
