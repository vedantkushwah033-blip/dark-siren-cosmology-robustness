# Five-arm matched comparison: preregistration and implementation contract

**Date:** 2026-11-13. **Status:** preregistered design contract; one engineering smoke replicate has been executed and independently audited, but no confirmatory five-arm experiment has been executed. **Scope:** synthetic low-redshift Euclidean dark-siren methodology, not LVK/GWTC analysis or observational H0 measurement.

> **Implementation status (2026-11-13):** A reflected weighted-KDE likelihood path is connected to the smoke runner. The runner now records each arm's realized redshift-space kernel width, and the independent raw-artifact auditor reconstructs those widths from the raw weights. The smoke weighted-Scott multiplier (0.7) is an engineering fixture only; it has not been selected by the preregistered parent-density training/holdout procedure. No arm-ranking or scientific result may be inferred until exact-head CI and the bandwidth-selection gate pass.

## 1. Research question

Under what conditions does estimating galaxy-catalogue inclusion probabilities improve dark-siren inference, and under what conditions do model misspecification, finite validation samples, and inverse-probability-weight instability eliminate or reverse that benefit?

## 2. Five distinct comparison arms

The informal labels “oracle inclusion weighting” and “known-probability IPW” can refer to the same calculation. They must not be counted as two methods. This protocol makes the five arms mathematically distinct:

| Arm | Host-population input | Interpretation |
|---|---|---|
| A. Catalogue-only | Observed thinned catalogue, weighted by the declared base intensity/merger-rate function b(z) | Uncorrected practical baseline |
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

**A — Catalogue-only.** Estimate host-redshift density from the observed analysis catalogue using the declared base-population/merger-rate weights (equal weights only when the declared base measure is uniform).

**B — Parent-population oracle.** Estimate host-redshift density from all galaxies in the simulated parent population before thinning, using the same base-population/merger-rate weighting rule as every other arm. This is a finite-parent reference, not the exact analytic population density unless the latter is explicitly used as a separate benchmark.

**C — Known-probability IPW.** For each observed catalogue galaxy use raw weight w_i=b(z_i)/pi_true(z_i), where b(z) is the same declared base-population/merger-rate weight used in A and B. Do not clip weights in the primary analysis. Report effective sample size and concentration diagnostics.

**D — Plug-in estimated IPW.** Fit the same degree-one logistic model to the validation data and use the posterior mode theta_MAP: w_i=b(z_i)/pi(z_i;theta_MAP). This defines the point-estimate arm. Do not silently replace it with posterior-mean probabilities; that would be a distinct sensitivity estimator.

**E — Uncertainty-integrated estimated IPW.** Use the existing proper validation posterior quadrature. At each coefficient node theta_b, calculate node-specific weights w_i(theta_b)=b(z_i)/pi(z_i;theta_b), evaluate every event numerator and the GW detection-efficiency normalization using that same theta_b, and marginalize the complete joint likelihood once:
L_E(H0) = sum_b omega_b exp{ell(H0,theta_b)}.
Normalize the H0 posterior once. Do not average separately normalized posterior curves or average event likelihood and selection efficiency independently. This is the existing modular/cut target: validation data determine p(theta|V), and GW data do not feed back into theta.

## 5. Common numerical choices and fairness

- Same H0 prior support/density, H0 grid, redshift grid, merger-rate/base-population function, event likelihood, and GW selection convention in all arms.
- Same preselected weighted-Scott multiplier in A-E, chosen in a separate parent-density training/holdout study without looking at H0 outcomes. The corrected reflected-KDE implementation computes each arm's Scott factor from that arm's original weighted sample, so the realized physical kernel widths can differ across arms as their weighted variance and effective sample size differ. This is a common bandwidth-selection rule, not an identical absolute redshift bandwidth; report each realized width and add an absolute-width sensitivity if conclusions depend on this distinction. Bandwidth sensitivity is a separate factorial analysis.
- **Current implementation status:** PR #44 now contains a separate common reflected-KDE likelihood path. The smoke runner records realized physical bandwidths for A-D and every E posterior node; the independent auditor reconstructs these widths. The smoke weighted-Scott multiplier is not a scientifically selected bandwidth; the training/holdout selection gate remains open.
- Apply the same declared base-population/merger-rate weighting rule b(z) in every arm: A uses b(z) on the observed catalogue; B uses b(z) on the full parent sample; C, D and E use b(z)/pi(z) with true, MAP and node-specific probabilities respectively. The parent and catalogue base-weight arrays must be evaluated from the same function and parameterization.
- Use the same reflected weighted-KDE estimator and preselected weighted-Scott multiplier on all five populations; record realized physical bandwidths for each arm and each E node. No arm may silently omit the merger-rate/base-population factor.
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

**The single engineering smoke is a pipeline check only. No estimator ranking, calibration result, or scientific conclusion is claimed by this document.**


## 11. Numerical-gate implementation update (2026-11-13)

The branch now includes a deterministic downstream H0-grid refinement regression test for the reflected-KDE five-arm likelihood. It evaluates the same fixed synthetic replicate at 81, 161, and 321 H0 points, normalizes each arm under the same uniform [50, 90] prior, and checks 161-to-321 stability for all five arms using maximum posterior-density difference <= 0.003 and selected posterior-quantile difference <= 0.20 km/s. This is a numerical regression gate, not a substitute for a representative multi-replicate grid-sensitivity study. The new test is wired into the Python-version CI matrix; its result remains pending until the exact-head workflow completes.


## 12. Multi-replicate orchestration pilot (2026-11-13)

A reusable engineering-pilot orchestrator is now implemented in `scripts/run_five_arm_replicates.py`. The single-replicate runner persists the shared observable and truth inputs, plus every nested validation sample, before fitting inclusion probabilities or evaluating the H0 likelihood; a downstream numerical failure therefore does not erase the raw generated replicate. It derives separate analysis-population and validation-frame seeds from a master `SeedSequence`, runs each replicate through the existing single-replicate runner, and invokes the independent raw-artifact auditor in a separate process. Each replicate uses one shared parent population, thinned catalogue, detected-event set, and nested validation frame across all five arms and all validation sizes.

The orchestrator writes:
- `replicate_ledger.csv`: one explicit SUCCESS/FAILED/NON_EVALUABLE lifecycle row for every replicate × validation-size × arm cell, including both seeds and audit status;
- `replicate_posterior_summaries.csv`: successful replicate-level posterior summaries only;
- `replicate_audits.csv`: independent raw-artifact audit outcomes and reconstruction discrepancies;
- `paired_contrasts.csv`: within-replicate differences in median error, absolute error, squared error, and 68% interval coverage against A and C;
- `aggregate_metrics.csv`: success-only metrics plus failure-inclusive bounds for bias, absolute error, RMSE, coverage and interval width;
- `config.json`: master seed, deterministic child-seed schedule, experiment settings, and explicit failure-bound assumptions;
- `manifest.json`: SHA-256 and byte-size records for every pilot artifact, including per-replicate raw inputs and audit reports.

The failure-inclusive bounds use the declared H0 prior support [50, 90] km/s/Mpc. For a generating value H0_true inside that support, the maximum possible absolute posterior-median error is M = max(|50 - H0_true|, |90 - H0_true|), and the squared-error bound is M²; at the default H0_true=70 these reduce to 20 and 400. The 68% coverage indicator is bounded by [0,1], and equal-tail interval width by [0,40]. These are conservative worst-case bounds for missing/non-evaluable replicate outcomes, not imputed point estimates. The paired contrast file is descriptive at pilot scale; no normal-theory p-values or superiority claims are produced.

A separate `scripts/audit_five_arm_replicates.py` verifies the top-level artifact manifest, reconstructs the deterministic seed schedule, cross-checks each replicate configuration against the frozen top-level design, verifies lifecycle completeness and per-replicate audit references, and independently reconstructs paired contrasts and aggregate metrics from the saved replicate summaries and ledger without importing the orchestration runner. The CI workflow now includes a two-replicate orchestration smoke on Python 3.12, plus CLI checks on all matrix versions and a regression test that deliberately triggers an invalid quadrature order to verify raw-input preservation. The pilot smoke is intended to exercise the orchestration and audit plumbing; it is not the planned N=4000, 30-event confirmatory study and does not close the corrected bandwidth-selection gate.


## 13. Failure-inclusive error-bound correction (2026-11-13)

The multi-replicate aggregation now derives the worst-case posterior-median error bound from the actual generating value and the declared prior support, rather than assuming the generating value is the midpoint. For H0_true in [50, 90], M = max(|50-H0_true|, |90-H0_true|); failure-inclusive absolute-error and RMSE upper bounds use M and M² respectively, and bias bounds use ±M. The independent auditor independently reconstructs this quantity from config.json and verifies the declared bound contract. A regression test covers an off-centre generating value H0_true=60, for which M=30. This is a correctness fix to conservative missing-outcome accounting, not a scientific result.


## 14. Training/holdout bandwidth-selection implementation (2026-11-13)

A separate selector and independent auditor are now included:
- `scripts/select_bandwidth_train_holdout.py` selects the weighted-Scott multiplier using only parent-sample KDE integrated absolute error against the known analytic merger-weighted parent density.
- Replicate IDs are split deterministically into disjoint training and holdout sets. The selected multiplier is minimized on training IAE only; all candidate holdout scores are reported without re-selection.
- `scripts/audit_bandwidth_train_holdout.py` independently reconstructs the split, training scores, selected multiplier, holdout scores, and source hashes.
- The selector does not read H0 posterior or coverage outcomes. This avoids choosing smoothing to make the cosmological inference look favorable.
- The density-recovery runner now accepts an explicit comma-separated bandwidth grid. The default grid is 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30 and 0.40 weighted-Scott multipliers. It remains necessary to generate a fresh source artifact with this candidate grid and adequate independent replicates before selecting a project bandwidth.
- A unit test verifies that a deliberately worse holdout score does not cause the training-selected multiplier to be changed.

The selector implementation is an analysis workflow, not a completed bandwidth study. No multiplier has been selected from actual training/holdout simulation artifacts yet.


### Reproducible execution commands for the bandwidth gate

The lightweight parent-only source experiment can be run without generating the much larger finite-validation factorial artifact:

```bash
python scripts/run_parent_density_bandwidth_pilot.py \
  --replicates 200 --n-parent 4000 --seed 20261114 \
  --grid-points 801 \
  --bandwidth-multipliers 0.05,0.08,0.10,0.12,0.15,0.20,0.25,0.30,0.40 \
  --output-dir artifacts/parent_density_bandwidth_source

python scripts/audit_parent_density_bandwidth_pilot.py \
  artifacts/parent_density_bandwidth_source \
  --output artifacts/parent_density_bandwidth_source_audit.json

python scripts/select_bandwidth_train_holdout.py \
  artifacts/parent_density_bandwidth_source \
  --output artifacts/parent_density_bandwidth_selection.json \
  --split-seed 20261113 --training-fraction 0.70

python scripts/audit_bandwidth_train_holdout.py \
  artifacts/parent_density_bandwidth_source \
  artifacts/parent_density_bandwidth_selection.json \
  --output artifacts/parent_density_bandwidth_selection_audit.json
```

The 200-replicate run is a training/holdout bandwidth-selection pilot, not the final confirmatory five-arm campaign. Freeze the selected multiplier and source hashes before applying it to the five-arm H0 study. Do not reselect after viewing H0 posterior results.


## 15. Dedicated bandwidth-selection execution workflow (2026-11-14)

The exact-head Python 3.10/3.11/3.12 numerical-diagnostics matrix passed at commit `b0f3a94e9e992e3c29d75cadd210bfa5a01d042c`. A dedicated GitHub Actions workflow, `.github/workflows/bandwidth-selection.yml`, now runs the preregistered 200-replicate parent-density bandwidth source experiment on this research branch, independently reconstructs every candidate IAE, performs the deterministic 70/30 training/holdout selection, and independently audits the selected result. The raw source data, selection report, and both audit reports are uploaded as a versioned workflow artifact. This execution is deliberately separate from the five-arm H0 outcomes: no cosmological posterior or coverage metric enters the bandwidth selection. The selected multiplier is eligible for use only if the complete source and selection audits pass; holdout scores remain a report-only assessment and cannot trigger re-selection.
