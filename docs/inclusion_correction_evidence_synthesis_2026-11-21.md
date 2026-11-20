# Evidence synthesis and decision gate: estimated catalogue inclusion correction

**Date:** 2026-11-21  
**Scope:** synthetic low-redshift Euclidean dark-siren benchmark only. No LVK/GWTC/public event data, real galaxy survey, or observational H0 measurement is used.

## Executive conclusion

The current evidence supports treating estimated inverse-probability weighting (IPW) as a sensitivity method whose performance depends on inclusion-model specification and weight stability. It does **not** support promoting a higher-degree inclusion model, coefficient integration, or any single pre-posterior weight-tail diagnostic as a generally reliable solution.

Three evidence streams now constrain the next steps:

1. A 50-truth prior-predictive five-arm pilot shows substantial improvement over catalogue-only inference, but estimated-IPW arms do not consistently match known-probability IPW; the full-parent arm is an oracle reference.
2. A prospective held-out test of two predeclared top-5%-mass predictors found zero of 18 strata meeting the Holm-adjusted predictive-value criterion.
3. A 30-event degree-one/degree-two quadrature pilot shows good adjacent order-7/order-9 numerical agreement in one realization, but severe degree-two host-weight concentration and no consistent posterior shift in favor of degree two.

These are complementary pilots, not a single pooled calibration study. Do not combine their replicate counts as if they came from one exchangeable experiment.

## Evidence ledger

### A. Five-arm prior-predictive pilot (50 independent truths)

- Design: 50 independently drawn H0 truths; 1,200 parent galaxies; 30 detected events; curvature {-1,0,+1}; validation sizes {100,300,500}; five paired arms.
- Execution: 2,250/2,250 planned arm-cells SUCCESS; 2,250/2,250 audit-PASS; 150/150 child audits; aggregate audit PASS.
- Pooled descriptive RMSE (km/s/Mpc): A catalogue-only 17.38; B full-parent oracle 3.60; C known-π IPW 4.70; D degree-one MAP IPW 5.80; E degree-one integrated IPW 5.46.
- Pooled 68% coverage: A 0.7%; B 70.0%; C 52.7%; D 45.8%; E 54.4%.
- E improves on D in pooled RMSE but does not consistently beat C. At validation n=500, gamma=+1, E's paired squared-error disadvantage versus C was +10.37 [95% paired t interval +3.43,+17.32].
- Interpretation: strong evidence of catalogue-only failure in this particular synthetic DGP; weaker and curvature-dependent evidence for estimated-IPW benefit. B is a matched synthetic oracle, not an implementable estimator.

### B. H0-grid refinement for the 50-truth archive

- Exact selected-likelihood recomputation conditional on archived KDE densities compared 561 uniform, 1,121 geometric, and 2,241 geometric H0 grids.
- 561→1,121 maximum posterior-density L1 difference was 0.00191816, exceeding the strict 0.001 shape criterion.
- 1,121→2,241 maximum density L1 difference was 0.00012453; maximum median shift 0.000542 km/s/Mpc.
- Thus 1,121 is stable against 2,241 for this archived set of curves, while 561 is not shape-converged under the chosen criterion.
- Limitation: this is exact H0 integration conditional on frozen KDEs, not raw-galaxy re-generation or a fully grid-converged end-to-end rerun.

### C. Prospective held-out diagnostic prediction

- Frozen design: two independent 100-replicate training streams and one untouched 100-replicate holdout stream; fixed H0=70; gamma {-1,0,+1}; validation sizes {100,300,500}; arms D and E.
- Predictors were prespecified top-5%-weight mass (D) and maximum node top-5%-mass (E); each was evaluated against a training-mean absolute-error baseline.
- Current-head workflow completed all 4,500/4,500 holdout lifecycle cells and 18 predeclared tests with 3,000 nested bootstrap resamples per cell.
- 12/18 point estimates improved on the baseline and 5/18 unadjusted intervals excluded zero, but **0/18 passed the predeclared Holm-adjusted criterion**.
- Interpretation: these two scalar diagnostics did not demonstrate prospective predictive value in this one fixed-truth synthetic DGP. This does not invalidate reporting ESS, maximum weight, tail mass, support gaps, or fitted curvature as descriptive warnings; it rejects the stronger claim that these two predeclared predictors provide reliable out-of-sample error prediction under this design.

### D. 30-event degree-one/degree-two quadrature sensitivity

- Design: one fixed H0=70 realization; 1,200 parent galaxies; 30 events; gamma {-1,0,+1}; validation sizes 100/500; inclusion degrees 1/2; Gauss-Hermite orders 5/7/9.
- All 18 child audit reports passed; 342 manifest payload files matched their declared byte counts and SHA-256 hashes.
- For integrated-IPW E, maximum order-7→9 density L1 difference was 0.0006544 (degree 1) and 0.0003966 (degree 2). Maximum median shifts were 0.002071 and 0.001670 km/s/Mpc.
- Degree-two minus degree-one E posterior-median shifts ranged from -1.18 to +4.11 km/s/Mpc across the six curvature × validation-size cells; the direction was not consistent.
- Degree-two E had minimum node-level host-weight Kish ESS as low as 1.24 and maximum node-level top-5%-weight mass as high as 0.9999. Quadrature convergence therefore did not imply stable host weighting.
- Interpretation: this is a numerical convergence and single-realization sensitivity pilot, not repeated-sampling calibration or degree-selection evidence.

## Decisions for the next phase

1. **Keep all current research PRs draft/unmerged.** Do not promote degree two or change production defaults on the basis of the one-realization quadrature pilot.
2. **Retain the five-arm comparison.** A catalogue-only baseline, full-parent oracle, known-π IPW, estimated MAP IPW, and integrated estimated IPW answer different diagnostic questions and should remain paired.
3. **Do not invent a universal ESS or top-tail cutoff.** The prospective test did not validate the two prespecified scalar predictors. Continue to publish the full diagnostic panel descriptively, without calling it a calibrated failure classifier.
4. **Make the next confirmatory run prior-predictive and independently seeded.** Target at least 100 independent H0 truths per curvature × validation-size cell; preserve common random numbers across arms within each truth; do not count repeated validation subsamples as independent truths.
5. **Use a grid that has passed the archived-curve refinement gate.** Start at 1,121 geometric H0 points on [30,110] (or demonstrate a stricter exact-grid criterion for any alternative). Keep KDE-generation and H0-grid integration convergence as separate gates.
6. **Resolve quadrature separately.** For coefficient-integrated arms, predeclare an adjacent-order criterion (at least 7→9 in the present design) and retain node-level host-weight diagnostics. Numerical quadrature ESS and host-weight Kish ESS are different quantities.
7. **Fail closed.** Every planned replicate/arm cell must remain in the ledger as SUCCESS, NON_EVALUABLE, FAILED, or NOT_RUN. Report the planned denominator and all failures; do not silently drop failed cells.
8. **Report calibration and accuracy separately.** Include rank/PIT diagnostics, 68/90/95% coverage, bias, RMSE, paired contrasts with replicate-level uncertainty, and numerical convergence. Do not treat density IAE or ESS as substitutes for H0 calibration.

## What this evidence does and does not establish

It establishes that the current synthetic pipeline can execute several large audited pilots, that catalogue-only inference is badly mismatched in the chosen synthetic setting, that estimated IPW has material curvature-dependent limitations, and that the currently tested tail-mass predictors did not pass a prospective familywise predictive criterion.

It does **not** establish estimator superiority, universal diagnostic thresholds, robust performance under arbitrary selection mechanisms, observational H0 accuracy, or validity for real galaxy catalogues and gravitational-wave events.
