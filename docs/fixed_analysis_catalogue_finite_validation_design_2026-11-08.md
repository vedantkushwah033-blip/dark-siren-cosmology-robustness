# Fixed-analysis-catalogue design for finite-validation IPW

**Status:** pre-implementation protocol. No experiment, code path, or numerical result is claimed by this document.  
**Branch target:** `research/independent-posterior-audit-claim-fix`.  
**Scientific scope:** synthetic low-redshift Euclidean host-redshift density recovery. This is not an observational H0 analysis.

## 1. Motivation and estimand

The existing disjoint finite-parent design samples validation units from the same finite parent frame and excludes them from the analysis catalogue. It is a useful operational design, but changing validation size changes two things simultaneously:

1. the information available to estimate the inclusion function; and
2. the number and identities of retained galaxies available for density reconstruction.

Consequently, an observed trend across validation sizes cannot be attributed solely to improved inclusion-probability estimation.

This protocol defines a complementary **fixed-analysis-catalogue (FAC) design**. Its primary estimand is the change in density-recovery performance as the number of independently observed validation labels increases, conditional on exactly the same analysis parent frame and exactly the same retained analysis catalogue within each replicate and DGP scenario.

The FAC design is not a replacement for the existing disjoint finite-parent design. It isolates an information-size effect under independent validation sampling; the existing design measures a practical finite-frame validation/catalogue trade-off. Results from the two designs must be presented separately.

## 2. Sampling construction

For each replicate (r):

1. Draw an analysis parent frame (A_r={z_i^A}_{i=1}^{N_A}) iid from the declared parent-redshift law.
2. Independently draw one maximum-size validation frame (V_r^{max}={z_j^V}_{j=1}^{n_{max}}) from the same parent-redshift law. The validation frame is independent of the analysis frame.
3. For each DGP curvature/retention scenario, generate independent retention uniforms for the analysis frame and validation frame. Apply the same declared inclusion law (pi_gamma(z)) to both frames.
4. Define the analysis catalogue once:
   [
   C_r={z_i^A:R_i^A=1}.
   ]
   This catalogue is frozen across every validation-size condition in that replicate/scenario.
5. Define nested validation samples using deterministic prefixes of a pre-randomized validation-frame order:
   [
   V_{r,n}={(z_j^V,R_j^V):j=1,ldots,n},
   quad nin{100,300,1000}.
   ]
   The validation labels are observed for sampled validation units, irrespective of their retention outcome. The prefix order is randomized independently of redshift and retention.

The validation frame is an independent sample from the parent population, not a second slice removed from the analysis frame. No validation object is an analysis-catalogue object. For a given replicate and scenario, (A_r), its retention labels, and (C_r) are identical for n=100, 300 and 1,000. Only the amount of validation information changes.

The validation-frame sampling mechanism must be stated explicitly. This protocol assumes simple random sampling from the parent law with complete ascertainment of the binary retention label for every sampled unit. It does not represent a real survey's follow-up selection unless that selection is separately modelled.

## 3. Paired scenarios and random-number discipline

Use common random numbers across paired curvature scenarios where possible, but do not accidentally reuse one uniform array in a way that makes outcomes across distinct frames dependent in an undocumented fashion. Store separate arrays for:

- analysis-frame redshifts;
- validation-frame redshifts;
- analysis-frame retention uniforms;
- validation-frame retention uniforms;
- validation-frame randomized order.

The same arrays are reused across gamma/retention conditions to create paired comparisons, with the exact coupling declared in the manifest. The analysis catalogue is allowed to differ across DGP scenarios because the inclusion law differs; it must not differ across validation sizes within a scenario.

Store the generated frames and uniforms in raw records, or store sufficient deterministic seeds plus a fully specified generator/version and independently reconstructable ordering. Prefer raw records for the confirmatory artifact.

## 4. Comparison arms

For each frozen analysis catalogue and each validation-size condition, evaluate:

1. **Analytic target diagnostic:** declared normalized parent intensity multiplied by the declared merger-weight function. This is the target, not an estimator.
2. **Analysis-parent KDE diagnostic:** KDE using all (N_A) analysis-frame objects, before retention. This quantifies finite analysis-frame and smoothing error.
3. **Catalogue-only KDE:** KDE using the fixed retained analysis catalogue with unit weights.
4. **Oracle IPW KDE:** same catalogue weighted by the true (pi_gamma(z_i^A)^{-1}). Oracle probabilities are simulation diagnostics only and must never be passed to fitted estimators.
5. **Estimated degree-1 IPW KDE:** fit the linear logistic inclusion model to (V_{r,n}), then apply its predictions to the frozen analysis catalogue.
6. **Estimated degree-2 IPW KDE:** analogous quadratic fit.

Every arm uses the same redshift grid, merger-weight convention, boundary/reflection rule, and predeclared Scott-multiplier bandwidth. No H0 outcomes may be used to tune bandwidth. Bandwidth sensitivity is secondary and must not be collapsed into a single selected winner.

For the estimated arms, fit only from ((z_j^V,R_j^V)). The fitter must not receive analysis-frame truth probabilities, validation-frame truth probabilities, or any oracle quantities. The auditor may use truth probabilities to reconstruct the DGP and score the fitted predictions.

## 5. Outcomes and failure-aware reporting

### Primary density endpoint

For a successful density estimate (hat f_{r,n,a}(z)), report integrated absolute error against the declared analytic target:
[
mathrm{IAE}_{r,n,a}=
int_{z_{min}}^{z_{max}}
|hat f_{r,n,a}(z)-f_star(z)|,dz.
]

Also retain ISE, maximum absolute density error, pointwise density error, normalization residual, and weight-concentration diagnostics. Density uncertainty intervals, if added, must distinguish posterior/fit uncertainty within a replicate from repeated-sampling variation across independent replicates.

### Paired contrasts

For each n, degree, curvature and bandwidth condition, define
[
Delta_{r,n,d,gamma,h}
=mathrm{IAE}_{mathrm{estimated IPW},r,n,d,gamma,h}
-mathrm{IAE}_{mathrm{catalogue only},r,gamma,h}.
]
Negative values favour estimated IPW on this density endpoint; positive values favour catalogue-only. The pairing is by replicate and frozen catalogue, not by individual galaxy or validation-frame object.

Report the number of planned pairs, successful pairs, and the complete distribution of paired differences. Do not treat validation units, galaxies, KDE grid points, or coefficient draws as independent experimental replicates.

### Failure ledger

Each planned arm/replicate/condition must have exactly one status:
- `SUCCESS`: valid fit and density estimate;
- `NON_EVALUABLE`: the data/design do not identify the requested estimator under the declared rules (for example, one-class validation labels or rank deficiency);
- `FAILED`: an implementation/numerical failure after a valid design was supplied;
- `NOT_RUN`: planned attempt not executed.

Do not silently drop failed attempts. Report per-arm and per-condition success probabilities and failure counts separately from
(E[mathrm{IAE}midmathrm{SUCCESS}]). Do not assign an arbitrary numeric IAE penalty to failures unless an operational loss is justified and predeclared. A comparison conditioned on both arms succeeding must be labelled as such and accompanied by the unmatched failure rates.

The auditor must reject a missing factorial attempt, duplicate key, impossible status transition, density attached to an unsuccessful attempt, or successful estimate lacking its required raw inputs.

## 6. Statistical analysis plan

The confirmatory analysis must be frozen before running the large experiment.

- Declare one primary DGP contrast and endpoint before execution. Since earlier exploratory pilots have already examined positive curvature and several bandwidths, any choice motivated by those pilots must be explicitly labelled **pilot-informed** rather than independent confirmation.
- Treat the remaining curvature × validation-size × model-degree × bandwidth contrasts as secondary.
- Publish all planned estimates and intervals, not only favourable cells.
- Use replicate-level paired uncertainty. A paired bootstrap must resample replicate IDs, carrying all arms and conditions for a selected replicate together.
- For a family of secondary contrasts, either report simultaneous intervals or apply a declared family-wise/multiple-testing procedure (for example Holm adjustment). State whether intervals are pointwise or simultaneous.
- Separate effect-size uncertainty from Monte Carlo uncertainty in any numerical integration or bootstrap calculation.
- Report sensitivity to the number of independent replicates; do not mistake a narrow conditional standard error based on many galaxies for high precision across simulated universes.

## 7. Independent audit contract

The independent auditor must not import the benchmark runner, production inclusion fitter, or production KDE wrapper. It should independently reconstruct:

1. frame sizes, unique IDs, and independence/disjointness of analysis and validation frames;
2. nested validation subsets and the frozen analysis-frame identity across n;
3. retention labels from stored uniforms and reconstructed DGP probabilities;
4. the retained analysis catalogue and its exact equality across validation-size conditions;
5. validation-only logistic score, objective, fitted probabilities, observed information and covariance diagnostics for successful MLE fits;
6. analysis-catalogue fitted probabilities and inverse-weight diagnostics;
7. each density vector using an independently implemented declared KDE protocol;
8. density normalization, IAE, ISE, maximum error and paired replicate-level contrasts;
9. the full factorial attempt ledger and status counts; and
10. manifest configuration hash, software metadata, byte counts and SHA-256 checksums.

As with the existing auditor, sharing SciPy's `gaussian_kde` primitive is a known limitation, not full algorithmic independence. Record that limitation in every audit report.

Adversarial tests should modify a fitted coefficient and then refresh the artifact checksum; the auditor must still reject the artifact through score/objective/covariance or density reconstruction. Additional attacks should alter validation ordering, change a frozen analysis-catalogue membership across n, duplicate an attempt, remove an attempt, and tamper with a stored contrast while refreshing the manifest.

## 8. Relationship to the existing disjoint finite-parent design

The two designs answer different questions:

| Feature | Fixed-analysis-catalogue design | Disjoint finite-parent design |
|---|---|---|
| Analysis catalogue across n | Identical within replicate/scenario | Changes as validation units are removed |
| Validation sample | Independent frame | Sampled from same parent frame |
| Main estimand | Information-size effect conditional on fixed catalogue | Operational validation/catalogue trade-off |
| Main confounding | Does not confound n with analysis-catalogue size | Validation size and catalogue size change together |
| Operational interpretation | Requires a separate validation frame | Represents finite-frame partitioning |
| Appropriate synthesis | Separate result panel | Separate result panel |

Do not pool the two designs into one regression or estimator ranking without an explicit design indicator and interaction model. Even then, interpret the design interaction rather than pretending they are exchangeable replicates.

## 9. Gate before H0 propagation

Passing the FAC density benchmark establishes only that the declared density-recovery experiment and its numerical artifact are internally reproducible. It does not establish posterior calibration.

Only after the FAC and disjoint designs have been independently audited should the project proceed to coefficient-uncertainty propagation and H0 inference. The latter requires a proper coefficient posterior (not the inverse-information covariance relabelled as a posterior), per-draw density normalization, the selected fixed-count likelihood evaluated separately for each density draw, log-sum-exp likelihood marginalization, and a separate audit of GW detection selection. Grid convergence, density recovery, fixed-truth coverage and simulation-based calibration remain distinct validation gates.

## 10. Pre-implementation checklist

- [ ] Freeze the parent intensity, inclusion law, curvature center and retention calibration.
- [ ] Freeze analysis-frame size, validation-frame maximum size, seeds and nested-order rule.
- [ ] Prove analysis catalogue equality across validation sizes in tests and raw artifact.
- [ ] Add a separate runner/config/artifact namespace; do not alter the disjoint-design artifact schema in place.
- [ ] Add explicit status ledger and failure-aware summaries.
- [ ] Pre-register primary/secondary contrasts and multiplicity handling.
- [ ] Add independent auditor and coefficient/partition/attempt tampering tests.
- [ ] Run focused tests and a small smoke artifact first.
- [ ] Run the confirmatory replicate count only after smoke audit passes.
- [ ] Keep H0 inference out of this stage.
