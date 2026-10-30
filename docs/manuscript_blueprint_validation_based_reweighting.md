# Manuscript blueprint: validation-based catalogue reweighting

**Working title:** *Finite-Validation Uncertainty and Selection-Model Misspecification in Inverse-Probability Host Weighting: A Controlled Dark-Siren Simulation Benchmark*

**Status:** manuscript plan only. No publication-ready manuscript or confirmatory result is claimed.

## One-sentence contribution

A reproducible, paired simulation benchmark quantifying how finite parent-frame validation samples, selection-model misspecification, and weak positivity jointly affect inverse-probability host-density reconstruction and downstream synthetic dark-siren posterior behaviour.

This is a candidate contribution, not a novelty claim. The literature audit must be updated as the study develops.

## Abstract structure (write only after confirmatory results)

1. **Context:** incomplete galaxy catalogues complicate host-redshift inference in dark-siren cosmology.
2. **Gap, narrowly stated:** the paper studies a specific external parent-frame validation design and quantifies the interaction of validation size, selection-model misspecification, and positivity in a controlled one-dimensional simulation. Do not say that no prior work studies incompleteness or weighting.
3. **Methods:** paired synthetic parent/catalogue/GW realizations; validation-based logistic inclusion estimation; catalogue-only and oracle controls; prespecified misspecification and positivity surfaces; validation-resampling sensitivity; independent artifact reconstruction.
4. **Results:** report only values from completed, archived, independently audited confirmatory artifacts. Include effect sizes and uncertainty, not only p-values.
5. **Limitations:** redshift-only population, simplified distance/selection likelihood, assumed validation sampling design, logistic family, no realistic sky geometry, no real LVK inference.
6. **Conclusion:** identify empirically demonstrated operating/failure regions; do not generalize beyond the simulator.

## Research questions

**RQ1 — finite validation:** How does validation sample size affect inclusion-probability error, host-density error, weight stability, and downstream posterior calibration when the fitted selection family is correct?

**RQ2 — model misspecification:** How quickly do those quantities deteriorate when the true retention law departs from the fitted linear-logit family?

**RQ3 — positivity:** Which observable diagnostics flag regimes in which inverse-probability weighting becomes unstable, and how often do those alerts precede downstream inference failures?

**RQ4 — uncertainty:** How much variation in posterior summaries is attributable to validation-record resampling conditional on the fixed catalogue and GW observations, and how much error remains because the assumed selection family is wrong?

RQ4 must explicitly separate sampling variability from structural model bias. A bootstrap within a misspecified logistic family cannot quantify the latter.

## Hypotheses to preregister

- **H1:** Under correct specification and adequate positivity, increasing validation size reduces selection-probability estimation error on average.
- **H2:** Under misspecification, increasing validation size does not eliminate asymptotic model bias; it may reduce variance while leaving systematic error.
- **H3:** Poor positivity manifests in lower ESS and heavier inverse-weight tails, and can increase host-density and posterior instability.
- **H4:** Validation bootstrap variation decreases with larger validation samples under a stable correctly specified model, but does not by itself reveal structural misspecification.
- **H5:** Observable-only pre-posterior warnings can identify some unstable regimes on held-out simulation replicates; performance must be assessed against simple baselines and with false-alert rates.

These are hypotheses, not findings.

## Methods requirements

### Generative model
Document the parent-redshift distribution, merger-rate weighting, retention probability, validation sampling design, GW-distance model, detection rule, and all priors. State clearly that the current simulator is one-dimensional in redshift.

### Estimators
Define each arm mathematically and in code:
- catalogue-only;
- oracle inverse-inclusion weighting (diagnostic upper-bound control);
- estimated logistic inverse-inclusion weighting;
- any flexible selection model only after a prespecified implementation and independent tests;
- any stabilized/truncated estimator as a separate, named arm with fixed threshold and sensitivity analysis.

No clipping or probability floor may be silently introduced.

### Pairing
Use the same parent population, retention realization, validation-frame ordering, and GW events across comparable arms. Record deterministic seed derivation. Contrasts must be paired by replicate ID.

### Outcomes
Primary:
- paired difference in absolute posterior-mean error between estimated weighting and catalogue-only.

Key calibration outcome:
- 68% fixed-truth credible-interval coverage with Wilson interval and fit-failure accounting.

Secondary:
- posterior bias, RMSE, interval width;
- host-density integrated squared error;
- inclusion-probability MAE on parent truth (scoring only);
- ESS, maximum/p95/p99 inverse weights, low-probability mass;
- optimizer boundary and failure rates;
- validation-bootstrap sensitivity distributions;
- pre-posterior warning metrics on held-out replicates.

Coverage, precision, and accuracy are different estimands. A narrower interval is not automatically better.

## Required figures

1. **Generative and information-flow diagram:** separate observed catalogue, validation records, inference payload, and scoring-only truth.
2. **Selection curves:** true and fitted inclusion probabilities across redshift, with validation observations.
3. **Positivity diagnostics:** ESS fraction and inverse-weight tails against selection severity.
4. **Host-density reconstruction:** true density shown only for evaluation, alongside catalogue-only, estimated, and oracle estimates.
5. **Posterior calibration:** coverage with Wilson intervals and posterior error/width distributions.
6. **Misspecification surface:** bias/RMSE/coverage across validation size and curvature.
7. **Bootstrap sensitivity:** within-replicate posterior-summary spread versus validation size, clearly distinguished from posterior credible intervals.
8. **Warning evaluation:** held-out ROC/precision-recall and calibration plots, plus false-alert rates.

Every figure should include denominators, failed-run counts, and uncertainty where relevant.

## Required tables

- Complete factor grid and simulation parameters.
- Attempted/successful/failed fits by arm and cell.
- Primary paired effects and uncertainty intervals.
- Coverage and Wilson intervals by arm and cell.
- Positivity and weight-tail diagnostics.
- Bootstrap success/failure accounting.
- Reproducibility manifest: commit, Python/dependency versions, seeds, artifact hashes, audit status.

## Threats to validity

1. **External validity:** one-dimensional redshift simulation is not a 3D sky catalogue.
2. **Likelihood realism:** the toy GW distance and detection model do not reproduce full LVK event-level likelihoods or detector selection.
3. **Validation design:** simple random parent-frame validation with retention labels may not match actual survey audit data.
4. **Selection-model risk:** a logistic family can be wrong; within-family bootstrap cannot remove this bias.
5. **Positivity:** inverse weights can be unstable where inclusion probabilities are small.
6. **Host population:** merger-rate weighting is specified rather than inferred jointly with cosmology.
7. **Monte Carlo uncertainty:** small replicate counts produce imprecise coverage and contrast estimates.
8. **Researcher degrees of freedom:** factor selection, thresholds, and primary outcomes must be frozen before confirmatory runs.

## Release gates before submission

- All tests pass on Python 3.10, 3.11, and 3.12 at the exact release commit.
- Confirmatory runs have fixed seeds, frozen configuration, and a complete runtime manifest.
- Independent auditors reconstruct every headline statistic from raw replicate records.
- All failure rows and exclusions are reported.
- Literature and novelty audit is updated and independently checked.
- A domain expert reviews the astrophysical assumptions and interpretation.
- A statistician reviews the estimands, calibration analysis, bootstrap, and power/precision plan.
- The paper and repository make clear that the results are synthetic and not a new H0 measurement.

## Positioning and venue strategy

Do not choose a journal based on prestige alone. First complete a defensible methods paper and obtain independent expert feedback. A broad statistical-methods journal would require a contribution that generalizes beyond this one application; a gravitational-wave/cosmology venue would require stronger astrophysical realism and direct comparison with established pipelines. A student project can be valuable for research applications and academic conversations without being accepted by a top journal.

Admissions committees at highly selective universities do not award admission merely for a GitHub repository, a preprint, or a claimed publication. The credible signal is demonstrated intellectual ownership, sound methods, transparent limitations, reproducible work, and serious external feedback. No project can guarantee admission to Caltech, Princeton, Stanford, or any other institution.
