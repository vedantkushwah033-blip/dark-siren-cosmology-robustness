# Prior-support × hidden-host factorial protocol

## Question and estimands

This prespecified synthetic experiment tests whether the effect of prior support interacts with hidden-host population misspecification. It extends two previously separate lines of diagnostics: fixed-truth prior-boundary sensitivity and paired hidden-host misspecification.

Primary estimands are, within each experimental cell:
1. the paired change in RMSE and empirical 68%/95% coverage when a hidden host is modelled versus omitted;
2. the paired change in those same metrics under prior supports [40,100] and [30,120] relative to [50,90], for each inference arm;
3. whether the direction/magnitude of those changes differs across hidden-host fractions, noise levels, and injected H0 values.

This is a synthetic-methodology experiment. It does not estimate the observed Hubble constant and cannot validate astrophysical assumptions by itself.

## Frozen design

- Injected H0: 55, 70, 85 km s^-1 Mpc^-1.
- Fractional distance-noise scale: 0.05 and 0.20 times c × 0.03 / H0.
- Hidden-host fraction: 0, 0.30, 0.45.
- Hidden-host redshift: 0.025 and 0.035.
- Prior supports: [50,90], [40,100], [30,120].
- Inference arms: catalogue weights 0.7/0.3; reversed weights 0.3/0.7; equal weights 0.5/0.5; and a hidden-host mixture with the simulated hidden fraction/redshift.
- 36 shared-data cells, 200 replicates per base cell, four inference arms and three prior supports.
- 7,200 shared simulated observations; 86,400 posterior fits; 432 cell/prior/arm summaries.
- Posterior grid: 4,001 points spanning H0=25–125.
- Uniform prior density on each stated interval; likelihood terms conditioned on positive observed distance, matching the generator's redraw rule.
- Data seed: 20261030. Bootstrap seed: 20261031. 1,000 paired bootstrap resamples per base cell, with replicate indices shared across priors and inference arms.

The primary prior contrast is the equal-cell paired difference in a metric under a wider prior versus [50,90], stratified by the base cell. The primary host contrast compares omitted versus modelled hidden population within each prior. Bootstrap intervals describe finite Monte Carlo uncertainty conditional on this generator and design; they do not represent model uncertainty or uncertainty about the appropriate prior.

## Negative controls and audit gates

- At hidden fraction zero, the base-weight and hidden-modelled arms must yield identical posterior quantiles for every prior support.
- Each prior support and inference arm uses exactly the same observation within a replicate.
- The independent auditor regenerates observations from the fixed seed, reconstructs each posterior directly from the likelihood equation, and independently recomputes summaries and both families of paired-bootstrap intervals.
- CI checks exact row counts, finite numeric values, interval ordering, and the independent audit. An adversarial unit test corrupts a prior-contrast interval and requires audit failure.
- All reported findings must be checked against Monte Carlo uncertainty and cell heterogeneity; no universal direction is inferred from a marginal average alone.

## Interpretation guardrails

The [50,90] prior includes truths 55 and 85 near its lower and upper boundaries, respectively; it does not exclude either truth. Its support can nevertheless truncate posterior mass near those boundaries, while the wider priors alter posterior support and interval geometry. Wider priors are not automatically superior; results can change because the prior is part of the statistical model. The hidden-host model is a deliberately simplified discrete mixture and assumes the hidden fraction and redshift are known to the inference procedure. This is a hybrid oracle-style diagnostic, not a realistic catalogue-completeness model.

The experiment omits realistic GW likelihoods, detector selection, peculiar velocities, galaxy-catalogue completeness uncertainty, and full hierarchical population inference. Do not use it to claim an observational H0 result, publication readiness, or validation against real gravitational-wave data.
