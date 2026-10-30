# Finite-parent target estimand and audit boundary

**Date:** 30 October 2026  
**Status:** implemented in the finite-parent SRSWOR runner; exact-head CI verification pending.

## Why a second target is necessary

The generating value of $H_0$ and the finite-parent target answer different questions. The former measures performance against the simulation truth. The latter asks how closely an estimator based on a validation sample reproduces the result that the *same declared estimated-inclusion pipeline* would produce if all $N$ paired validation records in this fixed synthetic parent frame were observed.

For scenario $s$, let the actual finite parent be
$\mathcal{F}_s=\{(z_i,R_i):i=1,\ldots,N\}$, with the catalogue and GW-distance sample held fixed. Define

$$
T_{s,N}=\operatorname{E}[H_0\mid D_{\rm GW},\widehat p_{s,N},\mathcal C],
$$

where $\widehat p_{s,N}$ is the bounded logistic inclusion model fitted to all $N$ actual paired records, the inverse-inclusion weighted host density is rebuilt with the declared reflected KDE, and the same selection-aware toy likelihood is used.

For an outer SRSWOR validation sample $S_{s,r}$ of size $n$, the sample-based estimate is $\widehat T_{s,r}$. The finite-frame error is

$$
e_{s,r}=\widehat T_{s,r}-T_{s,N}.
$$

This is **not** the error relative to the generating $H_0$, and it is **not** a bootstrap pseudo-population reference. Both targets should be retained and interpreted separately.

## Implementation

The finite-parent runner now:
1. regenerates one fixed synthetic parent, catalogue thinning realization, and GW-distance sample per scenario;
2. fits the estimated-inclusion pipeline once to all $N$ actual paired parent records;
3. stores the resulting finite-parent posterior mean, median, and 68% interval on each successful validation-draw record;
4. stores the signed validation-draw posterior-mean error relative to the finite-parent posterior mean;
5. reports the mean and RMSE of that finite-frame error by validation size.

The independent auditor separately regenerates the parent frame and catalogue, independently refits the logistic model, reconstructs the reflected KDE and selection-aware posterior using its own numerical implementation, and checks the stored finite-parent posterior mean. It also independently re-derives the finite-frame bias and RMSE aggregates. A tampered target must fail the audit.

## Interpretation boundary

- $T_{s,N}$ is a finite-frame **plug-in reference**, not the generating truth and not a claim that the logistic/KDE model is correctly specified.
- Since the fitted inclusion model is linear-logit while the data generator can include a quadratic term, the finite-frame reference can retain model misspecification. That is intentional: the estimand measures validation-sampling variability conditional on the declared estimator and fixed finite frame.
- Catalogue and GW distances are held fixed within scenario. Scenario-level variation and validation-draw variation must not be conflated.
- The current benchmark is a smoke/pilot runner, not a confirmatory comparison of bootstrap methods. It does not yet establish whether ordinary paired bootstrap, a finite-population correction, or the design-aware pseudo-population bootstrap is calibrated for this nonlinear posterior functional.
- Conditional $H_0$ credible intervals remain distinct from repeated-sampling uncertainty in $\widehat T_{s,r}$.

## Next gate

Use the stored exact validation indices and scenario provenance to run a nested comparison: outer independent SRSWOR samples; within each outer sample, paired-record bootstrap and design-aware pseudo-population bootstrap; then compare each estimated bootstrap SE and interval with the empirical outer-sample distribution of $\widehat T_{s,r}-T_{s,N}$. Freeze the scenario hierarchy, seeds, failure rules, and Monte Carlo precision targets before confirmatory runs. Do not infer theoretical validity from a small smoke run.
