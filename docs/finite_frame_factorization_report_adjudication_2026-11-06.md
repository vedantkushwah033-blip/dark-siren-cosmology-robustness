# Adjudication of the finite-validation factorization report (2026-11-06)

## Status

This note audits the supplied *Structured Literature & Mathematical Report: Finite-Validation Logistic Inclusion Estimation and Its Consequences for Inverse-Probability Host Reconstruction* against the current synthetic pilot implementation. It does not approve a new likelihood or claim that the current dark-siren posterior has been validated.

**Decision: the report is useful as a list of candidate frameworks, but Part B is not yet a valid single joint factorization.** It combines incompatible validation-design statements, does not define the GW host assignment, and uses a finite-frame notation alongside a superpopulation host-density likelihood without specifying the transition between them.

## 1. What the current code actually generates

For each replicate, the pilot:

1. Generates a finite parent array with redshifts (z_1,\ldots,z_N).
2. Generates catalogue retention independently conditional on those redshifts:
   [
   R_i\mid z_i,\theta_* \sim \mathrm{Bernoulli}(\pi_*(z_i)),
   \qquad C=\{i:R_i=1\}.
   ]
3. Independently permutes the (N) parent-unit indices and selects the first (m) as the validation sample. This is SRSWOR, not independent Bernoulli sampling:
   [
   S\mid N,m \sim \text{uniform over all subsets of size }m.
   ]
4. Records ((z_i,R_i)) for selected validation units and persists parent-array indices for catalogue units, validation units, and their exact intersection.

The validation sample is therefore selected independently of (R) in the simulation design, but its observed outcomes are not independent of catalogue records: a unit in the intersection has the same single realized (R_i) and the same (z_i) in both data views. The implementation's identity metadata makes this overlap observable; it does not itself determine an inferential likelihood.

## 2. Correct finite-frame sampling statement

Conditional on a fixed, fully specified parent redshift frame (Z=(z_1,\ldots,z_N)), and assuming the validation design is independently randomized SRSWOR, the *complete-data* design law is
[
p(R,S\mid Z,\theta)
=
\left\{\prod_{i=1}^N
\pi_i(\theta)^{R_i}[1-\pi_i(\theta)]^{1-R_i}
\right\}
\frac{\mathbf 1\{|S|=m\}}{\binom Nm}.
]

This is the correct joint law under those assumptions. It is not
[
\prod_i p(R_i\mid z_i,\theta)\,p(S_i\mid R_i,z_i)
]
with independent Bernoulli (S_i), unless the validation design really is independent Bernoulli sampling. The report's Part B factorization and its later SRSWOR section therefore describe different designs and cannot both be used as written.

For observed records, let (U=C\cup V_{\rm ids}) be the union of catalogue and validation unit IDs. Conditional on the fixed known frame (Z), the observed retention indicators contribute one Bernoulli factor per distinct unit in (U), not one factor per appearance in a data table:
[
L_{R}(\theta;R_U,Z)
=
\prod_{i\in U}
\pi_i(\theta)^{R_i}[1-\pi_i(\theta)]^{1-R_i}.
]
The SRSWOR design contributes (1/\binom Nm), which is constant in (\theta) when (N,m) are fixed. Units outside (U) have unobserved (R_i); integrating their Bernoulli outcomes gives one.

This expression is a conditional finite-frame likelihood only when the analyst conditions on the full (Z) and knows the relevant unit frame. That is a simulation-oracle analysis, not automatically the likelihood available to an observer who sees catalogue redshifts only for retained units and validation redshifts only for sampled units. An observer-facing likelihood must also model or integrate over unobserved parent redshifts and the observation process.

## 3. Overlap: what is and is not double-counting

The report is right that the same unit's (R_i) must not be entered twice as if it were two independent Bernoulli outcomes. Deduplicating by parent-unit ID is necessary in a full joint likelihood.

However, the report overstates this as a blanket objection to every modular/cut construction. A cut posterior is deliberately a modular inferential target, not the factorization of the ordinary joint posterior. It may intentionally use a validation module to construct (p_{\rm cut}(\theta\mid V)), then pass draws into a catalogue/GW module while blocking feedback. That choice does not become a full joint Bayesian posterior merely by calling the modules conditionally independent.

The important distinction is:

- **Full joint posterior:** the shared unit IDs and shared (R_i,z_i) must appear once in a single observed-data likelihood. Dependence induced by overlap is modeled explicitly.
- **Cut/modular target:** define the validation module and downstream module as an explicit inferential procedure, including exactly what is conditioned on and what feedback is blocked. Do not describe its product form as the ordinary factorization of (p(C,V,E\mid\theta,H_0)).
- **Independent validation and catalogue samples:** this is a different design and cannot be assumed for the current pilot.

Also, the statement that catalogue information necessarily must update \(\theta\) is not correct for a cut target; blocking that feedback is the defining decision. Conversely, claiming that the validation-only posterior is the full conditional (p(\theta\mid C,V)) is incorrect when the cut intentionally omits catalogue feedback.

## 4. GW likelihood: the report's Part B expression is incomplete

The Part B generative list writes (d_j\mid z_j,H_0), but does not define how each latent GW host (z_j) is drawn from the finite parent frame, whether hosts are sampled with merger-rate weights, or how repeated hosts are handled. It then writes a selection-normalized expression whose numerator and denominator use a continuous (p_\theta(z)). Those are not yet a single finite-frame generative model.

For the existing *toy fixed-count* likelihood, the downstream conditional likelihood is instead explicitly defined at the density level:
[
L(H_0;\theta,C)
=
\frac{\prod_{j=1}^{n}
\int p_{\theta,C}(z)\,
p(d_j\mid z,H_0)\,dz}
{\beta_{\theta,C}(H_0)^n},
\quad
\beta_{\theta,C}(H_0)
=
\int p_{\theta,C}(z)\,
p_{\rm det}(z,H_0)\,dz.
]
Here (p_{\theta,C}) is the reconstructed host density supplied to the toy likelihood. This is a conditional-on-(n), selection-normalized density-level model. It is not an event-rate likelihood, and it does not by itself establish that the density estimator is a correct finite-frame posterior predictive distribution.

If the scientific target is a draw-wise modular average, the full conditional likelihood ratio must be evaluated for each inclusion-parameter draw before averaging:
[
\bar L(H_0)
=
\int
\frac{\prod_j I_j(H_0,\theta)}
{\beta(H_0,\theta)^n}
\,p_{\rm mod}(\theta\mid V)\,d\theta.
]
Averaging numerator and denominator separately is generally not equivalent. This mathematical point in Part D is sound.

## 5. Specific claims that need correction or qualification

1. **“SRSWOR likelihood is a combinatorial term, not a product of Bernoulli terms.”** Correct for the selection indicators (S). It does not remove the Bernoulli retention likelihood for (R). Keep the two mechanisms separate.
2. **“The same galaxy contributes information to both catalogue and validation.”** Correct. In a full joint likelihood, its observed retention outcome is represented once, with the joint observation map and shared ID.
3. **“The cut posterior is coherent under the stated assumptions.”** Too strong without specifying the modules, their conditioning sets, and the target distribution. It is a defined modular inferential target; it is not automatically the ordinary posterior from the stated Part B factorization.
4. **“The validation sample is representative (MAR).”** The current pilot's SRSWOR design is stronger and more specific than a generic MAR label. The design should be stated as SRSWOR independent of retention conditional on the generated parent frame. Nonresponse or missing redshift measurement is not currently represented by this statement.
5. **“The finite-frame factorization is exact.”** Only under the fully specified conditioning and observation assumptions. The displayed report equation is not exact as written because it mixes SRSWOR with per-unit conditional selection factors and leaves the GW host-generating mechanism unspecified.
6. **“The recommended primary framework is the cut posterior.”** This is a research choice, not a mathematical consequence. It needs a declared rationale and an empirical comparison with fixed-(\theta), plug-in, and design-based repeated-validation procedures. A full joint model should not be built until its observation model is identified.

## 6. Required next experiment (small, matched, no large campaign)

Before implementing a new sampler, create one small synthetic experiment with the same generated parent frame, retention vector, GW data, and validation ID draw across all compared procedures. Report the procedures as distinct estimands:

- **A. Current plug-in modular baseline:** existing estimated inclusion fit and host-density reconstruction, with its present conditional interpretation stated plainly.
- **B. Validation-design repeated-sampling procedure:** repeat SRSWOR draws from the same fixed parent frame, refit the inclusion model, rebuild the host density and recompute the H0 summary. Report bias/dispersion/coverage across outer replicates; do not call the bootstrap distribution a posterior.
- **C. Finite-frame oracle diagnostic:** use the known simulated parent frame and true inclusion law only as a diagnostic upper-reference, not as an implementable estimator.
- **D. Candidate Bayesian validation module:** only after selecting a proper prior and a likelihood that respects SRSWOR and the observed binary outcomes. Compare it as a new method, not as an already validated correction.

The first comparison should isolate overlap/design effects, not combine them with nonlinear misspecification, solver version, bandwidth changes, and prior sensitivity all at once. Then add those stress dimensions one at a time.

Predeclare the outputs: host-density integrated squared error, posterior mean error, interval coverage and width, fit failure rate, and Monte Carlo uncertainty. Retain failed fits in the denominator for operational reliability; separately report conditional-on-success calibration.

## 7. Decision and blockers

- The current artifact-level overlap IDs are a useful provenance improvement.
- The uploaded report does not yet supply a valid complete joint likelihood for the current code.
- The existing code's SRSWOR selection design is verified from its permutation-based sample construction; the report's independent-Bernoulli validation factorization must not be applied to that arm.
- No Bayesian prior for \(\theta\), draw-wise density reconstruction, or draw-wise selection-efficiency integration has been implemented merely by writing these equations.
- Do not claim that U1 is closed, that a full joint posterior exists, or that propagated inclusion uncertainty has been validated.

**Next gate:** specify the observer-facing versus simulation-oracle estimand and write the exact observation map for (C,V) before coding a Bayesian sampler. Until then, continue with the current repeated-SRSWOR synthetic benchmark and clearly label it as a design-based sensitivity analysis.
