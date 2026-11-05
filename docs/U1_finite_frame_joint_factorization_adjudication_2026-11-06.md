# U1 mathematical adjudication: finite-frame catalogue, validation, and GW likelihood

**Date:** 6 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Status:** derivation note and implementation gate; not an implemented posterior and not an empirical result.

## Executive finding

The existing uncertainty protocol correctly warns that its proposed integral is only usable after the conditional GW likelihood is made coherent. The key unresolved issue is more specific than simply choosing “Bayesian” versus “frequentist”:

> The validation sample is drawn without replacement from the same finite parent frame that generated the retained catalogue. Therefore the catalogue and validation observations are generally dependent through shared frame units, and a validation-only inclusion posterior is a deliberate modular construction—not the full posterior conditional on all observed information.

The project should **not** implement a Bayesian sampler until the data-collection experiment, unit identity/overlap handling, and meaning of the host density have been stated as a single probability model.

This note provides the factorization target and identifies what the current code does and does not establish.

## 1. Finite-frame data-generating experiment

Let the finite parent frame contain (N) units, with latent or generated redshifts (Z=(z_1,ldots,z_N)). Let (R_iin{0,1}) indicate catalogue retention. Under the working inclusion model,

[
R_imid z_i,	heta sim operatorname{Bernoulli}(pi_	heta(z_i)),
qquad
pi_	heta(z)=operatorname{logit}^{-1}{x(z)^	op	heta}.
]

Let (S=(S_1,ldots,S_N)) denote the validation-sample indicators. Under fixed-size SRSWOR,

[
P(S=smid N,m)=
inom{N}{m}^{-1}
mathbf 1!left(sum_i s_i=might).
]

The retained catalogue is a deterministic observation of the retained units and their measured attributes. Write (C=g(Z,R)). The validation data are the selected frame records, including redshift and retention outcome, (V=h(Z,R,S)). If stable unit identifiers are available, catalogue-validation overlap is observable and a unit in (Ccap V) must be represented once in a joint data likelihood—not treated as two independent Bernoulli outcomes.

A superpopulation version may place a model (p(Zmidpsi)) on parent redshifts. Conditional on a fixed finite frame, (Z) is instead held fixed. These are different inferential targets and must not be mixed.

## 2. Full finite-frame factorization target

A candidate generative factorization, conditional on the declared frame size, is

[
p(Z,R,S,Emid H_0,	heta,psi)
=
p(Zmidpsi)
left{prod_{i=1}^{N}p(R_imid z_i,	heta)ight}
p(Smid N,m)
p(Emid H_0,Z,R),
]

with the observed catalogue and validation records imposed through their deterministic observation maps:

[
C=g(Z,R),qquad V=h(Z,R,S).
]

Thus, schematically,

[
p(H_0,	heta,psi,Z,R,Smid C,V,E)
propto
p(H_0)p(	heta)p(psi),
p(Zmidpsi)
prod_i p(R_imid z_i,	heta),
p(Smid N,m),
p(Emid H_0,Z,R),
mathbf 1{C=g(Z,R)}
mathbf 1{V=h(Z,R,S)}.
]

This is a **model template**, not a finished likelihood. The observation maps must encode exactly which redshifts, outcomes, and unit identities are observed. The GW term must define how hosts are drawn from the parent frame, how merger-rate weights enter, and how detection is conditioned upon. If redshifts are measured with error, a measurement model is also required.

The validation design factor (p(Smid N,m)) is independent of (	heta) under SRSWOR, but the observed validation outcomes are not independent of catalogue records when units overlap. Under a finite-frame design, the sample is a subset of fixed units; it is not an IID sample of new galaxies.

## 3. What can be learned from the current validation data?

The current validation pilot constructs a random permutation of the parent frame and records redshift and retention outcome for the resulting audit ordering. The protocol describes the validation sample as SRSWOR from the complete parent frame. This is materially different from sampling only retained catalogue objects.

For a model-based inclusion regression, the Bernoulli likelihood contribution for sampled validation units is

[
L_V(	heta;V)=
prod_{i:S_i=1}
pi_	heta(z_i)^{R_i}
{1-pi_	heta(z_i)}^{1-R_i}.
]

This can define a validation-informed posterior (p(	hetamid V)propto L_V(	heta;V)p(	heta)) **under the conditional Bernoulli model and the declared sampling design**. It is not automatically the full posterior given ((C,V)). In particular:

- If catalogue and validation records overlap, deduplicate by unit ID.
- Do not append catalogue positives as independent validation observations.
- Do not assume the catalogue contributes no information in a full generative analysis; whether it does depends on what is conditioned on and on the parent-redshift model.
- If the target is the realized finite-frame inclusion vector rather than a superpopulation logistic law, the parameter (	heta) is a working model parameter and its interpretation changes.
- The parent redshift distribution and the catalogue ascertainment mechanism must be represented if inference uses the observed catalogue redshift distribution to learn population structure.

If the implementation does not preserve unit identifiers, the overlap structure cannot be audited exactly. The protocol must then state whether overlap is impossible by design, ignored as a modular approximation, or unobserved and modeled.

## 4. Modular Bayesian construction: what it means

The proposed construction is

[
widetilde p(H_0mid E,C,V)
propto
p(H_0)
int L_{mathrm{GW,sel}}(H_0;	heta,C,E),
p(	hetamid V),d	heta.
]

This is defensible as a **modular/cut analysis** if all of the following are explicit:

1. (p(	hetamid V)) is derived from a declared validation likelihood and proper prior.
2. The module is intentionally protected from feedback from the catalogue and GW data.
3. (L_{mathrm{GW,sel}}) is a coherent conditional likelihood given the reconstructed host population and does not re-use catalogue-retention information as a second selection correction.
4. The host density for each (	heta) draw is defined by the full map: predicted inclusion probabilities (ightarrow) inverse weights (ightarrow) weighted KDE (ightarrow) merger-rate weighting (ightarrow) normalization.
5. The GW event numerator and detection normalizer use the same draw-specific host density and selection model.
6. The result is described as modular, not as the fully joint posterior conditional on all catalogue and validation information.

This modular target deliberately prevents (C) and (E) from updating (	heta). That may be a defensible protection against feedback through a misspecified parent-population model, but it is a modeling decision—not a consequence of SRSWOR.

## 5. Selection-conditioned GW likelihood

For a normalized pre-detection host-redshift density (p_	heta(z)), distance model (p(d_jmid z,H_0)), and detection probability (s(z,H_0)), define

[
I_j(H_0,	heta)=int p_	heta(z)p(d_jmid z,H_0),dz,
]

[
eta(H_0,	heta)=int p_	heta(z)s(z,H_0),dz.
]

Conditional on a fixed number (n) of detected events, and under conditional independence of event draws,

[
L_{mathrm{GW,sel}}(H_0;	heta)
=
prod_{j=1}^{n}rac{I_j(H_0,	heta)}{eta(H_0,	heta)}.
]

Consequently the modular marginal likelihood is

[
ar L(H_0)=
int
left[
prod_{j=1}^{n} I_j(H_0,	heta)
ight]
eta(H_0,	heta)^{-n}
p(	hetamid V),d	heta.
]

The integral is over the **complete conditional likelihood**. In general,

[
E_	heta[N(H_0,	heta)eta(H_0,	heta)^{-n}]

e
rac{E_	heta[N(H_0,	heta)]}{E_	heta[eta(H_0,	heta)]^n}.
]

The latter is not the same hierarchical marginalization. Nor is an unweighted average of already-normalized conditional (H_0) posteriors equivalent to integrating the likelihood and then normalizing, because component marginal likelihoods differ.

This likelihood conditions on (n). It does not include the information in the event count or infer an event rate. The toy threshold rule is not a realistic detector-selection model.

## 6. Host-density semantics and the double-correction guard

The inspected estimated-inclusion path reconstructs a weighted catalogue density intended to estimate the parent host-redshift distribution, then passes that density into the event integrals and the detection-efficiency calculation. In that path, adding a separate missing-host mixture that represents the same uncatalogued parent mass would double-correct the population.

The correct rule is semantic, not a universal assertion that an (f) parameter must always equal zero:

> If the IPW density already targets the complete parent host population, do not add another component that reintroduces the same missing parent mass.

This must be checked separately for each likelihood arm and density definition.

## 7. Repeated-sampling estimand is separate

For a fixed finite parent frame, fixed retention realization, fixed observed catalogue, and fixed GW dataset, repeatedly draw (S^{(b)}) under the actual SRSWOR design, refit the inclusion model, and recompute the entire density and (H_0) posterior. The resulting distribution estimates sensitivity to the validation-sampling design, conditional on those fixed objects.

Across outer synthetic replicates, one may additionally estimate fixed-truth bias, interval coverage, and operational failure probability. These are repeated-sampling properties of the complete procedure. They are not the posterior distribution of (	heta) or (H_0).

Do not use an IID bootstrap or a scalar finite-population correction as a substitute for the actual SRSWOR design without a separate derivation.

## 8. Decisions required before U1 can pass

| Decision | Required specification | Current status |
|---|---|---|
| Parent frame | Fixed finite frame or superpopulation model | Must be explicit per estimand |
| Unit identity | Stable IDs across catalogue and validation | Must verify in artifacts |
| Validation design | SRSWOR, sample size, nonresponse/measurement assumptions | Protocol says SRSWOR; verify persisted implementation |
| Inclusion model | Working Bernoulli law and coefficient prior | Bayesian prior not yet selected/validated |
| Catalogue role | Conditioned upon, modeled jointly, or modularly cut | Modular option proposed; not full joint |
| Host-density target | Parent pre-detection density or selected density | Must be explicit in every arm |
| GW hosts | Parent-frame draw and merger-rate weighting | Requires source-level generative contract |
| GW selection | Conditional-on-detection likelihood or event-rate model | Current toy code is fixed-count conditional |
| Overlap | Exact overlap representation/deduplication | Requires ID audit |
| Numerical integration | Draw-wise density and (eta(H_0,	heta)) recomputation | Not implemented |

## 9. Recommended near-term decision

Do not begin a fully joint Bayesian sampler now. The necessary frame-level observation model, overlap representation, and host/GW generative factorization have not been established in the inspected implementation.

The most tractable defensible sequence is:

1. **First**, audit whether the persisted validation records retain parent-unit identifiers and establish the exact SRSWOR sampling contract.
2. **Second**, specify a validation-only Bayesian inclusion model with an explicit prior and posterior predictive checks. Treat it as a modular analysis only if the intentional cut from catalogue/GW feedback is justified.
3. **Third**, derive the draw-specific conditional GW likelihood and verify that the reconstructed density represents the parent population and that the same density enters both event numerators and selection normalization.
4. **In parallel**, use repeated SRSWOR validation resampling as a separate design-based sensitivity/calibration analysis, conditional on fixed catalogue and GW data.
5. If the unit identities or likelihood contract cannot be verified, limit the primary result to host-density reconstruction and do not claim complete propagated dark-siren posterior uncertainty.

This is a proposed research decision, not an approval to implement the sampler. U1 remains open until the listed data and likelihood contracts are verified.


## 8. Artifact-level finite-parent unit identity contract (2026-11-06)

The estimated-inclusion pilot now records a zero-based positional `object_id` for
each unit in the generated finite parent frame. These IDs are local to one
replicate's generated parent array; they are not persistent astronomical object
identifiers and do not imply an observed real catalogue.

For every replicate-arm row, the artifact stores:
- `catalogue_object_ids_json`: parent-frame indices with retention indicator 1;
- `validation_object_ids_json`: parent-frame indices selected into that arm's
  validation sample (empty for catalogue-only and oracle arms);
- `overlap_object_ids_json`: the sorted set intersection of those two lists;
- `n_catalogue_validation_overlap`: the corresponding overlap count;
- `unit_overlap_schema = finite_parent_unit_identity_v1`.

The successful inference payload repeats the catalogue, validation, and overlap
IDs. The independent auditor checks the IDs are unique and within the declared
parent frame, recomputes the set intersection, checks the counts, and verifies
that the inference payload agrees with the replicate CSV. Failed rows are also
required to preserve the CSV-level identities.

This closes the artifact-level identity/overlap traceability gap for newly
generated pilot artifacts. It does **not** retroactively add identities to
historical artifacts, prove that validation outcomes are independent of
catalogue observations, or resolve the statistical choice between a modular
validation-informed posterior and a full joint finite-frame likelihood.
The validation sample remains an SRSWOR subset of parent units, and the same
unit may appear in both the observed catalogue and validation sample. Any
likelihood treating these as independent records would need an explicit
dependence/overlap model; these identifiers make that dependence auditable but
do not by themselves specify the correct likelihood.
