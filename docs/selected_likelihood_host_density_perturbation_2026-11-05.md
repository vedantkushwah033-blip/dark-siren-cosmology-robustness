# First-order host-density perturbation in the selected toy H0 likelihood

**Status:** deterministic mathematical regression added; CI verification pending.  
**Scope:** synthetic toy likelihood only. This is not a detector model or an observational H0 result.

## 1. Production likelihood audited

The implementation in `src/dark_siren_robustness/selection_toy.py` computes, on the supplied redshift grid,

```
I_j(H) = integral p(z) Normal(d_j | c z/H, sigma) dz
beta(H) = integral p(z) Phi((c z/H - tau)/sigma) dz
ell(H) = sum_j log I_j(H) - n log beta(H).
```

The implementation uses the same supplied normalized host density in every event integral and in the detection-efficiency integral. The event count is fixed. This is not an event-rate likelihood.

## 2. Perturbation and derivative

Let a differentiable normalized density path satisfy

```
p_epsilon(z) = p(z) [1 - epsilon g(z)] + O(epsilon^2),
integral p(z) g(z) dz = 0.
```

For an event-specific normalized redshift distribution

```
r_j(z | H) = p(z) q_j(z | H) / I_j(H),
q_j(z | H) = Normal(d_j | c z/H, sigma),
```

and the detection-conditioned redshift distribution

```
r_sel(z | H) = p(z) s(z | H) / beta(H),
s(z | H) = Phi((c z/H - tau)/sigma),
```

the fixed-H first derivative is

```
d ell(H; epsilon) / d epsilon at epsilon=0
  = - sum_j E_{r_j(.|H)}[g(z)]
    + n E_{r_sel(.|H)}[g(z)].
```

The second term is required because the selection normalizer changes with the host density. Omitting it differentiates only the unconditional event-integral component, not the selected likelihood implemented by the repository.

For a raw relative inclusion-probability perturbation delta, normalized IPW induces, at first order, `g = delta - E_p[delta]` in the convention above. In a finite-bandwidth KDE, the corresponding perturbation must first pass through the kernel operator; the fixed-bandwidth reflected-KDE regression covers that separate stage.

## 3. Deterministic numerical validation

Added `tests/test_selected_likelihood_host_density_perturbation.py`. It:

- computes event and selection expectations independently using SciPy probability densities/CDFs and trapezoidal node masses;
- compares the analytic derivative against finite differences of the production `selection_aware_toy_likelihood`;
- uses decreasing perturbation magnitudes and checks the first-order remainder decreases;
- separately asserts the selection-normalizer contribution is nonzero in a deliberately chosen example.

An independent calculation of the first fixture produced maximum absolute derivative discrepancies:

| epsilon | max absolute discrepancy |
|---:|---:|
| 0.0020 | 0.0006087 |
| 0.0010 | 0.0003042 |
| 0.0005 | 0.0001521 |

The approximately halving discrepancy is consistent with an O(epsilon) error after dividing the finite difference by epsilon. These values are a deterministic check of the stated discrete-grid calculation, not repeated-sampling performance estimates.

## 4. Limitations

- The derivative is conditional on fixed observed distances, fixed H0 values, fixed redshift grid and quadrature rule.
- It does not propagate uncertainty in the validation sample or fitted inclusion coefficients.
- It does not differentiate the adaptive KDE bandwidth, which depends on estimated weights and weighted moments.
- It does not establish calibration or robustness of the overall estimator under model misspecification.
- The detection rule in this toy model is a declared synthetic threshold, not a realistic gravitational-wave detection pipeline.
- CI for the newly added tests must pass on the exact commit before treating them as repository-verified tests.

## 5. Reproducible test commands

```bash
python -m pytest tests/test_ipw_kde_first_order_perturbation.py -q
python -m pytest tests/test_selected_likelihood_host_density_perturbation.py -q
```
