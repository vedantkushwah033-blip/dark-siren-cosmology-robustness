# Adaptive-bandwidth IPW-KDE first-order perturbation

**Status:** analytic derivation and finite-difference regression specification.  
**Scope:** the bounded reflected Gaussian KDE in `catalogue_density.py`; synthetic methodology only.

## 1. Production bandwidth

Let catalogue locations be x_i and positive baseline inverse-inclusion weights w_i. Write normalized weights a_i=w_i/Σw_i. The production implementation uses

```
mu = sum_i a_i x_i
v  = sum_i a_i (x_i-mu)^2
ESS = 1 / sum_i a_i^2
h = sqrt(v) ESS^(-1/5)
```

for its default Scott-style physical bandwidth. A supplied `bandwidth_factor` replaces `ESS^(-1/5)` with a fixed multiplier; that case is not adaptive through ESS.

Consider a local inclusion-probability perturbation whose unnormalized inverse weights follow
```
w_i(epsilon) = w_i / (1 + epsilon delta_i).
```
At epsilon zero, with `bar_delta = Σ a_i delta_i`,
```
a_i' = -a_i (delta_i - bar_delta).
```
The derivative is invariant to a common additive shift in delta, as expected after weight normalization.

## 2. Derivative of the adaptive bandwidth

The weighted mean derivative is
```
mu' = -Cov_a(x, delta).
```
Differentiating the centered weighted variance gives
```
v' = -Cov_a((x-mu)^2, delta).
```
Since ESS=(Σa_i²)^(-1),
```
ESS'/ESS = 2 { Σ a_i² delta_i / Σ a_i² - Σ a_i delta_i }.
```
Consequently,
```
h'/h = (1/2) v'/v - (1/5) ESS'/ESS.
```
This term is absent from the fixed-bandwidth calculation. It can be material when the perturbation changes both weighted spread and weight concentration.

## 3. Reflected KDE derivative

For bounds L,U, define the three reflected locations for each observation:
```
r_i0=x_i, r_iL=2L-x_i, r_iU=2U-x_i.
```
Let
```
K_h(t) = exp(-t^2/(2h^2)) / (sqrt(2 pi) h)
G(z) = (1/3) sum_i a_i sum_{r in {r_i0,r_iL,r_iU}} K_h(z-r).
```
Using `d log K_h(t)/dh = -1/h + t^2/h^3`,
```
G'(z) = (1/3) sum_i sum_r [
    a_i' K_h(z-r)
    + a_i K_h(z-r) (-1/h + (z-r)^2/h^3) h'
].
```
For fixed positive merger-rate multiplier m(z), the final normalized density is
```
p(z) = m(z)G(z) / integral m(u)G(u)du.
```
Its first derivative is
```
p'(z) = p(z) [q(z) - E_p q],  q(z)=G'(z)/G(z).
```
The expectation and normalization must use the same finite-grid quadrature rule as the production density.

## 4. Interpretation and boundaries

- This is a derivative conditional on the observed catalogue locations and the supplied perturbation direction `delta_i`; it does not integrate over catalogue-generation or validation-sample randomness.
- It includes the production default adaptive bandwidth, both reflected boundaries, merger weighting, and final finite-grid normalization.
- It does not differentiate logistic coefficient estimation itself. A fitted inclusion-model perturbation must first be mapped to the observation-level `delta_i`.
- The formula assumes positive weights, nonzero weighted variance, and a locally differentiable bandwidth branch. It does not regularize positivity or cure extreme-weight concentration.
- The fixed-bandwidth regression remains a distinct limiting/controlled check.
- A passing finite-difference test validates this local implementation derivative only; it is not evidence of H0 calibration or estimator superiority.
