# Independent likelihood validation protocol

## Purpose

This protocol checks the numerical implementation of the simplified low-redshift host-mixture posterior against an independently written reference calculation. It does **not** validate the astrophysical assumptions or turn the toy likelihood into a complete dark-siren analysis.

## Model being checked

For candidate host redshift $z_j$, host weight $w_j$, observed distance $d$, uncertainty $\sigma$, and trial Hubble constant $H_0$, the toy likelihood is

$$
L(H_0) = \sum_j \tilde w_j\,\mathcal N\!\left(d\mid \frac{c z_j}{H_0},\sigma\right),
\qquad \tilde w_j=\frac{w_j}{\sum_k w_k}.
$$

The posterior is proportional to this likelihood inside the specified uniform-prior interval and zero outside it. The density is normalized numerically over the supplied $H_0$ grid.

## Independent checks added

`tests/test_independent_likelihood_validation.py` compares the production posterior with a direct implementation that loops over hosts and evaluates Gaussian densities in ordinary probability space. Other tests verify:
- permutation invariance when host redshifts and weights are reordered together;
- invariance to multiplying all host weights by a common positive constant;
- the single-host likelihood maximum at $H_0=cz/d$ for the low-redshift distance approximation;
- zero posterior density outside the prior support.

These tests use analytically interpretable or independently calculated cases rather than reusing the production function to generate expected values.

## Numerical limitations

- Direct Gaussian densities can underflow for extreme inputs; the comparison test intentionally uses a moderate case where this is not expected.
- The finite grid introduces discretization error. The single-host mode test uses a fine grid and a stated tolerance.
- A passing numerical comparison shows agreement for the tested cases, not proof for all inputs.
- The model still omits full hierarchical population inference, cosmological-distance integration beyond $d_L\approx cz/H_0$, gravitational-wave detection-selection normalization, and survey-calibrated catalogue selection.

## Interpretation guardrail

Do not use these unit tests as evidence that empirical credible-interval coverage should equal its nominal level. Coverage is a repeated-sampling property that must be estimated in a dedicated simulation study under clearly stated data-generating assumptions. Do not present synthetic results as measurements from real gravitational-wave events.
