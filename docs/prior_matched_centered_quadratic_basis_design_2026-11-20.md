# Prior-matched centered quadratic inclusion basis

## Motivation

The original degree-two pilot used the raw polynomial basis `[1, x, x^2]` with independent coefficient priors. The generator instead uses

`eta(x) = alpha + beta*x + gamma*(x^2 - c_g)`.

The raw-basis coefficients obey `theta0 = alpha - c_g*gamma`, `theta1 = beta`, and `theta2 = gamma`. Independent Gaussian priors on the raw coefficients therefore induce a correlated prior on the generator's centered coefficients. The earlier PR #131 pilot correctly disclosed this confounding.

## Matched parameterization

The follow-up branch introduces an explicit `quadratic_center` in the inclusion-posterior representation and uses the basis

`X_1(x) = [1, x]`

`X_2(x) = [1, x, x^2 - c_g]`.

The value (c_g) is taken from the generator's parent-measure calibration output (`quadratic["config"]["x2_center"]`), not copied as a hard-coded numerical constant. The default Gaussian prior remains independent with scales ([5,2.5]) for degree one and ([5,2.5,2.5]) for degree two. Consequently the degree-one prior on ((alpha,beta)) is exactly the marginal prior on the first two centered degree-two coefficients; adding the independent gamma prior does not shift or correlate alpha and beta.

The same basis center is stored in the posterior object and used for:
- fitting the validation Bernoulli posterior;
- posterior predictive inclusion probabilities;
- node-specific IPW host weights;
- shared-node likelihood propagation;
- the independent child-artifact auditor's Newton/Gauss-Hermite reconstruction.

The auditor additionally checks that the fitted center equals the generator's declared (c_g).

## Regression tests

The tests verify:
1. ([1,x,x^2-c_g](alpha,beta,gamma)^T) reconstructs the generator linear predictor exactly;
2. non-finite basis centers are rejected;
3. degree-two five-arm MAP probabilities use the centered basis;
4. shared-node likelihood weights use the center stored in the posterior object;
5. the independent smoke auditor reconstructs the same centered parameterization.

## Scientific limits

This removes one known prior-parameterization confound; it does not establish that degree two is better. The model still requires paired numerical convergence checks, weight-stability assessment, larger replicate counts, and prior-predictive calibration. The degree-two posterior has an additional coefficient and can have greater finite-validation uncertainty. All work remains a synthetic low-redshift Euclidean methodology benchmark, not observational H0 inference.
