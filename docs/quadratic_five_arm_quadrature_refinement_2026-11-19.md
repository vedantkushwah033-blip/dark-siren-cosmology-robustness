# Post-artifact inclusion-posterior quadrature refinement

**Status:** exploratory numerical sensitivity check; not a preregistered repeated-sampling result.

## Input and method

This check uses the successful [quadratic five-arm 30-event smoke workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44447663357), artifact `14493993418` (ZIP SHA-256 `7211a4fac5f49a49f8d045916f57e8e327a904f34ab0361dbdbb07746f579851`). It uses the three archived curvature conditions `gamma = -1, 0, +1`, the largest nested validation sample (`n=500`), and the same 30 detected events in each paired curvature condition.

The degree-one logistic inclusion model, Gaussian coefficient prior (mean `[0,0]`, standard deviations `[5,2.5]`), coordinate transform, reflected weighted KDE, selection-conditioned likelihood, H0 prior `[50,90]`, and 281-point H0 grid were held fixed. Only tensor Gauss-Hermite order varied: 3, 5, 7, 9, and 11 (9, 25, 49, 81, and 121 coefficient nodes).

The coefficient posterior mode and Laplace covariance were recomputed with a separate damped-Newton implementation. Gauss-Hermite nodes and target-corrected importance weights were reconstructed independently. For each order, the complete Arm E likelihood was integrated over one shared coefficient node across all 30 events; nodes were not marginalized event-by-event.

## Reconstructed Arm E posterior quantiles

H0 values are in km s^-1 Mpc^-1.

| Curvature | GH order | q16 | Median | q84 | Quadrature ESS |
|---:|---:|---:|---:|---:|---:|
| -1 | 3 | 60.944420 | 64.740598 | 68.496344 | 4.010 |
| -1 | 5 | 60.941375 | 64.739901 | 68.499256 | 6.844 |
| -1 | 7 | 60.941443 | 64.739928 | 68.499332 | 9.453 |
| -1 | 9 | 60.941448 | 64.739932 | 68.499337 | 12.019 |
| -1 | 11 | 60.941448 | 64.739932 | 68.499337 | 14.579 |
| 0 | 3 | 68.293866 | 72.208431 | 75.844962 | 4.010 |
| 0 | 5 | 68.288330 | 72.208086 | 75.847206 | 6.845 |
| 0 | 7 | 68.288209 | 72.207995 | 75.847188 | 9.455 |
| 0 | 9 | 68.288212 | 72.207993 | 75.847187 | 12.021 |
| 0 | 11 | 68.288212 | 72.207993 | 75.847187 | 14.581 |
| +1 | 3 | 70.528020 | 74.081542 | 77.523239 | 4.010 |
| +1 | 5 | 70.516895 | 74.074434 | 77.519907 | 6.846 |
| +1 | 7 | 70.516780 | 74.074198 | 77.519708 | 9.456 |
| +1 | 9 | 70.516794 | 74.074200 | 77.519706 | 12.022 |
| +1 | 11 | 70.516795 | 74.074201 | 77.519706 | 14.583 |

## Numerical interpretation

Across these three curvature conditions, the largest absolute order-5-to-11 change was q16 **0.000118**, median **0.000233**, and q84 **0.000201**. The H0 grid spacing is approximately 0.142857, so these changes are far below one grid interval for this particular artifact. Orders 9 and 11 agree to the displayed precision.

Order 3 is less stable: its largest order-3-to-11 median difference is about 0.00734. The order-5 quadrature ESS is about 6.84 of 25 nodes; at order 11 it is about 14.58 of 121 nodes. This is an importance-weight concentration diagnostic, not a count of independent posterior samples.

## Limits

This is a numerical refinement of one fixed 30-event realization, at validation size 500, for three curvature values. It does **not** establish quadrature convergence across other seeds, validation sizes, stronger misspecification, or failure cases. It does not establish repeated-sampling calibration, estimator superiority, or observational validity. The five-arm experiment remains an engineering smoke; a repeated-seed study with explicit quadrature and failure accounting is still required.

The branch's independent auditor update reconstructs the selected-order logistic posterior mode, Laplace covariance, Gauss-Hermite nodes, and target-corrected quadrature weights. CI validation of that auditor update is tracked separately from this offline order-refinement calculation.
