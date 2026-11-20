# Paired degree-one / degree-two pilot readout — 2026-11-20

## Execution and independent audit

Workflow `44683355394` completed successfully on the corrected PR #131 head. The paired orchestrator reported PASS. Both degree-specific child workflows passed:

- 90/90 lifecycle cells SUCCESS for degree one; 9/9 child runs passed.
- 90/90 lifecycle cells SUCCESS for degree two; 9/9 child runs passed.
- Outer independent audit PASS.
- 354 outer manifest entries and 574 nested manifest entries verified.
- Nine paired child blocks and 90 paired summary rows reconstructed.
- A/B/C common-arm likelihoods invariant within tolerance; maximum absolute difference (5.68\times10^{-14}).
- Artifact ZIP SHA-256: `bc8f5268b27576fd2b01d389f113a05a514ea74abbec4549255df128966bfc7`.

This establishes successful execution, matched inputs, artifact integrity, and aggregate arithmetic. It does not establish calibration.

## Frozen design

Three paired replicates; gamma (-1,0,+1); parent size 1,200; 30 detected events; validation frame 500; nested validation sizes 100 and 500; fixed (H_0=70), prior [30,110], 281-point grid, quadrature order 5. Both degrees use identical random seeds and synthetic inputs.

## Selected posterior-error summaries

RMSE is in km/s/Mpc. Each value is based on only three replicates.

| gamma | validation n | arm | Degree-1 RMSE | Degree-2 RMSE |
|---:|---:|---|---:|---:|
| -1 | 100 | D MAP plug-in IPW | 3.170 | 6.113 |
| -1 | 100 | E integrated IPW | 3.561 | 8.855 |
| -1 | 500 | D MAP plug-in IPW | 1.847 | 5.100 |
| -1 | 500 | E integrated IPW | 2.163 | 6.747 |
| 0 | 100 | D MAP plug-in IPW | 3.392 | 4.873 |
| 0 | 100 | E integrated IPW | 4.354 | 6.904 |
| 0 | 500 | D MAP plug-in IPW | 4.565 | 5.081 |
| 0 | 500 | E integrated IPW | 4.859 | 5.836 |
| +1 | 100 | D MAP plug-in IPW | 8.097 | 6.483 |
| +1 | 100 | E integrated IPW | 8.937 | 8.153 |
| +1 | 500 | D MAP plug-in IPW | 8.630 | 5.416 |
| +1 | 500 | E integrated IPW | 8.877 | 6.238 |

The degree-two procedure is better in some positive-curvature cases but worse in the negative- and zero-curvature examples. This is not evidence of a consistent improvement. Coverage estimates are especially unstable at (n_{rep}=3).

## Interpretation constraints

1. **Prior mismatch:** the degree-one and degree-two coefficient priors are not exactly matched under raw versus centered polynomial coordinates. The current comparison is between complete procedures, not a clean isolation of polynomial degree.
2. **Quadrature:** order 5 was used. The earlier single-realization 5/7/9 pilot is not a substitute for convergence checks in a larger paired experiment.
3. **Calibration:** three replicates cannot support meaningful coverage or SBC conclusions.
4. **Scope:** synthetic low-redshift Euclidean methodology only. No observational (H_0), LVK/GWTC, or real-survey claim.

## Next experiment

Before increasing the replicate count, implement a prior-matched centered quadratic basis ([1,x,x^2-c_g]), with the same declared priors on ((\alpha,\beta,\gamma)) in the degree-two model and corresponding ((\alpha,\beta)) priors in degree one. Add direct tests that the fitted basis matches the generator parameterization. Then run a modest paired pilot with quadrature orders 5/7/9 and inspect weight stability. Only after that should a larger prior-predictive calibration study be scheduled.
