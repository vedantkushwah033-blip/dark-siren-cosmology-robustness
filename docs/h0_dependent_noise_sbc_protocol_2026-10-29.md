# H0-dependent-noise selection model: simulation-based calibration protocol

Date: 2026-10-29. Status: preregistered software validation plan; not a real-data analysis.

## Model

For discrete host redshift (z) and Hubble constant (H_0), define
[
mu(z,H_0)=c z/H_0,qquad
sigma(z,H_0)=sqrt{sigma_{\rm floor}^2+[f_\sigma\mu(z,H_0)]^2}.
]
The event distance likelihood is a discrete host mixture. Detection is (d\ge d_{\min}>0), and its efficiency is the population-weighted sum of Normal survival probabilities using the same (H_0)-dependent noise scale as the event likelihood. The fixed-size selected likelihood is
[
\log L_{\rm sel}(H_0)=\sum_i\log L_i(H_0)-N_{\rm det}\log\beta(H_0).
]
The model conditions on the number of detected events; it does not infer event rates. It is a distance-only synthetic model, not a realistic GW detector likelihood.

## Locked baseline

- Prior: uniform (H_0\in[50,90]) km s(^{-1}) Mpc(^{-1}), H0 grid 161 evenly spaced points including both endpoints.
- Population: 12 fixed redshift nodes evenly spaced on [0.025, 0.16], positive linearly increasing raw discrete host weights.
- Distance noise: fractional scale 0.12, floor 8 Mpc; detection threshold 300 Mpc.
- Three detected events per simulated dataset.
- Four conditions, 300 replicates each: truth drawn from the prior and fixed (H_0=60,70,80).
- Seeds are non-overlapping by condition: start seed plus replicate index, with offsets 0, 100000, 200000, and 300000.
- Each raw record saves the seed, true H0, and exact detected-distance vector. The runner writes SHA-256 hashes for raw records and summary.
- The independent audit replays the RNG from each seed, independently generates detected distances, recomputes the Normal-mixture event likelihood and H0-dependent beta, rebuilds all posterior summaries, and compares them with the raw table.

## Primary diagnostics

For the prior-drawn mode, report the mean posterior-CDF rank, its Monte Carlo standard error, a 10-bin descriptive chi-square statistic, and nominal 68% and 90% coverage with Wilson intervals. The chi-square p-value is descriptive and is not a stand-alone proof of calibration. For each fixed-truth mode, report conditional bias, bias SE, RMSE, coverage with Wilson intervals, posterior interval width, and CDF-rank behavior. All modes use the same selected-likelihood normalization; the unit tests also include an explicit negative control that removes this normalization.

## Interpretation and stopping rule

A pass means that software, generator, selected likelihood, and the declared synthetic model are internally consistent to the test tolerances. It does not validate realistic detector selection, galaxy-catalogue incompleteness, event-rate inference, uncertain redshifts, cosmological distance-redshift relations, or real H0 measurements. Any adverse calibration or weak coverage is a reported result, not grounds to change the locked model or rerun until a favorable outcome appears.
