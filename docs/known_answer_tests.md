# Known-answer validation: first inference baseline

## What was added

`src/dark_siren_robustness/inference.py` implements a deliberately simple
low-redshift Gaussian-mixture posterior for (H_0). It is intended for
validation exercises, not as a full gravitational-wave cosmology analysis.

The tests check:
- recovery of a known (H_0) value in a single-host case with a narrow distance likelihood;
- posterior normalization on the supplied grid;
- invariance to multiplying all host weights by the same positive constant;
- rejection of malformed grids, invalid redshifts, invalid uncertainty, and invalid weights;
- exclusion of grid points outside the specified uniform-prior bounds.

## Why start with one host?

For one known host redshift (z), the toy relation is

[
d_L approx \frac{c z}{H_0}.
]

If the observed distance is generated at a chosen (H_0) and the distance
uncertainty is small, the posterior mode should lie close to that chosen value,
provided it is inside the prior range. This is a useful known-answer check
before adding a population of candidate hosts.

## Limitations

Passing these unit tests would establish only that these specific code paths
behave as expected. It would not establish realistic catalogue modelling,
correct GW detection-selection treatment, full hierarchical inference,
simulation-based calibration, or publication readiness.

## Run locally

```powershell
python -m pip install -e ".[dev]"
python -m pytest
```

The tests have been committed to the repository, but they must still run in
GitHub Actions or a local environment before being described as passing.
