# Source-level model-consistency audit

## Scope

This note records a source inspection of the original repository at commit
`b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. It is not a substitute for
running the pinned simulation and inspecting its saved artifacts.

## Findings

### 1. Distance relation is internally consistent in the baseline

The pinned `src/cosmology.py` defines the simulated distance as
`d_L = c z / H0`. The inference likelihood in `src/bayesian_inference.py`
uses the same relation, `d = c z / H0_grid`. Therefore, a mismatch between
the distance-redshift relation used by the generator and the likelihood is
not evident in this baseline implementation. This is a deliberately
simplified low-redshift model, not a full luminosity-distance calculation in
a specified cosmology.

### 2. The event-distance noise model is Gaussian

The pinned `simulate_gw_event` draws the observed distance from a normal
distribution centered on the simulated true distance, with standard
deviation equal to 15% of that distance by default. The likelihood uses a
Gaussian density with the supplied observed distance and sigma. The generator
and likelihood therefore appear aligned on the conditional Gaussian form,
although the model remains an idealized distance-only experiment.

### 3. The unknown-host inference is a catalogue mixture

The complete-catalogue experiment generates a host uniformly from the
simulated catalogue, then evaluates a likelihood mixture over the catalogue
redshifts with uniform inference weights. This is the intended matched
control for the host-unknown case. The known-host diagnostic instead passes
only the selected host's redshift to the likelihood.

### 4. What the coverage difference can and cannot establish

The source structure makes catalogue marginalization a plausible contributor
to wider unknown-host posteriors. But the source inspection alone does not
establish why empirical coverage differs between the controls. Coverage,
interval width, bias, RMSE, and posterior-pull diagnostics must be evaluated
from the actual replicated outputs.

In particular, high coverage of nominal 68% intervals is evidence of
conservative intervals relative to the chosen repeated-sampling experiment;
it is not, by itself, proof that the inference is more accurate or that the
model describes real gravitational-wave observations.

## Remaining checks after the workflow run

1. Confirm the source commit and dependency versions in the run record.
2. Compare the reproduced summary with the values reported in the original
   calibration document.
3. Inspect coverage at nominal 50%, 68%, 80%, 90%, and 95% for both controls.
4. Compare mean interval widths and estimator bias/RMSE across controls.
5. Repeat with another seed before making a strong claim about calibration.
6. Treat the model as a synthetic low-redshift validation study unless and
   until realistic selection effects, cosmology, and observational data are
   incorporated.
