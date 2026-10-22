# Selection-mechanism interpretation audit

## Scope

The factorial experiment uses the selection functions in the pinned original simulator at commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. Its "redshift-dependent" mask is a deliberately simple synthetic mechanism, not a calibrated survey selection function.

## What the code does

- The random mask retains exactly `round(N × completeness)` galaxies, selected uniformly without replacement.
- The redshift-dependent mask retains the same exact number of galaxies, but samples without replacement with selection probabilities proportional to `1 / (z + 1e-9)`. At the same catalogue size, this favors lower-redshift galaxies.
- Thus, the random-versus-redshift comparison isolates the effect of this particular redshift-dependent selection rule while holding the retained count fixed.

## Interpretation limits

1. This mask is not derived from an instrument sensitivity curve, flux limit, galaxy-survey completeness map, or measured selection function.
2. The model does not claim that real catalogues universally become incomplete according to `1/z`; actual selection depends on survey depth, galaxy properties, sky position, redshift measurement, and catalogue construction.
3. Exact fixed-size masks condition on a chosen retained count. They do not model uncertainty in the total number of detected galaxies.
4. The simulation uses the low-redshift approximation `d_L ≈ c z / H0` for both generating distances and evaluating the likelihood. Conclusions therefore apply only within this simplified, internally matched model.
5. The luminosity proxy is sampled synthetically and used as a host-generation weight in selected scenarios. It is not a validated astrophysical host-rate model.

## Reporting recommendation

Use wording such as **"inverse-redshift-weighted synthetic missingness"** instead of implying a realistic or universal redshift selection function. Report the mask definition alongside every result. Treat conclusions as conditional sensitivity results under this controlled design, not direct predictions for LVK analyses or real galaxy surveys.

## Reproducibility

The exact definitions are in the pinned source files `src/completeness.py`, `src/cosmology.py`, `src/simulations_v2.py`, and `src/experiments.py` of the original repository. The current experiment runner records the pinned commit in its result metadata.
