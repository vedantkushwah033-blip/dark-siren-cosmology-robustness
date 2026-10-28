# Exploratory host-density-to-H0 benchmark

**Status:** prototype; synthetic only, not a real Hubble-constant measurement.

The runner compares three host-redshift density arms using paired simulation seeds and the existing fixed-count selection-conditioned toy likelihood:

1. **Oracle:** KDE from parent-population truth; reference only, not an observable-data analysis.
2. **Catalogue-only:** KDE from the observed incomplete catalogue without correction.
3. **Known-completeness:** KDE from observed catalogue redshifts with inverse-inclusion weights evaluated from the declared selection law.

The same detected-event distances, merger-rate weighting, KDE bandwidth, H0 grid, Gaussian distance-noise scale, and detection threshold are used within each replicate. Output includes replicate-level posterior mean, median, equal-tail 68% interval, interval width, coverage indicator, and minimum detection efficiency, plus arm-level bias, RMSE, empirical coverage, and interval width.

Limitations: no realistic detector response, peculiar velocities, event-rate inference, calibration uncertainty, or real catalogue systematics. The known-completeness arm assumes the selection function is exactly known. Small-run coverage is noisy, and bandwidth/grid sensitivity remain mandatory. Independently audit replicate outputs and compare likelihood quadrature before interpreting results.
