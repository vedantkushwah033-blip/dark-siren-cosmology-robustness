# Statistical model notes

This document will be expanded after auditing the existing project. It deliberately does not yet specify a complete dark-siren likelihood.

## Required components

1. **Cosmology:** mapping from redshift and cosmological parameters to luminosity distance.
2. **Event data model:** the likelihood or posterior representation for each gravitational-wave event.
3. **Host population:** a model for the distribution of potential host galaxies and any host weights.
4. **Catalogue inclusion:** a probability model for whether a galaxy appears in the catalogue.
5. **Detection selection:** a separate model for the probability that a gravitational-wave event is detected.
6. **Hierarchical inference:** a normalized likelihood for cosmological and population parameters, including the relevant selection normalization.
7. **Calibration design:** the repeated-sampling process used to evaluate bias, coverage, and posterior calibration.

## Important distinction

Galaxy-catalogue inclusion and gravitational-wave detection are different selection processes. They must be represented separately unless a carefully derived joint model explicitly connects them.

## Pending audit

Before specifying the full likelihood, inspect the existing project and document its generative process, likelihood, priors, host weights, selection treatment, and interval construction. Do not fill these gaps with guesses.
