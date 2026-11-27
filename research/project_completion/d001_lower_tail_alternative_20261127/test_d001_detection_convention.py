"""D-001 regression tests: detection-acceptance convention.

Self-contained ORACLE tests (run anywhere) + PRODUCTION adapter tests (skipped until
the adapters below are filled in against the repository's real functions; they have
NOT been run against production).

Convention under test: accepted event A = {d_obs < tau}, d_obs ~ N(mu, sigma)
(additive Gaussian in distance). P(A)=Phi((tau-mu)/sigma).
Opposite tail P(d_obs >= tau) = Phi((mu-tau)/sigma).

Run: python -m unittest test_d001_detection_convention -v (or pytest).
"""
import unittest

import numpy as np
from scipy.special import ndtr, log_ndtr

C = 299792.458


def p_accept_lt(mu, sigma, tau):
    return ndtr((tau - mu) / sigma)  # oracle: P(d < tau)


def p_accept_ge(mu, sigma, tau):
    return ndtr((mu - tau) / sigma)  # oracle: P(d >= tau)


def log_p_accept_lt(mu, sigma, tau):
    return log_ndtr((tau - mu) / sigma)


# Fill these in with read-only imports from the repository.
# Signatures are placeholders; adapters are deliberately not wired to production.
ADAPTERS = {
    "generator_accepts": None,         # callable(d_obs, tau) -> bool array
    "likelihood_log_selection": None,  # callable(means, sigma, tau) -> log P(accept)
    "selection_efficiency": None,      # callable(means, sigma, tau) -> P(accept)
    "dataset_norm": None,              # callable(means, sigma, tau) -> P(accept)
}


class TestOracleAlgebra(unittest.TestCase):
    def test_tails_are_complements(self):
        a = np.linspace(-8, 8, 161)
        np.testing.assert_allclose(ndtr(a) + ndtr(-a), 1.0, atol=1e-15)

    def test_deterministic_discriminating_point(self):
        # mu = tau - 2 sigma separates the two tails unambiguously.
        mu, sg, tau = 100.0, 10.0, 120.0
        self.assertAlmostEqual(
            p_accept_lt(mu, sg, tau), 0.9772498680518208, places=15
        )
        self.assertAlmostEqual(
            p_accept_ge(mu, sg, tau), 0.022750131948179195, places=15
        )

    def test_limits(self):
        self.assertAlmostEqual(p_accept_lt(50.0, 1e-9, 60.0), 1.0, places=12)
        self.assertAlmostEqual(p_accept_lt(70.0, 1e-9, 60.0), 0.0, places=12)
        # At mu=tau the tails coincide: this is NOT a discriminating test.
        self.assertAlmostEqual(p_accept_lt(60.0, 5.0, 60.0), 0.5, places=15)
        self.assertAlmostEqual(p_accept_lt(60.0, 5.0, 1e12), 1.0, places=12)
        self.assertAlmostEqual(p_accept_lt(60.0, 5.0, -1e12), 0.0, places=12)

    def test_monotone_in_H0(self):
        # mu=c*z/H0 decreases with H0, so P(d<tau) increases with H0.
        H = np.linspace(50, 90, 401)
        mu = C * 0.05 / H
        s = 0.1 * mu
        a = p_accept_lt(mu, s, 330.0)
        b = p_accept_ge(mu, s, 330.0)
        self.assertTrue(np.all(np.diff(a) >= 0))
        self.assertTrue(np.all(np.diff(b) <= 0))

        # A single host can saturate at exactly 1 or 0 in floating point.
        live = (a[:-1] > 1e-12) & (a[1:] < 1 - 1e-12)
        self.assertTrue(live.any())
        self.assertTrue(np.all(np.diff(a)[live] > 0))

        # Host-marginal beta(H): avoids single-host saturation.
        zg = (np.arange(400) + 0.5) * 0.1 / 400
        w = zg**2
        w /= w.sum()
        M = C * zg[None, :] / H[:, None]
        beta = (w * p_accept_lt(M, 0.1 * M, 330.0)).sum(1)
        self.assertTrue(np.all(np.diff(beta) > 0))
        self.assertTrue(np.all(np.diff(1 - beta) < 0))

    def test_monte_carlo_acceptance_matches_analytic(self):
        rng = np.random.default_rng(1)
        mu, sg, tau, n = 300.0, 40.0, 330.0, 400_000
        x = mu + sg * rng.standard_normal(n)
        p = p_accept_lt(mu, sg, tau)
        se = np.sqrt(p * (1 - p) / n)
        self.assertLess(abs((x < tau).mean() - p), 4 * se)
        self.assertGreater(
            abs((x < tau).mean() - p_accept_ge(mu, sg, tau)), 0.3
        )

    def test_host_marginal_mass_detects_wrong_tail(self):
        """Accepted density on d<tau integrates to 1 after division by beta.

        Dividing by P(d>=tau) instead integrates to beta/(1-beta).
        """
        zg = (np.arange(400) + 0.5) * 0.1 / 400
        w = zg**2
        w /= w.sum()
        mu = C * zg / 70.0
        s = 0.1 * mu
        tau = 330.0
        d = np.linspace(0, tau, 40001)
        N = (
            np.exp(-0.5 * ((d[:, None] - mu) / s) ** 2)
            / (s * np.sqrt(2 * np.pi))
        )
        num = np.trapezoid((N * w).sum(1), d)
        beta = (w * p_accept_lt(mu, s, tau)).sum()
        self.assertAlmostEqual(num / beta, 1.0, places=6)
        self.assertAlmostEqual(num / (1 - beta), beta / (1 - beta), places=6)
        self.assertGreater(abs(num / (1 - beta) - 1), 0.05)

    def test_log_space_stability_in_the_tail(self):
        v = log_ndtr(-40.0)
        self.assertTrue(np.isfinite(v))
        self.assertAlmostEqual(v, -804.6084420137538, places=6)
        # Naive log(Phi) underflows.
        self.assertEqual(np.log(ndtr(-40.0)), -np.inf)

    def test_selection_power_is_N_not_average(self):
        # N events => beta(H)^-N once, not a per-event average.
        b = np.array([0.5, 0.25])
        N = 7
        self.assertAlmostEqual(
            float(-N * np.log(b[1]) + N * np.log(b[0])),
            7 * np.log(2.0),
            places=12,
        )


def _need(name):
    if ADAPTERS[name] is None:
        raise unittest.SkipTest(
            f"adapter '{name}' not wired to production yet "
            "(NOT RUN against production)"
        )
    return ADAPTERS[name]


class TestProductionAgreement(unittest.TestCase):
    MU = np.array([100.0, 119.0, 120.0, 150.0])
    SG = np.array([10.0, 10.0, 10.0, 10.0])
    TAU = 120.0

    def test_generator_accepts_strictly_below_threshold(self):
        f = _need("generator_accepts")
        d = np.array([119.999, 120.0, 120.001])
        np.testing.assert_array_equal(
            np.asarray(f(d, 120.0), bool), [True, False, False]
        )

    def test_likelihood_normalisation_is_the_generators_event(self):
        f = _need("likelihood_log_selection")
        np.testing.assert_allclose(
            f(self.MU, self.SG, self.TAU),
            log_p_accept_lt(self.MU, self.SG, self.TAU),
            rtol=1e-12,
            atol=1e-12,
        )

    def test_every_selection_implementation_agrees(self):
        ref = p_accept_lt(self.MU, self.SG, self.TAU)
        for k in ("selection_efficiency", "dataset_norm"):
            np.testing.assert_allclose(
                _need(k)(self.MU, self.SG, self.TAU),
                ref,
                rtol=1e-12,
                atol=1e-12,
                err_msg=k,
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
