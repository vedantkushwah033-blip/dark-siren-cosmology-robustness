from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_h0_known_inclusion_oracle_sensitivity import selected_posterior_quantiles


def inputs():
    z = np.array([0.04, 0.08, 0.12, 0.18, 0.25])
    weights = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    distances = np.array([450.0, 620.0])
    h0 = np.linspace(30.0, 150.0, 241)
    return z, weights, distances, h0


def test_oracle_posterior_quantiles_are_ordered_and_finite():
    z, weights, distances, h0 = inputs()
    q = selected_posterior_quantiles(z, weights, distances, h0, 100.0, 300.0)
    values = [q[k] for k in ("q025", "q050", "q160", "q500", "q840", "q950", "q975")]
    assert np.isfinite(values).all()
    assert np.all(np.diff(values) > 0)
    assert h0[0] <= values[0] < values[-1] <= h0[-1]


def test_oracle_posterior_is_invariant_to_common_weight_scaling():
    z, weights, distances, h0 = inputs()
    a = selected_posterior_quantiles(z, weights, distances, h0, 100.0, 300.0)
    b = selected_posterior_quantiles(z, 37.0 * weights, distances, h0, 100.0, 300.0)
    assert a == pytest.approx(b, abs=1e-12)


def test_oracle_rejects_events_outside_the_declared_detection_region():
    z, weights, distances, h0 = inputs()
    with pytest.raises(ValueError, match="detected-event"):
        selected_posterior_quantiles(z, weights, np.array([250.0]), h0, 100.0, 300.0)
