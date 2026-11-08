import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.validation_design import split_disjoint_validation


def _parent_frame(n=500):
    rng = np.random.default_rng(4)
    z = rng.uniform(0.01, 0.40, n)
    p = 0.1 + 0.8 / (1.0 + np.exp((z - 0.22) / 0.05))
    retained = rng.binomial(1, p)
    return pd.DataFrame({
        "object_id": np.arange(n),
        "redshift": z,
        "retained": retained,
        "inclusion_probability": p,
    })


def test_validation_is_parent_srs_and_disjoint_from_analysis_catalogue():
    parent = _parent_frame()
    split = split_disjoint_validation(parent, n_validation=100, seed=19)
    assert len(split.validation) == 100
    assert np.intersect1d(split.validation_object_ids, split.analysis_object_ids).size == 0
    assert set(split.analysis_catalogue.columns) == {"object_id", "redshift"}
    assert set(split.validation.columns) == {"object_id", "redshift", "retained"}
    assert len(split.validation_object_ids) == len(np.unique(split.validation_object_ids))
    expected = parent.loc[
        (parent["retained"] == 1)
        & (~parent["object_id"].isin(split.validation_object_ids)),
        "object_id",
    ].to_numpy()
    assert np.array_equal(np.sort(split.analysis_object_ids), np.sort(expected))


def test_split_is_reproducible_and_ignores_hidden_probability_column():
    parent = _parent_frame()
    a = split_disjoint_validation(parent, n_validation=80, seed=8)
    b = split_disjoint_validation(parent, n_validation=80, seed=8)
    assert np.array_equal(a.validation_object_ids, b.validation_object_ids)
    assert "inclusion_probability" not in a.validation.columns


def test_validation_size_and_input_contracts_are_checked():
    parent = _parent_frame()
    with pytest.raises(ValueError):
        split_disjoint_validation(parent, n_validation=4, seed=1)
    with pytest.raises(ValueError):
        split_disjoint_validation(parent.drop(columns=["retained"]), n_validation=20, seed=1)
    duplicate = pd.concat([parent, parent.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        split_disjoint_validation(duplicate, n_validation=20, seed=1)
