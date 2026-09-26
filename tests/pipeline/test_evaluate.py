import numpy as np
import pandas as pd

from lincs_processing.pipeline.evaluate import (
    baselines,
    deg_direction_precision,
    metrics,
)


def test_perfect_prediction_scores_perfectly():
    rng = np.random.default_rng(0)
    y = rng.normal(size=(20, 50))
    m = metrics(y, y)
    assert m["pearson_mean"] == 1.0
    assert m["r2"] == 1.0
    assert m["rmse"] == 0.0
    assert m["deg100_direction"] == 1.0


def test_sign_flip_is_anti_correlated():
    y = np.random.default_rng(1).normal(size=(5, 30))
    m = metrics(y, -y)
    np.testing.assert_allclose(m["pearson_mean"], -1.0)
    assert deg_direction_precision(y, -y, k=10) == 0.0


def test_baselines_use_training_fold_only():
    obs = pd.DataFrame(
        {"cell_id": ["a", "a", "b", "b"], "compound_key": ["X", "Y", "X", "Z"]}
    )
    Y = np.array([[1.0, 0.0], [3.0, 0.0], [5.0, 2.0], [100.0, 100.0]])
    idx = {"train": np.array([0, 1, 2]), "test": np.array([3])}
    out = baselines(obs, Y, idx)
    np.testing.assert_allclose(out["train_mean"][0], [3.0, 2 / 3])
    np.testing.assert_allclose(out["cell_mean"][0], [5.0, 2.0])  # cell b in train
    np.testing.assert_allclose(out["compound_mean"][0], [3.0, 2 / 3])  # Z unseen
