import json

import numpy as np
import pandas as pd

from lincs_processing.pipeline.evaluate import evaluate
from lincs_processing.pipeline.model import Bundle, Preprocessor
from lincs_processing.pipeline.splits import make_split, split_record
from lincs_processing.pipeline.tracking import make_tracker
from lincs_processing.pipeline.train import fold_indices, train


def test_preprocessor_sees_only_the_training_fold(qc_adata):
    obs, Y = qc_adata.obs, np.asarray(qc_adata.X)
    fold = make_split(obs, "cold_cell", seed=0)
    tr = fold_indices(obs, fold)["train"]

    pre = Preprocessor("cell_center").fit(obs.iloc[tr], Y[tr])
    Y_poisoned = Y.copy()
    Y_poisoned[np.setdiff1d(np.arange(len(Y)), tr)] += 100.0  # val/test only
    pre2 = Preprocessor("cell_center").fit(obs.iloc[tr], Y_poisoned[tr])

    np.testing.assert_array_equal(pre.global_mean, pre2.global_mean)
    np.testing.assert_array_equal(pre.y_std, pre2.y_std)
    assert set(pre.cells) == set(obs.iloc[tr]["cell_id"])
    np.testing.assert_allclose(pre.global_mean, Y[tr].mean(0), rtol=1e-5)


def test_unseen_cell_falls_back_to_global_mean(qc_adata):
    obs, Y = qc_adata.obs, np.asarray(qc_adata.X)
    pre = Preprocessor("cell_center").fit(obs, Y)
    unseen = pd.DataFrame({"cell_id": ["NEVER_SEEN"]})
    np.testing.assert_array_equal(pre._offsets(unseen)[0], pre.global_mean)
    assert pre.cell_matrix(unseen).sum() == 0


def test_train_evaluate_round_trip(qc_adata, tmp_path):
    obs = qc_adata.obs
    fold = make_split(obs, "cold_compound", seed=0)
    split = split_record(fold, "cold_compound", 0, qc_adata.uns["data_hash"], 0.1, 0.2)
    summary = train(
        qc_adata,
        fold,
        split,
        seed=0,
        outdir=str(tmp_path),
        tracker=make_tracker("none"),
        epochs=3,
        batch_size=64,
        device="cpu",
    )
    assert summary["n_test"] == (fold == "test").sum()
    assert json.loads((tmp_path / "train_summary_seed0.json").read_text())["seed"] == 0

    preds = dict(np.load(tmp_path / "predictions_seed0.npz"))
    result = evaluate(qc_adata, fold, preds, summary)
    assert set(result["methods"]) == {
        "model",
        "train_mean",
        "cell_mean",
        "compound_mean",
    }

    bundle = Bundle(str(tmp_path / "model_seed0"))
    test_obs = obs[fold == "test"]
    np.testing.assert_allclose(bundle.predict(test_obs), preds["y_pred"], atol=1e-5)
