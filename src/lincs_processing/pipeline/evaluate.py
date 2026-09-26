"""Evaluate: baselines first, then the model, on the change from control.

Level-5 signatures are already expressed relative to vehicle control, so every
metric here scores predicted vs measured differential expression. The baselines
use only the training fold of the same split.
"""

import json

import anndata as ad
import numpy as np
import pandas as pd

from lincs_processing.pipeline.data import obs_frame
from lincs_processing.pipeline.train import fold_indices

DEG_K = 100


def _rowwise_pearson(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    yc = y - y.mean(1, keepdims=True)
    pc = p - p.mean(1, keepdims=True)
    denom = np.linalg.norm(yc, axis=1) * np.linalg.norm(pc, axis=1)
    return np.divide((yc * pc).sum(1), denom, out=np.zeros(len(y)), where=denom > 0)


def _rowwise_cosine(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(y, axis=1) * np.linalg.norm(p, axis=1)
    return np.divide((y * p).sum(1), denom, out=np.zeros(len(y)), where=denom > 0)


def deg_direction_precision(y: np.ndarray, p: np.ndarray, k: int = DEG_K) -> float:
    """Share of each signature's top-k |z| genes whose sign the prediction gets."""
    k = min(k, y.shape[1])
    top = np.argpartition(-np.abs(y), k - 1, axis=1)[:, :k]
    rows = np.arange(len(y))[:, None]
    return float((np.sign(y[rows, top]) == np.sign(p[rows, top])).mean())


def metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    y = y.astype(np.float64)
    p = p.astype(np.float64)
    pearson = _rowwise_pearson(y, p)
    sst = ((y - y.mean(0)) ** 2).sum()
    return {
        "pearson_mean": float(pearson.mean()),
        "pearson_median": float(np.median(pearson)),
        "cosine_mean": float(_rowwise_cosine(y, p).mean()),
        "r2": float(1 - ((y - p) ** 2).sum() / sst) if sst > 0 else 0.0,
        "rmse": float(np.sqrt(((y - p) ** 2).mean())),
        f"deg{DEG_K}_direction": deg_direction_precision(y, p),
    }


def _group_mean_prediction(
    train_keys: pd.Series, Y_train: np.ndarray, test_keys: pd.Series
) -> np.ndarray:
    """Per-group training mean; groups unseen in training get the global mean."""
    global_mean = Y_train.mean(0)
    means = {
        key: Y_train[rows].mean(0)
        for key, rows in train_keys.groupby(train_keys, observed=True).indices.items()
    }
    return np.stack([means.get(k, global_mean) for k in test_keys])


def baselines(obs: pd.DataFrame, Y: np.ndarray, idx: dict) -> dict[str, np.ndarray]:
    tr, te = idx["train"], idx["test"]
    obs_tr, obs_te = obs.iloc[tr], obs.iloc[te]
    return {
        "train_mean": np.broadcast_to(Y[tr].mean(0), (len(te), Y.shape[1])),
        # Informative for cold-compound (cells seen); collapses to train_mean
        # for cold-cell, where the test cells were never trained on.
        "cell_mean": _group_mean_prediction(
            obs_tr["cell_id"].reset_index(drop=True), Y[tr], obs_te["cell_id"]
        ),
        # The mirror image: informative for cold-cell only.
        "compound_mean": _group_mean_prediction(
            obs_tr["compound_key"].reset_index(drop=True),
            Y[tr],
            obs_te["compound_key"],
        ),
    }


def evaluate(
    adata: ad.AnnData, fold: pd.Series, predictions: dict, summary: dict
) -> dict:
    obs = obs_frame(adata)
    Y = np.asarray(adata.X, dtype=np.float32)
    idx = fold_indices(obs, fold)
    y_test = Y[idx["test"]]

    pred = pd.DataFrame(predictions["y_pred"], index=predictions["sig_id"])
    y_model = pred.loc[obs.index[idx["test"]]].to_numpy()

    results = {name: metrics(y_test, p) for name, p in baselines(obs, Y, idx).items()}
    results["model"] = metrics(y_test, y_model)
    return {
        "split_id": summary["split_id"],
        "split_type": summary["split_type"],
        "seed": summary["seed"],
        "run_id": summary.get("run_id"),
        "data_hash": summary["data_hash"],
        "split_hash": summary["split_hash"],
        "best_val_loss": summary["best_val_loss"],
        "n_test": int(len(idx["test"])),
        "methods": results,
        "beats_train_mean": results["model"]["pearson_mean"]
        > results["train_mean"]["pearson_mean"],
    }


def write_metrics(result: dict, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(result, fh, indent=2)
