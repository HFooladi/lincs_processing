"""Fit & train: one task per (split, seed).

Scaling, the cell vocabulary and batch correction are fitted here, after the
split, on the training fold only. Early stopping watches the validation fold;
the test fold is only touched to write predictions for Evaluate.
"""

import copy
import json
import os
import random
from dataclasses import asdict

import anndata as ad
import numpy as np
import pandas as pd
import torch
from torch import nn

from lincs_processing.pipeline.data import obs_frame
from lincs_processing.pipeline.model import (
    ModelConfig,
    PerturbationMLP,
    Preprocessor,
    save_bundle,
)
from lincs_processing.pipeline.tracking import Tracker


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def resolve_device(device: str) -> str:
    if device == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return device


def fold_indices(obs: pd.DataFrame, fold: pd.Series) -> dict[str, np.ndarray]:
    labels = fold.reindex(obs.index).fillna("excluded").to_numpy()
    return {f: np.flatnonzero(labels == f) for f in ("train", "val", "test")}


@torch.no_grad()
def _predict_scaled(
    model: nn.Module, X: np.ndarray, device: str, batch_size: int
) -> np.ndarray:
    model.eval()
    out = [
        model(torch.from_numpy(X[i : i + batch_size]).to(device)).cpu().numpy()
        for i in range(0, len(X), batch_size)
    ]
    return np.concatenate(out, axis=0)


def train(
    adata: ad.AnnData,
    fold: pd.Series,
    split: dict,
    seed: int,
    outdir: str,
    tracker: Tracker,
    epochs: int = 50,
    batch_size: int = 512,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    patience: int = 5,
    model_cfg: ModelConfig | None = None,
    batch_correction: str = "cell_center",
    cell_features: pd.DataFrame | None = None,
    device: str = "auto",
) -> dict:
    seed_everything(seed)
    device = resolve_device(device)
    model_cfg = model_cfg or ModelConfig()
    obs = obs_frame(adata)
    Y = np.asarray(adata.X, dtype=np.float32)
    idx = fold_indices(obs, fold)
    tr, va, te = idx["train"], idx["val"], idx["test"]

    pre = Preprocessor(batch_correction).fit(obs.iloc[tr], Y[tr], cell_features)
    X_tr, Z_tr = pre.inputs(obs.iloc[tr]), pre.targets(obs.iloc[tr], Y[tr])
    X_va, Z_va = pre.inputs(obs.iloc[va]), pre.targets(obs.iloc[va], Y[va])

    model = PerturbationMLP(X_tr.shape[1], Y.shape[1], model_cfg).to(device)
    optim = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    loss_fn = nn.MSELoss()
    X_tr_t, Z_tr_t = torch.from_numpy(X_tr), torch.from_numpy(Z_tr)
    gen = torch.Generator().manual_seed(seed)

    params = {
        "split_id": split["split_id"],
        "split_type": split["split_type"],
        "seed": seed,
        "epochs": epochs,
        "batch_size": batch_size,
        "lr": lr,
        "weight_decay": weight_decay,
        "patience": patience,
        "batch_correction": batch_correction,
        "cell_features": cell_features is not None,
        **asdict(model_cfg),
    }
    tracker.log_params(params)
    tracker.set_tags(
        {
            "data_hash": split["data_hash"],
            "split_hash": split["split_hash"],
            "split_id": split["split_id"],
            "device": device,
            "container": os.environ.get("LINCS_CONTAINER", "none"),
        }
    )

    best_loss, best_epoch, best_state, stale = float("inf"), -1, None, 0
    epoch = 0
    for epoch in range(epochs):
        model.train()
        order = torch.randperm(len(X_tr_t), generator=gen)
        total = 0.0
        for start in range(0, len(order), batch_size):
            batch = order[start : start + batch_size]
            xb = X_tr_t[batch].to(device, non_blocking=True)
            zb = Z_tr_t[batch].to(device, non_blocking=True)
            optim.zero_grad(set_to_none=True)
            loss = loss_fn(model(xb), zb)
            loss.backward()
            optim.step()
            total += loss.item() * len(batch)
        train_loss = total / len(order)
        val_loss = float(
            np.mean((_predict_scaled(model, X_va, device, 4096) - Z_va) ** 2)
        )
        tracker.log_metrics({"train_loss": train_loss, "val_loss": val_loss}, epoch)
        print(f"epoch {epoch:3d} train {train_loss:.4f} val {val_loss:.4f}")

        if val_loss < best_loss - 1e-6:
            best_loss, best_epoch, stale = val_loss, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1
            if stale >= patience:
                break

    assert best_state is not None
    model.load_state_dict(best_state)

    model_dir = os.path.join(outdir, f"model_seed{seed}")
    extra = {
        "split_id": split["split_id"],
        "split_type": split["split_type"],
        "split_hash": split["split_hash"],
        "data_hash": split["data_hash"],
        "seed": seed,
        "train_params": params,
    }
    genes = list(adata.var["gene_symbol"]) if "gene_symbol" in adata.var else []
    save_bundle(
        model_dir,
        model.cpu(),
        pre,
        model_cfg,
        X_tr.shape[1],
        genes or list(adata.var_names),
        extra,
    )

    obs_te = obs.iloc[te]
    Z_te = _predict_scaled(model, pre.inputs(obs_te), "cpu", 4096)
    np.savez_compressed(
        os.path.join(outdir, f"predictions_seed{seed}.npz"),
        sig_id=obs_te.index.to_numpy(dtype=str),
        y_pred=pre.inverse_targets(obs_te, Z_te),
    )

    summary = {
        **extra,
        "run_id": tracker.run_id,
        "device": device,
        "best_epoch": best_epoch,
        "best_val_loss": best_loss,
        "epochs_run": epoch + 1,
        "n_train": int(len(tr)),
        "n_val": int(len(va)),
        "n_test": int(len(te)),
    }
    tracker.log_metrics({"best_val_loss": best_loss})
    tracker.log_artifacts(model_dir, "model")
    with open(os.path.join(outdir, f"train_summary_seed{seed}.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary
