"""Model and the training-fold-only preprocessing that goes with it.

``Preprocessor.fit`` is only ever called on the training fold inside the
training task, so val/test signatures cannot influence scaling, the cell
vocabulary or the batch correction.
"""

import json
import os
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import torch
from torch import nn

from lincs_processing.pipeline.features import condition_features, morgan_fingerprints

BATCH_CORRECTIONS = ("none", "cell_center")


@dataclass
class ModelConfig:
    hidden: int = 1024
    layers: int = 2
    dropout: float = 0.2


class Preprocessor:
    """Input/target scaling, cell encoding and per-cell centring."""

    def __init__(self, batch_correction: str = "cell_center") -> None:
        if batch_correction not in BATCH_CORRECTIONS:
            raise ValueError(f"batch_correction must be one of {BATCH_CORRECTIONS}")
        self.batch_correction = batch_correction
        self.cells: list[str] = []
        self.cell_table: pd.DataFrame | None = None
        self.cond_mean = np.zeros(2, dtype=np.float32)
        self.cond_std = np.ones(2, dtype=np.float32)
        self.cell_mean: dict[str, np.ndarray] = {}
        self.global_mean = np.zeros(0, dtype=np.float32)
        self.y_std = np.ones(0, dtype=np.float32)

    def fit(
        self,
        obs: pd.DataFrame,
        Y: np.ndarray,
        cell_features: pd.DataFrame | None = None,
    ) -> "Preprocessor":
        cond = condition_features(obs)
        self.cond_mean = cond.mean(0).astype(np.float32)
        self.cond_std = (cond.std(0) + 1e-6).astype(np.float32)
        self.cells = sorted(obs["cell_id"].unique())
        if cell_features is not None:
            table = cell_features.astype(np.float32)
            train_rows = table.loc[table.index.intersection(self.cells)]
            mean = train_rows.mean(axis=0)
            std = train_rows.std(axis=0).fillna(0) + 1e-6
            self.cell_table = (table - mean) / std

        self.global_mean = Y.mean(0).astype(np.float32)
        if self.batch_correction == "cell_center":
            for cell, rows in obs.groupby("cell_id", observed=True).indices.items():
                self.cell_mean[str(cell)] = Y[rows].mean(0).astype(np.float32)
        self.y_std = (self._center(obs, Y).std(0) + 1e-6).astype(np.float32)
        return self

    # Unseen cells (cold-cell) fall back to the global training mean.
    def _offsets(self, obs: pd.DataFrame) -> np.ndarray:
        if self.batch_correction == "none":
            return np.broadcast_to(self.global_mean, (len(obs), self.global_mean.size))
        return np.stack(
            [self.cell_mean.get(c, self.global_mean) for c in obs["cell_id"]]
        )

    def _center(self, obs: pd.DataFrame, Y: np.ndarray) -> np.ndarray:
        return Y - self._offsets(obs)

    def cell_matrix(self, obs: pd.DataFrame) -> np.ndarray:
        if self.cell_table is not None:
            table = self.cell_table.reindex(obs["cell_id"]).fillna(0.0)
            return table.to_numpy(dtype=np.float32)
        lookup = {c: i for i, c in enumerate(self.cells)}
        out = np.zeros((len(obs), len(self.cells)), dtype=np.float32)
        for row, cell in enumerate(obs["cell_id"]):
            if cell in lookup:
                out[row, lookup[cell]] = 1.0
        return out

    def inputs(self, obs: pd.DataFrame) -> np.ndarray:
        cond = (condition_features(obs) - self.cond_mean) / self.cond_std
        fps = morgan_fingerprints(obs["smiles"])
        return np.concatenate([fps, self.cell_matrix(obs), cond], axis=1).astype(
            np.float32
        )

    def targets(self, obs: pd.DataFrame, Y: np.ndarray) -> np.ndarray:
        return (self._center(obs, Y) / self.y_std).astype(np.float32)

    def inverse_targets(self, obs: pd.DataFrame, Z: np.ndarray) -> np.ndarray:
        return (Z * self.y_std + self._offsets(obs)).astype(np.float32)

    def save(self, directory: str) -> None:
        arrays = {
            "cond_mean": self.cond_mean,
            "cond_std": self.cond_std,
            "global_mean": self.global_mean,
            "y_std": self.y_std,
            "cell_mean_keys": np.array(list(self.cell_mean), dtype=str),
            "cell_mean_values": np.stack(list(self.cell_mean.values()))
            if self.cell_mean
            else np.zeros((0, self.global_mean.size), dtype=np.float32),
        }
        np.savez(os.path.join(directory, "preprocessor.npz"), **arrays)  # type: ignore[arg-type]
        meta = {"batch_correction": self.batch_correction, "cells": self.cells}
        with open(os.path.join(directory, "preprocessor.json"), "w") as fh:
            json.dump(meta, fh, indent=2)
        if self.cell_table is not None:
            self.cell_table.to_parquet(os.path.join(directory, "cell_features.parquet"))

    @classmethod
    def load(cls, directory: str) -> "Preprocessor":
        with open(os.path.join(directory, "preprocessor.json")) as fh:
            meta = json.load(fh)
        pre = cls(meta["batch_correction"])
        pre.cells = meta["cells"]
        arrays = np.load(os.path.join(directory, "preprocessor.npz"))
        pre.cond_mean, pre.cond_std = arrays["cond_mean"], arrays["cond_std"]
        pre.global_mean, pre.y_std = arrays["global_mean"], arrays["y_std"]
        pre.cell_mean = dict(
            zip(arrays["cell_mean_keys"], arrays["cell_mean_values"], strict=True)
        )
        table = os.path.join(directory, "cell_features.parquet")
        if os.path.exists(table):
            pre.cell_table = pd.read_parquet(table)
        return pre


class PerturbationMLP(nn.Module):
    def __init__(self, n_in: int, n_out: int, cfg: ModelConfig) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        width = n_in
        for _ in range(cfg.layers):
            layers += [
                nn.Linear(width, cfg.hidden),
                nn.LayerNorm(cfg.hidden),
                nn.GELU(),
                nn.Dropout(cfg.dropout),
            ]
            width = cfg.hidden
        layers.append(nn.Linear(width, n_out))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def save_bundle(
    directory: str,
    model: PerturbationMLP,
    pre: Preprocessor,
    cfg: ModelConfig,
    n_in: int,
    genes: list[str],
    extra: dict | None = None,
) -> None:
    """A self-contained model directory: weights, preprocessing, genes, config."""
    os.makedirs(directory, exist_ok=True)
    torch.save(model.state_dict(), os.path.join(directory, "model.pt"))
    pre.save(directory)
    config = {"model": asdict(cfg), "n_in": n_in, "genes": genes, **(extra or {})}
    with open(os.path.join(directory, "config.json"), "w") as fh:
        json.dump(config, fh, indent=2)


class Bundle:
    """A loaded model directory that predicts signatures from metadata."""

    def __init__(self, directory: str, device: str = "cpu") -> None:
        with open(os.path.join(directory, "config.json")) as fh:
            self.config = json.load(fh)
        self.genes: list[str] = self.config["genes"]
        self.pre = Preprocessor.load(directory)
        self.model = PerturbationMLP(
            self.config["n_in"], len(self.genes), ModelConfig(**self.config["model"])
        )
        state = torch.load(
            os.path.join(directory, "model.pt"), map_location=device, weights_only=True
        )
        self.model.load_state_dict(state)
        self.model.to(device).eval()
        self.device = device

    @torch.no_grad()
    def predict(self, obs: pd.DataFrame, batch_size: int = 4096) -> np.ndarray:
        X = self.pre.inputs(obs)
        out = [
            self.model(torch.from_numpy(X[i : i + batch_size]).to(self.device))
            .cpu()
            .numpy()
            for i in range(0, len(X), batch_size)
        ]
        return self.pre.inverse_targets(obs, np.concatenate(out, axis=0))
