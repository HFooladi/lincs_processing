"""QC & filter: rule-based filtering of signatures. Nothing is fitted here.

Every rule is a fixed threshold on per-signature metadata, so QC can run before
the split without leaking anything from a future test fold.
"""

import json
from dataclasses import asdict, dataclass

import anndata as ad
import numpy as np
import pandas as pd

from lincs_processing.pipeline.data import obs_frame
from lincs_processing.pipeline.features import compound_key, valid_smiles
from lincs_processing.pipeline.hashing import data_hash


@dataclass
class QCConfig:
    min_cc_q75: float = 0.2  # replicate agreement
    min_tas: float = 0.0  # signature strength, comparable across releases
    min_ss: float | None = None  # release-specific scale; off by default
    min_rep: int = 2
    min_dose_um: float | None = None
    max_dose_um: float | None = None
    times_h: list[float] | None = None
    require_smiles: bool = True


def qc_masks(obs: pd.DataFrame, cfg: QCConfig) -> dict[str, np.ndarray]:
    """One boolean keep-mask per rule, in the order they are reported."""
    masks = {
        "has_dose_time": obs["dose_um"].notna() & obs["time_h"].notna(),
        "replicate_agreement": obs["cc_q75"] >= cfg.min_cc_q75,
        "signature_strength_tas": obs["tas"] >= cfg.min_tas,
        "min_replicates": obs["n_rep"] >= cfg.min_rep,
    }
    if cfg.min_ss is not None:
        masks["signature_strength_ss"] = obs["ss"] >= cfg.min_ss
    if cfg.min_dose_um is not None:
        masks["min_dose"] = obs["dose_um"] >= cfg.min_dose_um
    if cfg.max_dose_um is not None:
        masks["max_dose"] = obs["dose_um"] <= cfg.max_dose_um
    if cfg.times_h:
        masks["time"] = obs["time_h"].isin(cfg.times_h)
    if cfg.require_smiles:
        masks["valid_smiles"] = valid_smiles(obs["smiles"])
    return {rule: np.asarray(mask, dtype=bool) for rule, mask in masks.items()}


def apply_qc(adata: ad.AnnData, cfg: QCConfig) -> tuple[ad.AnnData, dict]:
    keep = np.ones(adata.n_obs, dtype=bool)
    dropped = {}
    for rule, mask in qc_masks(obs_frame(adata), cfg).items():
        dropped[rule] = int((keep & ~mask).sum())
        keep &= mask

    out = adata[keep].copy()
    out.obs["compound_key"] = [
        compound_key(s, k)
        for s, k in zip(out.obs["smiles"], out.obs["inchikey"], strict=True)
    ]
    # A compound without a structure cannot be grouped for a cold-compound split.
    has_key = (out.obs["compound_key"] != "").to_numpy()
    dropped["no_compound_key"] = int((~has_key).sum())
    out = out[has_key].copy()
    out.uns["qc"] = {k: v for k, v in asdict(cfg).items() if v is not None}
    out.uns["source_data_hash"] = adata.uns.get("data_hash", "")
    out.uns["data_hash"] = data_hash(
        np.asarray(out.X), out.obs_names, out.uns.get("release", "")
    )

    report = {
        "source_data_hash": out.uns["source_data_hash"],
        "data_hash": out.uns["data_hash"],
        "config": asdict(cfg),
        "n_in": int(adata.n_obs),
        "n_out": int(out.n_obs),
        "dropped_by_rule": dropped,
        "n_compounds": int(out.obs["compound_key"].nunique()),
        "n_cells": int(out.obs["cell_id"].nunique()),
    }
    return out, report


def write_report(report: dict, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2)
