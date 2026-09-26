"""Ingest: level-5 ``.gctx`` plus metadata tables -> one AnnData on Zarr.

``X`` holds the landmark-gene signatures (rows are signatures, already expressed
as change from control), ``obs`` the joined signature/compound metadata in the
canonical layout of ``schema.OBS_COLUMNS``, ``var`` the landmark genes and
``uns`` the release plus a sha256 of every raw input file. Nothing is filtered
beyond the perturbation type; that is QC's job.
"""

import datetime as dt
import json
import os
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

from lincs_processing.drug_info import _read_drug_info
from lincs_processing.gctx_parser import _landmark_row_ids
from lincs_processing.helper import sig_info_augment
from lincs_processing.pipeline.hashing import data_hash, file_sha256
from lincs_processing.pipeline.schema import LAYOUTS, MISSING, OBS_COLUMNS, RELEASES

_DOSE_TO_UM = {"um": 1.0, "µm": 1.0, "μm": 1.0, "nm": 1e-3, "mm": 1e3}
_TIME_TO_H = {"h": 1.0, "hr": 1.0, "hrs": 1.0, "m": 1 / 60, "min": 1 / 60}


def _to_numeric(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values.replace(MISSING, np.nan), errors="coerce")


def _scale(value: pd.Series, unit: pd.Series, table: dict[str, float]) -> pd.Series:
    factor = unit.astype(str).str.strip().str.lower().map(table)
    return _to_numeric(value) * factor


def _clean_str(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).replace({MISSING: "", "nan": ""})


def read_landmarks(gene_info: str, release: str) -> pd.DataFrame:
    """Landmark genes, indexed by the row id used in the ``.gctx`` file."""
    layout = LAYOUTS[release]
    ids = _landmark_row_ids(
        gene_info,
        id_col=layout.gene_id,
        flag_col=layout.landmark_col,
        flag_value=layout.landmark_value,
    )
    table = pd.read_csv(gene_info, sep="\t", dtype=str).set_index(layout.gene_id)
    genes = table.loc[ids.to_numpy(), [layout.gene_symbol]]
    genes.columns = ["gene_symbol"]
    genes.index.name = "gene_id"
    return genes


def read_signatures(
    sig_info: str,
    release: str,
    sig_metrics: str | None = None,
    pert_type: str = "trt_cp",
) -> pd.DataFrame:
    """Signature metadata in the canonical layout (compound columns excluded)."""
    layout = LAYOUTS[release]
    if release == "gse70138":
        # GSE70138 packs dose and time into "10 um"/"24 h"; split them out.
        sigs = sig_info_augment(sig_info)
    else:
        sigs = pd.read_csv(sig_info, sep="\t", low_memory=False)
    sigs = sigs[sigs["pert_type"] == pert_type]

    if not layout.metrics_in_sig_info:
        if sig_metrics is None:
            raise ValueError(f"release {release!r} needs --sig_metrics")
        metrics = pd.read_csv(sig_metrics, sep="\t", low_memory=False)
        keep = ["sig_id", layout.cc_col, layout.ss_col, "tas", layout.nrep_col]
        sigs = sigs.merge(metrics[keep], on="sig_id", how="left")

    return pd.DataFrame(
        {
            "sig_id": sigs["sig_id"].astype(str),
            "pert_id": sigs["pert_id"].astype(str),
            "pert_iname": _clean_str(sigs[layout.iname_col]),
            "cell_id": sigs[layout.cell_col].astype(str),
            "dose_um": _scale(sigs["pert_dose"], sigs["pert_dose_unit"], _DOSE_TO_UM),
            "time_h": _scale(sigs["pert_time"], sigs["pert_time_unit"], _TIME_TO_H),
            "cc_q75": _to_numeric(sigs[layout.cc_col]),
            "ss": _to_numeric(sigs[layout.ss_col]),
            "tas": _to_numeric(sigs["tas"]),
            "n_rep": _to_numeric(sigs[layout.nrep_col]),
        }
    ).reset_index(drop=True)


def _join_unique(values: pd.Series) -> str:
    items = {
        part
        for value in _clean_str(values)
        for part in value.split("|")
        if part.strip()
    }
    return "|".join(sorted(items))


def read_compounds(
    pert_info: str, release: str, drug_info: str | None = None
) -> pd.DataFrame:
    """Per-compound SMILES, InChIKey, MOA and target, indexed by ``pert_id``."""
    info = pd.read_csv(pert_info, sep="\t", dtype=str)
    info = info.rename(columns={"canonical_smiles": "smiles", "inchi_key": "inchikey"})
    for column in ("moa", "target"):
        if column not in info:
            info[column] = ""

    # compoundinfo_beta has one row per (compound, target); collapse to one row.
    compounds = info.groupby("pert_id").agg(
        smiles=("smiles", lambda s: _clean_str(s).iloc[0]),
        inchikey=("inchikey", lambda s: _clean_str(s).iloc[0]),
        moa=("moa", _join_unique),
        target=("target", _join_unique),
    )

    if drug_info is not None:
        name_col = LAYOUTS[release].iname_col
        names = info.drop_duplicates("pert_id").set_index("pert_id")[name_col]
        hub = _read_drug_info(drug_info).drop_duplicates("pert_iname")
        hub = hub.set_index("pert_iname")
        for column in ("moa", "target"):
            fill = names.reindex(compounds.index).map(hub[column])
            empty = compounds[column] == ""
            compounds.loc[empty, column] = _clean_str(fill[empty])
    return compounds


def read_matrix(
    gctx: str, row_ids: list[str], col_ids: list[str], chunk: int = 20_000
) -> np.ndarray:
    """Signatures x genes, rows in ``col_ids`` order and columns in ``row_ids``."""
    blocks = []
    for start in range(0, len(col_ids), chunk):
        cids = col_ids[start : start + chunk]
        data = parse(gctx, rid=row_ids, cid=cids).data_df
        blocks.append(data.loc[row_ids, cids].to_numpy(dtype=np.float32).T)
    return np.concatenate(blocks, axis=0)


def ingest(
    release: str,
    gctx: str,
    sig_info: str,
    gene_info: str,
    pert_info: str,
    sig_metrics: str | None = None,
    drug_info: str | None = None,
    pert_type: str = "trt_cp",
) -> ad.AnnData:
    if release not in RELEASES:
        raise ValueError(f"release must be one of {RELEASES}, got {release!r}")

    genes = read_landmarks(gene_info, release)
    in_gctx_rows = set(parse(gctx, row_meta_only=True).index.astype(str))
    genes = genes[genes.index.isin(in_gctx_rows)]

    sigs = read_signatures(sig_info, release, sig_metrics, pert_type)
    in_gctx_cols = set(parse(gctx, col_meta_only=True).index.astype(str))
    sigs = sigs[sigs["sig_id"].isin(in_gctx_cols)].drop_duplicates("sig_id")

    compounds = read_compounds(pert_info, release, drug_info)
    obs = sigs.join(compounds, on="pert_id")
    for column in ("smiles", "inchikey", "moa", "target"):
        obs[column] = _clean_str(obs[column])
    obs = obs[OBS_COLUMNS].set_index("sig_id", drop=False)
    obs.index.name = None

    X = read_matrix(gctx, list(genes.index), list(obs.index))

    inputs = {
        "gctx": gctx,
        "sig_info": sig_info,
        "gene_info": gene_info,
        "pert_info": pert_info,
        "sig_metrics": sig_metrics,
        "drug_info": drug_info,
    }
    uns = {
        "release": release,
        "pert_type": pert_type,
        "created": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "lincs_processing": version("lincs-processing"),
        "source_files": {
            name: {"file": os.path.basename(path), "sha256": file_sha256(path)}
            for name, path in inputs.items()
            if path is not None
        },
        "data_hash": data_hash(X, obs.index, release),
    }
    var = genes.copy()
    var.index = var.index.astype(str)
    return ad.AnnData(X=X, obs=obs, var=var, uns=uns)


def manifest(adata: ad.AnnData) -> dict:
    """The provenance record written next to every dataset version."""
    return {
        "release": adata.uns["release"],
        "data_hash": adata.uns["data_hash"],
        "n_signatures": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_compounds": int(adata.obs["pert_id"].nunique()),
        "n_cells": int(adata.obs["cell_id"].nunique()),
        "source_files": adata.uns["source_files"],
        "created": adata.uns["created"],
    }


def write_dataset(adata: ad.AnnData, zarr_path: str, manifest_path: str) -> None:
    adata.write_zarr(zarr_path)
    with open(manifest_path, "w") as fh:
        json.dump(manifest(adata), fh, indent=2, default=str)
