"""Reading dataset versions back from Zarr."""

import anndata as ad
import pandas as pd


def load_dataset(path: str) -> ad.AnnData:
    """Read an AnnData store with plain string (not categorical) obs columns.

    anndata stores string columns as categoricals; grouping on those also
    yields empty groups for categories absent from a fold, so undo it here.
    """
    adata = ad.read_zarr(path)
    for column in adata.obs.columns:
        if isinstance(adata.obs[column].dtype, pd.CategoricalDtype):
            adata.obs[column] = adata.obs[column].astype(str)
    return adata


def obs_frame(adata: ad.AnnData) -> pd.DataFrame:
    """``adata.obs`` typed as a DataFrame (anndata also allows lazy Dataset2D)."""
    obs = adata.obs
    assert isinstance(obs, pd.DataFrame)
    return obs
