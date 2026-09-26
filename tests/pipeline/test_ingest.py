import json

import numpy as np
import pytest

from lincs_processing.pipeline.ingest import ingest, write_dataset
from lincs_processing.pipeline.schema import OBS_COLUMNS
from lincs_processing.pipeline.testdata import CELLS, COMPOUNDS, N_LANDMARK


def _ingest(release, paths):
    return ingest(
        release,
        paths["gctx"],
        paths["sig_info"],
        paths["gene_info"],
        paths["pert_info"],
        sig_metrics=paths.get("sig_metrics"),
    )


@pytest.mark.parametrize("release", ["beta", "gse92742", "gse70138"])
def test_every_release_maps_to_the_same_dataset(mini_release, raw_adata, release):
    adata = _ingest(release, mini_release[release])

    n_expected = len(COMPOUNDS) * len(CELLS) * 3 * 2  # vehicle controls dropped
    assert adata.shape == (n_expected, N_LANDMARK)
    assert list(adata.obs.columns) == OBS_COLUMNS
    assert set(adata.obs["pert_id"]) == set(raw_adata.obs["pert_id"])
    assert sorted(adata.obs["dose_um"].unique()) == [0.1, 1.0, 10.0]
    assert sorted(adata.obs["time_h"].unique()) == [6.0, 24.0]
    np.testing.assert_allclose(
        adata[raw_adata.obs_names].X, raw_adata.X, rtol=1e-6, atol=1e-6
    )


def test_compoundinfo_rows_collapse_to_one_compound(raw_adata):
    vorinostat = raw_adata.obs[raw_adata.obs["pert_iname"] == "vorinostat"]
    assert set(vorinostat["target"]) == {"HDAC1|HDAC2"}
    assert set(vorinostat["moa"]) == {"HDAC inhibitor"}


def test_geo_requires_sig_metrics(mini_release):
    paths = dict(mini_release["gse92742"], sig_metrics=None)
    with pytest.raises(ValueError, match="sig_metrics"):
        _ingest("gse92742", paths)


def test_write_dataset_records_provenance(raw_adata, tmp_path):
    zarr_path, manifest_path = tmp_path / "d.zarr", tmp_path / "m.json"
    write_dataset(raw_adata, str(zarr_path), str(manifest_path))

    manifest = json.loads(manifest_path.read_text())
    assert manifest["data_hash"] == raw_adata.uns["data_hash"]
    assert set(manifest["source_files"]) == {
        "gctx",
        "sig_info",
        "gene_info",
        "pert_info",
    }
    assert all(len(f["sha256"]) == 64 for f in manifest["source_files"].values())

    from lincs_processing.pipeline.data import load_dataset

    back = load_dataset(str(zarr_path))
    assert back.obs["cell_id"].dtype == object
    np.testing.assert_array_equal(back.X, raw_adata.X)
