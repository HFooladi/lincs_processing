"""Fixtures built from the synthetic mini release (``pipeline.testdata``)."""

import pytest

pytest.importorskip("anndata")
pytest.importorskip("rdkit")
pytest.importorskip("torch")


@pytest.fixture(scope="session")
def mini_release(tmp_path_factory):
    from lincs_processing.pipeline.testdata import write_all

    root = tmp_path_factory.mktemp("lincs_mini")
    return write_all(str(root), ("beta", "gse92742", "gse70138"))


@pytest.fixture(scope="session")
def raw_adata(mini_release):
    from lincs_processing.pipeline.ingest import ingest

    p = mini_release["beta"]
    return ingest("beta", p["gctx"], p["sig_info"], p["gene_info"], p["pert_info"])


@pytest.fixture(scope="session")
def qc_adata(raw_adata):
    from lincs_processing.pipeline.qc import QCConfig, apply_qc

    out, _ = apply_qc(raw_adata, QCConfig())
    return out
