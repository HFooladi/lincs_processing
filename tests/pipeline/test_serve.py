import numpy as np
import pytest

pytest.importorskip("fastapi")
pytest.importorskip("gradio")

from fastapi.testclient import TestClient  # noqa: E402

from lincs_processing.pipeline.splits import make_split, split_record  # noqa: E402
from lincs_processing.pipeline.tracking import make_tracker  # noqa: E402
from lincs_processing.pipeline.train import train  # noqa: E402


@pytest.fixture(scope="module")
def client(qc_adata, tmp_path_factory):
    from serve.app import create_app

    out = tmp_path_factory.mktemp("serve")
    fold = make_split(qc_adata.obs, "cold_cell", seed=0)
    split = split_record(fold, "cold_cell", 0, qc_adata.uns["data_hash"], 0.1, 0.2)
    train(
        qc_adata,
        fold,
        split,
        seed=0,
        outdir=str(out),
        tracker=make_tracker("none"),
        epochs=2,
        batch_size=64,
        device="cpu",
    )
    return TestClient(create_app(str(out / "model_seed0")))


def test_health_reports_provenance(client):
    body = client.get("/health").json()
    assert body["status"] == "ok"
    assert body["split_id"].startswith("cold_cell-")
    assert len(body["data_hash"]) == 64


def test_predict_returns_a_signature(client):
    resp = client.post(
        "/predict",
        json={
            "smiles": "CCO",
            "cell_id": "NEW_CELL",
            "dose_um": 10,
            "time_h": 24,
            "top_k": 5,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert len(body["signature"]) == len(body["genes"])
    assert body["cell_seen_in_training"] is False
    assert len(body["top_up"]) == 5
    assert np.isfinite(body["signature"]).all()


def test_bad_smiles_is_rejected(client):
    resp = client.post(
        "/predict",
        json={"smiles": "C1CC(", "cell_id": "A375", "dose_um": 1, "time_h": 6},
    )
    assert resp.status_code == 422


def test_ui_is_mounted(client):
    assert client.get("/ui/").status_code == 200
