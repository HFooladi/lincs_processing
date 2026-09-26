import json
import os

import pandas as pd

from lincs_processing.pipeline.collect import collect
from lincs_processing.pipeline.package import package
from lincs_processing.pipeline.registry import gate

METRICS = [
    "pearson_mean",
    "pearson_median",
    "cosine_mean",
    "r2",
    "rmse",
    "deg100_direction",
]


def _metrics_file(tmp_path, split_type, seed, model, baseline, val_loss):
    path = tmp_path / f"metrics_{split_type}_seed{seed}.json"
    methods = {
        "model": dict.fromkeys(METRICS, model),
        "train_mean": dict.fromkeys(METRICS, baseline),
    }
    path.write_text(
        json.dumps(
            {
                "split_id": f"{split_type}-abc",
                "split_type": split_type,
                "seed": seed,
                "best_val_loss": val_loss,
                "methods": methods,
            }
        )
    )
    return str(path)


def test_collect_and_gate(tmp_path):
    files = [
        _metrics_file(tmp_path, "cold_cell", 0, 0.5, 0.2, 1.0),
        _metrics_file(tmp_path, "cold_cell", 1, 0.3, 0.2, 0.9),
        _metrics_file(tmp_path, "cold_compound", 0, 0.1, 0.2, 1.0),
    ]
    long, summary = collect(files)
    assert len(long) == 3 * 2 * len(METRICS)
    row = summary.query(
        "split_type == 'cold_cell' and method == 'model' and metric == 'pearson_mean'"
    )
    assert row["mean"].item() == 0.4 and row["n_seeds"].item() == 2

    assert gate(summary)["status"] == "pass"
    assert gate(summary, margin=0.5)["status"] == "fail"
    assert gate(summary, split_type="cold_compound")["status"] == "fail"
    assert gate(summary, split_type="cold_both")["status"] == "skipped"


def test_package_selects_seed_by_validation_loss(tmp_path):
    split_json = tmp_path / "cold_cell.json"
    split_json.write_text(
        json.dumps(
            {
                "split_id": "cold_cell-abc",
                "split_type": "cold_cell",
                "split_hash": "h" * 64,
                "data_hash": "d" * 64,
                "counts": {"train": 8, "val": 1, "test": 1, "excluded": 0},
            }
        )
    )
    models = []
    for seed in (0, 1):
        d = tmp_path / f"model_seed{seed}"
        d.mkdir()
        (d / "config.json").write_text(
            json.dumps({"seed": seed, "train_params": {"batch_correction": "none"}})
        )
        models.append(str(d))
    metrics = [
        _metrics_file(tmp_path, "cold_cell", 0, 0.5, 0.2, 1.0),
        _metrics_file(tmp_path, "cold_cell", 1, 0.3, 0.2, 0.9),
    ]
    dest = package(str(split_json), models, metrics, str(tmp_path / "packages"))

    assert os.path.basename(dest) == "cold_cell-abc"
    assert json.load(open(os.path.join(dest, "config.json")))["seed"] == 1
    card = open(os.path.join(dest, "MODEL_CARD.md")).read()
    assert "beats the training-mean baseline" in card
    assert "selected seed: 1" in card
    assert (
        pd.read_json(os.path.join(dest, "metrics.json"), typ="series")["selected_seed"]
        == 1
    )
