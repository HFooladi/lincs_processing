"""Registry gate: promote the cold-cell model only if it beats the baseline.

The decision is always written to ``gate_decision.json``. Registration uses
the MLflow model registry at ``MLFLOW_TRACKING_URI``.
"""

import json
import os

import pandas as pd


def gate(
    summary: pd.DataFrame,
    split_type: str = "cold_cell",
    metric: str = "pearson_mean",
    margin: float = 0.0,
) -> dict:
    rows = summary[
        (summary["split_type"] == split_type) & (summary["metric"] == metric)
    ]
    if rows.empty:
        return {"status": "skipped", "reason": f"no {split_type} split in this run"}
    means = rows.set_index("method")["mean"]
    model, baseline = float(means["model"]), float(means["train_mean"])
    passed = model > baseline + margin
    return {
        "status": "pass" if passed else "fail",
        "split_type": split_type,
        "split_id": str(rows["split_id"].iloc[0]),
        "metric": metric,
        "model": model,
        "train_mean_baseline": baseline,
        "margin": margin,
    }


def register(
    package_dir: str, decision: dict, model_name: str, experiment: str
) -> dict:
    from mlflow import MlflowClient
    from mlflow.exceptions import MlflowException

    from lincs_processing.pipeline.tracking import make_tracker

    with open(os.path.join(package_dir, "split.json")) as fh:
        split = json.load(fh)
    tracker = make_tracker("mlflow", experiment, run_name="registry-gate")
    tracker.set_tags(
        {
            "data_hash": split["data_hash"],
            "split_hash": split["split_hash"],
            "split_id": split["split_id"],
            "gate": decision["status"],
        }
    )
    tracker.log_metrics(
        {
            "gate_model": decision["model"],
            "gate_baseline": decision["train_mean_baseline"],
        }
    )
    tracker.log_artifacts(package_dir, "package")
    # The package is a plain artefact directory rather than an MLflow "logged
    # model", so create the version from its artifact URI directly.
    client = MlflowClient()
    try:
        client.create_registered_model(model_name)
    except MlflowException:
        pass  # already registered
    assert tracker.run_id is not None
    source = f"{client.get_run(tracker.run_id).info.artifact_uri}/package"
    version = client.create_model_version(
        model_name,
        source,
        run_id=tracker.run_id,
        tags={k: split[k] for k in ("data_hash", "split_hash", "split_id")},
    )
    client.set_registered_model_alias(model_name, "cold-cell-gated", version.version)
    tracker.end()
    return {"model_name": model_name, "version": str(version.version)}
