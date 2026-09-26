"""Experiment tracking: MLflow by default, W&B optional, or nothing.

The tracking URI comes from the environment (``MLFLOW_TRACKING_URI``), not from
task inputs, so changing where runs are logged never invalidates Nextflow's
task cache.
"""

import os
import random
import time
from collections.abc import Callable
from typing import Any

TRACKERS = ("mlflow", "wandb", "none")


def _retry(fn: Callable[[], Any], attempts: int = 6) -> Any:
    # Parallel training tasks may race to create the SQLite schema or the
    # experiment; back off and try again instead of failing the task.
    for attempt in range(attempts):
        try:
            return fn()
        except Exception:
            if attempt == attempts - 1:
                raise
            time.sleep(0.5 * 2**attempt + random.random())


class Tracker:
    """Minimal logging surface shared by the training and evaluation steps."""

    run_id: str | None = None

    def log_params(self, params: dict) -> None: ...
    def log_metrics(self, metrics: dict, step: int | None = None) -> None: ...
    def set_tags(self, tags: dict) -> None: ...
    def log_artifacts(self, directory: str, path: str | None = None) -> None: ...
    def end(self) -> None: ...


class MlflowTracker(Tracker):
    def __init__(
        self, experiment: str, run_name: str | None = None, run_id: str | None = None
    ) -> None:
        import mlflow

        self.mlflow = mlflow

        def setup() -> None:
            if mlflow.get_experiment_by_name(experiment) is None:
                artifacts = os.environ.get("LINCS_MLFLOW_ARTIFACTS")
                try:
                    mlflow.create_experiment(experiment, artifact_location=artifacts)
                except Exception:
                    if mlflow.get_experiment_by_name(experiment) is None:
                        raise
            mlflow.set_experiment(experiment)

        _retry(setup)
        run = _retry(lambda: mlflow.start_run(run_id=run_id, run_name=run_name))
        self.run_id = run.info.run_id

    def log_params(self, params: dict) -> None:
        _retry(lambda: self.mlflow.log_params({k: str(v) for k, v in params.items()}))

    def log_metrics(self, metrics: dict, step: int | None = None) -> None:
        _retry(lambda: self.mlflow.log_metrics(metrics, step=step))

    def set_tags(self, tags: dict) -> None:
        _retry(lambda: self.mlflow.set_tags({k: str(v) for k, v in tags.items()}))

    def log_artifacts(self, directory: str, path: str | None = None) -> None:
        _retry(lambda: self.mlflow.log_artifacts(directory, artifact_path=path))

    def end(self) -> None:
        _retry(self.mlflow.end_run)


class WandbTracker(Tracker):
    def __init__(
        self, experiment: str, run_name: str | None = None, run_id: str | None = None
    ) -> None:
        import wandb  # not a declared dependency; install it into the image/venv

        self.wandb = wandb
        self.run = wandb.init(
            project=experiment, name=run_name, id=run_id, resume="allow"
        )
        self.run_id = self.run.id

    def log_params(self, params: dict) -> None:
        self.run.config.update(params, allow_val_change=True)

    def log_metrics(self, metrics: dict, step: int | None = None) -> None:
        self.run.log(metrics, step=step)

    def set_tags(self, tags: dict) -> None:
        self.run.summary.update(tags)

    def log_artifacts(self, directory: str, path: str | None = None) -> None:
        artifact = self.wandb.Artifact(path or "model", type="model")
        artifact.add_dir(directory)
        self.run.log_artifact(artifact)

    def end(self) -> None:
        self.run.finish()


def make_tracker(
    kind: str,
    experiment: str = "lincs-perturbation",
    run_name: str | None = None,
    run_id: str | None = None,
) -> Tracker:
    if kind == "mlflow":
        return MlflowTracker(experiment, run_name, run_id)
    if kind == "wandb":
        return WandbTracker(experiment, run_name, run_id)
    if kind == "none":
        return Tracker()
    raise ValueError(f"tracker must be one of {TRACKERS}")
