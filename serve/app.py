"""Serve a packaged model: FastAPI endpoint plus a Gradio demo.

This sits outside the workflow and only consumes a published package
directory (``results/packages/<split_id>``)::

    LINCS_MODEL_DIR=results/packages/cold_cell-xxxx \\
        uv run uvicorn serve.app:create_app --factory --port 8000

Then ``POST /predict`` or open http://localhost:8000/ui.
"""

import json
import os

import gradio as gr
import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from lincs_processing.pipeline.features import mol_from_smiles
from lincs_processing.pipeline.model import Bundle


class PredictRequest(BaseModel):
    smiles: str
    cell_id: str
    dose_um: float = Field(gt=0)
    time_h: float = Field(gt=0)
    top_k: int = Field(default=10, ge=1, le=100)


class GeneScore(BaseModel):
    gene: str
    z: float


class PredictResponse(BaseModel):
    split_id: str
    cell_seen_in_training: bool
    genes: list[str]
    signature: list[float]
    top_up: list[GeneScore]
    top_down: list[GeneScore]


def _predict(bundle: Bundle, req: PredictRequest) -> PredictResponse:
    if mol_from_smiles(req.smiles) is None:
        raise ValueError(f"cannot parse SMILES {req.smiles!r}")
    obs = pd.DataFrame(
        {
            "smiles": [req.smiles],
            "cell_id": [req.cell_id],
            "dose_um": [req.dose_um],
            "time_h": [req.time_h],
        }
    )
    signature = bundle.predict(obs)[0]
    order = np.argsort(signature)
    k = min(req.top_k, len(order))
    return PredictResponse(
        split_id=bundle.config["split_id"],
        cell_seen_in_training=req.cell_id in bundle.pre.cells,
        genes=bundle.genes,
        signature=signature.round(4).tolist(),
        top_up=[
            GeneScore(gene=bundle.genes[i], z=signature[i]) for i in order[::-1][:k]
        ],
        top_down=[GeneScore(gene=bundle.genes[i], z=signature[i]) for i in order[:k]],
    )


def _demo(bundle: Bundle) -> gr.Blocks:
    def run(smiles: str, cell_id: str, dose_um: float, time_h: float):
        try:
            out = _predict(
                bundle,
                PredictRequest(
                    smiles=smiles, cell_id=cell_id, dose_um=dose_um, time_h=time_h
                ),
            )
        except ValueError as err:
            raise gr.Error(str(err)) from err
        up = pd.DataFrame([g.model_dump() for g in out.top_up])
        down = pd.DataFrame([g.model_dump() for g in out.top_down])
        return up, down

    with gr.Blocks(title="LINCS perturbation response") as demo:
        gr.Markdown(
            f"### Predicted L1000 signature\nModel `{bundle.config['split_id']}`"
        )
        with gr.Row():
            smiles = gr.Textbox(label="SMILES", value="CC(=O)OC1=CC=CC=C1C(=O)O")
            cell = gr.Dropdown(
                choices=bundle.pre.cells,
                value=bundle.pre.cells[0] if bundle.pre.cells else None,
                allow_custom_value=True,
                label="Cell line",
            )
            dose = gr.Number(label="Dose (µM)", value=10.0)
            time = gr.Number(label="Time (h)", value=24.0)
        button = gr.Button("Predict")
        with gr.Row():
            up = gr.Dataframe(label="Most up-regulated")
            down = gr.Dataframe(label="Most down-regulated")
        button.click(run, [smiles, cell, dose, time], [up, down])
    return demo


def create_app(model_dir: str | None = None) -> FastAPI:
    model_dir = model_dir or os.environ.get("LINCS_MODEL_DIR")
    if not model_dir:
        raise RuntimeError("set LINCS_MODEL_DIR to a results/packages/<split_id> dir")
    bundle = Bundle(model_dir)
    metrics_path = os.path.join(model_dir, "metrics.json")
    metrics = {}
    if os.path.exists(metrics_path):
        with open(metrics_path) as fh:
            metrics = json.load(fh)["summary"]

    app = FastAPI(title="LINCS perturbation response", version="0.3.0")

    @app.get("/health")
    def health() -> dict:
        return {
            "status": "ok",
            "split_id": bundle.config["split_id"],
            "data_hash": bundle.config["data_hash"],
            "split_hash": bundle.config["split_hash"],
            "n_genes": len(bundle.genes),
            "held_out_metrics": metrics.get("model", {}),
        }

    @app.get("/model-card")
    def model_card() -> dict:
        path = os.path.join(model_dir, "MODEL_CARD.md")
        text = open(path).read() if os.path.exists(path) else ""
        return {"markdown": text}

    @app.post("/predict", response_model=PredictResponse)
    def predict(req: PredictRequest) -> PredictResponse:
        try:
            return _predict(bundle, req)
        except ValueError as err:
            raise HTTPException(status_code=422, detail=str(err)) from err

    return gr.mount_gradio_app(app, _demo(bundle), path="/ui")
