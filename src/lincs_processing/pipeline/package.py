"""Package: one self-contained, published model directory per split.

A package holds the best seed's model (by validation loss, never test), its
preprocessing, the config, the data and split hashes, the metrics of every
seed and a generated model card.
"""

import json
import os
import shutil

import numpy as np

from lincs_processing.pipeline.evaluate import DEG_K

_CARD_METRICS = ["pearson_mean", "cosine_mean", "r2", "rmse", f"deg{DEG_K}_direction"]


def _load(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def _summarise(per_seed: list[dict]) -> dict:
    out: dict[str, dict] = {}
    for method in per_seed[0]["methods"]:
        out[method] = {}
        for metric in per_seed[0]["methods"][method]:
            values = [m["methods"][method][metric] for m in per_seed]
            out[method][metric] = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
            }
    return out


def model_card(split: dict, selected: dict, summary: dict, n_seeds: int) -> str:
    header = "| method | " + " | ".join(_CARD_METRICS) + " |"
    lines = [header, "|" + "---|" * (len(_CARD_METRICS) + 1)]
    for method, values in summary.items():
        cells = [
            f"{values[m]['mean']:.3f} ± {values[m]['std']:.3f}" for m in _CARD_METRICS
        ]
        lines.append(f"| {method} | " + " | ".join(cells) + " |")
    model_p = summary["model"]["pearson_mean"]["mean"]
    base_p = summary["train_mean"]["pearson_mean"]["mean"]
    verdict = "beats" if model_p > base_p else "does NOT beat"
    counts = ", ".join(f"{k}={v}" for k, v in split["counts"].items())
    cfg = selected["config"]
    return f"""# Model card: LINCS perturbation response, {split["split_id"]}

Predicts a level-5 L1000 signature (landmark genes, change from vehicle
control) from compound structure (Morgan fingerprint), cell line, dose and time.

## Provenance
- data hash: `{split["data_hash"]}`
- split: `{split["split_id"]}` ({split["split_type"]}), hash `{split["split_hash"]}`
- fold sizes: {counts}
- selected seed: {cfg["seed"]} (lowest validation loss of {n_seeds}; test was
  not used for selection)
- batch correction: {cfg["train_params"]["batch_correction"]}

## Held-out performance (mean ± sd over {n_seeds} seeds)
{chr(10).join(lines)}

On pearson_mean the model {verdict} the training-mean baseline
({model_p:.3f} vs {base_p:.3f}).

## How it was made
- QC applied fixed thresholds only; nothing was fitted before the split.
- Scaling, cell vocabulary and batch correction were fitted on the training
  fold inside the training task.
- A held-out unit keeps all its doses, times and replicates in one fold;
  compounds are grouped by InChIKey connectivity block.

## Limitations
- Without external cell features, unseen cell lines are encoded as an
  all-zero vector: cold-cell predictions then fall back on the compound,
  dose and time effect plus the training-mean offset.
- Only the landmark genes are predicted.
"""


def package(
    split_json: str, model_dirs: list[str], metrics_paths: list[str], outdir: str
) -> str:
    split = _load(split_json)
    per_seed = sorted((_load(p) for p in metrics_paths), key=lambda m: m["seed"])
    by_seed = {_load(os.path.join(d, "config.json"))["seed"]: d for d in model_dirs}
    best = min(per_seed, key=lambda m: m["best_val_loss"])
    source = by_seed[best["seed"]]

    dest = os.path.join(outdir, split["split_id"])
    shutil.copytree(source, dest, dirs_exist_ok=True)
    summary = _summarise(per_seed)
    with open(os.path.join(dest, "metrics.json"), "w") as fh:
        json.dump(
            {"selected_seed": best["seed"], "summary": summary, "per_seed": per_seed},
            fh,
            indent=2,
        )
    shutil.copy(split_json, os.path.join(dest, "split.json"))
    card = model_card(
        split,
        {"config": _load(os.path.join(source, "config.json"))},
        summary,
        len(per_seed),
    )
    with open(os.path.join(dest, "MODEL_CARD.md"), "w") as fh:
        fh.write(card)
    return dest
