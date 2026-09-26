"""``lincs-pipeline``: one subcommand per workflow step.

Nextflow processes call these; they are also usable by hand. Heavy imports
happen inside each command so ``--help`` stays fast.

Example::

    lincs-pipeline ingest --release beta --gctx level5.gctx \\
        --sig_info siginfo_beta.txt --gene_info geneinfo_beta.txt \\
        --pert_info compoundinfo_beta.txt --out dataset.zarr
"""

import argparse
import json
import os
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version


def _ingest(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.ingest import ingest, write_dataset

    adata = ingest(
        a.release,
        a.gctx,
        a.sig_info,
        a.gene_info,
        a.pert_info,
        sig_metrics=a.sig_metrics,
        drug_info=a.drug_info,
        pert_type=a.pert_type,
    )
    write_dataset(adata, a.out, a.manifest)
    print(f"ingested {adata.n_obs} signatures x {adata.n_vars} genes -> {a.out}")


def _qc(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.data import load_dataset
    from lincs_processing.pipeline.qc import QCConfig, apply_qc, write_report

    cfg = QCConfig(
        min_cc_q75=a.min_cc_q75,
        min_tas=a.min_tas,
        min_ss=a.min_ss,
        min_rep=a.min_rep,
        min_dose_um=a.min_dose_um,
        max_dose_um=a.max_dose_um,
        times_h=a.times_h,
    )
    out, report = apply_qc(load_dataset(a.adata), cfg)
    out.write_zarr(a.out)
    write_report(report, a.report)
    print(json.dumps(report["dropped_by_rule"], indent=2))


def _split(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.data import load_dataset, obs_frame
    from lincs_processing.pipeline.splits import make_split, split_record, write_split

    adata = load_dataset(a.adata)
    fold = make_split(obs_frame(adata), a.split_type, a.seed, a.val_frac, a.test_frac)
    record = split_record(
        fold, a.split_type, a.seed, adata.uns["data_hash"], a.val_frac, a.test_frac
    )
    write_split(fold, record, a.out_parquet, a.out_json)
    print(json.dumps(record, indent=2))


def _train(a: argparse.Namespace) -> None:
    import pandas as pd

    from lincs_processing.pipeline.data import load_dataset
    from lincs_processing.pipeline.model import ModelConfig
    from lincs_processing.pipeline.splits import read_split
    from lincs_processing.pipeline.tracking import make_tracker
    from lincs_processing.pipeline.train import train

    with open(a.split_json) as fh:
        split = json.load(fh)
    cell_features = (
        pd.read_csv(a.cell_features, sep="\t", index_col=0) if a.cell_features else None
    )
    tracker = make_tracker(
        a.tracker, a.experiment, run_name=f"{split['split_id']}-seed{a.seed}"
    )
    try:
        train(
            load_dataset(a.adata),
            read_split(a.split),
            split,
            a.seed,
            a.outdir,
            tracker,
            epochs=a.epochs,
            batch_size=a.batch_size,
            lr=a.lr,
            weight_decay=a.weight_decay,
            patience=a.patience,
            model_cfg=ModelConfig(a.hidden, a.layers, a.dropout),
            batch_correction=a.batch_correction,
            cell_features=cell_features,
            device=a.device,
        )
    finally:
        tracker.end()


def _evaluate(a: argparse.Namespace) -> None:
    import numpy as np

    from lincs_processing.pipeline.data import load_dataset
    from lincs_processing.pipeline.evaluate import evaluate, write_metrics
    from lincs_processing.pipeline.splits import read_split
    from lincs_processing.pipeline.tracking import make_tracker

    with open(a.train_summary) as fh:
        summary = json.load(fh)
    predictions = dict(np.load(a.predictions))
    result = evaluate(load_dataset(a.adata), read_split(a.split), predictions, summary)
    write_metrics(result, a.out)

    if a.tracker != "none" and summary.get("run_id"):
        tracker = make_tracker(a.tracker, a.experiment, run_id=summary["run_id"])
        tracker.log_metrics(
            {
                f"test_{method}_{metric}": value
                for method, values in result["methods"].items()
                for metric, value in values.items()
            }
        )
        tracker.set_tags({"beats_train_mean": result["beats_train_mean"]})
        tracker.end()
    print(json.dumps(result["methods"]["model"], indent=2))


def _collect(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.collect import collect

    long, summary = collect(a.metrics)
    long.to_csv(a.out, sep="\t", index=False)
    summary.to_csv(a.summary, sep="\t", index=False)
    print(summary[summary["metric"] == "pearson_mean"].to_string(index=False))


def _package(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.package import package

    print(package(a.split_json, a.models, a.metrics, a.outdir))


def _gate(a: argparse.Namespace) -> None:
    import pandas as pd

    from lincs_processing.pipeline.registry import gate, register

    decision = gate(pd.read_csv(a.summary, sep="\t"), a.split_type, a.metric, a.margin)
    if decision["status"] == "pass" and a.tracker == "mlflow":
        package_dir = os.path.join(a.packages, decision["split_id"])
        decision["registered"] = register(
            package_dir, decision, a.model_name, a.experiment
        )
    with open(a.out, "w") as fh:
        json.dump(decision, fh, indent=2)
    print(json.dumps(decision, indent=2))
    if a.strict and decision["status"] == "fail":
        raise SystemExit("registry gate failed: model does not beat the baseline")


def _testdata(a: argparse.Namespace) -> None:
    from lincs_processing.pipeline.testdata import write_all

    for release, paths in write_all(a.outdir, tuple(a.releases)).items():
        print(release, json.dumps(paths, indent=2))


def _versions(a: argparse.Namespace) -> None:
    packages = ["lincs-processing", "anndata", "torch", "rdkit", "mlflow", "cmapPy"]
    lines = [f'"{a.process}":']
    for name in packages:
        try:
            lines.append(f"    {name}: {version(name)}")
        except PackageNotFoundError:
            pass
    print("\n".join(lines))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lincs-pipeline", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("ingest", help="gctx + metadata -> AnnData zarr")
    p.add_argument("--release", required=True, choices=["beta", "gse92742", "gse70138"])
    p.add_argument("--gctx", required=True)
    p.add_argument("--sig_info", required=True)
    p.add_argument("--gene_info", required=True)
    p.add_argument("--pert_info", required=True)
    p.add_argument("--sig_metrics")
    p.add_argument("--drug_info")
    p.add_argument("--pert_type", default="trt_cp")
    p.add_argument("--out", default="dataset.zarr")
    p.add_argument("--manifest", default="ingest_manifest.json")
    p.set_defaults(func=_ingest)

    p = sub.add_parser("qc", help="rule-based filtering, nothing fitted")
    p.add_argument("--adata", required=True)
    p.add_argument("--min_cc_q75", type=float, default=0.2)
    p.add_argument("--min_tas", type=float, default=0.0)
    p.add_argument("--min_ss", type=float)
    p.add_argument("--min_rep", type=int, default=2)
    p.add_argument("--min_dose_um", type=float)
    p.add_argument("--max_dose_um", type=float)
    p.add_argument("--times_h", type=float, nargs="*")
    p.add_argument("--out", default="filtered.zarr")
    p.add_argument("--report", default="qc_report.json")
    p.set_defaults(func=_qc)

    p = sub.add_parser("split", help="group-aware cold splits")
    p.add_argument("--adata", required=True)
    p.add_argument("--split_type", required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--val_frac", type=float, default=0.1)
    p.add_argument("--test_frac", type=float, default=0.2)
    p.add_argument("--out_parquet", default="split.parquet")
    p.add_argument("--out_json", default="split.json")
    p.set_defaults(func=_split)

    p = sub.add_parser("train", help="fit preprocessing + model on one (split, seed)")
    p.add_argument("--adata", required=True)
    p.add_argument("--split", required=True)
    p.add_argument("--split_json", required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--outdir", default=".")
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--batch_size", type=int, default=512)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--patience", type=int, default=5)
    p.add_argument("--hidden", type=int, default=1024)
    p.add_argument("--layers", type=int, default=2)
    p.add_argument("--dropout", type=float, default=0.2)
    p.add_argument(
        "--batch_correction", default="cell_center", choices=["none", "cell_center"]
    )
    p.add_argument("--cell_features")
    p.add_argument("--device", default="auto")
    p.add_argument("--tracker", default="mlflow", choices=["mlflow", "wandb", "none"])
    p.add_argument("--experiment", default="lincs-perturbation")
    p.set_defaults(func=_train)

    p = sub.add_parser("evaluate", help="baselines + model metrics on the test fold")
    p.add_argument("--adata", required=True)
    p.add_argument("--split", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--train_summary", required=True)
    p.add_argument("--out", default="metrics.json")
    p.add_argument("--tracker", default="mlflow", choices=["mlflow", "wandb", "none"])
    p.add_argument("--experiment", default="lincs-perturbation")
    p.set_defaults(func=_evaluate)

    p = sub.add_parser("collect", help="merge metrics files into one table")
    p.add_argument("--metrics", nargs="+", required=True)
    p.add_argument("--out", default="metrics.tsv")
    p.add_argument("--summary", default="summary.tsv")
    p.set_defaults(func=_collect)

    p = sub.add_parser("package", help="bundle the best seed of one split")
    p.add_argument("--split_json", required=True)
    p.add_argument("--models", nargs="+", required=True)
    p.add_argument("--metrics", nargs="+", required=True)
    p.add_argument("--outdir", default="packages")
    p.set_defaults(func=_package)

    p = sub.add_parser("gate", help="registry gate on the cold-cell number")
    p.add_argument("--summary", required=True)
    p.add_argument("--packages", required=True)
    p.add_argument("--split_type", default="cold_cell")
    p.add_argument("--metric", default="pearson_mean")
    p.add_argument("--margin", type=float, default=0.0)
    p.add_argument("--model_name", default="lincs-perturbation")
    p.add_argument("--tracker", default="mlflow", choices=["mlflow", "wandb", "none"])
    p.add_argument("--experiment", default="lincs-perturbation")
    p.add_argument("--strict", action="store_true")
    p.add_argument("--out", default="gate_decision.json")
    p.set_defaults(func=_gate)

    p = sub.add_parser("testdata", help="write the synthetic mini release")
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--releases",
        nargs="+",
        default=["beta"],
        choices=["beta", "gse92742", "gse70138"],
    )
    p.set_defaults(func=_testdata)

    p = sub.add_parser("versions", help="print a versions.yml block")
    p.add_argument("--process", required=True)
    p.set_defaults(func=_versions)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
