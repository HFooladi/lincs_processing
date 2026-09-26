"""Collect every (split, seed) metrics file into one table."""

import json

import pandas as pd

KEYS = ["split_id", "split_type", "seed"]


def collect(paths: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for path in paths:
        with open(path) as fh:
            result = json.load(fh)
        for method, values in result["methods"].items():
            for metric, value in values.items():
                rows.append(
                    {
                        **{k: result[k] for k in KEYS},
                        "method": method,
                        "metric": metric,
                        "value": value,
                    }
                )
    long = pd.DataFrame(rows).sort_values(KEYS + ["method", "metric"])
    summary = (
        long.groupby(["split_type", "split_id", "method", "metric"])["value"]
        .agg(mean="mean", std="std", n_seeds="count")
        .reset_index()
    )
    summary["std"] = summary["std"].fillna(0.0)
    return long.reset_index(drop=True), summary
