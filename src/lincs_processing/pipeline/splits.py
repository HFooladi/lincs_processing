"""Split: group-aware train/val/test assignment, written as a versioned artefact.

The held-out unit is a compound (by InChIKey connectivity block), a cell line,
or both. Every dose, time and replicate of a held-out unit goes to the same
fold, so the test set never contains "the same compound at another dose".
"""

import json

import numpy as np
import pandas as pd

from lincs_processing.pipeline.hashing import text_sha256

SPLIT_TYPES = ("cold_compound", "cold_cell", "cold_both")
FOLDS = ("train", "val", "test", "excluded")

# Which obs columns must not be shared between the training and held-out folds.
_HELD_OUT_KEYS = {
    "cold_compound": ["compound_key"],
    "cold_cell": ["cell_id"],
    "cold_both": ["compound_key", "cell_id"],
}


class LeakageError(AssertionError):
    pass


def _partition(
    groups: pd.Series, val_frac: float, test_frac: float, rng: np.random.Generator
) -> dict[str, str]:
    unique = sorted(groups.unique())
    n = len(unique)
    n_test = max(1, round(test_frac * n))
    n_val = max(1, round(val_frac * n)) if val_frac > 0 else 0
    if n_test + n_val >= n:
        raise ValueError(
            f"only {n} groups in {groups.name!r}; cannot hold out "
            f"{n_test} test + {n_val} val groups and keep any for training"
        )
    order = [unique[i] for i in rng.permutation(n)]
    fold = {g: "test" for g in order[:n_test]}
    fold.update({g: "val" for g in order[n_test : n_test + n_val]})
    fold.update({g: "train" for g in order[n_test + n_val :]})
    return fold


def make_split(
    obs: pd.DataFrame,
    split_type: str,
    seed: int = 0,
    val_frac: float = 0.1,
    test_frac: float = 0.2,
) -> pd.Series:
    """Fold label per signature, indexed like ``obs``."""
    if split_type not in SPLIT_TYPES:
        raise ValueError(f"split_type must be one of {SPLIT_TYPES}")
    rng = np.random.default_rng(seed)

    if split_type in ("cold_compound", "cold_cell"):
        key = _HELD_OUT_KEYS[split_type][0]
        fold = obs[key].map(_partition(obs[key], val_frac, test_frac, rng))
    else:
        c = obs["compound_key"].map(
            _partition(obs["compound_key"], val_frac, test_frac, rng)
        )
        cell = obs["cell_id"].map(_partition(obs["cell_id"], val_frac, test_frac, rng))
        # Train only on (seen compound, seen cell); test only on (unseen, unseen).
        # Pairs with exactly one test-side unit would leak, so they are excluded.
        fold = pd.Series("excluded", index=obs.index)
        fold[(c == "train") & (cell == "train")] = "train"
        touches_val = (c == "val") | (cell == "val")
        fold[touches_val & (c != "test") & (cell != "test")] = "val"
        fold[(c == "test") & (cell == "test")] = "test"

    fold = fold.astype(str).rename("fold")
    for required in ("train", "val", "test"):
        if not (fold == required).any():
            raise ValueError(f"{split_type} split produced an empty {required} fold")
    assert_no_leakage(obs, fold, split_type)
    return fold


def assert_no_leakage(obs: pd.DataFrame, fold: pd.Series, split_type: str) -> None:
    """Held-out units must never appear in the training fold."""
    if obs.index.duplicated().any():
        raise LeakageError("duplicate signature ids: a replicate would span folds")
    for key in _HELD_OUT_KEYS[split_type]:
        train = set(obs.loc[fold == "train", key])
        for held_out in ("val", "test"):
            if split_type == "cold_both" and held_out == "val":
                continue  # cold-both val shares seen cells/compounds by design
            shared = train & set(obs.loc[fold == held_out, key])
            if shared:
                raise LeakageError(
                    f"{split_type}: {len(shared)} {key} value(s) in both train and "
                    f"{held_out}, e.g. {sorted(shared)[:3]}"
                )


def split_hash(fold: pd.Series) -> str:
    lines = sorted(f"{sig}\t{f}" for sig, f in fold.items())
    return text_sha256("\n".join(lines))


def split_record(
    fold: pd.Series,
    split_type: str,
    seed: int,
    data_hash: str,
    val_frac: float,
    test_frac: float,
) -> dict:
    digest = split_hash(fold)
    counts = fold.value_counts()
    return {
        "split_id": f"{split_type}-{digest[:10]}",
        "split_type": split_type,
        "split_hash": digest,
        "data_hash": data_hash,
        "seed": seed,
        "val_frac": val_frac,
        "test_frac": test_frac,
        "counts": {f: int(counts.get(f, 0)) for f in FOLDS},
    }


def write_split(fold: pd.Series, record: dict, parquet_path: str, json_path: str):
    frame = fold.rename_axis("sig_id").reset_index()
    frame.to_parquet(parquet_path, index=False)
    with open(json_path, "w") as fh:
        json.dump(record, fh, indent=2)


def read_split(parquet_path: str) -> pd.Series:
    return pd.read_parquet(parquet_path).set_index("sig_id")["fold"]
