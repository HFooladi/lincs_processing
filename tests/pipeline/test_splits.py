import pandas as pd
import pytest

from lincs_processing.pipeline.splits import (
    SPLIT_TYPES,
    LeakageError,
    assert_no_leakage,
    make_split,
    split_hash,
    split_record,
)


@pytest.mark.parametrize("split_type", SPLIT_TYPES)
def test_splits_have_all_folds_and_no_leakage(qc_adata, split_type):
    fold = make_split(qc_adata.obs, split_type, seed=0)
    assert {"train", "val", "test"} <= set(fold)
    assert_no_leakage(qc_adata.obs, fold, split_type)


def test_cold_compound_keeps_every_dose_time_and_duplicate_id_together(qc_adata):
    obs = qc_adata.obs
    fold = make_split(obs, "cold_compound", seed=0)
    per_key = fold.groupby(obs["compound_key"], observed=True).nunique()
    assert (per_key == 1).all()
    # The two BRD ids of caffeine share one fold because they share a key.
    caffeine = fold[obs["pert_iname"].isin(["caffeine", "caffeine-dup"])]
    assert caffeine.nunique() == 1


def test_cold_cell_holds_out_whole_cell_lines(qc_adata):
    obs = qc_adata.obs
    fold = make_split(obs, "cold_cell", seed=0)
    assert (fold.groupby(obs["cell_id"], observed=True).nunique() == 1).all()


def test_cold_both_tests_only_unseen_pairs(qc_adata):
    obs = qc_adata.obs
    fold = make_split(obs, "cold_both", seed=0)
    train, test = obs[fold == "train"], obs[fold == "test"]
    assert not set(train["cell_id"]) & set(test["cell_id"])
    assert not set(train["compound_key"]) & set(test["compound_key"])
    assert (fold == "excluded").any()


def test_split_is_deterministic_and_seed_dependent(qc_adata):
    a = make_split(qc_adata.obs, "cold_compound", seed=0)
    b = make_split(qc_adata.obs, "cold_compound", seed=0)
    c = make_split(qc_adata.obs, "cold_compound", seed=1)
    assert split_hash(a) == split_hash(b)
    assert split_hash(a) != split_hash(c)
    record = split_record(a, "cold_compound", 0, "abc", 0.1, 0.2)
    assert record["split_id"] == f"cold_compound-{split_hash(a)[:10]}"
    assert sum(record["counts"].values()) == len(a)


def test_leak_is_detected():
    obs = pd.DataFrame(
        {"compound_key": ["A", "A", "B"], "cell_id": ["x", "y", "x"]},
        index=["s1", "s2", "s3"],
    )
    leaky = pd.Series(["train", "test", "val"], index=obs.index)
    with pytest.raises(LeakageError):
        assert_no_leakage(obs, leaky, "cold_compound")


def test_too_few_groups_is_an_error(qc_adata):
    obs = qc_adata.obs[qc_adata.obs["cell_id"].isin(["A375", "A549"])]
    with pytest.raises(ValueError, match="groups"):
        make_split(obs, "cold_cell", seed=0)
