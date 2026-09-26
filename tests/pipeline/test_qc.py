from lincs_processing.pipeline.qc import QCConfig, apply_qc


def test_rules_drop_low_quality_and_unparseable(raw_adata, qc_adata):
    obs = qc_adata.obs
    assert (obs["cc_q75"] >= 0.2).all()
    assert (obs["n_rep"] >= 2).all()
    assert "broken" not in set(obs["pert_iname"])
    assert qc_adata.n_obs < raw_adata.n_obs


def test_report_accounts_for_every_dropped_signature(raw_adata):
    out, report = apply_qc(raw_adata, QCConfig(times_h=[24]))
    assert report["n_in"] - report["n_out"] == sum(report["dropped_by_rule"].values())
    assert set(out.obs["time_h"]) == {24.0}
    assert report["data_hash"] != report["source_data_hash"]


def test_same_molecule_under_two_ids_shares_a_compound_key(qc_adata):
    obs = qc_adata.obs
    keys = obs.groupby("pert_iname", observed=True)["compound_key"].first()
    assert keys["caffeine"] == keys["caffeine-dup"]
    assert (
        obs.loc[obs["pert_iname"] == "caffeine", "pert_id"].iloc[0]
        != (obs.loc[obs["pert_iname"] == "caffeine-dup", "pert_id"].iloc[0])
    )
