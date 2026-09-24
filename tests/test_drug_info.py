from lincs_processing import drug_info


def test_print_drug_statistics(drug_info_file, capsys):
    drug_info.print_drug_statistics(drug_info_file)
    out = capsys.readouterr().out
    assert "Number of available drugs in the datasets: 3" in out
    assert "Number of Unique mechanism of actions: 3" in out
    # PTGS1, PTGS2, HDAC1, HDAC2 and the "nan" of the empty target cell
    assert "Number of Unique targets: 5" in out


def test_drug_pert_retrieval_keeps_touchstone_only(drug_info_file, pert_info_file):
    supp, names = drug_info.drug_pert_retrieval(drug_info_file, pert_info_file)
    # aspirin is touchstone (BRD-A / BRD-C); vorinostat is not; nothere is absent
    assert names == ["aspirin"]
    assert list(supp.moa) == ["cyclooxygenase inhibitor"]


def test_drug_pert_retrieval_without_touchstone_column(tmp_path, drug_info_file):
    path = tmp_path / "pert_info_gse70138.txt"
    path.write_text(
        "pert_id\tpert_iname\tpert_type\n"
        "BRD-A\taspirin\ttrt_cp\n"
        "BRD-B\tvorinostat\ttrt_cp\n"
    )
    _, names = drug_info.drug_pert_retrieval(drug_info_file, str(path))
    assert names == ["aspirin", "vorinostat"]
