import pytest

from lincs_processing import pert_info


def test_print_pert_statistics(pert_info_file, capsys):
    pert_info.print_pert_statistics(pert_info_file)
    out = capsys.readouterr().out
    assert "Number of available perturbations: 5" in out
    assert "Number of available trt_cp: 4" in out
    assert "Number of Touchstone of trt_cp perturbations: 2" in out


def test_print_pert_statistics_rejects_unknown_type(pert_info_file):
    with pytest.raises(AssertionError):
        pert_info.print_pert_statistics(pert_info_file, pert_type="trt_xpr")


def test_pert_touchstone(pert_info_file):
    by_id, by_name = pert_info.pert_touchstone(pert_info_file)
    assert by_id == {"BRD-A": 1, "BRD-B": 0, "BRD-C": 1, "BRD-D": 0}
    assert by_name["vorinostat"] == 0
    assert "BRD-E" not in by_id


def test_duplicate_pert_name(pert_info_file):
    assert pert_info.duplicate_pert_name(pert_info_file) == ["aspirin"]


def test_mapping_id_iname(pert_info_file):
    mapping = pert_info.mapping_id_iname(pert_info_file)
    assert mapping["BRD-C"] == "aspirin"
    assert len(mapping) == 5
