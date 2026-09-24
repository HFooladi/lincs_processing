import numpy as np
import pytest

from lincs_processing import utils


def _cells(records):
    return {line[0][0] for line in records}


def test_pickle_round_trip(tmp_path, records):
    path = str(tmp_path / "sample.pkl")
    utils.write_pickle(path, records)
    loaded = utils.load_pickle(path)
    assert len(loaded) == len(records)
    np.testing.assert_array_equal(loaded[3][1], records[3][1])


def test_functions_accept_a_path(tmp_path, records):
    path = str(tmp_path / "sample.pkl")
    utils.write_pickle(path, records)
    assert len(utils.cell_line_list(path, ["PC3"])) == 4


def test_print_helpers_do_not_crash(records, capsys):
    utils.print_statistics(records)
    utils.print_most_frequent(records, n=2)
    out = capsys.readouterr().out
    assert "Number of unique Cell Lines: 3" in out
    assert "Most frequent Cell Lines: [('MCF7', 6), ('PC3', 4)]" in out


def test_cell_line_list_default_is_mcf7(records):
    assert _cells(utils.cell_line_list(records)) == {"MCF7"}


def test_cell_line_list_keeps_only_requested(records):
    out = utils.cell_line_list(records, ["HL60", "PC3"])
    assert _cells(out) == {"HL60", "PC3"}
    assert len(out) == 6


def test_cell_line_frequent(records):
    out = utils.cell_line_frequent(records, n=2)
    assert _cells(out) == {"MCF7", "PC3"}


def test_cell_line_frequent_warns_when_n_too_large(records):
    with pytest.warns(UserWarning):
        out = utils.cell_line_frequent(records, n=10)
    assert len(out) == len(records)


@pytest.mark.parametrize(
    ("indicator", "query", "expected"),
    [
        (0, ["HL60"], 2),
        (1, ["BRD-A"], 5),
        (2, [0.1], 3),
        (3, [6], 4),
    ],
)
def test_parse_list_each_indicator(records, indicator, query, expected):
    assert len(utils.parse_list(records, indicator, query)) == expected


def test_parse_list_rejects_bad_indicator(records):
    with pytest.raises(AssertionError):
        utils.parse_list(records, 4, ["x"])


def test_parse_most_frequent_compounds(records):
    out = utils.parse_most_frequent(records, indicator=1, n=1)
    assert {line[0][1] for line in out} == {"BRD-A"}
    assert len(out) == 5


def test_parse_most_frequent_rejects_n_beyond_unique(records):
    with pytest.raises(AssertionError):
        utils.parse_most_frequent(records, indicator=0, n=4)


def test_parse_chunk_frequent_skips_top_item(records):
    out = utils.parse_chunk_frequent(records, indicator=0, start=1, end=3)
    assert _cells(out) == {"PC3", "HL60"}


def test_parse_chunk_frequent_allows_end_equal_to_unique_count(records):
    # The old implementation asserted end < n_unique, so this was impossible.
    out = utils.parse_chunk_frequent(records, indicator=0, start=0, end=3)
    assert len(out) == len(records)


def test_parse_chunk_frequent_accepts_time_indicator(records):
    out = utils.parse_chunk_frequent(records, indicator=3, start=0, end=1)
    assert {line[0][5] for line in out} == {24}


def test_parse_dose_range_is_exclusive_and_accepts_floats(records):
    out = utils.parse_dose_range(records, dose_min=0.5, dose_max=10.0)
    assert {line[0][3] for line in out} == {1.0}


def test_to_dataframe_layout(records):
    df = utils.to_dataframe(records)
    assert df.shape == (12, 5 + 4)
    assert list(df.columns[-4:]) == ["cell_lines", "compounds", "doses", "times"]
    assert df["cell_lines"].iloc[0] == "MCF7"
    assert df.iloc[0, 0] == 0.0


def test_parse_list_v2_with_data_and_no_path(extended_records):
    out = utils.parse_list_v2(indicator=4, query=[1], data=extended_records)
    assert {line[0][1] for line in out} == {"BRD-A", "BRD-C"}


def test_parse_list_v2_splits_multivalued_fields(extended_records):
    out = utils.parse_list_v2(
        indicator=6, query=["MEK inhibitor"], data=extended_records
    )
    assert {line[0][1] for line in out} == {"BRD-C"}

    out = utils.parse_list_v2(indicator=7, query=["ERBB2"], data=extended_records)
    assert {line[0][1] for line in out} == {"BRD-A"}


def test_parse_list_v2_reads_from_path(tmp_path, extended_records):
    path = str(tmp_path / "ext.pkl")
    utils.write_pickle(path, extended_records)
    out = utils.parse_list_v2(path, indicator=0, query=["HL60"])
    assert len(out) == 2


def test_parse_list_v2_requires_a_source():
    with pytest.raises(AssertionError):
        utils.parse_list_v2(indicator=0, query=["HL60"])
