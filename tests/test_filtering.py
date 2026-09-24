from lincs_processing import filtering
from lincs_processing.utils import load_pickle, write_pickle


def test_filter_records_applies_each_criterion(records):
    out = filtering.filter_records(records, cells=["MCF7"], times=[24])
    assert len(out) == 4
    assert all(line[0][0] == "MCF7" and line[0][5] == 24 for line in out)


def test_main_reads_filters_and_writes(tmp_path, records):
    src = tmp_path / "in.pkl"
    dst = tmp_path / "out.pkl"
    write_pickle(str(src), records)

    filtering.main(
        [
            "--dataset_dir",
            str(src),
            "--compounds",
            "BRD-A",
            "--doses",
            "10",
            "--output_dir",
            str(dst),
        ]
    )

    out = load_pickle(str(dst))
    assert len(out) == 3
    assert {line[0][1] for line in out} == {"BRD-A"}
    assert {line[0][3] for line in out} == {10.0}


def test_cli_help_lists_flags(capsys):
    parser = filtering.build_parser()
    parser.print_help()
    out = capsys.readouterr().out
    for flag in ("--dataset_dir", "--cells", "--compounds", "--doses", "--times"):
        assert flag in out
