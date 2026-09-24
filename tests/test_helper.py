from lincs_processing.helper import sig_info_augment


def test_sig_info_augment_splits_dose_and_time(sig_info_file):
    df = sig_info_augment(sig_info_file)

    assert df.shape == (3, 5 + 4)
    assert list(df.pert_dose) == [10.0, -666.0, 0.1]
    assert list(df.pert_dose_unit) == ["um", "-666", "um"]
    assert list(df.pert_time) == [24.0, 6.0, -666.0]
    assert list(df.pert_time_unit) == ["h", "h", "-666"]
