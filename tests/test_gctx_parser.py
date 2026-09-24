"""The gctx parser needs real LINCS files, so only its wiring is tested here.

Set ``LINCS_DATA_DIR`` to a folder holding ``Level3_INF_mlr12k_n1319138x12328.gctx``,
``inst_info.txt`` and ``gene_info.txt`` to run the slow end-to-end test.
"""

import os

import numpy as np
import pandas as pd
import pytest

from lincs_processing import gctx_parser


def test_module_imports_with_cmappy():
    assert callable(gctx_parser.parsing_level3_cp)
    assert callable(gctx_parser.parsing_level5_cp)


def test_records_from_gctoo_aligns_metadata_to_columns():
    meta = pd.DataFrame(
        {
            "inst_id": ["s1", "s2"],
            "cell_id": ["MCF7", "PC3"],
            "pert_id": ["BRD-A", "BRD-B"],
            "pert_type": ["trt_cp", "trt_cp"],
            "pert_dose": ["10", "0.1"],
            "pert_dose_unit": ["um", "um"],
            "pert_time": [24, 6],
            "pert_time_unit": ["h", "h"],
        }
    ).set_index("inst_id")
    # Columns deliberately in the opposite order to the metadata rows.
    data_df = pd.DataFrame(
        {"s2": [1.0, 2.0, 3.0], "s1": [4.0, 5.0, 6.0]}, index=["g1", "g2", "g3"]
    )

    out = gctx_parser._records_from_gctoo(meta, data_df)

    assert out[0][0] == ("PC3", "BRD-B", "trt_cp", 0.1, "um", 6, "h")
    np.testing.assert_array_equal(out[0][1], [1.0, 2.0, 3.0])
    assert out[1][0][0] == "MCF7"
    np.testing.assert_array_equal(out[1][1], [4.0, 5.0, 6.0])


@pytest.mark.slow
@pytest.mark.skipif("LINCS_DATA_DIR" not in os.environ, reason="needs LINCS files")
def test_parsing_level3_cp_end_to_end():
    root = os.environ["LINCS_DATA_DIR"]
    out = gctx_parser.parsing_level3_cp(
        os.path.join(root, "Level3_INF_mlr12k_n1319138x12328.gctx"),
        os.path.join(root, "inst_info.txt"),
        os.path.join(root, "gene_info.txt"),
    )
    assert out and out[0][1].shape == (978,)
