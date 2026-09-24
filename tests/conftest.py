"""Synthetic LINCS-like fixtures so the tests never need real data files."""

import numpy as np
import pytest

# (cell_line, pert_id, dose, time) for 12 records. Frequencies are deliberately
# uneven so "most frequent" behaviour is testable: MCF7 x6, PC3 x4, HL60 x2.
_META = [
    ("MCF7", "BRD-A", 10.0, 24),
    ("MCF7", "BRD-A", 1.0, 24),
    ("MCF7", "BRD-B", 10.0, 6),
    ("MCF7", "BRD-B", 0.1, 24),
    ("MCF7", "BRD-C", 10.0, 24),
    ("MCF7", "BRD-D", 1.0, 6),
    ("PC3", "BRD-A", 10.0, 24),
    ("PC3", "BRD-A", 0.1, 6),
    ("PC3", "BRD-B", 1.0, 24),
    ("PC3", "BRD-C", 10.0, 24),
    ("HL60", "BRD-A", 10.0, 24),
    ("HL60", "BRD-D", 0.1, 6),
]


def _vector(i: int) -> np.ndarray:
    return np.arange(5, dtype=float) + i


@pytest.fixture
def records() -> list:
    """Seven-field records as produced by ``gctx_parser``."""
    return [
        [(cell, pert, "trt_cp", dose, "um", time, "h"), _vector(i)]
        for i, (cell, pert, dose, time) in enumerate(_META)
    ]


@pytest.fixture
def extended_records() -> list:
    """Eleven-field records as consumed by ``parse_list_v2``.

    Fields 8-10 (clinical phase, MOA, target) are one-element tuples holding a
    ``|``-separated string, matching the shape ``parse_list_v2`` splits on.
    """
    extras = {
        "BRD-A": (1, ("Launched",), ("EGFR inhibitor",), ("EGFR|ERBB2",)),
        "BRD-B": (0, ("Phase 2",), ("HDAC inhibitor",), ("HDAC1",)),
        "BRD-C": (1, ("Launched",), ("EGFR inhibitor|MEK inhibitor",), ("MAP2K1",)),
        "BRD-D": (0, ("Preclinical",), ("unknown",), ("-666",)),
    }
    return [
        [(cell, pert, "trt_cp", dose, "um", time, "h", *extras[pert]), _vector(i)]
        for i, (cell, pert, dose, time) in enumerate(_META)
    ]


@pytest.fixture
def pert_info_file(tmp_path):
    """A tiny ``pert_info.txt`` with a duplicated pert_iname and a touchstone flag."""
    path = tmp_path / "pert_info.txt"
    rows = [
        "pert_id\tpert_iname\tpert_type\tis_touchstone",
        "BRD-A\taspirin\ttrt_cp\t1",
        "BRD-B\tvorinostat\ttrt_cp\t0",
        "BRD-C\taspirin\ttrt_cp\t1",
        "BRD-D\tunknown-cp\ttrt_cp\t0",
        "BRD-E\tEGFR\ttrt_sh\t1",
    ]
    path.write_text("\n".join(rows) + "\n")
    return str(path)


@pytest.fixture
def drug_info_file(tmp_path):
    """A tiny Repurposing Hub export: nine comment lines then a tab table."""
    path = tmp_path / "repurposing_drugs.txt"
    header = [f"! comment line {i}" for i in range(9)]
    rows = [
        "pert_iname\tclinical_phase\tmoa\ttarget",
        "aspirin\tLaunched\tcyclooxygenase inhibitor\tPTGS1|PTGS2",
        "vorinostat\tLaunched\tHDAC inhibitor\tHDAC1|HDAC2",
        "nothere\tPhase 1\tunknown\t",
    ]
    path.write_text("\n".join(header + rows) + "\n", encoding="latin-1")
    return str(path)


@pytest.fixture
def sig_info_file(tmp_path):
    """A GSE70138-style ``sig_info.txt`` with combined dose/time columns."""
    path = tmp_path / "sig_info.txt"
    rows = [
        "sig_id\tpert_id\tcell_id\tpert_idose\tpert_itime",
        "S1\tBRD-A\tMCF7\t10 um\t24 h",
        "S2\tBRD-B\tPC3\t-666\t6 h",
        "S3\tBRD-C\tHL60\t0.1 um\t-666",
    ]
    path.write_text("\n".join(rows) + "\n")
    return str(path)
