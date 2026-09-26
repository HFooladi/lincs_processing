"""Turn LINCS ``.gctx`` files plus their metadata tables into record lists.

Importing this module pulls in ``cmapPy`` and ``h5py``; the rest of the package
does not, so it is not re-exported from ``lincs_processing.__init__``.
"""

from collections import Counter

import numpy as np
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def _landmark_row_ids(
    gene_info_dir: str,
    id_col: str = "gene_id",
    flag_col: str = "is_lm",
    flag_value: str = "1",
) -> pd.Series:
    """Row ids of the landmark genes.

    The column names differ between releases: GEO uses ``pr_gene_id``/``pr_is_lm``
    and the 2020 CLUE beta release uses ``gene_id``/``feature_space == "landmark"``.
    """
    gene_info = pd.read_csv(gene_info_dir, sep="\t", dtype=str)
    print(f"Number of measured genes in the dataset: {gene_info.shape[0]}")

    landmark_gene_row_ids = gene_info[id_col][gene_info[flag_col] == flag_value]
    print(f"Number of landmark genes in the dataset: {landmark_gene_row_ids.shape[0]}")
    return landmark_gene_row_ids


def _records_from_gctoo(query_trt: pd.DataFrame, data_df: pd.DataFrame) -> list[list]:
    """Zip one metadata row and one expression column into a record per sample."""
    query_trt = query_trt.reindex(data_df.columns)
    expression = data_df.to_numpy()

    parse_list = []
    for i, row in enumerate(query_trt.itertuples(index=False)):
        parse_list.append(
            [
                (
                    row.cell_id,
                    row.pert_id,
                    row.pert_type,
                    float(str(row.pert_dose)),
                    row.pert_dose_unit,
                    row.pert_time,
                    row.pert_time_unit,
                ),
                np.array(expression[:, i]),
            ]
        )
    return parse_list


def parsing_level3_cp(
    dataset_dir: str,
    inst_info_dir: str,
    gene_info_dir: str,
    pert_type: str = "trt_cp",
    landmarks: bool = True,
) -> list[list]:
    """Parsing the data to keep desired sig_ids

    This function takes the directory of dataset, perturbation type, and
    whether we want to only keep landmark genes or not. It returns a list
    based on the inputs.


    Parameters
    ----------
    dataset_dir: str
      It must be string file that shows the directory of the dataset.
      dataset should be a gctx file. e.g., valid argument is something like this:
      './Data/Level3_INF_mlr12k_n1319138x12328.gctx'
    inst_info_dir: str
      directory of inst_info. It contains the information about the
      experiment, perturbation type and cell line. For example:
      './Data/inst_info.txt'
    gene_info_dir: str
      directory of gene_info. It contains the information about the genes.
      For example: './Data/gene_info.txt'
    pert_type: str (default= "trt_cp")
      String object that determine which perturbation type you want to parse.
      Default='trt_cp'
    landmarks: bool
      boolean which determines whether you want to just keep landmark genes
      after parsing or you want to keep all the genes. Default=True

    Returns
    ------
    parse_list: List
      Output list (Train, Validation, Test) Format:
      line[0]:(cell_line,
      drug,
      drug_type,
      does,
      does_type,
      time,
      time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
    """
    assert isinstance(pert_type, str), "pert_type must be a string object"
    assert isinstance(landmarks, bool), "landmarks must be a boolean object"

    landmark_gene_row_ids = _landmark_row_ids(gene_info_dir)

    inst_info = pd.read_csv(inst_info_dir, sep="\t")
    print(f"Number of availbale gene expression profiles: {inst_info.shape[0]}")
    print(f"Unique perturbation types: {inst_info.pert_type.unique()}")

    assert pert_type in inst_info.pert_type.unique(), "pert_type is not valid!!"

    query_trt = inst_info[inst_info["pert_type"] == pert_type]
    print(
        f"Number of availbale gene expression profiles of {pert_type}: "
        f"{query_trt.shape[0]}"
    )

    print((query_trt == "-666").sum())
    print(f"Number of different pert_dose_units: {Counter(query_trt.pert_dose_unit)}")

    # Compound treatments come with a handful of dose units; keep only the
    # micromolar ones so doses are comparable across samples.
    if pert_type == "trt_cp":
        query_trt = query_trt[query_trt.pert_dose_unit != "-666"]
        query_trt = query_trt[query_trt.pert_dose_unit == "um"]

    query_ids = query_trt.inst_id
    print(f"Number of samples at the end: {query_ids.shape[0]}")

    print("=================================================================")
    print("Please wait while we are parsing the data ...")

    if landmarks:
        query_gctoo = parse(dataset_dir, rid=landmark_gene_row_ids, cid=query_ids)
    else:
        query_gctoo = parse(dataset_dir, cid=query_ids)

    print("Parse Completed")
    print(f"Size of the data after parsing: {query_gctoo.data_df.shape}")

    query_trt = query_trt.set_index(query_trt.inst_id)
    return _records_from_gctoo(query_trt, query_gctoo.data_df)


def parsing_level5_cp(
    dataset_dir: str,
    sig_info_dir: str,
    gene_info_dir: str,
    pert_type: str = "trt_cp",
    landmarks: bool = True,
    cell_line: str | None = None,
) -> list[list]:
    """Parsing the data to keep desired sig_ids

    This function takes the directory of dataset, perturbation type, and
    whether we want to only keep landmark genes or not. It returns a list
    based on the inputs.


    Parameters
    ----------
    dataset_dir: str
      It must be string file that shows the directory of the dataset.
      dataset should be a gctx file. e.g., valid argument is something like this:
      './Data/Level3_INF_mlr12k_n1319138x12328.gctx'
    sig_info_dir: str
      directory of sig_info. It contains the information about the
      experiment, perturbation type and cell line. For example:
      './Data/sig_info.txt'
    gene_info_dir: str
      directory of gene_info. It contains the information about the genes.
      For example: './Data/gene_info.txt'
    pert_type: str (default="trt_cp")
      String object that determine which perturbation type you want to parse.
      Default='trt_cp'
    landmarks: bool (default=True)
      boolean which determines whether you want to just keep landmark genes
      after parsing or you want to keep all the genes. Default=True
    cell_line: str (default=None)
      Whether you want to select a particular cell_line and parse data just
      for that cell line or not. Default=None Which means parse information of
      all the cell lines.

    Returns
    -------
    parse_list: List
      Output list (Train, Validation, Test) Format:
      line[0]:(cell_line,
      drug,
      drug_type,
      does,
      does_type,
      time,
      time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
    """
    assert isinstance(pert_type, str), "pert_type must be a string object"
    assert isinstance(landmarks, bool), "landmarks must be a boolean object"

    landmark_gene_row_ids = _landmark_row_ids(gene_info_dir)

    sig_info = pd.read_csv(sig_info_dir, sep="\t")
    print(f"Number of availbale gene expression profiles: {sig_info.shape[0]}")
    print(f"Unique perturbation types: {sig_info.pert_type.unique()}")

    assert pert_type in sig_info.pert_type.unique(), "pert_type is not valid!!"

    query_trt = sig_info[sig_info["pert_type"] == pert_type]
    print(
        f"Number of availbale gene expression profiles of {pert_type}: "
        f"{query_trt.shape[0]}"
    )

    print((query_trt == "-666").sum())
    print(f"Number of different pert_dose_units: {Counter(query_trt.pert_dose_unit)}")

    if pert_type == "trt_cp":
        query_trt = query_trt[query_trt.pert_dose_unit != "-666"]

    if cell_line is not None:
        query_trt = query_trt[query_trt["cell_id"] == cell_line]

    query_ids = query_trt.sig_id
    print(f"Number of samples at the end: {query_ids.shape[0]}")

    print("=================================================================")
    print("Please wait while we are parsing the data ...")

    if landmarks:
        query_gctoo = parse(dataset_dir, rid=landmark_gene_row_ids, cid=query_ids)
    else:
        query_gctoo = parse(dataset_dir, cid=query_ids)

    print("Parse Completed")
    print(f"Size of the data after parsing: {query_gctoo.data_df.shape}")

    query_trt = query_trt.set_index(query_trt.sig_id)
    return _records_from_gctoo(query_trt, query_gctoo.data_df)
