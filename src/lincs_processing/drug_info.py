"""
Information about the compounds in the Drug Repurposing Hub.
Provided are annotations for 6,125 drug and tool compounds
(2,369 FDA-approved drugs, 1,619 drugs that reached phases 1-3 of clinical
development, 96 compounds that were previously approved but withdrawn from use,
and 2,041 preclinical or tool compounds). Annotations include compound name,
chemical structure, clinical trial status, mechanism of action, protein targets,
disease areas, approved indications (where applicable), purity of the purchased
sample, and vendor ID.
"""

from collections import Counter

import pandas as pd

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def _read_drug_info(drug_info_dir: str) -> pd.DataFrame:
    """Read a Repurposing Hub ``repurposing_drugs_*.txt`` export.

    The file starts with nine comment lines and is latin-1 encoded.
    """
    return pd.read_csv(drug_info_dir, sep="\t", skiprows=9, encoding="latin-1")


def print_drug_statistics(drug_info_dir: str) -> None:
    """Print basic statisctics about drug repurposing hub

    This function takes the directory of drug_info.txt
    and print some information about perturbations.

    Parameters
    ----------
    drug_info_dir: str
      The directory of drug_info file. E.g., './Data/repurposing_drugs_20180907.txt'
    """
    assert isinstance(drug_info_dir, str), "The dataset_dir must be a string object"

    drug_info = _read_drug_info(drug_info_dir)

    print("=================================================================")
    print("Data Statistics\n")
    print(f"Number of available drugs in the datasets: {drug_info.shape[0]}")
    print(f"Number of Unique mechanism of actions: {drug_info.moa.unique().shape[0]}")

    print(
        "Information about Number of drugs in different clinical phases: "
        f"{Counter(drug_info.clinical_phase)}"
    )

    targets = {item for targets in drug_info.target for item in str(targets).split("|")}
    print(f"Number of Unique targets: {len(targets)}")


def drug_pert_retrieval(
    drug_info_dir: str, pert_info_dir: str, pert_type: str = "trt_cp"
) -> tuple[pd.DataFrame, list]:
    """Mode of action and supplementary information of drugs.

    This function takes the directory of drug_info.txt, pert_info.txt,
    and pert_type; and return drug information correspondant to
    perturbation information.

    Parameters
    ----------
    drug_info_dir: str
      The directory of drug_info file. E.g., './Data/repurposing_drugs_20180907.txt'
    pert_info_dir: str
      The directory of pert_info file. E.g., './Data/pert_info.txt'
    pert_type: str, optional (default 'trt_cp')
      perturbation type that you want to extract information about it.
      Default='trt_cp'

    Returns
    -------
    pert_supp_info: pandas.DataFrame
      Pandas dataframe which contains supplementary information
      about touchstone compounds such as mode of actions, targets and ...
    pert_list: List[str]
      List of strings which indicates the compounds that we have additional
      information for them.
      Example: ['abiraterone', 'ABT-737', 'ABT-751', 'AC-55649',...]
    """
    assert isinstance(drug_info_dir, str), "The dataset_dir must be a string object"
    assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
    assert isinstance(pert_type, str), "The pert_type must be a string object"

    drug_info = _read_drug_info(drug_info_dir)
    pert_info = pd.read_csv(pert_info_dir, sep="\t")

    assert pert_type in pert_info.pert_type.unique(), (
        "pert_type should be in the list of available perturbations"
    )

    x = pert_info[pert_info.pert_type == pert_type]
    if "is_touchstone" in x.columns:
        x = x[x.is_touchstone.astype(bool)]

    print("=================================================================")
    print(f"Number of Touchstone of {pert_type}: {x.shape[0]}")

    pert_supp_info = drug_info[drug_info.pert_iname.isin(x.pert_iname)]
    pert_list = list(pert_supp_info.pert_iname)

    print(
        f"Number of {pert_type} that we have additional information about them: "
        f"{len(pert_list)}"
    )

    return pert_supp_info, pert_list
