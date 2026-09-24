from collections import Counter

import pandas as pd

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def print_pert_statistics(pert_info_dir: str, pert_type: str = "trt_cp") -> None:
    """Print pert_info statistics

    This function takes the directory of pert_info.txt
    file and perturbation type, and print some information about perturbations.

    Parameters
    ----------
    pert_info_dir: str
      The directory of pert_info file. E.g., './Data/pert_info.txt'
    pert_type: str, optional (default 'trt_cp')
      perturbation type that you want to extract information
      about it. Default='trt_cp'
    """
    assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
    assert isinstance(pert_type, str), "The pert_type must be a string object"

    pert_info = pd.read_csv(pert_info_dir, sep="\t")

    assert pert_type in pert_info.pert_type.unique(), (
        "pert_type should be in the list of available perturbations"
    )

    print("Data Statistics\n")
    print(f"Number of available perturbations: {pert_info.shape[0]}")
    print(f"Number of all Touchstone perturbations: {pert_info.is_touchstone.sum()}")

    x = pert_info[pert_info.pert_type == pert_type]
    print(f"Number of available {pert_type}: {x.shape[0]}")
    print(f"Number of Touchstone of {pert_type} perturbations: {x.is_touchstone.sum()}")


def pert_touchstone(pert_info_dir: str, pert_type: str = "trt_cp") -> tuple[dict, dict]:
    """Check the Touchstone

    This function takes the directory of pert_info.txt
    file and perturbation type, and return a dictionary that
    maps perturbations to 0 and 1 (whether they are touchstone or not)
    For example: {'BRD-A00100033':1, BRD-A00150179:0, ...}


    Parameters
    ----------
    pert_info_dir: str
      The directory of pert_info file. E.g., './Data/pert_info.txt'
    pert_type: str, optional (default 'trt_cp')
      perturbation type that you want to extract information
      about it. Default='trt_cp'

    Returns
    -------
    pert_dict_id: Dict
      A dictionary that maps perturbation ids
      to 0 and 1 (whether they are touchstone or not)
      Keys are pert_id (str) and values are 0 or 1.
      For example: {'BRD-A00100033':1, BRD-A00150179:0, ...}
    pert_dict_iname: Dict
      A dictionary that maps perturbation names
      to 0 and 1 (whether they are touchstone or not)
      Keys are pert_iname (str) and values are 0 or 1.
      For example:
      {'BRD-A00100033':1, BRD-A00150179:0, ...}
    """
    assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
    assert isinstance(pert_type, str), "The pert_type must be a string object"

    pert_info = pd.read_csv(pert_info_dir, sep="\t")

    assert pert_type in pert_info.pert_type.unique(), (
        "pert_type should be in the list of available perturbations"
    )

    x = pert_info[pert_info.pert_type == pert_type]

    pert_dict_id = dict(zip(x.pert_id, x.is_touchstone, strict=True))
    pert_dict_iname = dict(zip(x.pert_iname, x.is_touchstone, strict=True))

    return pert_dict_id, pert_dict_iname


def duplicate_pert_name(pert_info_dir: str) -> list[str]:
    """Checking duplicate pert_name

    There are some perturbation name that maps to more than one
    perturbation ID in Lincs dataset. Here, I am going to find those
    perturbations that have more than one IDs.
    You can find more information about this by following this link:
    https://clue.io/connectopedia/some_perts_have_over_one_brdid


    Parameters
    ----------
    pert_info_dir: str
      The directory of pert_info file. E.g., './Data/pert_info.txt'

    Returns
    -------
    duplicate_list: List
      A list of pert_inames that have multiple pert_ids.
    """
    assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"

    pert_info = pd.read_csv(pert_info_dir, sep="\t")

    duplicate_list: list[str] = [
        name for name, count in Counter(pert_info.pert_iname).items() if count > 1
    ]
    print(f"Number of pert_iname that have multiple pert_ids: {len(duplicate_list)}")

    return duplicate_list


def mapping_id_iname(pert_info_dir: str) -> dict:
    """Finding a mapping (dictionary) from pert_id to pert_iname.

    Parameters
    ----------
    pert_info_dir: str
      The directory of pert_info file. E.g., './Data/pert_info.txt'

    Returns
    -------
    mapping: Dict
      A dictionary that maps pert_id to pert_iname
    """
    assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"

    pert_info = pd.read_csv(pert_info_dir, sep="\t")

    return dict(zip(pert_info.pert_id.values, pert_info.pert_iname.values, strict=True))
