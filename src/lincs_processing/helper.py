import pandas as pd

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def sig_info_augment(sig_info_dir: str) -> pd.DataFrame:
    """Unification between GSE70138 and GSE92742

    This function has been written for working with GSE70138.
    Unfortunately, format of sig_info (number of columns) differs between
    GSE92742 and GSE70138. So, I have written this function for unification
    between these two dataset. Particularly, I am going to add 4 columns
    (pert_dose, pert_dose_unit, pert_time, pert_time_unit) to sig_info of
    GSE70138.

    Parameters
    ----------
    sig_info_dir: str
      The directory of sig_info file. E.g., './Data/sig_info.txt'

    Returns
    -------
    sig_info_v1: pd.DataFrame
      A dataframe. It is like the input file, except it has four more columns.
    """
    assert isinstance(sig_info_dir, str), "The dataset_dir must be a string object"

    sig_info = pd.read_csv(sig_info_dir, sep="\t")

    print("Data Statistics\n")
    print(f"Number of available gene expression signature: {sig_info.shape[0]}")
    print(f"Number of available columns: {sig_info.shape[1]}")

    def split(value: str) -> list[str]:
        value = str(value)
        return value.split() if value != "-666" else ["-666", "-666"]

    dose = [split(x) for x in sig_info.pert_idose]
    time = [split(x) for x in sig_info.pert_itime]

    adding = pd.DataFrame(
        {
            "pert_dose": [float(x[0]) for x in dose],
            "pert_dose_unit": [x[1] for x in dose],
            "pert_time": [float(x[0]) for x in time],
            "pert_time_unit": [x[1] for x in time],
        }
    )

    sig_info_v1 = pd.concat([sig_info, adding], axis=1)
    print(f"Number of available columns after augmentation: {sig_info_v1.shape[1]}")

    return sig_info_v1
