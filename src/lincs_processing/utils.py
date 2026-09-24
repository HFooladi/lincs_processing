"""Filtering and reshaping helpers for pickled LINCS record lists.

Every function accepts ``data`` either as a path to a pickle file or as an
already-loaded list of records. A record is a two-element list/tuple::

    record[0]: (cell_line, pert_id, pert_type, dose, dose_unit, time, time_unit)
    record[1]: numpy array with 978 (landmark) or 12328 gene values

``parse_list_v2`` expects the extended eleven-field metadata tuple that
additionally carries touchstone, clinical phase, MOA and target.
"""

import pickle
import warnings
from collections import Counter

import pandas as pd
from tqdm import tqdm

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"

# Position of each filterable field inside record[0].
_FIELD_INDEX = {0: 0, 1: 1, 2: 3, 3: 5, 4: 7, 5: 8, 6: 9, 7: 10}
_FIELD_NAME = {
    0: "cell_lines",
    1: "compounds",
    2: "doses",
    3: "time",
    4: "touchstone",
    5: "clinical_phase",
    6: "moa",
    7: "target",
}


def _load(data: str | list) -> list:
    """Return ``data`` as a list, reading it from a pickle file if it is a path."""
    assert isinstance(data, (str, list)), "The data should be string or list object"
    if isinstance(data, str):
        with open(data, "rb") as f:
            return pickle.load(f)
    return data


def _print_loading() -> None:
    print("=================================================================")
    print("Data Loading..")


def load_pickle(dataset_dir: str) -> list:
    """Loading (reading) a pickle file

    Parameters
    ----------
    dataset_dir: str
      It must be string file that shows the directory of the dataset.

    Returns
    -------
    List
    """
    assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"

    with open(dataset_dir, "rb") as fp:
        return pickle.load(fp)


def write_pickle(dataset_dir: str, data: list) -> None:
    """Writing a file (data) into a pickle file (dataset_dir)

    Parameters
    ----------
    dataset_dir: str
      It must be string file that shows the directory for writing.
    data: List
      the object that should be written into the pickle file.
    """
    assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"

    with open(dataset_dir, "wb") as fp:
        pickle.dump(data, fp)


def print_statistics(data: str | list) -> None:
    """Print data statistics

    This function takes the directory of dataset and
    returns some useful statistics about the data.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
    """
    _print_loading()
    train = _load(data)

    print("Data Statistics\n")
    print(f"Number of Train Data: {len(train)}")

    print("Please wait while we are retriving information ...")
    cell_lines = [line[0][0] for line in train]
    compounds = [line[0][1] for line in train]
    doses = [line[0][3] for line in train]
    times = [line[0][5] for line in train]

    print(f"Number of unique Cell Lines: {len(set(cell_lines))}")
    print(f"Number of unique Compounds: {len(set(compounds))}")
    print(f"Number of unique doses: {len(set(doses))}")
    print(f"Number of unique times: {len(set(times))}")


def print_most_frequent(data: str | list, n: int = 3) -> None:
    """Print most frequent cell line, compounds, and does.

    This function takes the directory of dataset (or a list object) and integer n
    and returns The n most frequent cell lines, compounds and
    doses in the dataset.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
    n: int, optional (default 3)
      An integer which determine number of frequent statistics we want
      to retrieve. Default=3.
    """
    _print_loading()

    assert isinstance(n, int), "The parameter n must be an integer"
    train = _load(data)

    print("Please wait while we are retriving information ...")
    cell_lines = []
    compounds = []
    doses = []
    for line in tqdm(train):
        cell_lines.append(line[0][0])
        compounds.append(line[0][1])
        doses.append(line[0][3])

    print("loop finished !!!")

    print(f"Most frequent Cell Lines: {Counter(cell_lines).most_common(n)}")
    print(f"Most frequent Compounds: {Counter(compounds).most_common(n)}")
    print(f"Most frequent Doses: {Counter(doses).most_common(n)}")


def cell_line_frequent(data: str | list, n: int = 3) -> list:
    """Returns list of data belongs to most frequent cell lines

    This function takes the directory of dataset (or a list object) and integer n,
    and parse the data to keep only the data that belongs to n
    most frequent cell lines.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    n: int, optional (default 3)
      An integer which determine number of frequent statistics we want
      to retrieve. Default=3.

    Returns
    -------
    parse_data: List
      A list containing data that belongs to n most frequent cell lines.
    """
    _print_loading()

    assert isinstance(n, int), "The parameter n must be an integer"
    train = _load(data)

    print("Please wait while we are retriving information ...")
    cell_lines = [line[0][0] for line in tqdm(train)]

    print(f"Number of unique Cell Lines: {len(set(cell_lines))}")
    print(f"Most frequent Cell Lines: {Counter(cell_lines).most_common(n)}")

    if n > len(set(cell_lines)):
        warnings.warn(
            "n is greater than number of unique cell lines available in the dataset",
            stacklevel=2,
        )

    # List of n most frequent cell lines
    x = [item for item, _ in Counter(cell_lines).most_common(n)]

    return [line for line in train if line[0][0] in x]


def cell_line_list(data: str | list, cells: list[str] | None = None) -> list:
    """Filter data based on desired cell line list

    This function takes the directory of dataset (or alist object) and a list cells,
    and parse the data to keep only the data that belongs to cells list.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    cells: List[str]
      list of cell lines that we want to keep their data to retrieve.
      Default=['MCF7']

    Returns
    -------
    parse_data: list
      A list containing data that belongs to desired list.
    """
    if cells is None:
        cells = ["MCF7"]
    assert isinstance(cells, list), "The parameter cells must be a list"

    _print_loading()
    train = _load(data)

    print(f"Number of Train Data: {len(train)}")

    parse_data = [line for line in train if line[0][0] in cells]

    print(f"Number of Data after parsing: {len(parse_data)}")
    return parse_data


def parse_list(data: str | list, indicator: int = 0, query: list | None = None) -> list:
    """Filter the data based on compound, cell line, dose or time

    This function takes the directory of dataset, indicator that indicates
    whether you want to subset the data based on cell line, compound, dose, or time
    and a list which shows what part of the data you want to keep.
    The output will be a list of desired parsed dataset.


    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    indicator: int
      it must be an integer from 0 1 2 and 3 that shows whether
      we want to retrieve the data based on cells, compound or dose.
      0: cell_lines
      1:compounds
      2:doses
      3:time
      Default=0 (cell_lines)

    query: List
      list of cells or compounds or doses that we want to retrieve.
      The list depends on the indicator. If the indicator is 0, you should enter
      the list of desired cell lines and so on. Default=['MCF7']

    Returns
    -------
    parse_data: List
      A list containing data that belongs to desired list.
    """
    if query is None:
        query = ["MCF7"]
    assert isinstance(indicator, int), "The indicator must be an int object"
    assert indicator in [0, 1, 2, 3], "You should choose indicator from 0, 1, 2, 3"
    assert isinstance(query, list), "The parameter query must be a list"

    _print_loading()
    train = _load(data)

    k = _FIELD_INDEX[indicator]

    print(f"Number of Train Data: {len(train)}")
    print(f"You are parsing the data base on {_FIELD_NAME[indicator]}")

    parse_data = [line for line in train if line[0][k] in query]

    print(f"Number of Data after parsing: {len(parse_data)}")
    return parse_data


def parse_most_frequent(data: str | list, indicator: int = 0, n: int = 3) -> list:
    """Returns most frequent data (based on cell line, compound, ...)

    This function takes the directory of dataset, indicator that indicates
    whether you want to subset the data based on cell line, compound, dose, or time
    and a n which how much frequent items you want to keep.
    The output will be a list of desired parsed dataset.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    indicator: int, optional (default n=0)
      It must be an integer from 0 1 2 and 3 that shows whether
      we want to retrieve the data based on cells, compound or dose.
      0: cell_lines
      1:compounds
      2:doses
      3:time
      Default=0

    n: int, optional (default n=3)
      number of most frequent cells or compounds or doses that we want to retrieve.
      The list depends on the indicator. If the indicator is 0, you should enter
      the number of desired cell lines and so on. Default=3

    Returns
    -------
    parse_data: List
      A list containing data that belongs to desired list.
    """
    assert isinstance(indicator, int), "The indicator must be an int object"
    assert indicator in [0, 1, 2, 3], "You should choose indicator from 0, 1, 2, 3"
    assert isinstance(n, int), "The parameter n must be an integer"

    _print_loading()
    train = _load(data)

    k = _FIELD_INDEX[indicator]
    name = _FIELD_NAME[indicator]

    mylist = [line[0][k] for line in tqdm(train)]

    print(f"Number of unique {name}: {len(set(mylist))}")
    print(f"Most frequent {name}: {Counter(mylist).most_common(n)}")

    assert n <= len(set(mylist)), "n is out of valid range!"

    y = [item for item, _ in Counter(mylist).most_common(n)]

    return [line for line in train if line[0][k] in y]


def parse_chunk_frequent(
    data: str | list, indicator: int = 0, start: int = 0, end: int = 3
) -> list:
    """Keep the records whose field ranks between ``start`` and ``end`` by frequency.

    This function takes the directory of dataset, indicator that indicates
    whether you want to subset the data based on cell line, compound, dose, or time
    and a start and end which shows what chunk of data is desirable.
    E.g., if start=0 and end=3, you are subsetting 3 most frequent data.
    The output will be a list of desired parsed dataset.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    indicator: int, optional (default n=0)
      It must be an integer from 0 1 2 and 3 that shows whether
      we want to retrieve the data based on cells, compound or dose.
      0: cell_lines
      1:compounds
      2:doses
      3:time
      Default=0

    start: int
      indicates the start of the list you want to subset. Default=0
    end: int
      indicates the end of the list you want to subset (exclusive, like a
      Python slice). Default=3

    Returns
    -------
    parse_data: List
      A list containing data that belongs to desired list.
    """
    assert isinstance(indicator, int), "The indicator must be an int object"
    assert indicator in [0, 1, 2, 3], "You should choose indicator from 0, 1, 2, 3"
    assert isinstance(start, int), "The parameter start must be an integer"
    assert isinstance(end, int), "The parameter end must be an integer"
    assert start <= end, "The start should be less than the end!!"

    _print_loading()
    train = _load(data)

    k = _FIELD_INDEX[indicator]
    name = _FIELD_NAME[indicator]

    mylist = [line[0][k] for line in train]

    print(f"Number of unique {name}: {len(set(mylist))}")

    assert end <= len(set(mylist)), "end is out of valid range!"

    y = [item for item, _ in Counter(mylist).most_common()][start:end]

    print(f"Desired {name}: {y}")

    return [line for line in train if line[0][k] in y]


def parse_dose_range(
    data: str | list, dose_min: float = 0, dose_max: float = 5
) -> list:
    """Keep the records whose dose lies strictly between ``dose_min`` and ``dose_max``.

    This function takes the directory of dataset minimum and maximum dose
    and return a list of data that are within the desired range.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    dose_min: float, optional (default dose_min=0)
      minimum dose (exclusive). Default=0
    dose_max: float, optional (default dose_max=5)
      maximum dose (exclusive). Default=5

    Returns
    --------
    parse_data: List
      A list containing data that belongs to desired list (
      Desired range of doses).
    """
    assert isinstance(dose_min, (int, float)), "dose_min must be a number"
    assert isinstance(dose_max, (int, float)), "dose_max must be a number"
    assert dose_min < dose_max, "The minimum dose must be less than the maximum dose !!"

    _print_loading()
    train = _load(data)

    print(f"Number of Train Data: {len(train)}")

    parse_data = [line for line in train if dose_min < line[0][3] < dose_max]

    print(f"Number of Data after parsing: {len(parse_data)}")
    return parse_data


def to_dataframe(data: str | list) -> pd.DataFrame:
    """This takes a list and produce a pandas datframe of data

    The input to this function is a list which contains metadata
    (such as cell lines, compounds, ..) and gene expression. this
    function returns a pandas dataframe where the first columns
    belongs to gene expression and last four columns contain metaddata
    cell line, compound, dose, and time in this order.

    Parameters
    ----------
    data: Union[str, List]
      the data can be a string which is the directory of the dataset.
      dataset should be a pickle file. e.g., valid argument is something like this:
      './Data/level3_trt_cp_landmark.pkl'
      or it can be a list which contains the gene expression and metadata.
      It must be a list of tuples with the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

    Returns
    -------
    pd.DataFrame
      This is a pandas dataframe where the first columns contains
      gene expression (978 or 12328-dimension) and the last four columns
      contains cell line, pert_id, dose, and time
    """
    train = _load(data)

    genes = [line[1] for line in train]
    metadata = {
        "cell_lines": [line[0][0] for line in train],
        "compounds": [line[0][1] for line in train],
        "doses": [line[0][3] for line in train],
        "times": [line[0][5] for line in train],
    }

    return pd.concat([pd.DataFrame(genes), pd.DataFrame(metadata)], axis=1)


def parse_list_v2(
    dataset_dir: str | None = None,
    indicator: int = 0,
    query: list | None = None,
    data: list | None = None,
) -> list:
    """Filter extended records by cell line, compound, dose, time, touchstone, ...

    This function takes the directory of dataset, indicator that indicates
    whether you want to subset the data based on cell line, compound, dose, time,
    touchstone, clinical phase, MOA or target. Moreover, it takes a list which
    shows what part of the data you want to keep.
    The output will be a list of desired parsed dataset.

    Parameters
    ----------
    dataset_dir: str, optional
      Path of the pickle file, e.g. './Data/level3_trt_cp_landmark_allinfo.pkl'.
      Ignored when ``data`` is given. The pickle file should be a list with
      the following format:
      line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type,
      touchstone, clinical phase, moa, target)
      line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
    indicator: int
      it must be an integer from 0 1 2 3 4 5 6 7 that shows whether
      we want to retrieve the data based on cells, compound, dose, touchstone,
      clinical phase, moa or target.
      0: cell_lines
      1: compounds
      2: doses
      3: time
      4: touchstone
      5: clinical phase
      6: moa
      7: target
      Default=0 (cell_lines)
    query: list
      list of cells or compounds or doses or time or touchstone or clinical
      phase or MOA or target that we want to retrieve. The list depends on the
      indicator. If the indicator is 0, you should enter the list of desired
      cell lines and so on. Default=['MCF7']
    data: list, optional
      Already-loaded records in the format described above. When provided,
      ``dataset_dir`` is ignored.

    Returns
    -------
    parse_data: list
      A list containing data that belongs to desired list.
    """
    if query is None:
        query = ["MCF7"]
    assert isinstance(indicator, int), "The indicator must be an int object"
    assert indicator in range(8), "You should choose indicator from 0 to 7"
    assert isinstance(query, list), "The parameter query must be a list"

    _print_loading()
    if data is None:
        assert isinstance(dataset_dir, str), (
            "The dataset_dir must be a string object when data is not given"
        )
        train = _load(dataset_dir)
    else:
        assert isinstance(data, list), "The data must be a list object"
        train = data

    k = _FIELD_INDEX[indicator]

    print(f"Number of Train Data: {len(train)}")
    print(f"You are parsing the data base on {_FIELD_NAME[indicator]}")

    parse_data = []
    if indicator in [0, 1, 2, 3, 4]:
        parse_data = [line for line in train if line[0][k] in query]
    else:
        # clinical phase, moa and target are "|"-separated multi-valued fields
        for line in train:
            values = line[0][k][0].split("|")
            if any(value in query for value in values):
                parse_data.append(line)

    print(f"Number of Data after parsing: {len(parse_data)}")
    return parse_data
