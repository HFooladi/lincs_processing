"""Command-line filter for pickled LINCS record lists.

Example::

    lincs-filter --dataset_dir Data/level3_trt_cp_landmark.pkl \\
        --cells MCF7 PC3 --times 24 --output_dir Data/after_parsing.pkl
"""

import argparse
from collections.abc import Sequence

from lincs_processing.utils import load_pickle, write_pickle

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parsing LINCS")
    parser.add_argument(
        "--dataset_dir", type=str, default="Data/level3_trt_cp_landmark.pkl"
    )
    parser.add_argument("--cells", type=str, nargs="+", default=None)
    parser.add_argument("--compounds", type=str, nargs="+", default=None)
    parser.add_argument("--doses", type=float, nargs="+", default=None)
    parser.add_argument("--times", type=int, nargs="+", default=None)
    parser.add_argument("--output_dir", type=str, default="Data/after_parsing.pkl")
    return parser


def filter_records(
    train: list,
    cells: Sequence[str] | None = None,
    compounds: Sequence[str] | None = None,
    doses: Sequence[float] | None = None,
    times: Sequence[int] | None = None,
) -> list:
    """Apply each non-``None`` criterion in turn and return the surviving records."""
    if cells is not None:
        train = [line for line in train if line[0][0] in cells]
        print(
            f"Number of training data after parsing based on cell lines: {len(train)}"
        )

    if compounds is not None:
        train = [line for line in train if line[0][1] in compounds]
        print(f"Number of training data after parsing based on compounds: {len(train)}")

    if doses is not None:
        train = [line for line in train if line[0][3] in doses]
        print(f"Number of training data after parsing based on doses: {len(train)}")

    if times is not None:
        train = [line for line in train if line[0][5] in times]
        print(f"Number of training data after parsing based on times: {len(train)}")

    return train


def main(argv: Sequence[str] | None = None) -> None:
    flags = build_parser().parse_args(argv)

    print("=================================================================")
    print("Data Loading..")
    train = load_pickle(flags.dataset_dir)
    print(f"Number of Train Data: {len(train)}")

    train = filter_records(
        train,
        cells=flags.cells,
        compounds=flags.compounds,
        doses=flags.doses,
        times=flags.times,
    )

    print(f"Number of final training data after parsing: {len(train)}")
    write_pickle(flags.output_dir, train)


if __name__ == "__main__":
    main()
