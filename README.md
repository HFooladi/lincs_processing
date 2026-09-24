# LINCS processing

Helpers for filtering and reshaping data from the
[Library of Integrated Network-Based Cellular Signatures (LINCS)](https://lincsproject.org/)
L1000 project. You can turn the raw `.gctx` matrices plus their metadata into a
compact list of records, then slice that list by cell line, compound, dose,
time, or Repurposing Hub annotations.

## Install

The project is managed with [uv](https://docs.astral.sh/uv/) and needs
Python 3.10 or newer.

```bash
uv sync --all-groups   # creates .venv with runtime + dev dependencies
uv run pytest          # run the test suite
```

Or install it into any environment with `pip install .`.

`cmapPy` (used to read `.gctx` files) has not had a release since 2019 and does
not work with pandas 3, so pandas is pinned to `<3`.

## Data files

Nothing is bundled. Download what you need from:

- GEO [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)
  (Phase I) and [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)
  (Phase II): the level 3/5 `.gctx` matrices and the `inst_info`, `sig_info`,
  `gene_info`, and `pert_info` tables.
- The [Drug Repurposing Hub](https://clue.io/repurposing) `repurposing_drugs_*.txt`
  annotation file for mechanism of action and target information.

## Record format

Every filtering function works on a list of records. Each record is a
two-element list:

```python
record[0] == (cell_line, pert_id, pert_type, dose, dose_unit, time, time_unit)
record[1] == np.ndarray  # 978 landmark genes, or 12328 with landmarks=False
```

`parse_list_v2` uses an extended eleven-field tuple that also carries
`touchstone, clinical_phase, moa, target`.

## Usage

Build the record list from a level 3 `.gctx` file (this needs `cmapPy` and a
lot of memory for the full matrix):

```python
from lincs_processing.gctx_parser import parsing_level3_cp
from lincs_processing import write_pickle

records = parsing_level3_cp(
    "Data/Level3_INF_mlr12k_n1319138x12328.gctx",
    "Data/inst_info.txt",
    "Data/gene_info.txt",
    pert_type="trt_cp",
    landmarks=True,
)
write_pickle("Data/level3_trt_cp_landmark.pkl", records)
```

Filter records from Python. Every function accepts either a pickle path or an
already-loaded list:

```python
from lincs_processing import parse_list, parse_most_frequent, to_dataframe

mcf7 = parse_list("Data/level3_trt_cp_landmark.pkl", indicator=0, query=["MCF7"])
top3_compounds = parse_most_frequent(mcf7, indicator=1, n=3)
df = to_dataframe(top3_compounds)  # genes first, then cell_lines/compounds/doses/times
```

`indicator` selects the field: `0` cell line, `1` compound, `2` dose, `3` time.

Or filter from the command line:

```bash
uv run lincs-filter --dataset_dir Data/level3_trt_cp_landmark.pkl \
    --cells MCF7 PC3 --times 24 --output_dir Data/after_parsing.pkl
```

Metadata helpers:

```python
from lincs_processing import (
    pert_touchstone,
    duplicate_pert_name,
    drug_pert_retrieval,
    sig_info_augment,
)

by_id, by_name = pert_touchstone("Data/pert_info.txt")
dupes = duplicate_pert_name("Data/pert_info.txt")
supp, names = drug_pert_retrieval(
    "Data/repurposing_drugs_20180907.txt", "Data/pert_info.txt"
)
# adds split pert_dose / pert_time columns to a GSE70138-style sig_info
sig_info = sig_info_augment("Data/GSE70138_sig_info.txt")
```

## Development

```bash
uv run ruff check .          # lint
uv run ruff format .         # format
uv run mypy src              # type check
uv run pytest -q             # tests (set LINCS_DATA_DIR to run the slow gctx test)
```

CI runs the same commands on Python 3.11 and 3.13.
