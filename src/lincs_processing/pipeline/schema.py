"""Canonical metadata layout shared by every pipeline step.

The three supported releases name the same things differently. Ingest maps
them onto ``OBS_COLUMNS`` so nothing downstream needs to know which release a
dataset came from.
"""

from dataclasses import dataclass

RELEASES = ("beta", "gse92742", "gse70138")

# One row per level-5 signature (consensus of replicates, change from control).
OBS_COLUMNS = [
    "sig_id",
    "pert_id",
    "pert_iname",
    "cell_id",
    "dose_um",
    "time_h",
    "cc_q75",  # replicate agreement: 75th pct of pairwise replicate correlation
    "ss",  # signature strength (scale differs between releases)
    "tas",  # transcriptional activity score, comparable across releases
    "n_rep",
    "smiles",
    "inchikey",
    "moa",
    "target",
]

MISSING = "-666"


@dataclass(frozen=True)
class ReleaseLayout:
    """Column names of one release's metadata tables."""

    gene_id: str
    gene_symbol: str
    landmark_col: str
    landmark_value: str
    cell_col: str
    iname_col: str
    cc_col: str
    ss_col: str
    nrep_col: str
    # beta and GSE92742 put everything in one table; GEO splits metrics out.
    metrics_in_sig_info: bool


LAYOUTS = {
    "beta": ReleaseLayout(
        gene_id="gene_id",
        gene_symbol="gene_symbol",
        landmark_col="feature_space",
        landmark_value="landmark",
        cell_col="cell_iname",
        iname_col="cmap_name",
        cc_col="cc_q75",
        ss_col="ss_ngene",
        nrep_col="nsample",
        metrics_in_sig_info=True,
    ),
    "gse92742": ReleaseLayout(
        gene_id="pr_gene_id",
        gene_symbol="pr_gene_symbol",
        landmark_col="pr_is_lm",
        landmark_value="1",
        cell_col="cell_id",
        iname_col="pert_iname",
        cc_col="distil_cc_q75",
        ss_col="distil_ss",
        nrep_col="distil_nsample",
        metrics_in_sig_info=False,
    ),
}
LAYOUTS["gse70138"] = LAYOUTS["gse92742"]
