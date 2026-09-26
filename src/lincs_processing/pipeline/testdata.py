"""A tiny synthetic LINCS release for tests and ``-profile test``.

Signatures are generated from a known function of compound structure, cell,
dose and time plus noise, so a working pipeline has something real to learn.
The fixture deliberately contains the traps the pipeline must handle: one
molecule under two BRD ids, a compound with an unparseable SMILES, a compound
with two targets (two rows in compoundinfo), low-quality signatures, and
vehicle-control signatures that ingest must drop.
"""

import os

import h5py
import numpy as np
import pandas as pd
from rdkit import Chem

from lincs_processing.pipeline.features import morgan_fingerprints

COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "vorinostat": "ONC(=O)CCCCCCC(=O)NC1=CC=CC=C1",
    "imatinib": (
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
    ),
    "gefitinib": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
    "erlotinib": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
    "dasatinib": "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "metformin": "CN(C)C(=N)N=C(N)N",
    "tamoxifen": "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3",
    "bortezomib": "CC(C)CC(B(O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)C2=NC=CN=C2",
    "sorafenib": "CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F",
    "trichostatin-a": "CC(C=C(C)C=CC(=O)NO)C(=O)C1=CC=C(C=C1)N(C)C",
    "methotrexate": (
        "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(C(=O)O)CCC(=O)O"
    ),
    "camptothecin": "CCC1(C2=C(COC1=O)C(=O)N3CC4=CC5=CC=CC=C5N=C4C3=C2)O",
    "caffeine-dup": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # same molecule, second BRD id
    "broken": "C1CC(",  # unparseable: QC must drop it
}
CELLS = ["A375", "A549", "HA1E", "HT29", "MCF7", "PC3", "VCAP"]
DOSES_UM = [0.1, 1.0, 10.0]
TIMES_H = [6, 24]
N_GENES = 120
N_LANDMARK = 80


def _frames(seed: int = 7) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    names = list(COMPOUNDS)
    pert_ids = {n: f"BRD-K{1000 + i:05d}" for i, n in enumerate(names)}
    smiles = pd.Series([COMPOUNDS[n] for n in names], index=names)
    fps = morgan_fingerprints(smiles)

    genes = [str(10_000 + i) for i in range(N_GENES)]
    W = rng.normal(size=(fps.shape[1], N_GENES)) * (rng.random((fps.shape[1], 1)) < 0.2)
    raw = np.tanh(fps @ W / np.sqrt(fps.sum(1, keepdims=True).clip(min=1))) * 2
    effect = pd.DataFrame(raw, index=names)
    cell_offset = {c: rng.normal(scale=0.4, size=N_GENES) for c in CELLS}
    cell_gain = {c: rng.uniform(0.7, 1.3) for c in CELLS}

    rows, X = [], []
    for name in names:
        for cell in CELLS:
            for dose in DOSES_UM:
                for time in TIMES_H:
                    scale = (np.log10(dose) + 2) / 3 * (1.0 if time == 24 else 0.6)
                    X.append(
                        cell_gain[cell] * scale * effect.loc[name].to_numpy()
                        + cell_offset[cell]
                        + rng.normal(scale=0.3, size=N_GENES)
                    )
                    rows.append(
                        {
                            "sig_id": f"MINI_{cell}_{pert_ids[name]}_{dose}um_{time}h",
                            "pert_id": pert_ids[name],
                            "pert_iname": name,
                            "pert_type": "trt_cp",
                            "cell_id": cell,
                            "dose": dose,
                            "time": time,
                            "cc_q75": rng.uniform(0.1, 0.9),
                            "ss": rng.uniform(1, 8),
                            "tas": rng.uniform(0.1, 0.7),
                            "n_rep": int(rng.choice([1, 2, 3, 3, 3])),
                        }
                    )
    for cell in CELLS:  # vehicle controls live in the same matrix
        X.append(rng.normal(scale=0.1, size=N_GENES))
        rows.append(
            {
                "sig_id": f"MINI_{cell}_DMSO_24h",
                "pert_id": "DMSO",
                "pert_iname": "DMSO",
                "pert_type": "ctl_vehicle",
                "cell_id": cell,
                "dose": np.nan,
                "time": 24,
                "cc_q75": 0.5,
                "ss": 1.0,
                "tas": 0.1,
                "n_rep": 3,
            }
        )
    sigs = pd.DataFrame(rows)
    matrix = pd.DataFrame(np.array(X, dtype=np.float32).T, index=genes)
    matrix.columns = sigs["sig_id"]

    compounds = []
    for name in names:
        mol = Chem.MolFromSmiles(COMPOUNDS[name])
        compounds.append(
            {
                "pert_id": pert_ids[name],
                "pert_iname": name,
                "smiles": COMPOUNDS[name],
                "inchikey": Chem.MolToInchiKey(mol) if mol is not None else "",
                "moa": "HDAC inhibitor" if "vorinostat" in name else "",
                "target": "HDAC1" if "vorinostat" in name else "",
            }
        )
    # Two targets -> two rows in compoundinfo_beta, as in the real file.
    compounds.append({**compounds[1], "target": "HDAC2"})
    gene_table = pd.DataFrame(
        {
            "gene_id": genes,
            "gene_symbol": [f"G{i}" for i in range(N_GENES)],
            "is_landmark": [i < N_LANDMARK for i in range(N_GENES)],
        }
    )
    return {
        "sigs": sigs,
        "matrix": matrix,
        "compounds": pd.DataFrame(compounds),
        "genes": gene_table,
    }


def _write_gctx(matrix: pd.DataFrame, path: str) -> None:
    """Minimal GCTX 1.0 writer (cmapPy's writer uses the removed ``np.string_``).

    Layout as read by ``cmapPy.pandasGEXpress.parse``: the matrix is stored
    samples x genes under ``/0/DATA/0/matrix`` with ids under ``/0/META``.
    """
    with h5py.File(path, "w") as h5:
        h5.attrs["version"] = np.bytes_("GCTX1.0")
        h5.create_dataset(
            "0/DATA/0/matrix", data=matrix.to_numpy(np.float32).T, compression="gzip"
        )
        for node, ids in (("ROW", matrix.index), ("COL", matrix.columns)):
            h5.create_dataset(f"0/META/{node}/id", data=np.array(ids, dtype="S"))


def _tsv(frame: pd.DataFrame, path: str) -> None:
    frame.to_csv(path, sep="\t", index=False)


def write_beta(outdir: str, frames: dict) -> dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    s, c, g = frames["sigs"], frames["compounds"], frames["genes"]
    paths = {
        "gctx": os.path.join(outdir, "level5_beta_trt_cp_mini.gctx"),
        "sig_info": os.path.join(outdir, "siginfo_beta.txt"),
        "pert_info": os.path.join(outdir, "compoundinfo_beta.txt"),
        "gene_info": os.path.join(outdir, "geneinfo_beta.txt"),
    }
    _write_gctx(frames["matrix"].copy(), paths["gctx"])
    _tsv(
        pd.DataFrame(
            {
                "sig_id": s.sig_id,
                "pert_id": s.pert_id,
                "cmap_name": s.pert_iname,
                "pert_type": s.pert_type,
                "cell_iname": s.cell_id,
                "pert_dose": s.dose,
                "pert_dose_unit": np.where(s.dose.notna(), "uM", ""),
                "pert_time": s.time,
                "pert_time_unit": "h",
                "cc_q75": s.cc_q75,
                "ss_ngene": (s.ss * 10).round(),
                "tas": s.tas,
                "nsample": s.n_rep,
            }
        ),
        paths["sig_info"],
    )
    _tsv(
        c.rename(
            columns={
                "pert_iname": "cmap_name",
                "smiles": "canonical_smiles",
                "inchikey": "inchi_key",
            }
        ),
        paths["pert_info"],
    )
    _tsv(
        pd.DataFrame(
            {
                "gene_id": g.gene_id,
                "gene_symbol": g.gene_symbol,
                "feature_space": np.where(g.is_landmark, "landmark", "inferred"),
            }
        ),
        paths["gene_info"],
    )
    return paths


def write_geo(outdir: str, frames: dict, release: str) -> dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    s, c, g = frames["sigs"], frames["compounds"], frames["genes"]
    paths = {
        "gctx": os.path.join(outdir, "level5_mini.gctx"),
        "sig_info": os.path.join(outdir, "sig_info.txt"),
        "sig_metrics": os.path.join(outdir, "sig_metrics.txt"),
        "pert_info": os.path.join(outdir, "pert_info.txt"),
        "gene_info": os.path.join(outdir, "gene_info.txt"),
    }
    _write_gctx(frames["matrix"].copy(), paths["gctx"])
    idose = [f"{d} um" if pd.notna(d) else "-666" for d in s.dose]
    itime = [f"{t} h" for t in s.time]
    base = {
        "sig_id": s.sig_id,
        "pert_id": s.pert_id,
        "pert_iname": s.pert_iname,
        "pert_type": s.pert_type,
        "cell_id": s.cell_id,
    }
    if release == "gse92742":
        sig_info = pd.DataFrame(
            {
                **base,
                "pert_dose": s.dose.fillna(-666),
                "pert_dose_unit": np.where(s.dose.notna(), "µM", "-666"),
                "pert_idose": idose,
                "pert_time": s.time,
                "pert_time_unit": "h",
                "pert_itime": itime,
            }
        )
    else:  # GSE70138 only has the combined "10 um" / "24 h" columns
        sig_info = pd.DataFrame({**base, "pert_idose": idose, "pert_itime": itime})
    _tsv(sig_info, paths["sig_info"])
    _tsv(
        pd.DataFrame(
            {
                "sig_id": s.sig_id,
                "pert_id": s.pert_id,
                "distil_cc_q75": s.cc_q75,
                "distil_ss": s.ss,
                "tas": s.tas,
                "distil_nsample": s.n_rep,
            }
        ),
        paths["sig_metrics"],
    )
    pert = c.drop_duplicates("pert_id")
    _tsv(
        pd.DataFrame(
            {
                "pert_id": pert.pert_id,
                "pert_iname": pert.pert_iname,
                "pert_type": "trt_cp",
                "canonical_smiles": pert.smiles.replace("", "-666"),
                "inchi_key": pert.inchikey.replace("", "-666"),
            }
        ),
        paths["pert_info"],
    )
    _tsv(
        pd.DataFrame(
            {
                "pr_gene_id": g.gene_id,
                "pr_gene_symbol": g.gene_symbol,
                "pr_is_lm": g.is_landmark.astype(int),
            }
        ),
        paths["gene_info"],
    )
    return paths


def write_all(outdir: str, releases: tuple[str, ...] = ("beta",)) -> dict[str, dict]:
    frames = _frames()
    out = {}
    for release in releases:
        target = os.path.join(outdir, release)
        if release == "beta":
            out[release] = write_beta(target, frames)
        else:
            out[release] = write_geo(target, frames, release)
    return out
