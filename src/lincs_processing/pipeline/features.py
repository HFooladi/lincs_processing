"""Deterministic featurisation: nothing here is fitted to data.

Anything that learns from data (scalers, cell vocabularies, batch correction)
lives in ``model.Preprocessor`` and is fitted inside the training task.
"""

from functools import lru_cache

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]

FP_RADIUS = 2
FP_BITS = 2048


@lru_cache(maxsize=1)
def _generator():  # type: ignore[no-untyped-def]
    return rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def valid_smiles(smiles: pd.Series) -> pd.Series:
    unique = {s: mol_from_smiles(s) is not None for s in smiles.unique()}
    return smiles.map(unique).astype(bool)


def compound_key(smiles: str, inchikey: str) -> str:
    """InChIKey connectivity block: the same molecule under different BRD ids,
    salts or stereo-annotations maps to one key, so it lands in one fold."""
    if not inchikey:
        mol = mol_from_smiles(smiles)
        inchikey = Chem.MolToInchiKey(mol) if mol is not None else ""
    return inchikey.split("-")[0] if inchikey else ""


def morgan_fingerprints(smiles: pd.Series) -> np.ndarray:
    """Morgan bit vectors, one row per entry; computed once per unique SMILES."""
    unique = smiles.unique()
    table = np.zeros((len(unique), FP_BITS), dtype=np.float32)
    for i, smi in enumerate(unique):
        mol = mol_from_smiles(smi)
        if mol is not None:
            table[i] = _generator().GetFingerprintAsNumPy(mol)
    index = pd.Series(np.arange(len(unique)), index=unique)
    return table[index.loc[smiles.to_numpy()].to_numpy()]


def condition_features(obs: pd.DataFrame) -> np.ndarray:
    """log10 dose (µM) and log2 time (h); scaled later on the training fold."""
    dose = np.log10(obs["dose_um"].to_numpy(dtype=np.float64).clip(min=1e-4))
    time = np.log2(obs["time_h"].to_numpy(dtype=np.float64).clip(min=1e-2))
    return np.stack([dose, time], axis=1).astype(np.float32)
