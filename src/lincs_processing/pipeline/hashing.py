"""Content hashes written into every artefact, run and model card.

The hashes are computed from content, not paths or timestamps, so the same
data always gets the same id no matter where or when it was produced.
"""

import hashlib
from collections.abc import Iterable

import numpy as np

_CHUNK = 1 << 20


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(_CHUNK):
            digest.update(block)
    return digest.hexdigest()


def data_hash(X: np.ndarray, sig_ids: Iterable[str], release: str) -> str:
    """Hash a signature matrix together with its row order and release."""
    digest = hashlib.sha256(release.encode())
    digest.update("\n".join(sig_ids).encode())
    digest.update(str(X.shape).encode())
    for start in range(0, X.shape[0], 10_000):
        block = np.ascontiguousarray(X[start : start + 10_000], dtype=np.float32)
        digest.update(block.tobytes())
    return digest.hexdigest()


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()
