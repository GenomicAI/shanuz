from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import is_sparse


def calc_n(matrix, margin: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Compute nCount (column sums) and nFeature (non-zero count per column).

    margin=2 operates on columns (cells), matching R's default.
    Returns (nCount, nFeature) as 1-D float arrays.
    """
    if sp.issparse(matrix):
        ncount = np.asarray(matrix.sum(axis=0)).flatten()
        nfeature = np.diff(matrix.tocsc().indptr).astype(float)
    else:
        arr = np.asarray(matrix)
        ncount = arr.sum(axis=0).flatten()
        nfeature = (arr != 0).sum(axis=0).flatten()
    return ncount.astype(float), nfeature.astype(float)


def match_cells(
    query: list[str],
    target: list[str],
    error_missing: bool = False,
) -> np.ndarray:
    """Return integer indices of query items in target (like R match()).

    Returns -1 for items not found. Raises if error_missing=True and any are absent.
    """
    target_index = {v: i for i, v in enumerate(target)}
    idx = np.array([target_index.get(q, -1) for q in query], dtype=int)
    if error_missing and (idx == -1).any():
        missing = [q for q, i in zip(query, idx) if i == -1]
        raise KeyError(f"Cells not found in target: {missing}")
    return idx


def intersect_names(a: list[str], b: list[str]) -> list[str]:
    b_set = set(b)
    return [x for x in a if x in b_set]


def unique_names(names: list[str]) -> bool:
    return len(names) == len(set(names))


def validate_cell_names(names) -> list[str]:
    names = list(names)
    if not unique_names(names):
        raise ValueError("Cell names must be unique.")
    return names


def validate_feature_names(names) -> list[str]:
    names = list(names)
    if not unique_names(names):
        raise ValueError("Feature names must be unique.")
    return names
