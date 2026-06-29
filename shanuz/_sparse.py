from __future__ import annotations

from typing import Union

import numpy as np
import scipy.sparse as sp


def is_sparse(x) -> bool:
    return sp.issparse(x)


def as_sparse(x, format: str = "csc") -> sp.spmatrix:
    """Convert a dense matrix or another sparse format to a scipy sparse matrix."""
    if sp.issparse(x):
        return x.asformat(format)
    arr = np.asarray(x)
    converter = {"csc": sp.csc_matrix, "csr": sp.csr_matrix, "coo": sp.coo_matrix}
    return converter.get(format, sp.csc_matrix)(arr)


def as_dense(x) -> np.ndarray:
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def is_matrix_empty(x) -> bool:
    """Mirrors R IsMatrixEmpty — True when any dimension is 0."""
    if x is None:
        return True
    if sp.issparse(x):
        return x.shape[0] == 0 or x.shape[1] == 0
    arr = np.asarray(x)
    return arr.size == 0 or 0 in arr.shape


def empty_sparse(nrow: int = 0, ncol: int = 0) -> sp.csc_matrix:
    return sp.csc_matrix((nrow, ncol), dtype=np.float64)


def empty_dense(nrow: int = 0, ncol: int = 0) -> np.ndarray:
    return np.empty((nrow, ncol), dtype=np.float64)


def check_matrix(x, rows: list[str] | None = None, cols: list[str] | None = None) -> None:
    """Validate matrix dimensions against expected row/col names."""
    if rows is not None and x.shape[0] != len(rows):
        raise ValueError(
            f"Matrix has {x.shape[0]} rows but {len(rows)} row names were supplied."
        )
    if cols is not None and x.shape[1] != len(cols):
        raise ValueError(
            f"Matrix has {x.shape[1]} columns but {len(cols)} col names were supplied."
        )


def stitch_matrix(
    matrices: list,
    row_names: list[str],
    col_names: list[str],
) -> sp.csc_matrix:
    """Combine a list of sub-matrices into one sparse matrix using given index mappings.

    Mirrors R StitchMatrix — used when merging assays across samples.
    """
    blocks = [sp.csc_matrix(m) if not sp.issparse(m) else m for m in matrices]
    return sp.hstack(blocks, format="csc")
