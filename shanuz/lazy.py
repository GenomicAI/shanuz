"""BPCells-style lazy, on-disk matrices for out-of-core analysis (v0.8.0 Scale).

The problem this solves
=======================
A single dense ``float64`` matrix of a million cells by twenty-thousand genes is
160 GB — it never fits in memory, and even the sparse counts (a few GB) leave no
headroom for the copies every analysis step makes. Seurat's answer is the
``BPCells`` package: it keeps the matrix **on disk** in a compressed sparse
format, memory-maps it, and streams operations over it a block of cells at a time
so that peak RAM is a function of the *block* size, not the dataset size. Nothing
is loaded until it is needed, and only the slice that is needed is loaded.

``LazyMatrix`` is the Python analogue, built on the one dependency shanuz already
has — NumPy. A matrix is written to a directory as three memory-mapped arrays
(the ``data`` / ``indices`` / ``indptr`` triple of a compressed-sparse-**column**
matrix, exactly scipy's ``csc_matrix`` layout) plus a small JSON header. Opening
it maps those arrays without reading them; a slice reads only the touched
columns' non-zeros off disk and returns an ordinary in-memory ``scipy.sparse``
matrix, so every downstream call site that already accepts a sparse layer keeps
working unchanged.

Why compressed-sparse-*column*
==============================
shanuz matrices are ``features × cells``, and the operations that dominate at
scale — sketching, cell subsetting, per-cell normalisation — select or stream
over **cells**, i.e. columns. CSC stores each column contiguously, so reading an
arbitrary set of columns costs only their own non-zeros; ``col_blocks`` walks the
whole matrix in cell-blocks for a streaming reduction. Selecting a subset of
*features* (rows) from a CSC store still requires scanning the touched columns,
which is why a subset-of-both index (``m[np.ix_(genes, cells)]``) first narrows to
the cell columns and only then to the gene rows — the feature filter then runs
over a small in-memory block. (BPCells keeps both orientations on disk to make
row-only slicing cheap too; that is a natural future extension here.)

How it slots into an assay
==========================
An :class:`~shanuz.assay5.Assay5` layer may hold *any* array-like object, so a
``LazyMatrix`` can be dropped in wherever a ``counts`` / ``data`` layer would sit::

    lazy = write_lazy_matrix(assay.layers["counts"], "counts.mat")
    assay.set_layer_data("counts", lazy)

``Assay5.layer_data`` indexes the layer with ``m[np.ix_(rows, cols)]`` and gets a
concrete sparse block back; the hot paths that then call ``.toarray()`` on that
block are untouched. Code that densifies the whole layer — ``as_dense`` /
``np.asarray`` — routes through :meth:`LazyMatrix.__array__`, the deliberate
"materialise everything" escape hatch you avoid on the million-cell path but keep
for the small ones.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Iterator, Optional, Tuple, Union

import numpy as np
import scipy.sparse as sp

from ._sparse import as_sparse

_FORMAT = "shanuz-lazy-csc-v1"
_META = "meta.json"
_DATA = "data.npy"
_INDICES = "indices.npy"
_INDPTR = "indptr.npy"


class LazyMatrix:
    """A memory-mapped, on-disk compressed-sparse-column matrix.

    Instances are created by :func:`write_lazy_matrix` (persist an in-memory
    matrix) or :func:`open_lazy_matrix` (map an existing store); they are not
    constructed directly. The three CSC arrays are memory-mapped read-only, so
    the object is cheap to hold and its footprint is the slices you touch — not
    the whole matrix.

    Supports the slicing idioms the assay layer accessors use — ``m[rows, cols]``,
    ``m[np.ix_(rows, cols)]``, ``m[idx, :]``, ``m[:, idx]``, contiguous slices —
    returning a ``scipy.sparse.csc_matrix`` block. Tuple indexing is always an
    **outer** (cross-product) selection, matching ``np.ix_`` and how layers are
    block-subset throughout shanuz; element-wise pair indexing is not supported.
    """

    __slots__ = ("_path", "_shape", "_dtype", "_data", "_indices", "_indptr")

    def __init__(
        self,
        path: Union[str, Path],
        shape: Tuple[int, int],
        data: np.ndarray,
        indices: np.ndarray,
        indptr: np.ndarray,
    ) -> None:
        self._path = Path(path)
        self._shape = (int(shape[0]), int(shape[1]))
        self._data = data
        self._indices = indices
        self._indptr = indptr

    # ------------------------------------------------------------------
    # Array-like metadata
    # ------------------------------------------------------------------

    @property
    def shape(self) -> Tuple[int, int]:
        return self._shape

    @property
    def nrow(self) -> int:
        return self._shape[0]

    @property
    def ncol(self) -> int:
        return self._shape[1]

    @property
    def ndim(self) -> int:
        return 2

    @property
    def dtype(self) -> np.dtype:
        return self._data.dtype

    @property
    def nnz(self) -> int:
        return int(self._data.shape[0])

    @property
    def path(self) -> Path:
        return self._path

    def __len__(self) -> int:
        return self._shape[0]

    def __repr__(self) -> str:
        return (
            f"LazyMatrix(shape={self._shape}, nnz={self.nnz}, "
            f"dtype={self.dtype}, path='{self._path}')"
        )

    # ------------------------------------------------------------------
    # Lazy column reads — the heart of the out-of-core story
    # ------------------------------------------------------------------

    def _read_columns(self, cols: Union[slice, np.ndarray]) -> sp.csc_matrix:
        """Materialise the requested columns as an in-memory CSC block.

        Only the touched columns' ``data`` / ``indices`` ranges are read off the
        memory-mapped arrays; the rest of the matrix is never faulted in.
        """
        nrow = self.nrow
        indptr = self._indptr

        # Fast path: a contiguous column slice is one contiguous disk read.
        if isinstance(cols, slice):
            start, stop, step = cols.indices(self.ncol)
            if step == 1:
                s, e = int(indptr[start]), int(indptr[stop])
                data = np.array(self._data[s:e])
                indices = np.array(self._indices[s:e])
                new_indptr = np.array(indptr[start:stop + 1]) - int(indptr[start])
                return sp.csc_matrix(
                    (data, indices, new_indptr), shape=(nrow, stop - start)
                )
            cols = np.arange(start, stop, step, dtype=np.intp)

        cols = np.asarray(cols, dtype=np.intp)
        starts = np.asarray(indptr[cols])
        ends = np.asarray(indptr[cols + 1])
        sizes = ends - starts

        new_indptr = np.empty(cols.shape[0] + 1, dtype=np.int64)
        new_indptr[0] = 0
        np.cumsum(sizes, out=new_indptr[1:])
        total = int(new_indptr[-1])

        new_data = np.empty(total, dtype=self.dtype)
        new_indices = np.empty(total, dtype=self._indices.dtype)
        for k in range(cols.shape[0]):
            s, e = int(starts[k]), int(ends[k])
            a, b = int(new_indptr[k]), int(new_indptr[k + 1])
            new_data[a:b] = self._data[s:e]
            new_indices[a:b] = self._indices[s:e]
        return sp.csc_matrix(
            (new_data, new_indices, new_indptr), shape=(nrow, cols.shape[0])
        )

    @staticmethod
    def _resolve(sel, n: int) -> Union[slice, np.ndarray]:
        """Normalise one axis selector to a plain slice or an integer index array."""
        if isinstance(sel, slice):
            return sel
        arr = np.asarray(sel)
        if arr.dtype == bool:
            return np.flatnonzero(arr)
        return arr.ravel().astype(np.intp)

    def __getitem__(self, key) -> sp.csc_matrix:
        if isinstance(key, tuple):
            if len(key) != 2:
                raise IndexError("LazyMatrix supports 2-D indexing only.")
            row_sel, col_sel = key
        else:
            row_sel, col_sel = key, slice(None)

        rows = self._resolve(row_sel, self.nrow)
        cols = self._resolve(col_sel, self.ncol)

        sub = self._read_columns(cols)
        if not (isinstance(rows, slice) and rows == slice(None, None, None)):
            sub = sub[rows, :]
        return sp.csc_matrix(sub)

    def col_blocks(
        self, block_size: int = 10_000
    ) -> Iterator[Tuple[int, int, sp.csc_matrix]]:
        """Stream the matrix in blocks of ``block_size`` columns (cells).

        Yields ``(start, stop, block)`` where ``block`` is an in-memory
        ``csc_matrix`` of columns ``[start:stop)``. This is the primitive for an
        out-of-core reduction: process a million cells at bounded peak memory.
        """
        if block_size <= 0:
            raise ValueError("block_size must be positive.")
        for start in range(0, self.ncol, block_size):
            stop = min(start + block_size, self.ncol)
            yield start, stop, self._read_columns(slice(start, stop))

    # ------------------------------------------------------------------
    # Streaming reductions (single pass over the on-disk arrays)
    # ------------------------------------------------------------------

    def sum(self, axis: Optional[int] = None):
        """Sum over ``axis`` (0 → per-cell, 1 → per-feature, None → scalar)."""
        data = np.asarray(self._data)
        if axis is None:
            return float(data.sum())
        if data.size == 0:
            return np.zeros(self.ncol if axis == 0 else self.nrow, dtype=float)
        if axis == 0:
            counts = np.diff(np.asarray(self._indptr))
            col_ids = np.repeat(np.arange(self.ncol), counts)
            return np.bincount(col_ids, weights=data, minlength=self.ncol).astype(float)
        if axis == 1:
            return np.bincount(
                np.asarray(self._indices), weights=data, minlength=self.nrow
            ).astype(float)
        raise ValueError("axis must be 0, 1, or None.")

    def mean(self, axis: Optional[int] = None):
        """Mean over ``axis``, dividing the streamed sums by the matrix extent."""
        total = self.sum(axis)
        if axis is None:
            return total / (self.nrow * self.ncol)
        denom = self.nrow if axis == 0 else self.ncol
        return total / denom

    def nnz_per_col(self) -> np.ndarray:
        """Non-zeros per column (cell) — the ``nFeature`` count, read from indptr."""
        return np.diff(np.asarray(self._indptr)).astype(np.int64)

    # ------------------------------------------------------------------
    # Materialisation (the escape hatches — avoid on the huge path)
    # ------------------------------------------------------------------

    def to_scipy(self) -> sp.csc_matrix:
        """Read the whole matrix into an in-memory ``csc_matrix``."""
        return sp.csc_matrix(
            (np.array(self._data), np.array(self._indices), np.array(self._indptr)),
            shape=self._shape,
        )

    def toarray(self) -> np.ndarray:
        """Read the whole matrix into a dense ``ndarray``."""
        return self.to_scipy().toarray()

    def __array__(self, dtype=None) -> np.ndarray:
        arr = self.toarray()
        return arr.astype(dtype) if dtype is not None else arr

    # ------------------------------------------------------------------
    # Resource management
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Release the memory-mapped arrays."""
        for name in ("_data", "_indices", "_indptr"):
            arr = getattr(self, name, None)
            mm = getattr(arr, "_mmap", None) if arr is not None else None
            if mm is not None:
                mm.close()
            setattr(self, name, None)

    def __enter__(self) -> "LazyMatrix":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


# ----------------------------------------------------------------------
# Persist / open
# ----------------------------------------------------------------------


def write_lazy_matrix(
    matrix, path: Union[str, Path], *, overwrite: bool = False
) -> LazyMatrix:
    """Write ``matrix`` to ``path`` as an on-disk CSC store and open it lazily.

    ``matrix`` may be a scipy sparse matrix, a dense array-like, or another
    :class:`LazyMatrix`; it is canonicalised to sorted, duplicate-summed CSC
    before being saved as three ``.npy`` arrays plus a JSON header. Returns a
    :class:`LazyMatrix` mapping the freshly written store.
    """
    path = Path(path)
    if isinstance(matrix, LazyMatrix):
        matrix = matrix.to_scipy()
    m = as_sparse(matrix, "csc")
    m.sum_duplicates()
    m.sort_indices()

    if path.exists():
        looks_like_store = (path / _META).exists()
        is_empty_dir = path.is_dir() and not any(path.iterdir())
        if not overwrite:
            raise FileExistsError(
                f"'{path}' already exists; pass overwrite=True to replace it."
            )
        if not (looks_like_store or is_empty_dir):
            raise ValueError(
                f"Refusing to overwrite '{path}': it is not a lazy-matrix store."
            )
        shutil.rmtree(path)
    path.mkdir(parents=True)

    np.save(path / _DATA, m.data)
    np.save(path / _INDICES, m.indices)
    np.save(path / _INDPTR, m.indptr)
    meta = {
        "format": _FORMAT,
        "shape": [int(m.shape[0]), int(m.shape[1])],
        "dtype": str(m.data.dtype),
        "nnz": int(m.nnz),
    }
    (path / _META).write_text(json.dumps(meta, indent=2))
    return open_lazy_matrix(path)


def open_lazy_matrix(path: Union[str, Path]) -> LazyMatrix:
    """Memory-map an on-disk CSC store written by :func:`write_lazy_matrix`."""
    path = Path(path)
    meta_path = path / _META
    if not meta_path.exists():
        raise FileNotFoundError(f"No lazy-matrix store at '{path}'.")
    meta = json.loads(meta_path.read_text())
    if meta.get("format") != _FORMAT:
        raise ValueError(f"Unknown lazy-matrix format: {meta.get('format')!r}.")

    data = np.load(path / _DATA, mmap_mode="r")
    indices = np.load(path / _INDICES, mmap_mode="r")
    indptr = np.load(path / _INDPTR, mmap_mode="r")
    return LazyMatrix(path, tuple(meta["shape"]), data, indices, indptr)


def is_lazy(x) -> bool:
    """True if ``x`` is a :class:`LazyMatrix` (an on-disk, memory-mapped layer)."""
    return isinstance(x, LazyMatrix)
