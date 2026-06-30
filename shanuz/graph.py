from __future__ import annotations

from typing import Optional, Union

import numpy as np
import scipy.sparse as sp


class Graph:
    """Sparse graph object for cell-cell relationships (e.g. SNN graph).

    Mirrors R's Graph class from graph.R, which extends dgCMatrix.
    In Python we wrap (not inherit) a scipy CSC matrix to avoid
    scipy subclassing pitfalls.

    Slots
    -----
    _matrix    : scipy.sparse.csc_matrix   underlying adjacency matrix
    assay_used : Optional[str]             assay that generated this graph
    _cell_names: list[str]                 row/col names (cells)
    """

    __slots__ = ("_matrix", "assay_used", "_cell_names")

    def __init__(
        self,
        matrix: sp.spmatrix,
        cell_names: list[str],
        assay_used: Optional[str] = None,
    ) -> None:
        if not sp.issparse(matrix):
            matrix = sp.csc_matrix(matrix)
        else:
            matrix = matrix.tocsc()

        n = len(cell_names)
        if matrix.shape != (n, n):
            raise ValueError(
                f"matrix shape {matrix.shape} does not match "
                f"cell_names length {n} × {n}."
            )
        if len(cell_names) != len(set(cell_names)):
            raise ValueError("cell_names must be unique.")

        self._matrix = matrix
        self._cell_names = list(cell_names)
        self.assay_used = assay_used
        self._validate()

    def _validate(self) -> None:
        if self.assay_used is not None and len(self.assay_used) == 0:
            raise ValueError("assay_used must be None or a non-empty string.")

    # ------------------------------------------------------------------
    # Delegate matrix operations to _matrix
    # ------------------------------------------------------------------

    def __getattr__(self, name: str):
        if name.startswith("_") or name == "assay_used":
            raise AttributeError(name)
        return getattr(self._matrix, name)

    @property
    def shape(self) -> tuple[int, int]:
        return self._matrix.shape

    def toarray(self) -> np.ndarray:
        return self._matrix.toarray()

    def tocsr(self) -> sp.csr_matrix:
        return self._matrix.tocsr()

    def tocsc(self) -> sp.csc_matrix:
        return self._matrix

    def __getitem__(self, key):
        return self._matrix[key]

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        return list(self._cell_names)

    def subset(self, cells: list[str]) -> "Graph":
        """Return a new Graph restricted to ``cells`` (cell×cell submatrix)."""
        idx_map = {c: i for i, c in enumerate(self._cell_names)}
        keep = [c for c in cells if c in idx_map]
        idx = [idx_map[c] for c in keep]
        sub = self._matrix[np.ix_(idx, idx)]
        return Graph(matrix=sub, cell_names=keep, assay_used=self.assay_used)

    def default_assay(self) -> Optional[str]:
        return self.assay_used

    def set_default_assay(self, value: Optional[str]) -> None:
        if value is not None and len(value) == 0:
            raise ValueError("assay_used must be None or a non-empty string.")
        self.assay_used = value

    # ------------------------------------------------------------------
    # Conversion to Neighbor
    # ------------------------------------------------------------------

    def as_neighbor(self) -> "Neighbor":
        from .neighbor import Neighbor

        mat = self._matrix.tocoo()
        n = self._matrix.shape[0]
        rows, cols, data = mat.row, mat.col, mat.data

        # Build dense nn_idx / nn_dist from COO entries
        k_per_row = np.bincount(rows, minlength=n)
        k = int(k_per_row.max()) if n > 0 else 0

        nn_idx = np.zeros((n, k), dtype=int)
        nn_dist = np.zeros((n, k), dtype=float)
        counters = np.zeros(n, dtype=int)

        for r, c, d in zip(rows, cols, data):
            pos = counters[r]
            if pos < k:
                nn_idx[r, pos] = c
                nn_dist[r, pos] = d
                counters[r] += 1

        return Neighbor(
            nn_idx=nn_idx,
            nn_dist=nn_dist,
            cell_names=self._cell_names,
        )

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n = len(self._cell_names)
        assay = f", assay={self.assay_used!r}" if self.assay_used else ""
        return f"Graph({n}×{n} sparse{assay})"


# ------------------------------------------------------------------
# Factory: as_graph
# ------------------------------------------------------------------

def as_graph(
    x: Union[np.ndarray, sp.spmatrix, "Neighbor"],
    cell_names: Optional[list[str]] = None,
    assay_used: Optional[str] = None,
    weighted: bool = True,
) -> Graph:
    """Convert a matrix or Neighbor to a Graph.  Mirrors R as.Graph()."""
    from .neighbor import Neighbor

    if isinstance(x, Graph):
        return x

    if isinstance(x, Neighbor):
        g = x.as_graph(weighted=weighted)
        if assay_used is not None:
            g.assay_used = assay_used
        return g

    if sp.issparse(x) or isinstance(x, np.ndarray):
        if cell_names is None:
            n = x.shape[0]
            cell_names = [str(i) for i in range(n)]
        return Graph(matrix=x, cell_names=cell_names, assay_used=assay_used)

    raise TypeError(f"Cannot convert {type(x).__name__} to Graph.")
