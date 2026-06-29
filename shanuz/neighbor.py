from __future__ import annotations

from typing import Any, Optional

import numpy as np
import scipy.sparse as sp


class Neighbor:
    """Nearest-neighbor results for a set of cells.

    Mirrors R's Neighbor class from neighbor.R.

    Slots
    -----
    nn_idx       : int matrix  (n_cells × k), neighbor indices (1-based in R; 0-based here)
    nn_dist      : float matrix (n_cells × k), corresponding distances
    alg_idx      : Any         algorithm index object (e.g. annoy index)
    alg_info     : dict        metadata about the algorithm used
    cell_names   : list[str]   cell barcodes, length n_cells
    """

    __slots__ = ("nn_idx", "nn_dist", "alg_idx", "alg_info", "cell_names")

    def __init__(
        self,
        nn_idx: np.ndarray,
        nn_dist: np.ndarray,
        cell_names: list[str],
        alg_idx: Any = None,
        alg_info: Optional[dict] = None,
    ) -> None:
        self.nn_idx = np.asarray(nn_idx, dtype=int)
        self.nn_dist = np.asarray(nn_dist, dtype=float)
        self.cell_names = list(cell_names)
        self.alg_idx = alg_idx
        self.alg_info = alg_info or {}
        self._validate()

    def _validate(self) -> None:
        if self.nn_idx.shape != self.nn_dist.shape:
            raise ValueError("nn_idx and nn_dist must have the same shape.")
        if self.nn_idx.shape[0] != len(self.cell_names):
            raise ValueError(
                f"nn_idx has {self.nn_idx.shape[0]} rows but "
                f"{len(self.cell_names)} cell names were provided."
            )

    # ------------------------------------------------------------------
    # Accessors — mirrors R generics
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        return list(self.cell_names)

    def indices(self) -> np.ndarray:
        mat = self.nn_idx.copy()
        return mat

    def distances(self) -> np.ndarray:
        mat = self.nn_dist.copy()
        return mat

    def index(self) -> Any:
        return self.alg_idx

    def dim(self) -> tuple[int, int]:
        return self.nn_idx.shape

    # ------------------------------------------------------------------
    # Cell renaming
    # ------------------------------------------------------------------

    def rename_cells(
        self,
        new_names: Optional[list[str]] = None,
        old_names: Optional[list[str]] = None,
    ) -> "Neighbor":
        if new_names is not None and old_names is None:
            if len(new_names) != len(self.cell_names):
                raise ValueError("new_names must have the same length as cell_names.")
            return Neighbor(
                nn_idx=self.nn_idx.copy(),
                nn_dist=self.nn_dist.copy(),
                cell_names=list(new_names),
                alg_idx=self.alg_idx,
                alg_info=dict(self.alg_info),
            )
        if old_names is not None and new_names is not None:
            mapping = dict(zip(old_names, new_names))
            updated = [mapping.get(c, c) for c in self.cell_names]
            return Neighbor(
                nn_idx=self.nn_idx.copy(),
                nn_dist=self.nn_dist.copy(),
                cell_names=updated,
                alg_idx=self.alg_idx,
                alg_info=dict(self.alg_info),
            )
        raise ValueError("Provide either new_names alone or both old_names and new_names.")

    # ------------------------------------------------------------------
    # Conversion to Graph
    # ------------------------------------------------------------------

    def as_graph(self, weighted: bool = True) -> "Graph":
        from .graph import Graph

        n = len(self.cell_names)
        k = self.nn_idx.shape[1]

        row_idx = np.repeat(np.arange(n), k)
        col_idx = self.nn_idx.flatten()
        data = self.nn_dist.flatten() if weighted else np.ones(n * k, dtype=float)

        mat = sp.coo_matrix((data, (row_idx, col_idx)), shape=(n, n))
        return Graph(matrix=mat.tocsc(), cell_names=self.cell_names)

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n, k = self.nn_idx.shape
        return f"Neighbor for {n} cells, {k} neighbors"
