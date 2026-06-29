from __future__ import annotations

from typing import Optional

import numpy as np

from ._sparse import empty_dense
from .jackstraw import JackStrawData
from .mixins import KeyMixin


class DimReduc(KeyMixin):
    """Stores a dimensionality reduction (PCA, UMAP, tSNE, …).

    Mirrors R's DimReduc class from dimreduc.R.

    Slots
    -----
    cell_embeddings              : np.ndarray  (n_cells × n_dims)   required
    feature_loadings             : np.ndarray  (n_features × n_dims) optional
    feature_loadings_projected   : np.ndarray  projected loadings    optional
    assay_used                   : str         source assay name
    global_                      : bool        if True, persists when assay is removed
    stdev                        : np.ndarray  per-dimension std devs
    jackstraw                    : JackStrawData
    misc                         : dict
    _key                         : str         prefix, e.g. "PC_"
    """

    __slots__ = (
        "cell_embeddings",
        "feature_loadings",
        "feature_loadings_projected",
        "assay_used",
        "global_",
        "stdev",
        "jackstraw",
        "misc",
        "_key",
        "_cell_names",
        "_feature_names",
    )

    def __init__(
        self,
        cell_embeddings: np.ndarray,
        cell_names: list[str],
        feature_loadings: Optional[np.ndarray] = None,
        feature_names: Optional[list[str]] = None,
        feature_loadings_projected: Optional[np.ndarray] = None,
        assay_used: str = "",
        global_: bool = False,
        stdev: Optional[np.ndarray] = None,
        jackstraw: Optional[JackStrawData] = None,
        misc: Optional[dict] = None,
        key: str = "PC_",
    ) -> None:
        self._key = key
        self.cell_embeddings = np.asarray(cell_embeddings)
        self._cell_names = list(cell_names)

        if len(self._cell_names) != self.cell_embeddings.shape[0]:
            raise ValueError(
                f"cell_names length {len(self._cell_names)} != "
                f"cell_embeddings rows {self.cell_embeddings.shape[0]}."
            )

        n_dims = self.cell_embeddings.shape[1]

        if feature_loadings is not None:
            self.feature_loadings = np.asarray(feature_loadings)
        else:
            self.feature_loadings = empty_dense(0, n_dims)

        self._feature_names = list(feature_names) if feature_names else []

        if feature_loadings_projected is not None:
            self.feature_loadings_projected = np.asarray(feature_loadings_projected)
        else:
            self.feature_loadings_projected = empty_dense(0, n_dims)

        self.assay_used = assay_used
        self.global_ = global_
        self.stdev = np.asarray(stdev) if stdev is not None else np.array([])
        self.jackstraw = jackstraw if jackstraw is not None else JackStrawData()
        self.misc = misc or {}

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def embeddings(self) -> np.ndarray:
        return self.cell_embeddings

    def loadings(self, projected: bool = False) -> np.ndarray:
        return self.feature_loadings_projected if projected else self.feature_loadings

    def set_loadings(self, value: np.ndarray, projected: bool = False) -> None:
        if projected:
            self.feature_loadings_projected = np.asarray(value)
        else:
            self.feature_loadings = np.asarray(value)

    def cells(self) -> list[str]:
        return list(self._cell_names)

    def features(self, projected: bool = False) -> list[str]:
        return list(self._feature_names)

    def default_assay(self) -> str:
        return self.assay_used

    def set_default_assay(self, value: str) -> None:
        self.assay_used = value

    def is_global(self) -> bool:
        return self.global_

    # ------------------------------------------------------------------
    # Cell renaming
    # ------------------------------------------------------------------

    def rename_cells(self, new_names: list[str]) -> "DimReduc":
        if len(new_names) != len(self._cell_names):
            raise ValueError("new_names must match the number of cells.")
        obj = DimReduc(
            cell_embeddings=self.cell_embeddings.copy(),
            cell_names=list(new_names),
            feature_loadings=self.feature_loadings.copy(),
            feature_names=list(self._feature_names),
            feature_loadings_projected=self.feature_loadings_projected.copy(),
            assay_used=self.assay_used,
            global_=self.global_,
            stdev=self.stdev.copy(),
            jackstraw=self.jackstraw,
            misc=dict(self.misc),
            key=self._key,
        )
        return obj

    # ------------------------------------------------------------------
    # Subsetting
    # ------------------------------------------------------------------

    def subset(
        self,
        cells: Optional[list[str]] = None,
        dims: Optional[list[int]] = None,
    ) -> "DimReduc":
        cell_idx = (
            [self._cell_names.index(c) for c in cells]
            if cells is not None
            else list(range(len(self._cell_names)))
        )
        dim_idx = dims if dims is not None else list(range(self.cell_embeddings.shape[1]))

        new_emb = self.cell_embeddings[np.ix_(cell_idx, dim_idx)]
        new_load = (
            self.feature_loadings[:, dim_idx]
            if self.feature_loadings.shape[0] > 0
            else self.feature_loadings
        )
        new_stdev = self.stdev[dim_idx] if self.stdev.size > 0 else self.stdev

        return DimReduc(
            cell_embeddings=new_emb,
            cell_names=[self._cell_names[i] for i in cell_idx],
            feature_loadings=new_load,
            feature_names=list(self._feature_names),
            assay_used=self.assay_used,
            global_=self.global_,
            stdev=new_stdev,
            misc=dict(self.misc),
            key=self._key,
        )

    # ------------------------------------------------------------------
    # Operators
    # ------------------------------------------------------------------

    def __getitem__(self, key):
        if isinstance(key, (list, np.ndarray)):
            return self.subset(cells=key)
        raise TypeError(f"Unsupported index type: {type(key).__name__}")

    def __repr__(self) -> str:
        n_cells = len(self._cell_names)
        n_dims = self.cell_embeddings.shape[1] if self.cell_embeddings.ndim == 2 else 0
        return (
            f"DimReduc: {n_cells} cells × {n_dims} dimensions\n"
            f"  Key: {self._key!r}  Assay: {self.assay_used!r}  Global: {self.global_}"
        )
