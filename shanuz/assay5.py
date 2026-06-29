from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import Optional, Type, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import as_sparse, empty_sparse, is_matrix_empty
from ._utils import validate_cell_names, validate_feature_names
from .logmap import LogMap
from .mixins import KeyMixin


class StdAssay(KeyMixin, ABC):
    """Abstract base for layered assays (v5 architecture).

    Mirrors R's StdAssay virtual class from assay5.R.
    Unlike the legacy Assay (v3), StdAssay stores arbitrary named layers
    and uses LogMap to track which cells/features belong to each layer.

    Slots
    -----
    layers      : dict[str, AnyMatrix]   named expression matrices (features × cells)
    cells       : LogMap                 per-layer boolean cell membership
    features    : LogMap                 per-layer boolean feature membership
    default     : int                    index of the default layer
    assay_orig  : Optional[str]
    meta_data   : pd.DataFrame           per-feature metadata
    misc        : dict
    _key        : str                    (from KeyMixin)
    """

    __slots__ = (
        "layers",
        "_cells",
        "_features",
        "default",
        "assay_orig",
        "meta_data",
        "misc",
        "_key",
        "_all_cell_names",
        "_all_feature_names",
        "_scaled_features",
        "_var_features",
    )

    def __init__(
        self,
        layers: dict[str, Union[np.ndarray, sp.spmatrix]],
        feature_names: list[str],
        cell_names: list[str],
        assay_orig: Optional[str] = None,
        meta_data: Optional[pd.DataFrame] = None,
        misc: Optional[dict] = None,
        key: str = "rna_",
        default: int = 0,
    ) -> None:
        self._key = key
        self.assay_orig = assay_orig
        self.misc = misc or {}
        self.default = default

        validate_feature_names(feature_names)
        validate_cell_names(cell_names)
        self._all_feature_names = list(feature_names)
        self._all_cell_names = list(cell_names)

        self.layers: dict[str, Union[np.ndarray, sp.spmatrix]] = {}
        self._cells = LogMap()
        self._features = LogMap()
        self._scaled_features: list[str] = []
        self._var_features: list[str] = []

        for name, mat in layers.items():
            self._add_layer(name, mat)

        self.meta_data = (
            meta_data
            if meta_data is not None
            else pd.DataFrame(index=self._all_feature_names)
        )

    # ------------------------------------------------------------------
    # Internal layer management
    # ------------------------------------------------------------------

    def _add_layer(self, name: str, mat: Union[np.ndarray, sp.spmatrix]) -> None:
        n_feat = len(self._all_feature_names)
        n_cells = len(self._all_cell_names)

        if mat.shape[0] != n_feat:
            raise ValueError(
                f"Layer '{name}' has {mat.shape[0]} rows but "
                f"{n_feat} features are registered."
            )
        if mat.shape[1] != n_cells:
            raise ValueError(
                f"Layer '{name}' has {mat.shape[1]} columns but "
                f"{n_cells} cells are registered."
            )

        self.layers[name] = mat
        self._cells[name] = np.ones(n_cells, dtype=bool)
        self._features[name] = np.ones(n_feat, dtype=bool)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def cells(self, layer: Optional[str] = None) -> list[str]:
        if layer is not None and layer in self._cells:
            mask = self._cells[layer]
            return [c for c, m in zip(self._all_cell_names, mask) if m]
        return list(self._all_cell_names)

    def features(self, layer: Optional[str] = None) -> list[str]:
        if layer is not None and layer in self._features:
            mask = self._features[layer]
            return [f for f, m in zip(self._all_feature_names, mask) if m]
        return list(self._all_feature_names)

    def layers_list(self, pattern: Optional[str] = None) -> list[str]:
        names = list(self.layers)
        if pattern:
            names = [n for n in names if re.search(pattern, n)]
        return names

    # ------------------------------------------------------------------
    # Default layer
    # ------------------------------------------------------------------

    @property
    def default_layer(self) -> Optional[str]:
        names = list(self.layers)
        if not names:
            return None
        idx = min(self.default, len(names) - 1)
        return names[idx]

    @default_layer.setter
    def default_layer(self, value: str) -> None:
        names = list(self.layers)
        if value not in names:
            raise KeyError(f"Layer '{value}' not found.")
        self.default = names.index(value)

    # ------------------------------------------------------------------
    # Layer data access
    # ------------------------------------------------------------------

    def layer_data(
        self,
        layer: Optional[str] = None,
        cells: Optional[list[str]] = None,
        features: Optional[list[str]] = None,
    ):
        if layer is None:
            layer = self.default_layer
        if layer is None:
            raise ValueError("No layers available.")
        if layer not in self.layers:
            raise KeyError(f"Layer '{layer}' not found.")

        mat = self.layers[layer]
        all_feat = self._all_feature_names
        all_cells = self._all_cell_names

        row_idx = (
            [all_feat.index(f) for f in features] if features is not None else slice(None)
        )
        col_idx = (
            [all_cells.index(c) for c in cells] if cells is not None else slice(None)
        )

        if isinstance(row_idx, list) or isinstance(col_idx, list):
            r = row_idx if isinstance(row_idx, list) else list(range(mat.shape[0]))
            c = col_idx if isinstance(col_idx, list) else list(range(mat.shape[1]))
            return mat[np.ix_(r, c)]
        return mat

    def set_layer_data(
        self,
        layer: str,
        value: Union[np.ndarray, sp.spmatrix],
        cell_names: Optional[list[str]] = None,
        feature_names: Optional[list[str]] = None,
    ) -> None:
        if layer in self.layers:
            self.layers[layer] = value
        else:
            self._add_layer(layer, value)

    # ------------------------------------------------------------------
    # Variable features
    # ------------------------------------------------------------------

    @property
    def variable_features(self) -> list[str]:
        # Return the ordered list (selection-rank order, not index order)
        if self._var_features:
            return list(self._var_features)
        # Fallback: highly_variable column in meta_data (unordered)
        col = "highly_variable"
        if col in self.meta_data.columns:
            return list(self.meta_data.index[self.meta_data[col].astype(bool)])
        return []

    @variable_features.setter
    def variable_features(self, value: list[str]) -> None:
        self._var_features = list(value)
        # Also keep a boolean column for downstream use
        self.meta_data["highly_variable"] = self.meta_data.index.isin(value)

    # ------------------------------------------------------------------
    # Join / split layers
    # ------------------------------------------------------------------

    def join_layers(self, layers: Optional[list[str]] = None) -> "StdAssay":
        names = layers if layers is not None else list(self.layers)
        if not names:
            return self
        mats = [self.layers[n] for n in names if n in self.layers]
        if not mats:
            return self
        if sp.issparse(mats[0]):
            joined = sp.hstack(mats, format="csc")
        else:
            joined = np.hstack(mats)
        # For simplicity: joined layer keeps all cells from all joined layers
        new_obj = self._copy()
        for n in names:
            del new_obj.layers[n]
        new_obj.layers["joined"] = joined
        return new_obj

    def split_layers(self, f: list[str], layer: Optional[str] = None) -> "StdAssay":
        if layer is None:
            layer = self.default_layer
        if layer is None:
            raise ValueError("No layer to split.")
        mat = self.layers[layer]
        if len(f) != mat.shape[1]:
            raise ValueError("f must have one entry per cell.")

        groups = {}
        for i, g in enumerate(f):
            groups.setdefault(g, []).append(i)

        new_obj = self._copy()
        del new_obj.layers[layer]
        for g, idxs in groups.items():
            sub = mat[:, idxs]
            new_obj.layers[f"{layer}_{g}"] = sub
        return new_obj

    # ------------------------------------------------------------------
    # Cast assay layer types
    # ------------------------------------------------------------------

    def cast_assay(self, to_sparse: bool = True) -> "StdAssay":
        new_obj = self._copy()
        for name, mat in new_obj.layers.items():
            if to_sparse and not sp.issparse(mat):
                new_obj.layers[name] = sp.csc_matrix(mat)
            elif not to_sparse and sp.issparse(mat):
                new_obj.layers[name] = mat.toarray()
        return new_obj

    # ------------------------------------------------------------------
    # Subsetting
    # ------------------------------------------------------------------

    def subset(
        self,
        cells: Optional[list[str]] = None,
        features: Optional[list[str]] = None,
    ) -> "StdAssay":
        all_feat = self._all_feature_names
        all_cells = self._all_cell_names

        row_idx = (
            [all_feat.index(f) for f in features]
            if features is not None
            else list(range(len(all_feat)))
        )
        col_idx = (
            [all_cells.index(c) for c in cells]
            if cells is not None
            else list(range(len(all_cells)))
        )

        new_layers = {}
        for name, mat in self.layers.items():
            new_layers[name] = mat[np.ix_(row_idx, col_idx)]

        new_features = [all_feat[i] for i in row_idx]
        new_cells = [all_cells[i] for i in col_idx]
        new_meta = self.meta_data.iloc[row_idx].copy()

        return self.__class__(
            layers=new_layers,
            feature_names=new_features,
            cell_names=new_cells,
            assay_orig=self.assay_orig,
            meta_data=new_meta,
            misc=dict(self.misc),
            key=self._key,
            default=self.default,
        )

    # ------------------------------------------------------------------
    # Cell renaming
    # ------------------------------------------------------------------

    def rename_cells(self, new_names: list[str]) -> "StdAssay":
        if len(new_names) != len(self._all_cell_names):
            raise ValueError("new_names must match the number of cells.")
        obj = self._copy()
        obj._all_cell_names = list(new_names)
        return obj

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def merge(
        self,
        y: Union["StdAssay", list["StdAssay"]],
        add_cell_ids: Optional[list[str]] = None,
    ) -> "StdAssay":
        others = [y] if isinstance(y, StdAssay) else y
        all_assays = [self] + others

        if add_cell_ids is not None and len(add_cell_ids) != len(all_assays):
            raise ValueError("add_cell_ids must have one entry per assay.")

        new_cell_names: list[str] = []
        for idx, a in enumerate(all_assays):
            prefix = add_cell_ids[idx] if add_cell_ids else None
            for c in a._all_cell_names:
                new_cell_names.append(f"{prefix}_{c}" if prefix else c)

        shared = set(self._all_feature_names)
        for a in others:
            shared &= set(a._all_feature_names)
        shared_features = [f for f in self._all_feature_names if f in shared]

        shared_layers = set(self.layers)
        for a in others:
            shared_layers &= set(a.layers)

        new_layers: dict = {}
        for layer_name in shared_layers:
            mats = []
            for a in all_assays:
                feat_idx = [a._all_feature_names.index(f) for f in shared_features]
                mat = a.layers[layer_name][feat_idx, :]
                mats.append(mat)
            if sp.issparse(mats[0]):
                new_layers[layer_name] = sp.hstack(mats, format="csc")
            else:
                new_layers[layer_name] = np.hstack(mats)

        return self.__class__(
            layers=new_layers,
            feature_names=shared_features,
            cell_names=new_cell_names,
            assay_orig=self.assay_orig,
            misc=dict(self.misc),
            key=self._key,
        )

    # ------------------------------------------------------------------
    # Calc nCount / nFeature
    # ------------------------------------------------------------------

    def calc_n(self) -> "pd.DataFrame":
        from ._utils import calc_n as _calc_n
        default = self.default_layer
        if default is None:
            n = len(self._all_cell_names)
            return pd.DataFrame(
                {"nCount": np.zeros(n), "nFeature": np.zeros(n)},
                index=self._all_cell_names,
            )
        mat = self.layers[default]
        ncount, nfeature = _calc_n(mat)
        return pd.DataFrame(
            {"nCount": ncount, "nFeature": nfeature},
            index=self._all_cell_names,
        )

    # ------------------------------------------------------------------
    # Copy helper
    # ------------------------------------------------------------------

    def _copy(self) -> "StdAssay":
        new = self.__class__(
            layers={k: v.copy() if hasattr(v, "copy") else v for k, v in self.layers.items()},
            feature_names=list(self._all_feature_names),
            cell_names=list(self._all_cell_names),
            assay_orig=self.assay_orig,
            meta_data=self.meta_data.copy(),
            misc=dict(self.misc),
            key=self._key,
            default=self.default,
        )
        new._var_features = list(self._var_features)
        new._scaled_features = list(self._scaled_features)
        return new

    # ------------------------------------------------------------------
    # Operators
    # ------------------------------------------------------------------

    def __getitem__(self, key):
        if isinstance(key, tuple):
            features, cells = key
        else:
            features, cells = key, None
        return self.subset(
            cells=cells if isinstance(cells, list) else None,
            features=features if isinstance(features, list) else None,
        )

    def __repr__(self) -> str:
        n_feat = len(self._all_feature_names)
        n_cells = len(self._all_cell_names)
        layer_names = list(self.layers)
        return (
            f"{self.__class__.__name__} with {n_feat} features and {n_cells} cells\n"
            f"  Key: {self._key!r}\n"
            f"  Layers: {layer_names}\n"
            f"  Default layer: {self.default_layer!r}"
        )


class Assay5(StdAssay):
    """Modern layered assay (v5).

    Mirrors R's Assay5 class from assay5.R.
    Extends StdAssay with no additional slots.
    """

    __slots__ = ()


# ------------------------------------------------------------------
# Factory
# ------------------------------------------------------------------

def create_assay5_object(
    counts=None,
    data=None,
    min_cells: int = 0,
    min_features: int = 0,
    feature_names: Optional[list[str]] = None,
    cell_names: Optional[list[str]] = None,
    key: str = "rna_",
) -> Assay5:
    matrix = counts if counts is not None else data
    if matrix is None:
        raise ValueError("Provide at least one of 'counts' or 'data'.")

    if sp.issparse(matrix):
        mat = matrix.tocsc()
    else:
        mat = sp.csc_matrix(np.asarray(matrix))

    n_features, n_cells = mat.shape

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(n_features)]
    if cell_names is None:
        cell_names = [f"cell_{i}" for i in range(n_cells)]

    # Filter features by min_cells
    if min_cells > 0:
        nnz_per_feat = np.diff(mat.T.tocsc().indptr)
        keep_feat = np.where(nnz_per_feat >= min_cells)[0]
        mat = mat[keep_feat, :]
        feature_names = [feature_names[i] for i in keep_feat]

    # Filter cells by min_features
    if min_features > 0:
        nnz_per_cell = np.diff(mat.tocsc().indptr)
        keep_cell = np.where(nnz_per_cell >= min_features)[0]
        mat = mat[:, keep_cell]
        cell_names = [cell_names[i] for i in keep_cell]

    layer_name = "counts" if counts is not None else "data"
    return Assay5(
        layers={layer_name: mat},
        feature_names=list(feature_names),
        cell_names=list(cell_names),
        key=key,
    )
