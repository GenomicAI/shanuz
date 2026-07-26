from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import (
    as_sparse,
    empty_dense,
    empty_sparse,
    is_matrix_empty,
)
from ._utils import calc_n, validate_cell_names, validate_feature_names
from .mixins import KeyMixin


class Assay(KeyMixin):
    """Legacy (v3) Assay object.

    Mirrors R's Assay class from assay.R.

    Slots
    -----
    counts        : raw counts / TPMs  (features × cells)
    data          : normalised expression (features × cells)
    scale_data    : scaled expression  (features × cells, dense) — a *subset*
                    of the features, since ScaleData defaults to the variable
                    ones. R's slot is a matrix and carries its own rownames;
                    a bare ndarray does not, so the labels live alongside it in
                    `_scaled_features` and every read of the layer goes through
                    `features("scale_data")`.
    assay_orig    : name of original assay this was derived from
    var_features  : list of highly variable feature names
    meta_features : per-feature metadata DataFrame (features × cols)
    misc          : dict for miscellaneous storage
    _key          : string key prefix (inherited from KeyMixin)
    """

    __slots__ = (
        "counts",
        "data",
        "scale_data",
        "assay_orig",
        "var_features",
        "meta_features",
        "misc",
        "_key",
        "_feature_names",
        "_cell_names",
        "_scaled_features",
    )

    def __init__(
        self,
        counts: Optional[Union[np.ndarray, sp.spmatrix]] = None,
        data: Optional[Union[np.ndarray, sp.spmatrix]] = None,
        scale_data: Optional[np.ndarray] = None,
        scaled_features: Optional[list[str]] = None,
        feature_names: Optional[list[str]] = None,
        cell_names: Optional[list[str]] = None,
        assay_orig: Optional[str] = None,
        var_features: Optional[list[str]] = None,
        meta_features: Optional[pd.DataFrame] = None,
        misc: Optional[dict] = None,
        key: str = "rna_",
    ) -> None:
        self.key = key  # validated by KeyMixin setter

        # Resolve which matrix to use
        if counts is not None and data is None:
            matrix = counts
        elif data is not None and counts is None:
            matrix = data
        elif counts is not None and data is not None:
            matrix = counts  # counts takes precedence for shape
        else:
            raise ValueError("Provide at least one of 'counts' or 'data'.")

        n_features = matrix.shape[0]
        n_cells = matrix.shape[1]

        if feature_names is None:
            feature_names = [f"feature_{i}" for i in range(n_features)]
        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(n_cells)]

        validate_feature_names(feature_names)
        validate_cell_names(cell_names)

        if len(feature_names) != n_features:
            raise ValueError(
                f"feature_names length {len(feature_names)} != matrix rows {n_features}."
            )
        if len(cell_names) != n_cells:
            raise ValueError(
                f"cell_names length {len(cell_names)} != matrix cols {n_cells}."
            )

        self._feature_names = list(feature_names)
        self._cell_names = list(cell_names)

        self.counts = counts if counts is not None else empty_sparse(n_features, n_cells)
        self.data = data if data is not None else self.counts
        self.scale_data = scale_data if scale_data is not None else empty_dense(0, n_cells)
        self._scaled_features = self._resolve_scaled_features(scaled_features)
        self.assay_orig = assay_orig
        self.var_features = list(var_features) if var_features else []
        self.meta_features = (
            meta_features
            if meta_features is not None
            else pd.DataFrame(index=self._feature_names)
        )
        self.misc = misc or {}

    # ------------------------------------------------------------------
    # Cell / feature accessors
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        return list(self._cell_names)

    def _resolve_scaled_features(
        self, scaled_features: Optional[list[str]]
    ) -> list[str]:
        """Labels for `scale_data`'s rows, checked against its height.

        Only a full-height (or empty) `scale_data` can be labelled without being
        told: any other height is a subset, and which subset is not recoverable
        from the matrix. Refusing to guess is the whole point — the guess that
        used to stand here took the *first* n features, and ScaleData's subset is
        the variable ones, which are scattered through the assay.
        """
        n_rows = self.scale_data.shape[0]
        if scaled_features is not None:
            if len(scaled_features) != n_rows:
                raise ValueError(
                    f"scaled_features has {len(scaled_features)} names for "
                    f"{n_rows} rows of scale_data."
                )
            return list(scaled_features)
        if n_rows == 0:
            return []
        if n_rows == len(self._feature_names):
            return list(self._feature_names)
        raise ValueError(
            f"scale_data has {n_rows} rows but the assay has "
            f"{len(self._feature_names)} features, so it holds a subset — pass "
            f"`scaled_features=` naming those rows. They cannot be inferred."
        )

    def features(self, layer: Optional[str] = None) -> list[str]:
        if layer in ("scale_data", "scale.data"):
            return list(self._scaled_features)
        return list(self._feature_names)

    # ------------------------------------------------------------------
    # Layer accessors  (mirrors GetAssayData / SetAssayData / LayerData)
    # ------------------------------------------------------------------

    def get_assay_data(self, layer: str = "data"):
        return self.layer_data(layer)

    def set_assay_data(self, layer: str, new_data, features: Optional[list[str]] = None) -> None:
        self._set_layer(layer, new_data, features)

    def layer_data(
        self,
        layer: str = "data",
        cells: Optional[list[str]] = None,
        features: Optional[list[str]] = None,
    ):
        mat = self._get_layer(layer)
        if features is not None or cells is not None:
            row_idx = (
                [self._feature_names.index(f) for f in features]
                if features is not None
                else slice(None)
            )
            col_idx = (
                [self._cell_names.index(c) for c in cells]
                if cells is not None
                else slice(None)
            )
            mat = mat[row_idx, :][:, col_idx] if isinstance(row_idx, list) else mat[row_idx, col_idx]
        return mat

    def _get_layer(self, layer: str):
        layers = {"counts": self.counts, "data": self.data, "scale_data": self.scale_data}
        if layer not in layers:
            raise KeyError(f"Layer '{layer}' not found. Choose from {list(layers)}.")
        return layers[layer]

    def _set_layer(self, layer: str, value, features: Optional[list[str]] = None) -> None:
        if layer == "counts":
            self.counts = value
        elif layer == "data":
            self.data = value
        elif layer in ("scale_data", "scale.data"):
            self.scale_data = np.asarray(value)
            # Relabel with the new matrix, never leave the old names behind: a
            # stale label list is worse than none, because every reader trusts it.
            self._scaled_features = self._resolve_scaled_features(features)
        else:
            raise KeyError(f"Unknown layer '{layer}'. Must be counts, data, or scale_data.")

    # ------------------------------------------------------------------
    # Variable features
    # ------------------------------------------------------------------

    @property
    def variable_features(self) -> list[str]:
        return list(self.var_features)

    @variable_features.setter
    def variable_features(self, value: list[str]) -> None:
        self.var_features = list(value)

    # ------------------------------------------------------------------
    # Calc nCount / nFeature
    # ------------------------------------------------------------------

    def calc_n(self) -> pd.DataFrame:
        mat = self.counts if not is_matrix_empty(self.counts) else self.data
        ncount, nfeature = calc_n(mat)
        return pd.DataFrame(
            {"nCount": ncount, "nFeature": nfeature},
            index=self._cell_names,
        )

    # ------------------------------------------------------------------
    # Cell renaming
    # ------------------------------------------------------------------

    def rename_cells(self, new_names: list[str]) -> "Assay":
        if len(new_names) != len(self._cell_names):
            raise ValueError("new_names must match the number of cells.")
        obj = self._copy()
        obj._cell_names = list(new_names)
        obj.meta_features = obj.meta_features.copy()
        return obj

    def _copy(self) -> "Assay":
        return Assay(
            counts=self.counts,
            data=self.data,
            scale_data=self.scale_data.copy() if self.scale_data is not None else None,
            scaled_features=list(self._scaled_features),
            feature_names=list(self._feature_names),
            cell_names=list(self._cell_names),
            assay_orig=self.assay_orig,
            var_features=list(self.var_features),
            meta_features=self.meta_features.copy(),
            misc=dict(self.misc),
            key=self._key,
        )

    # ------------------------------------------------------------------
    # Subsetting
    # ------------------------------------------------------------------

    def subset(
        self,
        cells: Optional[list[str]] = None,
        features: Optional[list[str]] = None,
    ) -> "Assay":
        col_idx = (
            [self._cell_names.index(c) for c in cells]
            if cells is not None
            else list(range(len(self._cell_names)))
        )
        row_idx = (
            [self._feature_names.index(f) for f in features]
            if features is not None
            else list(range(len(self._feature_names)))
        )

        def _sub(mat):
            if is_matrix_empty(mat):
                return mat
            return mat[row_idx, :][:, col_idx]

        new_features = [self._feature_names[i] for i in row_idx]
        new_cells = [self._cell_names[i] for i in col_idx]
        new_meta = self.meta_features.iloc[row_idx].copy()
        new_var = [f for f in self.var_features if f in set(new_features)]

        # Subset scale_data by *name*, not by the assay's row indices: it holds
        # its own subset of features, so `row_idx` addresses a different matrix.
        # Dropping the layer whenever it was not full-height — which is what
        # stood here, and which is every assay ScaleData has touched — silently
        # cost the subset its scaled values.
        sd = self.scale_data
        if is_matrix_empty(sd):
            new_sd = empty_dense(0, len(new_cells))
            new_scaled = []
        else:
            keep = set(new_features)
            sd_rows = [i for i, f in enumerate(self._scaled_features) if f in keep]
            new_sd = sd[sd_rows, :][:, col_idx]
            new_scaled = [self._scaled_features[i] for i in sd_rows]

        return Assay(
            counts=_sub(self.counts),
            data=_sub(self.data),
            scale_data=new_sd,
            scaled_features=new_scaled,
            feature_names=new_features,
            cell_names=new_cells,
            assay_orig=self.assay_orig,
            var_features=new_var,
            meta_features=new_meta,
            misc=dict(self.misc),
            key=self._key,
        )

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def merge(
        self,
        y: Union["Assay", list["Assay"]],
        add_cell_ids: Optional[list[str]] = None,
    ) -> "Assay":
        others = [y] if isinstance(y, Assay) else y
        all_assays = [self] + others

        if add_cell_ids is not None:
            if len(add_cell_ids) != len(all_assays):
                raise ValueError("add_cell_ids must have one entry per assay.")
            new_cell_names = []
            for prefix, assay in zip(add_cell_ids, all_assays):
                new_cell_names.extend([f"{prefix}_{c}" for c in assay._cell_names])
        else:
            new_cell_names = []
            for assay in all_assays:
                new_cell_names.extend(assay._cell_names)

        # Use shared features (intersection)
        shared = set(self._feature_names)
        for a in others:
            shared &= set(a._feature_names)
        shared_features = [f for f in self._feature_names if f in shared]

        def _stack_layer(layer_name):
            mats = []
            for a in all_assays:
                try:
                    m = a.layer_data(layer_name, features=shared_features)
                    mats.append(m)
                except Exception:
                    return None
            if any(m is None for m in mats):
                return None
            if sp.issparse(mats[0]):
                return sp.hstack(mats, format="csc")
            return np.hstack(mats)

        merged_counts = _stack_layer("counts")
        merged_data = _stack_layer("data")

        return Assay(
            counts=merged_counts,
            data=merged_data,
            feature_names=shared_features,
            cell_names=new_cell_names,
            key=self._key,
        )

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
        n_feat = len(self._feature_names)
        n_cells = len(self._cell_names)
        var = len(self.var_features)
        return (
            f"Assay (v3) with {n_feat} features and {n_cells} cells\n"
            f"  Key: {self._key!r}\n"
            f"  Variable features: {var}\n"
            f"  Layers: counts, data, scale_data"
        )


# ------------------------------------------------------------------
# Factory
# ------------------------------------------------------------------

def create_assay_object(
    counts=None,
    data=None,
    min_cells: int = 0,
    min_features: int = 0,
    feature_names: Optional[list[str]] = None,
    cell_names: Optional[list[str]] = None,
    key: str = "rna_",
) -> Assay:
    """Create an Assay, optionally filtering by min_cells / min_features."""
    matrix = counts if counts is not None else data

    if matrix is None:
        raise ValueError("Provide at least one of 'counts' or 'data'.")

    if sp.issparse(matrix):
        mat_csc = matrix.tocsc()
    else:
        mat_csc = sp.csc_matrix(np.asarray(matrix))

    n_features, n_cells = mat_csc.shape

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(n_features)]
    if cell_names is None:
        cell_names = [f"cell_{i}" for i in range(n_cells)]

    # Filter features by min_cells
    if min_cells > 0:
        cell_counts_per_feature = np.diff(mat_csc.T.tocsc().indptr)
        keep_feat = np.where(cell_counts_per_feature >= min_cells)[0]
        mat_csc = mat_csc[keep_feat, :]
        feature_names = [feature_names[i] for i in keep_feat]

    # Filter cells by min_features
    if min_features > 0:
        feat_counts_per_cell = np.diff(mat_csc.tocsc().indptr)
        keep_cell = np.where(feat_counts_per_cell >= min_features)[0]
        mat_csc = mat_csc[:, keep_cell]
        cell_names = [cell_names[i] for i in keep_cell]

    if counts is not None:
        return Assay(counts=mat_csc, feature_names=feature_names, cell_names=cell_names, key=key)
    return Assay(data=mat_csc, feature_names=feature_names, cell_names=cell_names, key=key)
