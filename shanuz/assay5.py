from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from typing import Optional, Type, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import as_sparse, empty_sparse, is_matrix_empty
from ._utils import validate_cell_names, validate_feature_names
from .lazy import is_lazy
from .logmap import LogMap
from .mixins import KeyMixin

#: Separator between a layer's stem and its split group, as in ``counts.batch1``.
#: Seurat's spelling, and not configurable — it is what ``Layers()`` returns and
#: what users write patterns against.
_SPLIT_SEP = "."


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
        "_layer_features",
        "_layer_cells",
        "default",
        "assay_orig",
        "meta_data",
        "misc",
        "_key",
        "_all_cell_names",
        "_all_feature_names",
        "_scaled_features",
        "_var_features",
        "_split_stems",
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
        layer_features: Optional[dict[str, list[str]]] = None,
        layer_cells: Optional[dict[str, list[str]]] = None,
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
        # Per-layer *ordered* feature / cell names. A layer may legitimately
        # span only a subset of the assay's features (e.g. scale.data holds
        # only the variable features, as in Seurat) or cells.
        self._layer_features: dict[str, list[str]] = {}
        self._layer_cells: dict[str, list[str]] = {}
        self._scaled_features: list[str] = []
        self._var_features: list[str] = []
        # Layer name -> the layer it was split from. See `_stem_of`.
        self._split_stems: dict[str, str] = {}

        layer_features = layer_features or {}
        layer_cells = layer_cells or {}
        for name, mat in layers.items():
            self._add_layer(
                name, mat,
                feature_names=layer_features.get(name),
                cell_names=layer_cells.get(name),
            )

        self.meta_data = (
            meta_data
            if meta_data is not None
            else pd.DataFrame(index=self._all_feature_names)
        )

    # ------------------------------------------------------------------
    # Internal layer management
    # ------------------------------------------------------------------

    def _add_layer(
        self,
        name: str,
        mat: Union[np.ndarray, sp.spmatrix],
        feature_names: Optional[list[str]] = None,
        cell_names: Optional[list[str]] = None,
    ) -> None:
        """Register a layer, optionally spanning a feature / cell subset.

        ``feature_names`` / ``cell_names`` give the *ordered* names for the
        layer's rows / columns. When omitted the layer is assumed to span all
        of the assay's features / cells (the common case for counts / data).
        """
        lf = list(feature_names) if feature_names is not None else list(self._all_feature_names)
        lc = list(cell_names) if cell_names is not None else list(self._all_cell_names)

        if mat.shape[0] != len(lf):
            raise ValueError(
                f"Layer '{name}' has {mat.shape[0]} rows but "
                f"{len(lf)} feature names were supplied."
            )
        if mat.shape[1] != len(lc):
            raise ValueError(
                f"Layer '{name}' has {mat.shape[1]} columns but "
                f"{len(lc)} cell names were supplied."
            )

        self.layers[name] = mat
        self._layer_features[name] = lf
        self._layer_cells[name] = lc

        feat_set = set(lf)
        cell_set = set(lc)
        self._features[name] = np.array(
            [f in feat_set for f in self._all_feature_names], dtype=bool
        )
        self._cells[name] = np.array(
            [c in cell_set for c in self._all_cell_names], dtype=bool
        )

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def cells(self, layer: Optional[str] = None) -> list[str]:
        if layer is not None and layer in self._layer_cells:
            return list(self._layer_cells[layer])
        return list(self._all_cell_names)

    def features(self, layer: Optional[str] = None) -> list[str]:
        if layer is not None and layer in self._layer_features:
            return list(self._layer_features[layer])
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
        # Index against the layer's *own* feature / cell names, which may be a
        # subset of the assay (e.g. scale.data holds only variable features).
        layer_feat = self._layer_features.get(layer, self._all_feature_names)
        layer_cell = self._layer_cells.get(layer, self._all_cell_names)

        row_idx = (
            [layer_feat.index(f) for f in features] if features is not None else slice(None)
        )
        col_idx = (
            [layer_cell.index(c) for c in cells] if cells is not None else slice(None)
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
        """Store (or replace) a layer.

        ``feature_names`` / ``cell_names`` declare which features / cells the
        matrix spans, so a layer may legitimately cover only a subset (e.g.
        scale.data over the variable features). When replacing an existing
        layer without supplying names, the previous span is reused.
        """
        if layer in self.layers:
            if feature_names is None:
                feature_names = self._layer_features.get(layer)
            if cell_names is None:
                cell_names = self._layer_cells.get(layer)
            del self.layers[layer]
        self._add_layer(layer, value, feature_names=feature_names, cell_names=cell_names)

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

    def _drop_layer(self, name: str) -> None:
        """Forget a layer and every index that mentions it."""
        self.layers.pop(name, None)
        self._layer_features.pop(name, None)
        self._layer_cells.pop(name, None)
        self._split_stems.pop(name, None)
        if name in self._features:
            del self._features[name]
        if name in self._cells:
            del self._cells[name]

    def _stem_of(self, name: str) -> Optional[str]:
        """The layer this one was split from, or ``None`` if it was not.

        Provenance is *recorded* by :meth:`split_layers` rather than parsed back
        out of the name, because the name cannot be parsed reliably: Seurat's own
        ``scale.data`` contains the separator, so splitting on the last ``.``
        would read it as layer ``scale`` of group ``data`` and rejoin it into
        something that never existed.
        """
        return self._split_stems.get(name)

    def join_layers(self, layers: Optional[list[str]] = None) -> "StdAssay":
        """Rejoin split layers, restoring the name, order and contents.

        Mirrors R's ``JoinLayers``. Each split *stem* is rejoined separately —
        ``counts.batch1`` and ``counts.batch2`` become ``counts`` again — and
        layers that were never split are left alone, which is what makes the
        no-argument call safe on a prepared assay that also holds ``data`` and a
        variable-features-only ``scale.data``.

        The rejoined columns come back in the **assay's** cell order, not in the
        order the split happened to produce. The assay's own cell vector never
        moved during the split, so anything else would leave the matrix silently
        transposed against the metadata that indexes it.
        """
        if layers is not None:
            # Explicit list: caller names the parts, so take the stem from the
            # recorded provenance and fall back to the shared prefix.
            parts = [n for n in layers if n in self.layers]
            if not parts:
                return self
            stems = {self._stem_of(n) for n in parts}
            stem = stems.pop() if len(stems) == 1 and None not in stems else None
            if stem is None:
                stem = os.path.commonprefix(parts).rstrip(_SPLIT_SEP)
            groups = {stem or parts[0]: parts}
        else:
            groups = {}
            for name in self.layers:
                stem = self._stem_of(name)
                if stem is not None:
                    groups.setdefault(stem, []).append(name)
        if not groups:
            return self

        new_obj = self._copy()
        for stem, parts in groups.items():
            features = list(self._layer_features.get(
                parts[0], self._all_feature_names))
            for part in parts[1:]:
                if list(self._layer_features.get(part, self._all_feature_names)) != features:
                    raise ValueError(
                        f"Cannot join layers {parts!r}: they do not span the "
                        f"same features."
                    )

            mats = [self.layers[p] for p in parts]
            combined = (sp.hstack(mats, format="csc") if sp.issparse(mats[0])
                        else np.hstack(mats))

            # Column j of `combined` belongs to concat_cells[j]; put them back
            # into the assay's order.
            concat_cells = [c for p in parts
                            for c in self._layer_cells.get(p, self._all_cell_names)]
            position = {c: j for j, c in enumerate(concat_cells)}
            ordered = [c for c in self._all_cell_names if c in position]
            combined = combined[:, [position[c] for c in ordered]]

            for part in parts:
                new_obj._drop_layer(part)
            new_obj._add_layer(stem, combined,
                               feature_names=features, cell_names=ordered)
        return new_obj

    def split_layers(self, f: list[str], layer: Optional[str] = None) -> "StdAssay":
        """Split one layer into per-group layers, as R's ``split()`` does.

        The parts are named ``<layer>.<group>`` — Seurat's spelling, which users
        match on with ``Layers(obj, pattern = "counts")`` — and each records the
        layer it came from so :meth:`join_layers` can put it back.
        """
        if layer is None:
            layer = self.default_layer
        if layer is None:
            raise ValueError("No layer to split.")
        mat = self.layers[layer]
        if len(f) != mat.shape[1]:
            raise ValueError("f must have one entry per cell.")

        groups: dict[str, list[int]] = {}
        for i, g in enumerate(f):
            groups.setdefault(str(g), []).append(i)

        new_obj = self._copy()
        base_feats = list(new_obj._layer_features.get(layer, new_obj._all_feature_names))
        base_cells = list(new_obj._layer_cells.get(layer, new_obj._all_cell_names))
        new_obj._drop_layer(layer)
        for g, idxs in groups.items():
            key = f"{layer}{_SPLIT_SEP}{g}"
            new_obj._add_layer(
                key, mat[:, idxs],
                feature_names=base_feats,
                cell_names=[base_cells[i] for i in idxs],
            )
            new_obj._split_stems[key] = layer
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
        feat_set = set(all_feat)
        cell_set = set(all_cells)

        # New global axes (preserve assay order, keep only requested members).
        new_features = (
            [f for f in features if f in feat_set] if features is not None else list(all_feat)
        )
        new_cells = (
            [c for c in cells if c in cell_set] if cells is not None else list(all_cells)
        )

        # Subset each layer against its *own* feature / cell span.
        new_layers: dict = {}
        new_layer_features: dict = {}
        new_layer_cells: dict = {}
        keep_feat = set(new_features)
        keep_cell = set(new_cells)
        for name, mat in self.layers.items():
            lf = self._layer_features.get(name, all_feat)
            lc = self._layer_cells.get(name, all_cells)
            fsel = [f for f in lf if f in keep_feat]
            csel = [c for c in lc if c in keep_cell]
            ridx = [lf.index(f) for f in fsel]
            cidx = [lc.index(c) for c in csel]
            new_layers[name] = mat[np.ix_(ridx, cidx)]
            new_layer_features[name] = fsel
            new_layer_cells[name] = csel

        new_meta = self.meta_data.reindex(new_features).copy()

        new = self.__class__(
            layers=new_layers,
            feature_names=new_features,
            cell_names=new_cells,
            assay_orig=self.assay_orig,
            meta_data=new_meta,
            misc=dict(self.misc),
            key=self._key,
            default=self.default,
            layer_features=new_layer_features,
            layer_cells=new_layer_cells,
        )
        new._var_features = [f for f in self._var_features if f in keep_feat]
        new._scaled_features = [f for f in self._scaled_features if f in keep_feat]
        # Carry the split provenance across: without it a subset of a split
        # assay can never be rejoined, and `join_layers` would report success
        # having done nothing, because it would find no stems to group.
        new._split_stems = {k: v for k, v in self._split_stems.items()
                            if k in new_layers}
        return new

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
        new_layer_features: dict = {}
        for layer_name in shared_layers:
            # A layer may span only a subset of the assay's features (e.g.
            # scale.data holds the variable features), so index against each
            # assay's *layer* feature names, not the assay-wide list.
            per_assay_feats = [
                a._layer_features.get(layer_name, a._all_feature_names)
                for a in all_assays
            ]
            common = set(per_assay_feats[0])
            for lf in per_assay_feats[1:]:
                common &= set(lf)
            layer_features = [f for f in per_assay_feats[0] if f in common]

            mats = []
            for a, lf in zip(all_assays, per_assay_feats):
                pos = {f: i for i, f in enumerate(lf)}
                feat_idx = [pos[f] for f in layer_features]
                mats.append(a.layers[layer_name][feat_idx, :])
            if sp.issparse(mats[0]):
                new_layers[layer_name] = sp.hstack(mats, format="csc")
            else:
                new_layers[layer_name] = np.hstack(mats)
            new_layer_features[layer_name] = layer_features

        merged = self.__class__(
            layers=new_layers,
            feature_names=shared_features,
            cell_names=new_cell_names,
            assay_orig=self.assay_orig,
            misc=dict(self.misc),
            key=self._key,
            layer_features=new_layer_features,
        )
        # Preserve the scaled-feature bookkeeping when every input agrees.
        scaled = new_layer_features.get("scale.data")
        if scaled is not None:
            merged._scaled_features = list(scaled)
        return merged

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
            layer_features={k: list(v) for k, v in self._layer_features.items()},
            layer_cells={k: list(v) for k, v in self._layer_cells.items()},
        )
        new._var_features = list(self._var_features)
        new._scaled_features = list(self._scaled_features)
        new._split_stems = dict(self._split_stems)
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

    if is_lazy(matrix):
        # Keep the on-disk layer on disk. `np.asarray` below would read the
        # whole store into a dense array here in the constructor, so the
        # obvious way to use the feature -- open a store, build an object
        # around it -- would end the laziness before any analysis began.
        mat = matrix
    elif sp.issparse(matrix):
        mat = matrix.tocsc()
    else:
        mat = sp.csc_matrix(np.asarray(matrix))

    n_features, n_cells = mat.shape

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(n_features)]
    if cell_names is None:
        cell_names = [f"cell_{i}" for i in range(n_cells)]

    # Filter features by min_cells. Either filter subsets the lazy store
    # through its own indexer, which yields an in-memory sparse block of just
    # the kept entries -- a filtered assay cannot stay on disk without writing
    # a second store, but it need never be dense to get there.
    if min_cells > 0:
        nnz_per_feat = (mat.nnz_per_row() if is_lazy(mat)
                        else np.diff(mat.T.tocsc().indptr))
        keep_feat = np.where(nnz_per_feat >= min_cells)[0]
        mat = mat[keep_feat, :]
        feature_names = [feature_names[i] for i in keep_feat]

    # Filter cells by min_features
    if min_features > 0:
        nnz_per_cell = (mat.nnz_per_col() if is_lazy(mat)
                        else np.diff(mat.tocsc().indptr))
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
