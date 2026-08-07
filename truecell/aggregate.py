"""Pseudobulk aggregation.

Mirrors Seurat's AggregateExpression(), which sums raw counts within groups of
cells (e.g. per cell type × donor) to produce a pseudobulk matrix — the standard
input for sample-level differential expression (DESeq2 / edgeR-style testing).
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .markers import _get_expression_matrix

# Metadata keys that resolve to the object's active identities rather than a
# meta_data column, matching Seurat's group.by = "ident".
_IDENT_KEYS = ("ident", "idents", "seurat_annotations")


def _group_labels(seurat, group_by: list[str]) -> pd.Series:
    """Per-cell group label: the group_by columns joined with '_' (Seurat's rule)."""
    cells = seurat.cell_names()
    meta = seurat.meta_data
    parts = pd.DataFrame(index=cells)
    for col in group_by:
        if col in _IDENT_KEYS and col not in meta.columns:
            parts[col] = [str(i) for i in seurat.idents]
        elif col in meta.columns:
            parts[col] = meta.loc[cells, col].astype(str).to_numpy()
        else:
            raise KeyError(
                f"group_by column {col!r} not found in meta_data "
                f"(columns: {list(meta.columns)})."
            )
    return parts.astype(str).agg("_".join, axis=1)


def aggregate_expression(
    seurat,
    group_by: Union[str, list[str]] = "ident",
    assays: Optional[Union[str, list[str]]] = None,
    features: Optional[list[str]] = None,
    layer: str = "counts",
    return_object: bool = False,
    normalization_method: str = "LogNormalize",
    scale_factor: float = 10000.0,
):
    """Sum counts within cell groups to form a pseudobulk profile.

    Mirrors R's ``AggregateExpression(obj, group.by = c("celltype", "donor"))``.

    Parameters
    ----------
    group_by      : metadata column(s) defining the groups. Multiple columns are
                    combined into a single label joined by ``"_"`` (as Seurat
                    does). ``"ident"`` uses the object's active identities.
    assays        : assay name(s) to aggregate (default: the active assay).
    features      : restrict to these features (default: all).
    layer         : layer to aggregate (default ``"counts"`` — pseudobulk is
                    defined on raw counts). Seurat's ``AggregateExpression`` has
                    no such argument and always sums ``counts``; this is a
                    superset with a matching default.
    return_object : if True, return a new :class:`Truecell` object with one "cell"
                    per group; if False (default), return a ``pd.DataFrame``
                    (features × groups), or a ``dict`` of them when several
                    assays are requested.
    normalization_method : how to fill the returned object's ``data`` layer when
                    ``return_object=True``. ``"LogNormalize"`` (Seurat's default)
                    or ``None`` to leave ``data`` as the raw sums.
    scale_factor  : the scale factor for that normalization (Seurat: 10000).

    Notes
    -----
    ``return_object=True`` **normalizes**, which is easy to miss and was wrong
    here until it was checked against R: Seurat's ``return.seurat = TRUE`` runs
    ``NormalizeData`` over the pseudobulk, so ``data`` holds
    ``log1p(sums / colSums × 10000)`` and not the sums. Leaving the sums in
    ``data`` — which is what this did — hands every downstream function that
    reads that layer un-normalized library-size-confounded values.

    This is the one place ``AggregateExpression`` and ``AverageExpression``
    diverge on their object output: :func:`average_expression` writes plain
    ``log1p`` of the averages, with no library-size step. Verified against
    Seurat 5.5.1 for both.

    Returns
    -------
    ``pd.DataFrame`` | ``dict[str, pd.DataFrame]`` | ``Truecell``
        A single DataFrame when one assay is aggregated, a dict keyed by assay
        name when several are, or a Truecell object when ``return_object=True``.
    """
    group_cols = [group_by] if isinstance(group_by, str) else list(group_by)
    if assays is None:
        assay_names = [seurat.active_assay]
    else:
        assay_names = [assays] if isinstance(assays, str) else list(assays)

    labels = _group_labels(seurat, group_cols)
    groups = sorted(labels.unique())
    n_cells = len(labels)
    n_groups = len(groups)

    # One-hot cells × groups indicator; counts(features×cells) @ indicator
    # sums each group's columns in a single sparse matmul.
    group_index = {g: j for j, g in enumerate(groups)}
    col_idx = np.fromiter((group_index[g] for g in labels), dtype=int, count=n_cells)
    indicator = sp.csr_matrix(
        (np.ones(n_cells), (np.arange(n_cells), col_idx)),
        shape=(n_cells, n_groups),
    )

    agg_frames: dict[str, pd.DataFrame] = {}
    for name in assay_names:
        assay_obj = seurat.assays[name]
        data, feature_names = _get_expression_matrix(assay_obj, layer)
        if features is not None:
            feat_set = set(features)
            keep = np.array([f in feat_set for f in feature_names])
            feature_names = [f for f, k in zip(feature_names, keep) if k]
            data = data[keep, :]
        summed = data @ indicator  # features × groups
        if sp.issparse(summed):
            summed = summed.toarray()
        agg_frames[name] = pd.DataFrame(
            np.asarray(summed), index=feature_names, columns=groups
        )

    if return_object:
        obj = _to_truecell(seurat, agg_frames, labels, groups, group_cols)
        if normalization_method is not None:
            from .preprocessing import normalize_data
            for name in assay_names:
                normalize_data(obj, assay=name,
                               normalization_method=normalization_method,
                               scale_factor=scale_factor)
        return obj

    if len(assay_names) == 1:
        return agg_frames[assay_names[0]]
    return agg_frames


def average_expression(
    seurat,
    group_by: Union[str, list[str]] = "ident",
    assays: Optional[Union[str, list[str]]] = None,
    features: Optional[list[str]] = None,
    layer: str = "data",
    return_object: bool = False,
):
    """Mean expression within cell groups.

    Mirrors R's ``AverageExpression(obj, group.by = "celltype")``, which is a
    different function from :func:`aggregate_expression` and not a rescaling of
    it. Two things differ, and both matter.

    **It averages rather than sums**, and **on the ``data`` layer it averages the
    back-transformed values**: ``mean(expm1(x))``, not ``mean(x)`` and not
    ``expm1(mean(x))``. The ``data`` layer holds log1p-normalized expression, so
    a mean taken on it is a geometric-ish mean of nothing in particular; Seurat
    undoes the log first, averages on the linear scale, and returns that. The
    difference is not cosmetic — on a small Poisson object the first gene reads
    **332.84** under ``AverageExpression`` against a count mean of **3.17**.

    The back-transform applies to the ``data`` layer **only**. ``counts`` and
    ``scale.data`` are not log-normalized, so those are averaged as they stand.
    Verified against Seurat 5.5.1 for all three layers.

    Parameters
    ----------
    group_by      : metadata column(s) defining the groups; ``"ident"`` uses the
                    object's active identities. Several are joined by ``"_"``.
    assays        : assay name(s) to average (default: the active assay).
    features      : restrict to these features (default: all).
    layer         : layer to average (default ``"data"``, as Seurat's is).
    return_object : if True, return a :class:`Truecell` with one "cell" per group.
                    Seurat puts the averaged matrix in that object's ``counts``
                    layer and ``log1p`` of it in ``data``; this does the same.

    Returns
    -------
    ``pd.DataFrame`` | ``dict[str, pd.DataFrame]`` | ``Truecell``

    See Also
    --------
    aggregate_expression : sums raw counts, which is what pseudobulk DE wants.
    """
    from .markers import expm1_keeping_sparsity

    group_cols = [group_by] if isinstance(group_by, str) else list(group_by)
    if assays is None:
        assay_names = [seurat.active_assay]
    else:
        assay_names = [assays] if isinstance(assays, str) else list(assays)

    labels = _group_labels(seurat, group_cols)
    groups = sorted(labels.unique())
    n_cells = len(labels)

    group_index = {g: j for j, g in enumerate(groups)}
    col_idx = np.fromiter((group_index[g] for g in labels), dtype=int, count=n_cells)
    indicator = sp.csr_matrix(
        (np.ones(n_cells), (np.arange(n_cells), col_idx)),
        shape=(n_cells, len(groups)),
    )
    # Divide by the group's own size, not by the mean size: the groups are
    # rarely balanced and Seurat averages within each.
    sizes = np.asarray(indicator.sum(axis=0)).ravel()

    avg_frames: dict[str, pd.DataFrame] = {}
    for name in assay_names:
        assay_obj = seurat.assays[name]
        data, feature_names = _get_expression_matrix(assay_obj, layer)
        if features is not None:
            feat_set = set(features)
            keep = np.array([f in feat_set for f in feature_names])
            feature_names = [f for f, k in zip(feature_names, keep) if k]
            data = data[keep, :]
        if layer == "data":
            data = expm1_keeping_sparsity(data)
        summed = data @ indicator
        if sp.issparse(summed):
            summed = summed.toarray()
        avg_frames[name] = pd.DataFrame(
            np.asarray(summed) / sizes, index=feature_names, columns=groups
        )

    if return_object:
        obj = _to_truecell(seurat, avg_frames, labels, groups, group_cols)
        # Seurat leaves the averaged values in `counts` and writes log1p of them
        # to `data`, rather than running NormalizeData over a matrix that is
        # already a per-group average.
        from .preprocessing import _set_layer
        for name in assay_names:
            _set_layer(obj.assays[name], "data",
                       np.log1p(avg_frames[name].to_numpy()))
        return obj

    if len(assay_names) == 1:
        return avg_frames[assay_names[0]]
    return avg_frames


def _to_truecell(seurat, agg_frames, labels, groups, group_cols):
    """Assemble aggregated frames into a Truecell object (one 'cell' per group)."""
    from .truecell import create_truecell_object

    # Decode one representative row of group_by values per group.
    decoded = pd.DataFrame(index=seurat.cell_names())
    for col in group_cols:
        if col in _IDENT_KEYS and col not in seurat.meta_data.columns:
            decoded[col] = [str(i) for i in seurat.idents]
        else:
            decoded[col] = seurat.meta_data.loc[decoded.index, col].astype(str).to_numpy()
    decoded["__group__"] = labels.to_numpy()
    group_meta = (
        decoded.drop_duplicates("__group__")
        .set_index("__group__")
        .loc[groups, group_cols]
    )

    primary = seurat.active_assay if seurat.active_assay in agg_frames else next(iter(agg_frames))
    base = agg_frames[primary]
    obj = create_truecell_object(
        counts=sp.csc_matrix(base.to_numpy()),
        assay=primary,
        feature_names=list(base.index),
        cell_names=groups,
        meta_data=group_meta,
        project=seurat.project_name,
    )
    # Attach any additional aggregated assays as extra Assay5 layers-of-record.
    for name, frame in agg_frames.items():
        if name == primary:
            continue
        from .assay5 import create_assay5_object

        obj.assays[name] = create_assay5_object(
            counts=sp.csc_matrix(frame.to_numpy()),
            feature_names=list(frame.index),
            cell_names=groups,
            key=f"{name.lower()}_",
        )
    return obj
