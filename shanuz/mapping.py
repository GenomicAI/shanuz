"""Reference mapping, part two — Seurat's ``ProjectUMAP`` / ``MapQuery``.

``transfer.py`` gets you *anchors* between a fixed reference and a query and
carries labels or expression across them. This module does the last thing you
usually want from a reference: **place the query in the reference's own UMAP**,
so the new cells land on the atlas you already know how to read.

The trick is that a UMAP is not an equation you can apply to a new point — it is
a stochastic optimisation. ``umap-learn`` solves this by keeping the fitted
model around: :func:`shanuz.run_umap` stashes it in
``reduction.misc["umap_model"]``, and its ``.transform`` places new points into
the *existing* embedding by pinning them near the training points they resemble
and running a short, held-fixed optimisation. So projection is a two-step map:

  1. **query → reference PCA.** The reference's PCA is a plain linear map from
     genes to components (``embeddings = scaledᵀ @ loadings``). Push the query's
     scaled expression through the *reference's* loadings and the query lands in
     the same principal-component space the reference UMAP was trained on — the
     same "project into a space the query never helped define" logic that makes
     :func:`shanuz.find_transfer_anchors`'s ``pcaproject`` robust.
  2. **reference PCA → reference UMAP.** Run the reference's fitted UMAP model in
     transform-only mode on those projected coordinates.

:func:`project_umap` is that two-step map on its own. :func:`map_query` is the
Seurat ``MapQuery`` convenience that composes the whole reference-mapping
workflow: take the anchors, :func:`~shanuz.transfer_data` the labels onto the
query's metadata, and :func:`project_umap` the query into the reference UMAP —
one call from anchors to an annotated, atlas-placed query.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .dimreduc import DimReduc
from .transfer import TransferAnchors, transfer_data


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def project_umap(
    query,
    reference,
    reduction: str = "pca",
    umap_reduction: str = "umap",
    dims: Optional[Union[list[int], range]] = None,
    reduction_name: str = "ref.umap",
    reduction_key: str = "refUMAP_",
    layer: str = "scale.data",
) -> DimReduc:
    """Project a query into a reference's UMAP (Seurat's ``ProjectUMAP``).

    Mirrors ``ProjectUMAP(query, reference, reduction.model = "umap")``. The
    query is first projected into the reference's PCA (through the reference's
    loadings), then run through the reference's *fitted* UMAP model in
    transform-only mode — so the query cells land in the reference's existing
    embedding rather than in a fresh, unrelated one. Stores the result as
    ``query.reductions[reduction_name]`` and returns it.

    The reference must already carry both a PCA reduction (``reduction``, for the
    loadings) and a UMAP reduction (``umap_reduction``) fitted with a returnable
    model — i.e. ``run_umap(reference)``, which stashes the ``umap-learn`` model
    in ``reduction.misc["umap_model"]``.

    Parameters
    ----------
    query          : query Shanuz object (normalized + scaled on the reference's
                     PCA features).
    reference      : reference Shanuz object carrying the fitted PCA + UMAP.
    reduction      : reference reduction whose loadings project the query
                     (default ``"pca"``).
    umap_reduction : reference reduction holding the fitted UMAP model (default
                     ``"umap"``).
    dims           : which PCA dimensions to feed the UMAP model (0-indexed).
                     Defaults to the dimensions the model was trained on.
    reduction_name : storage key for the projection in ``query.reductions``.
    reduction_key  : prefix for the projected dimension names.
    layer          : layer to draw the query's expression from.

    Returns
    -------
    DimReduc
        The query's cells in the reference UMAP; also stored on ``query``.
    """
    if reduction not in reference.reductions:
        raise KeyError(
            f"Reference reduction {reduction!r} not found; run run_pca(reference) "
            "first."
        )
    if umap_reduction not in reference.reductions:
        raise KeyError(
            f"Reference reduction {umap_reduction!r} not found; run "
            "run_umap(reference) first."
        )

    umap_dr = reference.reductions[umap_reduction]
    model = umap_dr.misc.get("umap_model")
    if model is None:
        raise ValueError(
            f"Reference reduction {umap_reduction!r} has no fitted UMAP model in "
            "misc['umap_model']; run_umap(reference) must be run from a reduction "
            "(not a graph) so the model is stored for transform-only projection."
        )

    query_pca = _project_into_reference_pca(query, reference, reduction, layer)
    query_pca = _select_model_dims(query_pca, model, dims)

    umap_coords = np.asarray(model.transform(query_pca))

    dim_names = [f"{reduction_key}{i + 1}" for i in range(umap_coords.shape[1])]
    projected = DimReduc(
        cell_embeddings=umap_coords,
        cell_names=query.cell_names(),
        feature_names=dim_names,
        assay_used=query.active_assay,
        key=reduction_key,
        misc={
            "projected_from": umap_reduction,
            "reference_reduction": reduction,
            "query_pca": query_pca,
        },
    )
    query.reductions[reduction_name] = projected
    return projected


def map_query(
    anchors: TransferAnchors,
    refdata: Optional[Union[str, np.ndarray, list, pd.Series]] = None,
    reference_reduction: str = "pca",
    reduction_model: str = "umap",
    reduction_name: str = "ref.umap",
    reduction_key: str = "refUMAP_",
    k_weight: int = 50,
    sd_weight: float = 1.0,
    refdata_features: Optional[list[str]] = None,
    layer: str = "scale.data",
) -> Optional[pd.DataFrame]:
    """Annotate and place a query in a reference (Seurat's ``MapQuery``).

    Mirrors ``MapQuery(anchorset, query, reference, refdata = "celltype")``. The
    single call that turns transfer anchors into a mapped query:

    1. :func:`~shanuz.transfer_data` carries ``refdata`` across the anchors, and
       — for a categorical label — writes ``predicted.id`` /
       ``prediction.score.*`` straight onto ``query.meta_data``.
    2. :func:`project_umap` projects the query into the reference's UMAP, stored
       as ``query.reductions[reduction_name]``.

    Both steps mutate the query object (the anchors' ``query``) in place.

    Parameters
    ----------
    anchors             : a :class:`~shanuz.TransferAnchors` from
                          :func:`~shanuz.find_transfer_anchors`.
    refdata             : reference labels to transfer (metadata column name or a
                          per-reference-cell array) or a 2-D
                          ``features × reference-cells`` matrix to impute. Pass
                          ``None`` to skip transfer and only project the UMAP.
    reference_reduction : reference reduction whose loadings project the query
                          into the UMAP's input space (default ``"pca"``).
    reduction_model     : reference reduction holding the fitted UMAP model.
    reduction_name      : storage key for the projection in ``query.reductions``.
    reduction_key       : prefix for the projected dimension names.
    k_weight, sd_weight : anchor-weighting knobs passed to
                          :func:`~shanuz.transfer_data`.
    refdata_features    : row names for imputation output (2-D ``refdata``).
    layer               : layer to draw the query's expression from.

    Returns
    -------
    pandas.DataFrame or None
        The transferred predictions (classification) or imputed expression
        (imputation); ``None`` when ``refdata`` is not given.
    """
    query = anchors.query

    predictions: Optional[pd.DataFrame] = None
    if refdata is not None:
        predictions = transfer_data(
            anchors,
            refdata,
            k_weight=k_weight,
            sd_weight=sd_weight,
            refdata_features=refdata_features,
        )
        classification = isinstance(refdata, str) or np.ndim(refdata) == 1
        if classification:
            # Write predicted.id / prediction.score.* onto the query metadata,
            # as Seurat's MapQuery does. Imputation output is returned, not stored.
            for col in predictions.columns:
                query.meta_data[col] = (
                    predictions[col].reindex(query.meta_data.index).to_numpy()
                )

    project_umap(
        query,
        anchors.reference,
        reduction=reference_reduction,
        umap_reduction=reduction_model,
        reduction_name=reduction_name,
        reduction_key=reduction_key,
        layer=layer,
    )
    return predictions


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _project_into_reference_pca(query, reference, reduction: str, layer: str) -> np.ndarray:
    """Project the query into the reference's PCA (``n_query × n_pcs``).

    The reference PCA is a linear map ``scaledᵀ @ loadings``; pushing the query's
    scaled expression through the *reference's* loadings lands it in the same
    principal-component space. Only the reference PCA features the query actually
    carries (scaled) are used, and the loadings are subset to match, so a query
    missing a few of the reference's variable features still projects cleanly.
    """
    ref_pca = reference.reductions[reduction]
    ref_feats = ref_pca.features()
    loadings = np.asarray(ref_pca.feature_loadings)  # (n_ref_feats × n_pcs)
    if loadings.shape[0] != len(ref_feats):
        raise ValueError(
            f"Reference reduction {reduction!r} carries no feature loadings; it "
            "cannot project a query. Use run_pca(reference)."
        )

    assay = query.get_assay()
    scaled_here = set(_scaled_feature_names(assay, layer))
    feat_pos = {f: i for i, f in enumerate(ref_feats)}
    keep = [f for f in ref_feats if f in scaled_here]
    if not keep:
        raise ValueError(
            "The query shares none of the reference PCA's features (scaled); "
            "ensure both were scaled on a common feature set."
        )

    loadings_kept = loadings[[feat_pos[f] for f in keep], :]  # (n_keep × n_pcs)
    query_scaled = assay.layer_data(layer, features=keep)     # (n_keep × n_query)
    if sp.issparse(query_scaled):
        query_scaled = query_scaled.toarray()
    query_scaled = np.asarray(query_scaled, dtype=float)

    return query_scaled.T @ loadings_kept                     # (n_query × n_pcs)


def _select_model_dims(
    query_pca: np.ndarray, model, dims: Optional[Union[list[int], range]]
) -> np.ndarray:
    """Take the PCA dimensions the fitted UMAP model expects as input.

    With ``dims`` given, use exactly those columns. Otherwise default to the
    dimensions the model was trained on (``run_umap``'s default is every PCA
    dimension, so this is usually a no-op), inferred from the fitted model when
    ``umap-learn`` exposes the training data.
    """
    if dims is not None:
        return query_pca[:, list(dims)]

    trained = getattr(model, "_raw_data", None)
    if trained is not None and getattr(trained, "ndim", 0) == 2:
        n_in = trained.shape[1]
        if n_in <= query_pca.shape[1]:
            return query_pca[:, :n_in]
    return query_pca


def _scaled_feature_names(assay, layer: str) -> list[str]:
    """Which features the assay actually carries in its scaled layer."""
    from .assay5 import Assay5

    if isinstance(assay, Assay5):
        return list(assay._layer_features.get(layer, assay._all_feature_names))
    # Classic Assay v3 stores scale_data over the assay's feature list.
    return list(assay._feature_names)
