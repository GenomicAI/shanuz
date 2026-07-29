"""Reference mapping — Seurat's ``FindTransferAnchors`` / ``TransferData``.

Integration (``anchors.py``) *corrects* several datasets onto each other so they
share a coordinate system. Reference mapping is the asymmetric cousin: the
reference is an annotated atlas you trust and never touch, and the query is a
new, unlabelled dataset you want to *annotate* by borrowing from the reference.

The anchor machinery is the same as integration — build a shared space, find
mutual nearest neighbours, score and filter them — but two things differ:

  1. **The space is directional.** By default (``reduction="pcaproject"``) the
     query is *projected into the reference's PCA* rather than a jointly-learned
     CCA space: the reference's principal axes are computed once on the shared
     anchor features and the query is pushed through the same loadings. Because
     those axes were learned without ever seeing the query, batch-specific
     structure the reference doesn't have simply lands nowhere — which is what
     makes projection robust for annotation. ``reduction="cca"`` is also
     supported for the harder cross-modality / cross-species cases.
  2. **The reference is fixed.** Anchors run *reference → query* and are used by
     :func:`transfer_data` to carry labels (or expression) across, not to move
     any cells.

:func:`transfer_data` turns the anchors into a per-query-cell weight over the
anchors (the same distance-weighted, score-scaled Gaussian kernel that
``IntegrateData`` uses) and then either

  * **classifies** — a categorical reference label becomes a one-hot matrix over
    the anchor reference cells; ``weights @ onehot`` gives each query cell a
    probability per class (``predicted.id`` = argmax, ``prediction.score.*`` =
    the probabilities), or
  * **imputes** — a continuous ``features × reference-cells`` matrix is carried
    across the same weights to predict query expression.

This reuses :class:`truecell.anchors` end to end; only the projection and the
label/expression transfer are new. ``MapQuery`` / ``ProjectUMAP`` (placing the
query in the reference UMAP) build on top of this and land separately.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd

from .anchors import (
    _anchor_feature_matrix,
    _anchor_weights,
    _cca,
    _filter_anchors,
    _integration_features,
    _l2_normalize_cols,
    _l2_normalize_rows,
    _mutual_nn_from,
    _neighbor_sets,
    _pca_loadings,
    _score_anchors,
)

REDUCTIONS = ("pcaproject", "cca")


# ----------------------------------------------------------------------
# Result container
# ----------------------------------------------------------------------


class TransferAnchors:
    """Anchors linking a fixed reference to a query, for label/data transfer.

    Slots
    -----
    - ``anchors`` — DataFrame with columns ``cell1, cell2, score``.
      ``cell1`` is a *within-reference* 0-based cell index,
      ``cell2`` a *within-query* one; the reference is never
      moved.
    - ``reference`` — the reference Truecell object (annotated atlas).
    - ``query`` — the query Truecell object (to be annotated).
    - ``reduction`` — ``"pcaproject"`` or ``"cca"`` — how the shared space was
      built.
    - ``anchor_features`` — the features the anchors run on.
    - ``dims`` — number of shared dimensions used.
    - ``query_embedding`` — ``(n_query_cells × dims)`` — the query cells in the shared
      space, used by :func:`transfer_data` to weight anchors.
    """

    __slots__ = (
        "anchors",
        "reference",
        "query",
        "reduction",
        "anchor_features",
        "dims",
        "query_embedding",
    )

    def __init__(
        self,
        anchors: pd.DataFrame,
        reference,
        query,
        reduction: str,
        anchor_features: list[str],
        dims: int,
        query_embedding: np.ndarray,
    ) -> None:
        self.anchors = anchors
        self.reference = reference
        self.query = query
        self.reduction = reduction
        self.anchor_features = list(anchor_features)
        self.dims = dims
        self.query_embedding = query_embedding

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"TransferAnchors: {len(self.anchors)} anchors "
            f"(reference {len(self.reference)} cells → query "
            f"{len(self.query)} cells)\n"
            f"  reduction={self.reduction!r}  dims={self.dims}  "
            f"features={len(self.anchor_features)}"
        )


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def find_transfer_anchors(
    reference,
    query,
    anchor_features: Optional[list[str]] = None,
    reduction: str = "pcaproject",
    dims: int = 30,
    k_anchor: int = 5,
    k_filter: int = 200,
    k_score: int = 30,
    layer: str = "scale.data",
    seed: int = 42,
) -> TransferAnchors:
    """Find transfer anchors from a reference to a query (Seurat's ``FindTransferAnchors``).

    Mirrors ``FindTransferAnchors(reference, query, reduction = "pcaproject")``.
    Anchors are mutual nearest neighbours between the reference and the query in
    a shared space — by default the reference's own PCA, into which the query is
    projected — scored by neighbourhood consistency and filtered against the
    original expression space.

    Parameters
    ----------
    reference       : annotated reference Truecell object (normalized + with
                      variable features / scaled data).
    query           : query Truecell object to annotate (same preprocessing).
    anchor_features : features to anchor on (default: variable features shared
                      by reference and query).
    reduction       : ``"pcaproject"`` (project the query into the reference's
                      PCA; the default and most robust for annotation) or
                      ``"cca"`` (a jointly-learned space, for harder
                      cross-modality/species cases). The reference PCA is
                      computed on the shared anchor features, so the reference
                      need not already carry a ``pca`` reduction.
    dims            : number of shared dimensions to use.
    k_anchor        : neighbours for the mutual-nearest-neighbour search.
    k_filter        : neighbourhood size for the feature-space anchor filter
                      (set to 0 or None to skip filtering).
    k_score         : neighbourhood size for anchor scoring.
    layer           : layer to draw the shared-space expression from.
    seed            : random seed for the PCA/neighbour steps.

    Returns
    -------
    TransferAnchors
    """
    reduction = reduction.lower()
    if reduction not in REDUCTIONS:
        raise ValueError(
            f"Unknown reduction {reduction!r}. Supported: {REDUCTIONS}."
        )

    features = _integration_features([reference, query], anchor_features)

    ref_scaled = _l2_normalize_cols(
        _anchor_feature_matrix(reference, features, layer)
    )
    query_scaled = _l2_normalize_cols(
        _anchor_feature_matrix(query, features, layer)
    )
    if ref_scaled.shape[0] != query_scaled.shape[0]:
        raise ValueError(
            "Reference and query disagree on the number of anchor features "
            f"({ref_scaled.shape[0]} vs {query_scaled.shape[0]}); ensure the "
            "anchor features are scaled in both objects."
        )

    used = min(
        dims,
        ref_scaled.shape[1] - 1,
        query_scaled.shape[1] - 1,
        ref_scaled.shape[0],
    )

    if reduction == "pcaproject":
        # Project the query through the reference's PCA loadings — both datasets
        # end up in the reference's principal-component space.
        load_ref = _pca_loadings(ref_scaled, used, seed=seed)
        ref_emb = _l2_normalize_rows(ref_scaled.T @ load_ref)
        query_emb = _l2_normalize_rows(query_scaled.T @ load_ref)
    else:  # cca — a jointly-learned shared space
        ref_emb, query_emb, _loadings = _cca(ref_scaled, query_scaled, used)

    # FindTransferAnchors goes through the same FindAnchors internals as
    # FindIntegrationAnchors, so it gets the same four-way neighbour tables.
    nbrs = _neighbor_sets(ref_emb, query_emb, max(k_anchor, k_score))
    pairs = _mutual_nn_from(nbrs, k_anchor)
    n_ref = ref_emb.shape[0]
    # Filter first: the score is rescaled against the percentiles of the set it
    # is handed, so scoring before filtering shifts every surviving score.
    if k_filter:
        pairs = _filter_anchors(pairs, ref_scaled, query_scaled, k_filter)
    scores = _score_anchors(pairs, nbrs, n_ref, k_score)

    rows = [(int(i), int(j), float(s)) for (i, j), s in zip(pairs, scores)]
    anchors = pd.DataFrame(rows, columns=["cell1", "cell2", "score"])

    return TransferAnchors(
        anchors=anchors,
        reference=reference,
        query=query,
        reduction=reduction,
        anchor_features=features,
        dims=used,
        query_embedding=query_emb,
    )


def transfer_data(
    anchors: TransferAnchors,
    refdata: Union[str, np.ndarray, list, pd.Series],
    k_weight: int = 50,
    sd_weight: float = 1.0,
    refdata_features: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Transfer labels or expression from reference to query (Seurat's ``TransferData``).

    Mirrors ``TransferData(anchorset, refdata = "celltype")``. Each query cell
    gets a weight over the anchors — a distance-weighted, anchor-score-scaled
    Gaussian kernel in the shared space (the same weighting ``IntegrateData``
    uses) — and the reference information is carried across those weights.

    ``refdata`` selects the mode:

    * **classification** — a metadata column name (``str``) or a 1-D array of
      per-reference-cell labels. Returns a DataFrame indexed by query cell with
      ``predicted.id``, one ``prediction.score.<class>`` column per reference
      class (each query cell's rows sum to 1), and ``prediction.score.max``.
    * **imputation** — a 2-D ``features × reference-cells`` matrix. Returns a
      DataFrame of predicted query expression (features × query cells); name the
      rows with ``refdata_features``.

    Parameters
    ----------
    anchors          : a :class:`TransferAnchors` from
                       :func:`find_transfer_anchors`.
    refdata          : reference labels (str column / 1-D array) or a 2-D
                       ``features × reference-cells`` matrix to impute.
    k_weight         : anchors used to weight each query cell.
    sd_weight        : bandwidth multiplier for the Gaussian anchor kernel.
    refdata_features : row names for the imputation output (2-D ``refdata``).

    Returns
    -------
    pandas.DataFrame
    """
    anchor_df = anchors.anchors
    if len(anchor_df) == 0:
        raise ValueError(
            "No anchors between reference and query; cannot transfer. Try a "
            "larger k_anchor or k_filter=0 in find_transfer_anchors."
        )

    query_cells = anchors.query.cell_names()
    ref_cell1 = anchor_df["cell1"].to_numpy()
    query_cell2 = anchor_df["cell2"].to_numpy()
    anchor_scores = anchor_df["score"].to_numpy()

    weights = _anchor_weight_matrix(
        anchors.query_embedding, query_cell2, anchor_scores, k_weight, sd_weight
    )  # (n_query × n_anchor)

    classification = isinstance(refdata, str) or np.ndim(refdata) == 1
    if classification:
        return _transfer_labels(
            anchors.reference, refdata, ref_cell1, weights, query_cells
        )
    return _transfer_expression(
        refdata, refdata_features, ref_cell1, weights, query_cells
    )


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _anchor_weight_matrix(
    query_emb: np.ndarray,
    anchor_query_idx: np.ndarray,
    anchor_scores: np.ndarray,
    k_weight: int,
    sd_weight: float,
) -> np.ndarray:
    """Per-query-cell weight over the anchors (``n_query × n_anchor``).

    ``TransferData`` weights through the same ``FindWeights`` as
    ``IntegrateData``; only the space differs (the query projected into the
    reference, rather than a PCA of the merged pair). So this is
    :func:`~truecell.anchors._anchor_weights` transposed.
    """
    return _anchor_weights(
        query_emb, anchor_query_idx, anchor_scores, k_weight, sd_weight
    ).T


def _transfer_labels(
    reference,
    refdata: Union[str, np.ndarray, list, pd.Series],
    ref_cell1: np.ndarray,
    weights: np.ndarray,
    query_cells: list[str],
) -> pd.DataFrame:
    """Categorical label transfer → per-class prediction scores + argmax id."""
    if isinstance(refdata, str):
        if refdata not in reference.meta_data.columns:
            raise KeyError(
                f"refdata column {refdata!r} not found in reference meta_data."
            )
        labels = reference.meta_data.loc[reference.cell_names(), refdata].to_numpy()
    else:
        labels = np.asarray(refdata)
        if labels.shape[0] != len(reference):
            raise ValueError(
                f"refdata has {labels.shape[0]} labels but the reference has "
                f"{len(reference)} cells."
            )

    labels = labels.astype(object)
    classes = sorted({str(x) for x in labels})
    class_index = {c: i for i, c in enumerate(classes)}

    # One-hot the anchor reference cells' labels: (n_anchor × n_classes).
    anchor_labels = labels[ref_cell1]
    onehot = np.zeros((len(anchor_labels), len(classes)))
    onehot[np.arange(len(anchor_labels)), [class_index[str(x)] for x in anchor_labels]] = 1.0

    scores = weights @ onehot  # (n_query × n_classes), rows sum to 1
    predicted = [classes[i] for i in scores.argmax(axis=1)]
    score_max = scores.max(axis=1)

    out = pd.DataFrame(
        {"predicted.id": predicted}, index=list(query_cells)
    )
    for c in classes:
        out[f"prediction.score.{c}"] = scores[:, class_index[c]]
    out["prediction.score.max"] = score_max
    return out


def _transfer_expression(
    refdata,
    refdata_features: Optional[list[str]],
    ref_cell1: np.ndarray,
    weights: np.ndarray,
    query_cells: list[str],
) -> pd.DataFrame:
    """Continuous imputation → predicted query expression (features × cells)."""
    refmat = np.asarray(refdata, dtype=float)
    if refmat.ndim != 2:
        raise ValueError("Continuous refdata must be a 2-D features × cells matrix.")

    n_features = refmat.shape[0]
    anchor_ref = refmat[:, ref_cell1]          # features × n_anchor
    imputed = anchor_ref @ weights.T           # features × n_query

    if refdata_features is None:
        refdata_features = [f"feature_{i}" for i in range(n_features)]
    elif len(refdata_features) != n_features:
        raise ValueError(
            f"refdata_features has {len(refdata_features)} names but refdata has "
            f"{n_features} rows."
        )

    return pd.DataFrame(
        imputed, index=list(refdata_features), columns=list(query_cells)
    )
