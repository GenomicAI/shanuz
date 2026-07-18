"""Anchor-based dataset integration (Seurat's CCA / RPCA ``IntegrateData``).

Harmony (``integration.py``) corrects an *existing* joint embedding. The
anchor approach of Stuart et al. (2019) works the other way round: it never
assumes the datasets share a coordinate system to begin with. Instead it

  1. builds a *shared* low-dimensional space for a pair of datasets — either by
     **canonical correlation analysis** (``reduction="cca"``: the SVD of the
     cross-covariance ``AᵀB``, whose singular vectors are the directions along
     which the two datasets co-vary most strongly) or by **reciprocal PCA**
     (``reduction="rpca"``: project each dataset into the *other's* PCA space);
  2. finds **mutual nearest neighbours** in that space — cell *i* of dataset A
     and cell *j* of dataset B are an *anchor* only if each is among the other's
     k nearest neighbours. A mutual pair is evidence the two cells are the same
     biological state seen in two batches;
  3. **scores** every anchor by how much of its neighbourhood the two members
     share (a consistent anchor sits in a coherent local structure) and
     **filters** ones that are not even near each other in the original
     expression space;
  4. **corrects** each query dataset onto the reference by adding, to every
     query cell, a distance-weighted average of the anchor *correction vectors*
     ``expr_ref − expr_query`` — pulling matched populations on top of each
     other while leaving genuinely reference-only structure alone.

The output of :func:`integrate_data` is a merged object carrying an
``"integrated"`` assay whose ``data`` is the batch-corrected expression of the
anchor features — exactly what you then ``scale_data`` + ``run_pca`` on to get
an embedding that clusters by cell type rather than by batch.

Only the anchor pairs and the reference-facing bookkeeping are Seurat-specific;
the same :class:`IntegrationAnchors` object is what v0.3.0's reference mapping
(``FindTransferAnchors`` / ``TransferData``) is built to reuse.

This is a *reference-based* implementation: anchors are found between the
reference (``reference=0`` by default) and each other dataset, and every other
dataset is corrected onto the reference. That is one of Seurat's supported
integration modes and keeps the guide-tree bookkeeping out of the first cut.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

REDUCTIONS = ("cca", "rpca")


# ----------------------------------------------------------------------
# Result container
# ----------------------------------------------------------------------


class IntegrationAnchors:
    """Anchors linking a reference dataset to one or more query datasets.

    Slots
    -----
    anchors         : DataFrame with columns ``dataset1, cell1, dataset2,
                      cell2, score``. ``dataset1`` is always the reference; the
                      cell columns hold *within-dataset* 0-based row indices.
    objects         : the list of Shanuz objects passed to
                      :func:`find_integration_anchors` (order preserved).
    reference       : index into ``objects`` of the reference dataset.
    reduction       : ``"cca"`` or ``"rpca"`` — how the shared space was built.
    anchor_features : the features the anchors (and correction) run on.
    dims            : number of shared dimensions used.
    weight_embeddings : ``{query_index: (n_query_cells × dims) array}`` — each
                      query dataset's cells in the shared space, used by
                      :func:`integrate_data` to weight the correction.
    """

    __slots__ = (
        "anchors",
        "objects",
        "reference",
        "reduction",
        "anchor_features",
        "dims",
        "weight_embeddings",
    )

    def __init__(
        self,
        anchors: pd.DataFrame,
        objects: list,
        reference: int,
        reduction: str,
        anchor_features: list[str],
        dims: int,
        weight_embeddings: dict[int, np.ndarray],
    ) -> None:
        self.anchors = anchors
        self.objects = objects
        self.reference = reference
        self.reduction = reduction
        self.anchor_features = list(anchor_features)
        self.dims = dims
        self.weight_embeddings = weight_embeddings

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"IntegrationAnchors: {len(self.anchors)} anchors across "
            f"{len(self.objects)} datasets\n"
            f"  reduction={self.reduction!r}  reference={self.reference}  "
            f"dims={self.dims}  features={len(self.anchor_features)}"
        )


# ----------------------------------------------------------------------
# Linear algebra helpers
# ----------------------------------------------------------------------


def _l2_normalize_rows(mat: np.ndarray) -> np.ndarray:
    """L2-normalize each row of ``mat`` (a per-cell embedding)."""
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return mat / norms


def _l2_normalize_cols(mat: np.ndarray) -> np.ndarray:
    """L2-normalize each column of ``mat`` (features × cells → per-cell unit)."""
    norms = np.linalg.norm(mat, axis=0, keepdims=True)
    norms[norms == 0] = 1.0
    return mat / norms


def _cca(A: np.ndarray, B: np.ndarray, dims: int) -> tuple[np.ndarray, np.ndarray]:
    """Shared CCA embedding of two datasets.

    ``A`` (features × cells_a) and ``B`` (features × cells_b) are L2-normalized
    per cell. The SVD of the cross-covariance ``AᵀB`` (cells_a × cells_b) yields
    left/right singular vectors that place the two datasets' cells in one space;
    each cell's coordinates are then L2-normalized (Seurat's ``L2Dim``).
    """
    cross = A.T @ B
    U, _s, Vt = np.linalg.svd(cross, full_matrices=False)
    d = min(dims, U.shape[1])
    emb_a = _l2_normalize_rows(U[:, :d])
    emb_b = _l2_normalize_rows(Vt[:d, :].T)
    return emb_a, emb_b


def _pca_loadings(mat: np.ndarray, dims: int, seed: int = 42) -> np.ndarray:
    """Top-``dims`` PCA loadings (features × dims) of a features × cells matrix."""
    from sklearn.decomposition import PCA

    d = min(dims, min(mat.shape) - 1)
    pca = PCA(n_components=max(d, 1), random_state=seed)
    pca.fit(mat.T)  # samples (cells) × features
    return pca.components_.T  # features × dims


# ----------------------------------------------------------------------
# Neighbours / MNN
# ----------------------------------------------------------------------


def _nearest(query: np.ndarray, reference: np.ndarray, k: int) -> np.ndarray:
    """For each row of ``query``, the indices of its ``k`` nearest ``reference`` rows."""
    from sklearn.neighbors import NearestNeighbors

    k = min(k, reference.shape[0])
    nn = NearestNeighbors(n_neighbors=k).fit(reference)
    return nn.kneighbors(query, return_distance=False)


def _mutual_nn(
    emb_a_for_b: np.ndarray,
    emb_b_for_b: np.ndarray,
    emb_b_for_a: np.ndarray,
    emb_a_for_a: np.ndarray,
    k: int,
) -> list[tuple[int, int]]:
    """Mutual nearest neighbours between datasets A and B.

    ``emb_*_for_b`` are the coordinates used to find *B*-cells around each
    *A*-cell (i.e. both datasets expressed in the space where B lives), and
    ``emb_*_for_a`` the reverse. For CCA the two spaces coincide; for RPCA they
    are the two reciprocal PCA projections.
    """
    a_to_b = _nearest(emb_a_for_b, emb_b_for_b, k)  # A-cell → nearby B-cells
    b_to_a = _nearest(emb_b_for_a, emb_a_for_a, k)  # B-cell → nearby A-cells
    b_sets = [set(row) for row in b_to_a]
    pairs = []
    for i, neigh in enumerate(a_to_b):
        for j in neigh:
            if i in b_sets[j]:
                pairs.append((int(i), int(j)))
    return pairs


# ----------------------------------------------------------------------
# Scoring & filtering
# ----------------------------------------------------------------------


def _score_anchors(
    pairs: list[tuple[int, int]],
    combined: np.ndarray,
    n_a: int,
    k_score: int,
) -> np.ndarray:
    """Shared-neighbourhood score in [0, 1] for each anchor.

    ``combined`` stacks the two datasets ([A; B]) in a common space. An anchor
    ``(i, j)`` scores high when the reference cell ``i`` and query cell
    ``j = n_a + j`` share many of their ``k_score`` nearest neighbours — i.e.
    the anchor sits in locally consistent structure. Scores are rescaled to
    [0, 1] with Seurat's 1st/90th-percentile clamp.
    """
    if not pairs:
        return np.array([])
    knn = _nearest(combined, combined, min(k_score, combined.shape[0]))
    neigh_sets = [set(row) for row in knn]
    raw = np.empty(len(pairs))
    for m, (i, j) in enumerate(pairs):
        raw[m] = len(neigh_sets[i] & neigh_sets[n_a + j])
    lo, hi = np.quantile(raw, 0.01), np.quantile(raw, 0.90)
    if hi <= lo:
        return np.ones_like(raw)
    return np.clip((raw - lo) / (hi - lo), 0.0, 1.0)


def _filter_anchors(
    pairs: list[tuple[int, int]],
    scores: np.ndarray,
    A_feat: np.ndarray,
    B_feat: np.ndarray,
    k_filter: int,
) -> tuple[list[tuple[int, int]], np.ndarray]:
    """Drop anchors whose cells are not neighbours in the original feature space.

    ``A_feat`` / ``B_feat`` are the L2-normalized (features × cells) anchor
    matrices. An anchor ``(i, j)`` survives only if query cell ``j`` is among
    reference cell ``i``'s ``k_filter`` nearest neighbours in expression space.
    """
    if not pairs:
        return pairs, scores
    keep_neigh = _nearest(A_feat.T, B_feat.T, min(k_filter, B_feat.shape[1]))
    allowed = [set(row) for row in keep_neigh]
    kept, kept_scores = [], []
    for m, (i, j) in enumerate(pairs):
        if j in allowed[i]:
            kept.append((i, j))
            kept_scores.append(scores[m])
    return kept, np.asarray(kept_scores)


# ----------------------------------------------------------------------
# Data extraction
# ----------------------------------------------------------------------


def _anchor_feature_matrix(obj, features: list[str], layer: str) -> np.ndarray:
    """Scaled (features × cells) matrix for the shared anchor features."""
    from .reduction import _get_scaled_data

    assay = obj.get_assay()
    return _get_scaled_data(assay, features, layer)


def _data_matrix(obj, features: list[str]) -> np.ndarray:
    """Log-normalized ``data`` (features × cells) for the given features.

    Uses ``layer_data`` so the rows come back in ``features`` order regardless
    of how the assay stores its layer (a v5 layer may hold its own subset).
    """
    import scipy.sparse as sp

    mat = obj.get_assay().layer_data("data", features=list(features))
    if sp.issparse(mat):
        return mat.toarray().astype(float)
    return np.asarray(mat).astype(float)


def _integration_features(objects, anchor_features: Optional[list[str]]) -> list[str]:
    """The features anchors run on: the caller's, else shared variable features."""
    from .reduction import _default_features

    if anchor_features is not None:
        common = anchor_features
    else:
        per_object = [set(_default_features(obj.get_assay(), None)) for obj in objects]
        shared = set.intersection(*per_object) if per_object else set()
        # Preserve the first object's ordering for determinism.
        first = _default_features(objects[0].get_assay(), None)
        common = [f for f in first if f in shared]
        if not common:  # fall back to the shared raw feature set
            feat_sets = [set(obj.get_assay().features()) for obj in objects]
            shared_all = set.intersection(*feat_sets)
            common = [f for f in objects[0].get_assay().features() if f in shared_all]
    # Keep only features present (with scale/data) in every object.
    for obj in objects:
        present = set(obj.get_assay().features())
        common = [f for f in common if f in present]
    if not common:
        raise ValueError("No shared anchor features across the objects.")
    return list(common)


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def find_integration_anchors(
    objects: list,
    anchor_features: Optional[list[str]] = None,
    reduction: str = "cca",
    dims: int = 30,
    k_anchor: int = 5,
    k_filter: int = 200,
    k_score: int = 30,
    reference: int = 0,
    layer: str = "scale.data",
    seed: int = 42,
) -> IntegrationAnchors:
    """Find anchors linking each dataset to the reference (Seurat's ``FindIntegrationAnchors``).

    Mirrors ``FindIntegrationAnchors(object.list, reduction = "cca")``. Anchors
    are mutual nearest neighbours in a shared CCA (or reciprocal-PCA) space,
    scored by neighbourhood consistency and filtered against the original
    expression space.

    Parameters
    ----------
    objects         : list of Shanuz objects (each normalized + with variable
                      features / scaled data). ``objects[reference]`` is treated
                      as the reference every other dataset is anchored to.
    anchor_features : features to integrate on (default: variable features
                      shared across all objects).
    reduction       : ``"cca"`` or ``"rpca"``.
    dims            : number of shared dimensions to use.
    k_anchor        : neighbours for the mutual-nearest-neighbour search.
    k_filter        : neighbourhood size for the feature-space anchor filter
                      (set to 0 or None to skip filtering).
    k_score         : neighbourhood size for anchor scoring.
    reference       : index of the reference dataset in ``objects``.
    layer           : layer to draw the shared-space expression from.
    seed            : random seed for the PCA/neighbour steps.

    Returns
    -------
    IntegrationAnchors
    """
    reduction = reduction.lower()
    if reduction not in REDUCTIONS:
        raise ValueError(
            f"Unknown reduction {reduction!r}. Supported: {REDUCTIONS}."
        )
    if len(objects) < 2:
        raise ValueError("find_integration_anchors needs at least two objects.")
    if not 0 <= reference < len(objects):
        raise IndexError(f"reference index {reference} out of range.")

    features = _integration_features(objects, anchor_features)

    ref_scaled = _l2_normalize_cols(
        _anchor_feature_matrix(objects[reference], features, layer)
    )

    rows = []
    weight_embeddings: dict[int, np.ndarray] = {}
    used_dims = min(dims, ref_scaled.shape[1] - 1, ref_scaled.shape[0])

    for d in range(len(objects)):
        if d == reference:
            continue
        query_scaled = _l2_normalize_cols(
            _anchor_feature_matrix(objects[d], features, layer)
        )
        used = min(used_dims, query_scaled.shape[1] - 1)

        if reduction == "cca":
            emb_ref, emb_query = _cca(ref_scaled, query_scaled, used)
            pairs = _mutual_nn(
                emb_ref, emb_query, emb_query, emb_ref, k_anchor
            )
            combined = np.vstack([emb_ref, emb_query])
            weight_emb = emb_query
        else:  # rpca — reciprocal PCA projections
            load_ref = _pca_loadings(ref_scaled, used, seed=seed)
            load_query = _pca_loadings(query_scaled, used, seed=seed)
            ref_in_ref = ref_scaled.T @ load_ref
            query_in_ref = query_scaled.T @ load_ref
            ref_in_query = ref_scaled.T @ load_query
            query_in_query = query_scaled.T @ load_query
            # Reciprocal search: find B(query)-neighbours of A(ref) in the query's
            # PCA space, and A-neighbours of B in the ref's space. _mutual_nn wants
            # (A-in-B, B-in-B, B-in-A, A-in-A) — args 1 & 4 are the reference
            # (n_ref rows), args 2 & 3 the query (n_query rows). Passing them in a
            # different order silently mismatches the index spaces and, when
            # n_query > n_ref, indexes past the b_sets list (IndexError).
            pairs = _mutual_nn(
                ref_in_query, query_in_query, query_in_ref, ref_in_ref, k_anchor
            )
            combined = np.vstack([ref_in_ref, query_in_ref])
            weight_emb = query_in_ref

        n_ref = ref_scaled.shape[1]
        scores = _score_anchors(pairs, combined, n_ref, k_score)
        if k_filter:
            pairs, scores = _filter_anchors(
                pairs, scores, ref_scaled, query_scaled, k_filter
            )

        weight_embeddings[d] = weight_emb
        for (i, j), s in zip(pairs, scores):
            rows.append((reference, int(i), d, int(j), float(s)))

    anchors = pd.DataFrame(
        rows, columns=["dataset1", "cell1", "dataset2", "cell2", "score"]
    )
    return IntegrationAnchors(
        anchors=anchors,
        objects=objects,
        reference=reference,
        reduction=reduction,
        anchor_features=features,
        dims=used_dims,
        weight_embeddings=weight_embeddings,
    )


def integrate_data(
    anchors: IntegrationAnchors,
    new_assay: str = "integrated",
    k_weight: int = 100,
    sd_weight: float = 1.0,
    add_cell_ids: Optional[list[str]] = None,
) -> "object":
    """Batch-correct query datasets onto the reference (Seurat's ``IntegrateData``).

    Mirrors ``IntegrateData(anchors)``. For every query dataset, each cell is
    corrected by a distance-weighted sum of anchor correction vectors
    (``expr_ref − expr_query``); the reference is left unchanged. The corrected
    expression of the anchor features is stored as the ``data`` layer of a new
    ``"integrated"`` assay on a merged object, which becomes the active assay.

    Downstream: ``scale_data`` + ``run_pca`` on the integrated assay yields an
    embedding that clusters by cell type rather than by batch.

    Parameters
    ----------
    anchors      : an :class:`IntegrationAnchors` from
                   :func:`find_integration_anchors`.
    new_assay    : name for the corrected assay (default ``"integrated"``).
    k_weight     : anchors used to weight each query cell's correction.
    sd_weight    : bandwidth multiplier for the Gaussian anchor kernel.
    add_cell_ids : optional per-object prefixes for the merged cell names.

    Returns
    -------
    Shanuz
        A merged object carrying the ``new_assay`` assay (active) alongside the
        original assay.
    """
    from .assay import Assay

    objects = anchors.objects
    ref = anchors.reference
    features = anchors.anchor_features

    # Merge order: reference first, then the remaining datasets in list order.
    order = [ref] + [d for d in range(len(objects)) if d != ref]
    ref_obj = objects[ref]
    others = [objects[d] for d in order[1:]]

    if add_cell_ids is not None:
        ordered_ids = [add_cell_ids[d] for d in order]
    else:
        ordered_ids = None

    merged = ref_obj.merge(others, add_cell_ids=ordered_ids)

    ref_data = _data_matrix(ref_obj, features)  # features × cells_ref (unchanged)
    corrected_blocks = [ref_data]

    for d in order[1:]:
        query_data = _data_matrix(objects[d], features)  # features × cells_q
        pair = anchors.anchors[anchors.anchors["dataset2"] == d]
        weight_emb = anchors.weight_embeddings[d]

        if len(pair) == 0:
            # No anchors to this dataset — leave it uncorrected.
            corrected_blocks.append(query_data)
            continue

        i_idx = pair["cell1"].to_numpy()
        j_idx = pair["cell2"].to_numpy()
        anchor_scores = pair["score"].to_numpy()

        # Correction vectors in feature space: reference minus query at anchors.
        bv = ref_data[:, i_idx] - query_data[:, j_idx]  # features × n_anchor

        anchor_pos = weight_emb[j_idx]  # n_anchor × dims (query anchor cells)
        k = min(k_weight, anchor_pos.shape[0])
        dist = _pairwise_to_anchors(weight_emb, anchor_pos, k)
        nn_idx = dist[1]        # cells_q × k  → indices into the anchor list
        nn_dist = dist[0]       # cells_q × k

        # Gaussian kernel over the k nearest anchors, scaled by anchor score.
        bandwidth = nn_dist[:, -1:].copy()
        bandwidth[bandwidth == 0] = 1.0
        weights = np.exp(-((nn_dist / (bandwidth / sd_weight)) ** 2))
        weights = weights * anchor_scores[nn_idx]
        wsum = weights.sum(axis=1, keepdims=True)
        wsum[wsum == 0] = 1.0
        weights = weights / wsum

        # correction[:, c] = Σ_k bv[:, nn_idx[c, k]] * weights[c, k]
        contrib = bv[:, nn_idx]  # features × cells_q × k
        correction = np.einsum("fck,ck->fc", contrib, weights)
        corrected_blocks.append(query_data + correction)

    integrated = np.hstack(corrected_blocks)  # features × total_cells
    cell_names = merged.cell_names()

    integrated_assay = Assay(
        data=integrated,
        feature_names=list(features),
        cell_names=list(cell_names),
        var_features=list(features),
        key=f"{new_assay.lower()}_",
    )
    merged.assays[new_assay] = integrated_assay
    merged.active_assay = new_assay
    return merged


def _pairwise_to_anchors(
    query_emb: np.ndarray, anchor_pos: np.ndarray, k: int
) -> tuple[np.ndarray, np.ndarray]:
    """Distances and indices from each query cell to its ``k`` nearest anchors."""
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k).fit(anchor_pos)
    dist, idx = nn.kneighbors(query_emb)
    return dist, idx
