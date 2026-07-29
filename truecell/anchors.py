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
  3. **filters** anchors whose members are not even near each other in the
     original log-normalized expression space, then **scores** the survivors by
     how much of their neighbourhood the two members share (a consistent anchor
     sits in coherent local structure). The order matters: the score is rescaled
     against the percentiles of whatever set it is given;
  4. **corrects** each query dataset onto the reference by adding, to every
     query cell, a distance-weighted average of the anchor *correction vectors*
     ``expr_ref − expr_query`` — pulling matched populations on top of each
     other while leaving genuinely reference-only structure alone. The weights
     are found in a fresh PCA of the merged pair, not in the anchor space.

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

import warnings
from typing import Optional

import numpy as np
import pandas as pd

from .dimreduc import DimReduc

REDUCTIONS = ("cca", "rpca")


# ----------------------------------------------------------------------
# Result container
# ----------------------------------------------------------------------


class IntegrationAnchors:
    """Anchors linking a reference dataset to one or more query datasets.

    Slots
    -----
    - ``anchors`` — DataFrame with columns ``dataset1, cell1, dataset2,
      cell2, score``. ``dataset1`` is always the reference; the
      cell columns hold *within-dataset* 0-based row indices.
    - ``objects`` — the list of Truecell objects passed to
      :func:`find_integration_anchors` (order preserved).
    - ``reference`` — index into ``objects`` of the reference dataset.
    - ``reduction`` — ``"cca"`` or ``"rpca"`` — how the shared space was built.
    - ``anchor_features`` — the features the anchors (and correction) run on.
    - ``dims`` — number of shared dimensions used.
    - ``weight_embeddings`` — ``{query_index: (n_query_cells × dims) array}`` — each
      query dataset's cells in the shared anchor space. Kept for
      inspection; :func:`integrate_data` does *not* weight with
      it, because Seurat weights in a fresh PCA of the merged
      pair instead.
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


def _standardize_and_l2(emb: np.ndarray) -> np.ndarray:
    """Normalize a reciprocal-PCA embedding the way Seurat's ``ReciprocalProject`` does.

    With ``l2.norm = TRUE`` (Seurat's default) each reciprocal space — the
    stacked ``[reference; query]`` cells in *one* dataset's PCA — is normalized
    in two steps before the neighbour search: every dimension (column) is divided
    by its standard deviation across the whole stack, then every cell (row) is
    L2-normalized. The per-dimension scaling is the load-bearing step: a
    projection's leading PC carries orders of magnitude more variance than its
    trailing ones, so without it the nearest-neighbour search collapses onto PC1
    and mates cells by library depth instead of type. (CCA needs no such step —
    its singular vectors already share a scale — which is why only RPCA was hit.)
    """
    sd = emb.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    return _l2_normalize_rows(emb / sd)


def _standardize_cols(mat: np.ndarray) -> np.ndarray:
    """Seurat's ``Standardize``: z-score every column (cell) of features × cells.

    This is *not* the per-cell L2 normalization it is easily mistaken for. The
    cross-covariance of two standardized matrices is a correlation matrix
    between cells; of two L2-normalized ones, a cosine-similarity matrix. They
    have different singular vectors, so using the wrong one moves every anchor:
    against Seurat 5.5.1 on ifnb, L2 recovers 70% of the anchors where
    ``Standardize`` recovers 100%.
    """
    mu = mat.mean(axis=0, keepdims=True)
    sd = mat.std(axis=0, ddof=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (mat - mu) / sd


def _check_features(mats: list[np.ndarray]) -> np.ndarray:
    """Indices of features with non-zero variance in *every* matrix.

    ``RunCCA`` runs ``CheckFeatures`` on both objects' scale.data and silently
    drops anything constant in either — 83 of the 2,000 anchor features on the
    ifnb pair. Left in, those columns standardize to NaN and contaminate the
    whole cross-covariance.
    """
    ok = np.ones(mats[0].shape[0], dtype=bool)
    for mat in mats:
        ok &= mat.std(axis=1, ddof=1) > 0
    return np.where(ok)[0]


def _cca(
    A: np.ndarray, B: np.ndarray, dims: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Shared CCA embedding of two datasets (Seurat's ``RunCCA``).

    ``A`` (features × cells_a) and ``B`` (features × cells_b) are standardized
    per cell. The SVD of the cross-covariance ``AᵀB`` (cells_a × cells_b) yields
    left/right singular vectors that place the two datasets' cells in one space;
    the stacked embedding is sign-fixed the way ``RunCCA.default`` does it, then
    each cell's coordinates are L2-normalized (Seurat's ``L2Dim``).

    Returns the two L2-normalized halves plus the feature loadings
    (``ProjectDim(do.center = FALSE)``: scale.data ᐧ the *pre*-L2 embedding),
    which is what picks the genes the anchor filter runs on.
    """
    d1, d2 = _standardize_cols(A), _standardize_cols(B)
    U, _s, Vt = np.linalg.svd(d1.T @ d2, full_matrices=False)
    d = min(dims, U.shape[1])
    ccv = np.vstack([U[:, :d], Vt[:d, :].T])
    # RunCCA forces each canonical vector's first entry positive. A whole-column
    # flip is distance-preserving, so this only matters for reproducing R's
    # signs, not for the anchors.
    flip = np.sign(ccv[0, :])
    flip[flip == 0] = 1.0
    ccv = ccv * flip
    loadings = np.hstack([A, B]) @ ccv
    n_a = A.shape[1]
    emb = _l2_normalize_rows(ccv)
    return emb[:n_a], emb[n_a:], loadings


def _top_dim_features(
    loadings: np.ndarray, features_per_dim: int = 100, max_features: int = 200
) -> list[int]:
    """Seurat's ``TopDimFeatures``: the genes ``FilterAnchors`` actually runs on.

    For each candidate count *y*, take the ``round(y/2)`` highest- and
    lowest-loading features of every dimension; keep the largest *y* whose union
    still fits under ``max_features``. The filter therefore runs in a ≤200-gene
    subspace picked from the CCA loadings, not in the full anchor-feature space.
    """
    n_dims = loadings.shape[1]
    max_features = max(n_dims * 2, max_features)
    order = np.argsort(-loadings, axis=0, kind="stable")

    def take(y: int) -> set[int]:
        half = 0 if y == 1 else int(np.round(y / 2.0))
        if half == 0:
            return set()
        got: set[int] = set()
        for d in range(n_dims):
            got.update(order[:half, d].tolist())
            got.update(order[-half:, d].tolist())
        return got

    counts = np.array([len(take(y)) for y in range(1, features_per_dim + 1)])
    ok = np.where(counts < max_features)[0]
    if len(ok) == 0:
        return sorted(take(1))
    return sorted(take(int(ok[np.argmax(counts[ok])]) + 1))


def _pca_loadings(mat: np.ndarray, dims: int, seed: int = 42) -> np.ndarray:
    """Top-``dims`` PCA loadings (features × dims) of a features × cells matrix.

    Uses an **exact** economy SVD, not sklearn's ``PCA``. For a matrix this
    shape sklearn defaults to *randomized* SVD, which is accurate for the
    leading components but visibly wrong in the trailing ones — on the ifnb
    batches only 12–14 of 30 PCs matched Seurat's irlba above 0.99. Ordinarily
    that trailing noise is harmless, but reciprocal PCA standardizes each
    projected dimension by its own SD before the neighbour search, which is not
    rotation-invariant, so a rotated trailing axis becomes a *different*
    reciprocal space and a different anchor. Exact SVD matches irlba to 1.0000
    on every PC, which takes RPCA anchor recall from 45% to full agreement.
    ``seed`` is accepted for signature compatibility but unused: the SVD is
    deterministic.

    Seurat's ``RunPCA`` runs irlba on ``t(scale.data)`` without re-centring
    (scale.data is already gene-centred), so this does the same.
    """
    d = max(min(dims, min(mat.shape) - 1), 1)
    # SVD of cells × features; the right singular vectors are the loadings.
    _u, _s, vt = np.linalg.svd(mat.T, full_matrices=False)
    return vt[:d, :].T  # features × dims


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
    return _mutual_nn_from({"ab": a_to_b, "ba": b_to_a}, k)


def _mutual_nn_from(nbrs: dict[str, np.ndarray], k_anchor: int) -> list[tuple[int, int]]:
    """``FindAnchorPairs``: mutual pairs within the first ``k_anchor`` columns.

    ``FindNN`` searches at ``max(k_anchor, k_score)`` so the same neighbour
    tables serve both the pairing and the scoring; only the leading
    ``k_anchor`` of them decide an anchor.
    """
    ab = nbrs["ab"][:, :k_anchor]
    ba_sets = [set(row) for row in nbrs["ba"][:, :k_anchor]]
    return [
        (int(i), int(j))
        for i, neigh in enumerate(ab)
        for j in neigh
        if i in ba_sets[j]
    ]


# ----------------------------------------------------------------------
# Scoring & filtering
# ----------------------------------------------------------------------


def _neighbor_sets(
    emb_a: np.ndarray, emb_b: np.ndarray, k: int
) -> dict[str, np.ndarray]:
    """Seurat's ``FindNN``: four searches, not one over the pooled stack.

    ``aa``/``bb`` are within-dataset (computed at ``k + 1`` so the self hit at
    position 0 is not one of the ``k`` kept), ``ab``/``ba`` across.
    """
    return {
        "aa": _nearest(emb_a, emb_a, k + 1),
        "bb": _nearest(emb_b, emb_b, k + 1),
        "ab": _nearest(emb_a, emb_b, k),
        "ba": _nearest(emb_b, emb_a, k),
    }


def _score_anchors(
    pairs: list[tuple[int, int]],
    nbrs: dict[str, np.ndarray],
    n_a: int,
    k_score: int,
) -> np.ndarray:
    """Shared-neighbourhood score in [0, 1] for each anchor (``ScoreAnchors``).

    Each anchor member's neighbourhood is ``k_score`` neighbours *within its own
    dataset* plus ``k_score`` *in the other* — 2·k_score cells drawn from four
    separate searches. A single kNN over the pooled stack is not the same thing
    and biases the scores down: with a batch effect present a cell's pooled
    neighbours are nearly all same-batch, so the two members of an anchor share
    almost nothing. Scores are rescaled with Seurat's 1st/90th-percentile clamp.
    """
    if not pairs:
        return np.array([])
    set_a = [
        set(nbrs["aa"][i, :k_score]) | {int(x) + n_a for x in nbrs["ab"][i, :k_score]}
        for i in range(nbrs["aa"].shape[0])
    ]
    set_b = [
        set(nbrs["ba"][j, :k_score]) | {int(x) + n_a for x in nbrs["bb"][j, :k_score]}
        for j in range(nbrs["bb"].shape[0])
    ]
    raw = np.array([len(set_a[i] & set_b[j]) for i, j in pairs], dtype=float)
    lo, hi = np.quantile(raw, 0.01), np.quantile(raw, 0.90)
    if hi <= lo:
        return np.ones_like(raw)
    return np.clip((raw - lo) / (hi - lo), 0.0, 1.0)


def _filter_anchors(
    pairs: list[tuple[int, int]],
    A_feat: np.ndarray,
    B_feat: np.ndarray,
    k_filter: int,
) -> list[tuple[int, int]]:
    """Drop anchors whose cells are not neighbours in the original feature space.

    ``A_feat`` / ``B_feat`` are (features × cells) **log-normalized** matrices
    already restricted to the top CCA-loading genes — ``FilterAnchors`` runs on
    the ``data`` layer, not on scale.data, and on ``TopDimFeatures`` rather than
    the whole anchor-feature set. An anchor ``(i, j)`` survives only if query
    cell ``j`` is among reference cell ``i``'s ``k_filter`` nearest neighbours
    there.
    """
    if not pairs:
        return pairs
    # Seurat keeps every anchor rather than shrinking k when a dataset is
    # smaller than k.filter: with too few cells the filter carries no
    # information, and quietly narrowing it would drop good anchors instead.
    if min(A_feat.shape[1], B_feat.shape[1]) < k_filter:
        return pairs
    cn1 = _l2_normalize_cols(A_feat).T
    cn2 = _l2_normalize_cols(B_feat).T
    allowed = [set(row) for row in _nearest(cn1, cn2, k_filter)]
    return [(i, j) for i, j in pairs if j in allowed[i]]


# ----------------------------------------------------------------------
# Data extraction
# ----------------------------------------------------------------------


def _anchor_feature_matrix(obj, features: list[str], layer: str) -> np.ndarray:
    """Scaled (features × cells) matrix for the shared anchor features.

    Every object must return the *same* rows in the same order — the reference
    and query matrices are compared row-for-row — so a per-object drop here is
    an error rather than a warning. `_integration_features` intersects against
    the layer precisely so this cannot fire.
    """
    from .reduction import _prep_dr

    assay = obj.get_assay()
    mat, used = _prep_dr(assay, features, layer)
    if used != list(features):
        missing = [f for f in features if f not in set(used)]
        raise ValueError(
            f"Object is missing {len(missing)} of the {len(features)} anchor "
            f"features in layer {layer!r} ({', '.join(missing[:10])}"
            f"{', ...' if len(missing) > 10 else ''}); the per-object matrices "
            f"would no longer describe the same genes row-for-row."
        )
    return mat


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


def _integration_features(
    objects, anchor_features: Optional[list[str]], layer: str = "scale.data"
) -> list[str]:
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
    # Keep only features every object carries *in the layer the anchors are
    # built from*. Checking `features()` instead — the assay's full list, which
    # is what stood here despite the comment — lets a feature through that one
    # object never scaled. `_anchor_feature_matrix` then drops it from that
    # object alone, and the reference and query matrices, which are multiplied
    # together row-for-row, come back describing different genes.
    from .reduction import _layer_feature_names

    for obj in objects:
        present = set(_layer_feature_names(obj.get_assay(), layer))
        common = [f for f in common if f in present]
    if not common:
        raise ValueError(
            f"No shared anchor features across the objects in layer {layer!r}. "
            f"Run scale_data() over a common feature set first."
        )
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
    objects         : list of Truecell objects (each normalized + with variable
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

    features = _integration_features(objects, anchor_features, layer)

    # Both reductions start from the per-object scale.data. CCA standardizes it
    # per cell inside _cca; reciprocal PCA runs straight on it (Seurat's
    # per-object RunPCA never column-normalizes for pca-based reductions).
    ref_raw = _anchor_feature_matrix(objects[reference], features, layer)

    rows = []
    weight_embeddings: dict[int, np.ndarray] = {}
    used_dims = min(dims, ref_raw.shape[1] - 1, ref_raw.shape[0])

    for d in range(len(objects)):
        if d == reference:
            continue
        query_raw = _anchor_feature_matrix(objects[d], features, layer)
        used = min(used_dims, query_raw.shape[1] - 1)

        if reduction == "cca":
            # CheckFeatures: RunCCA drops anchor features that are constant in
            # either object before standardizing, so the CCA and its loadings
            # live on this subset.
            keep = _check_features([ref_raw, query_raw])
            emb_ref, emb_query, loadings = _cca(ref_raw[keep], query_raw[keep], used)
            nbrs = _neighbor_sets(emb_ref, emb_query, max(k_anchor, k_score))
            pairs = _mutual_nn_from(nbrs, k_anchor)
            filter_feats = [features[keep[t]] for t in _top_dim_features(loadings)]
            weight_emb = emb_query
        else:  # rpca — reciprocal PCA projections (Seurat's ReciprocalProject)
            # The reciprocal spaces below are only as good as these loadings,
            # and reciprocal PCA standardizes each dimension before the
            # neighbour search (see _pca_loadings), so the trailing PCs have to
            # be right, not just the leading ones. Feeding this the exact SVD
            # loadings takes RPCA anchor agreement with Seurat from 45% to full.
            load_ref = _pca_loadings(ref_raw, used, seed=seed)
            load_query = _pca_loadings(query_raw, used, seed=seed)
            ref_in_ref = ref_raw.T @ load_ref        # ref cells, ref PCA
            query_in_ref = query_raw.T @ load_ref     # query projected into ref PCA
            ref_in_query = ref_raw.T @ load_query     # ref projected into query PCA
            query_in_query = query_raw.T @ load_query  # query cells, query PCA
            # Each object in its own unnormalized PCA — Seurat's plain "pca"
            # reduction, which is what the within-dataset neighbour tables use.
            own_ref, own_query = ref_in_ref, query_in_query

            # Seurat forms each reciprocal space as the stacked [ref; query] and,
            # with l2.norm=TRUE (its default), standardizes every dimension by its
            # SD over the whole stack before L2-normalizing each cell. Skipping
            # that let PC1's dominant variance swamp the neighbour search, so the
            # mutual pairs were wrong — RPCA under-integrated ifnb 4x (batch-mix
            # 0.22 vs Seurat's 0.91). Normalize each space, then split it back
            # into its ref/query halves for the reciprocal MNN.
            n_ref_cells = ref_in_ref.shape[0]
            ref_space = _standardize_and_l2(np.vstack([ref_in_ref, query_in_ref]))
            query_space = _standardize_and_l2(np.vstack([ref_in_query, query_in_query]))
            ref_in_ref, query_in_ref = ref_space[:n_ref_cells], ref_space[n_ref_cells:]
            ref_in_query, query_in_query = query_space[:n_ref_cells], query_space[n_ref_cells:]

            # Reciprocal search: find B(query)-neighbours of A(ref) in the query's
            # PCA space, and A-neighbours of B in the ref's space. The four
            # neighbour tables ScoreAnchors needs are asymmetric here — the
            # within-dataset ones come from each object's *own* PCA (Seurat's
            # nn.reduction stays "pca" on the rpca branch, so it never sees a
            # reciprocal space), the across-dataset ones from the reciprocal
            # projections. Getting the index spaces the wrong way round silently
            # mismatches them and, when n_query > n_ref, runs off the end.
            k_neighbor = max(k_anchor, k_score)
            nbrs = {
                "aa": _nearest(own_ref, own_ref, min(k_neighbor + 1, own_ref.shape[0])),
                "bb": _nearest(own_query, own_query,
                               min(k_neighbor + 1, own_query.shape[0])),
                "ab": _nearest(ref_in_query, query_in_query, k_neighbor),
                "ba": _nearest(query_in_ref, ref_in_ref, k_neighbor),
            }
            pairs = _mutual_nn_from(nbrs, k_anchor)
            filter_feats = None
            weight_emb = query_in_ref

        # FindAnchors filters BEFORE it scores, and the score is rescaled
        # against the 1st/90th percentiles of whatever set it is handed. Scoring
        # first sets those percentiles from anchors that are about to be thrown
        # away, which shifts every surviving score — same ranking, wrong values.
        #
        # Seurat forces k.filter <- NA for reciprocal-PCA (every pca-based
        # nn.reduction), skipping the expression-space filter entirely: the
        # reciprocal projection is itself an expression-space check, and
        # filtering on the shared anchor features here drops good anchors and
        # leaves RPCA under-integrating. Only CCA keeps the filter.
        if k_filter and filter_feats:
            pairs = _filter_anchors(
                pairs,
                _data_matrix(objects[reference], filter_feats),
                _data_matrix(objects[d], filter_feats),
                k_filter,
            )
        n_ref = ref_raw.shape[1]
        scores = _score_anchors(pairs, nbrs, n_ref, k_score)

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
    seed: int = 42,
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
    k_weight     : anchors used to weight each query cell's correction. Counts
                   anchors, not anchor cells — see :func:`_anchor_weights`.
    sd_weight    : bandwidth of the anchor kernel; enters as ``(2/sd_weight)²``.
    add_cell_ids : optional per-object prefixes for the merged cell names.
    seed         : random seed for the per-pair weight PCA.

    Returns
    -------
    Truecell
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

        if len(pair) == 0:
            # No anchors to this dataset — leave it uncorrected.
            corrected_blocks.append(query_data)
            continue

        i_idx = pair["cell1"].to_numpy()
        j_idx = pair["cell2"].to_numpy()
        anchor_scores = pair["score"].to_numpy()

        # Correction vectors in feature space: reference minus query at anchors.
        bv = ref_data[:, i_idx] - query_data[:, j_idx]  # features × n_anchor

        # The weights live in a PCA of *this pair*, not in the anchor space.
        # RunIntegration merges the reference and query, re-scales on the anchor
        # features and runs a fresh PCA, and FindWeights searches there. Reusing
        # the CCA embedding instead measures distances in a space built to make
        # the batches overlap, which is not the same neighbourhood.
        query_emb = _pair_weight_embedding(
            ref_obj, objects[d], features, anchors.dims, seed
        )
        weights = _anchor_weights(
            query_emb, j_idx, anchor_scores, k_weight, sd_weight
        )
        correction = bv @ weights           # features × cells_q
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


def integrate_embeddings(
    anchors: IntegrationAnchors,
    reduction,
    new_reduction: str = "integrated_dr",
    dims_to_integrate: Optional[list[int]] = None,
    k_weight: int = 100,
    sd_weight: float = 1.0,
) -> DimReduc:
    """Batch-correct an existing reduction (Seurat's ``IntegrateEmbeddings``).

    The v5 counterpart to :func:`integrate_data`, and a genuinely different
    algorithm rather than a wrapper over it. ``IntegrateData`` corrects
    *expression* and leaves you to re-scale and re-run PCA; ``IntegrateEmbeddings``
    corrects the **embedding itself**, so the output lives in the input
    reduction's basis and keeps its loadings.

    Seurat implements it by transposing the embedding into a fake assay whose
    "features" are the dimensions (``drtointegrate-1 …``) and pushing that
    through the very same anchor machinery, which is why this shares
    :func:`_anchor_weights` with :func:`integrate_data`. The one substantive
    difference is the weight space: ``RunIntegration``'s ``dims = NULL`` branch
    hands ``FindWeights`` the ``drtointegrate`` matrix itself, so neighbours are
    measured in the *uncorrected embedding* — not in the fresh per-pair PCA
    that the expression path builds.

    Parameters
    ----------
    anchors           : an :class:`IntegrationAnchors` from
                        :func:`find_integration_anchors`.
    reduction         : a :class:`DimReduc` covering every cell in
                        ``anchors.objects`` — the reduction to correct.
    new_reduction     : key for the returned reduction.
    dims_to_integrate : which dimensions to correct (0-indexed; default all).
    k_weight          : anchors used to weight each query cell's correction.
    sd_weight         : bandwidth of the anchor kernel; enters as ``(2/sd)²``.

    Returns
    -------
    DimReduc
        The corrected embedding, cells in merged order (reference first).
    """
    objects = anchors.objects
    ref = anchors.reference

    emb = np.asarray(reduction.cell_embeddings)
    pos = {c: i for i, c in enumerate(reduction.cells())}
    dims = list(range(emb.shape[1])) if dims_to_integrate is None \
        else list(dims_to_integrate)

    missing = [c for o in objects for c in o.cell_names() if c not in pos]
    if missing:
        raise ValueError(
            f"{len(missing)} cell(s) in the anchor objects are absent from "
            f"reduction {reduction.key!r} (first: {missing[0]!r}). "
            "IntegrateEmbeddings needs a reduction spanning every dataset."
        )

    def block(d):
        """The dataset's embedding as dims × cells — Seurat's drtointegrate."""
        cells = objects[d].cell_names()
        return emb[[pos[c] for c in cells]][:, dims].T

    order = [ref] + [d for d in range(len(objects)) if d != ref]
    ref_block = block(ref)
    blocks = {ref: ref_block}

    for d in order[1:]:
        query_block = block(d)
        pair = anchors.anchors[anchors.anchors["dataset2"] == d]
        if len(pair) == 0:
            blocks[d] = query_block
            continue
        i_idx = pair["cell1"].to_numpy()
        j_idx = pair["cell2"].to_numpy()
        bv = ref_block[:, i_idx] - query_block[:, j_idx]  # dims × n_anchor
        weights = _anchor_weights(
            query_block.T, j_idx, pair["score"].to_numpy(), k_weight, sd_weight
        )
        blocks[d] = query_block + bv @ weights

    corrected = np.hstack([blocks[d] for d in order]).T  # cells × dims
    cells = [c for d in order for c in objects[d].cell_names()]
    loadings = None
    if reduction.feature_loadings is not None and len(reduction.feature_loadings):
        loadings = np.asarray(reduction.feature_loadings)[:, dims]
    key = f"{new_reduction.replace('_', '').replace('.', '')}_"
    return DimReduc(
        cell_embeddings=corrected,
        feature_loadings=loadings,
        assay_used=reduction.assay_used,
        key=key,
        cell_names=cells,
        feature_names=(reduction.features() if loadings is not None else None),
    )


def _pair_weight_embedding(
    ref_obj, query_obj, features: list[str], dims: int, seed: int
) -> np.ndarray:
    """The query half of ``RunIntegration``'s weight PCA.

    Seurat merges the reference and query, re-runs ``ScaleData`` over the pair
    on the anchor features, runs a PCA of ``max(dims)`` components and hands
    *that* to ``FindWeights``. Returns the merged embedding restricted to the
    query cells, in the query object's own cell order.
    """
    from .preprocessing import scale_data
    from .reduction import run_pca

    merged = ref_obj.merge([query_obj])
    scale_data(merged, features=features)
    run_pca(merged, n_pcs=dims, features=features, reduction_name="_weights", seed=seed)
    emb = merged.reductions["_weights"].cell_embeddings
    order = {c: i for i, c in enumerate(merged.cell_names())}
    return emb[[order[c] for c in query_obj.cell_names()], :]


def _anchor_weights(
    query_emb: np.ndarray,
    anchor_cells: np.ndarray,
    scores: np.ndarray,
    k_weight: int,
    sd_weight: float,
) -> np.ndarray:
    """``FindWeights`` / ``FindWeightsC``: an (n_anchor × n_query) column-stochastic matrix.

    Three things here are easy to get wrong and all three move the answer:

    * the neighbour search runs over the **unique** query anchor cells, not over
      the anchor list — a cell that anchors five times is one candidate, not
      five;
    * ``k_weight`` then caps the number of **anchors** written, not the number
      of anchor cells. FindWeightsC walks those neighbours outwards, expands
      each into all of its anchor rows and stops at ``k_weight`` entries, so on
      a typical pair (~2.7 anchors per cell) only the nearest ~37 cells
      contribute;
    * the kernel is ``1 − exp(−d̃ · score / (2/sd)²)`` over ``d̃ = 1 − d/dₖ``,
      which *rises* with proximity and folds the anchor score into the exponent.
      A Gaussian in the raw distance multiplied by the score is a different
      curve with a different ranking.
    """
    from sklearn.neighbors import NearestNeighbors

    # RunIntegration lowers k.weight to the anchor count first, and only then
    # does FindWeights refuse when there are fewer distinct anchor cells left.
    if len(anchor_cells) < k_weight:
        warnings.warn(
            f"Number of anchors ({len(anchor_cells)}) is less than k_weight "
            f"({k_weight}); lowering k_weight to {len(anchor_cells)}.",
            stacklevel=2,
        )
        k_weight = len(anchor_cells)
    uniq, src = np.unique(anchor_cells, return_inverse=True)
    if len(uniq) < k_weight:
        raise ValueError(
            f"Number of anchor cells ({len(uniq)}) is less than k_weight "
            f"({k_weight}). Lower k_weight below {len(uniq)}, or raise k_anchor "
            "to find more anchors."
        )
    pos = query_emb[uniq, :]
    dist, idx = NearestNeighbors(n_neighbors=k_weight).fit(pos).kneighbors(query_emb)
    far = dist[:, -1:].copy()
    far[far == 0] = 1.0
    dscaled = 1.0 - dist / far

    rows_for: dict[int, list[int]] = {}
    for a, p in enumerate(src):
        rows_for.setdefault(int(p), []).append(a)

    weights = np.zeros((len(anchor_cells), query_emb.shape[0]))
    scale = (2.0 / sd_weight) ** 2
    for c in range(query_emb.shape[0]):
        filled = 0
        for slot in range(idx.shape[1]):
            if filled >= k_weight:
                break
            for a in rows_for.get(int(idx[c, slot]), ()):
                if filled >= k_weight:
                    break
                weights[a, c] = 1.0 - np.exp(-dscaled[c, slot] * scores[a] / scale)
                filled += 1
    total = weights.sum(axis=0, keepdims=True)
    total[total == 0] = 1.0
    return weights / total
