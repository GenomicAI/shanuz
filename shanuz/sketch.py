"""Leverage-score sketching — Seurat's ``SketchData`` / ``ProjectData``.

A million-cell atlas is mostly redundant: thousands of near-identical cells sit
on top of each other in the common states, while the rare states — the ones you
usually care about — are a handful of points each. Clustering, UMAP and marker
tests on the whole thing are slow *and* dominated by the redundant majority.
Sketching fixes both: pick a small, information-dense subset ("the sketch"),
do the expensive analysis on that, then :func:`project_data` the answers back
onto every cell.

The subset is not a uniform random sample — that would just reproduce the
majority. Cells are drawn with probability proportional to their **statistical
leverage**. Leverage measures how much a cell influences the geometry of the
data: a cell in a dense, redundant cloud has low leverage (drop it and the
subspace is unchanged), while a cell in a sparse, distinctive corner has high
leverage (it is the only evidence that corner exists). Sampling by leverage
therefore keeps the rare states that uniform sampling throws away.

**Computing leverage without an SVD of the whole matrix.** The exact leverage of
row *i* of an *n × d* matrix ``A`` is ``ℓ_i = ‖U_i‖²``, where ``U`` are the left
singular vectors — an *O(n d²)* SVD, the very cost we are trying to avoid. The
classical fix (Drineas et al. 2012, and what Seurat's ``LeverageScore`` does) is
a **subspace embedding**: hit the *n* rows with a sparse random ``CountSketch``
``S`` (one signed entry per column, so ``SA`` costs a single pass over the
non-zeros) to get a small ``nsketch × d`` matrix ``B`` with ``BᵀB ≈ AᵀA``.
Whitening ``A`` with ``B``'s right singular vectors, ``Z = A · V Σ⁻¹``, makes
``Z`` approximately orthonormal, and ``ℓ_i ≈ ‖Z_i‖²`` — leverage for every cell
from one small factorisation. When ``nsketch ≥ n`` (small data) no sketch is
taken and the scores are exact.

:func:`leverage_score` is that computation on its own; :func:`sketch_data` draws
the weighted subset and returns it as a standalone object; :func:`project_data`
is the inverse map — it pushes every full-dataset cell through the *sketch's*
PCA loadings (and, if present, its fitted UMAP model, reusing
:func:`shanuz.project_umap`), and optionally carries cluster labels from the
sketch to the full data via the :mod:`shanuz.transfer` anchors.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import scipy.sparse as sp

from .dimreduc import DimReduc
from .reduction import _default_features, _get_scaled_data


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def leverage_score(
    obj,
    nsketch: int = 5000,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    layer: str = "scale.data",
    var_name: Optional[str] = "leverage.score",
    eps: float = 1e-8,
    seed: int = 123,
) -> np.ndarray:
    """Per-cell statistical leverage (Seurat's ``LeverageScore``).

    A cell's leverage is how much it influences the column space of the data — low
    in a dense, redundant cloud, high in a sparse, distinctive corner. Leverage is
    approximated with a ``CountSketch`` subspace embedding so no full SVD of the
    data is needed; when ``nsketch`` is at least the number of cells the scores are
    computed exactly. The scores are non-negative, each is at most 1, and they sum
    to the dimension of the subspace (the rank).

    Parameters
    ----------
    obj      : a :class:`~shanuz.Shanuz` object (normalized + scaled).
    nsketch  : rows of the random sketch. Larger is more accurate and slower; the
               approximation kicks in only when ``nsketch`` is below the cell count.
    features : features to score on (default: the assay's variable features).
    assay    : assay to use (default: active assay).
    layer    : layer to draw the data from (default ``"scale.data"``).
    var_name : if given, the scores are also written to ``obj.meta_data[var_name]``.
    eps      : relative tolerance for dropping near-zero singular directions.
    seed     : random seed for the sketch.

    Returns
    -------
    numpy.ndarray
        One leverage score per cell, in ``obj.cell_names()`` order.
    """
    assay_name = assay or obj.active_assay
    assay_obj = obj.assays[assay_name]
    feats = _default_features(assay_obj, features)

    mat = _get_scaled_data(assay_obj, feats, layer)   # (n_features × n_cells)
    A = np.asarray(mat, dtype=float).T                # (n_cells × n_features)

    rng = np.random.default_rng(seed)
    scores = _leverage_from_matrix(A, nsketch, rng, eps)

    if var_name is not None:
        obj.meta_data[var_name] = scores
    return scores


def sketch_data(
    obj,
    ncells: int = 5000,
    method: str = "LeverageScore",
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    layer: str = "scale.data",
    nsketch: int = 5000,
    sketched_assay: str = "sketch",
    var_name: Optional[str] = "leverage.score",
    seed: int = 123,
):
    """Draw a leverage-weighted subset of cells (Seurat's ``SketchData``).

    Mirrors ``SketchData(object, ncells = 5000, method = "LeverageScore")``. Each
    cell is sampled without replacement with probability proportional to its
    :func:`leverage_score`, so the rare states a uniform sample would drop are kept
    (indeed over-represented). The leverage scores are written back onto ``obj``'s
    metadata, and the drawn subset is returned as a **standalone**
    :class:`~shanuz.Shanuz` object — run the expensive analysis (PCA, clustering,
    UMAP) on it and use :func:`project_data` to extend the results to every cell.

    This departs from Seurat, which stores the sketch as an extra assay on the same
    object; here it is a separate object, matching the roadmap and shanuz's
    ``subset`` model. Its active assay is renamed to ``sketched_assay`` so the
    provenance is visible, and ``obj.misc["sketch"]`` records how it was drawn.

    Parameters
    ----------
    obj            : a :class:`~shanuz.Shanuz` object (normalized + scaled).
    ncells         : cells to keep (capped at the number available).
    method         : only ``"LeverageScore"`` is implemented.
    features       : features to score on (default: variable features).
    assay          : assay to use (default: active assay).
    layer          : layer to draw the data from.
    nsketch        : sketch size passed to :func:`leverage_score`.
    sketched_assay : name the returned object's active assay is renamed to.
    var_name       : metadata column the scores are written to on ``obj``.
    seed           : random seed for scoring and sampling.

    Returns
    -------
    Shanuz
        The sketched subset (a new object).
    """
    if method != "LeverageScore":
        raise ValueError(
            f"Unknown sketch method {method!r}; only 'LeverageScore' is implemented."
        )

    scores = leverage_score(
        obj, nsketch=nsketch, features=features, assay=assay, layer=layer,
        var_name=var_name, seed=seed,
    )

    n = len(obj)
    k = min(ncells, n)
    total = float(scores.sum())
    probs = scores / total if total > 0 else np.full(n, 1.0 / n)

    rng = np.random.default_rng(seed + 1)
    idx = np.sort(rng.choice(n, size=k, replace=False, p=probs))
    cell_names = obj.cell_names()
    cells = [cell_names[i] for i in idx]

    sketched = obj.subset(cells=cells)
    if sketched_assay and sketched_assay != sketched.active_assay:
        src = sketched.active_assay
        sketched.assays[sketched_assay] = sketched.assays.pop(src)
        sketched.active_assay = sketched_assay
    sketched.misc["sketch"] = {
        "method": method,
        "ncells": k,
        "from_cells": n,
        "source_assay": obj.active_assay,
        "leverage_var": var_name,
    }
    return sketched


def project_data(
    full,
    sketch,
    reduction: str = "pca",
    full_reduction: str = "pca.full",
    umap_reduction: str = "umap",
    full_umap_reduction: str = "ref.umap",
    refdata: Optional[Union[str, dict, np.ndarray, list]] = None,
    project_umap: bool = True,
    dims: Optional[Union[list[int], range]] = None,
    k_weight: int = 50,
    sd_weight: float = 1.0,
    layer: str = "scale.data",
    seed: int = 42,
):
    """Extend a sketch's analysis to the full dataset (Seurat's ``ProjectData``).

    The inverse of :func:`sketch_data`: once the sketch has been reduced and
    (optionally) clustered, every full-dataset cell is placed into the sketch's
    coordinate system and, if asked, given the sketch's labels.

    1. **PCA.** Each full cell is pushed through the sketch's PCA loadings — the
       same "project into a space this cell never helped define" linear map that
       :func:`shanuz.project_umap` uses — and stored as
       ``full.reductions[full_reduction]``.
    2. **UMAP** (when ``project_umap`` and the sketch carries a fitted UMAP model):
       the projected cells are run through the sketch's UMAP via
       :func:`shanuz.project_umap`, stored as ``full.reductions[full_umap_reduction]``.
    3. **Labels** (when ``refdata`` is given): :func:`shanuz.find_transfer_anchors`
       links the sketch (reference) to the full data (query) and
       :func:`shanuz.transfer_data` carries the labels across, written onto
       ``full.meta_data``.

    ``full`` is mutated in place and returned.

    Parameters
    ----------
    full                : the full :class:`~shanuz.Shanuz` object (normalized +
                          scaled on the sketch's PCA features).
    sketch              : the sketched object from :func:`sketch_data`, already
                          carrying a PCA (and optionally a fitted UMAP).
    reduction           : sketch reduction whose loadings project the full data.
    full_reduction      : storage key for the projected PCA on ``full``.
    umap_reduction      : sketch reduction holding the fitted UMAP model.
    full_umap_reduction : storage key for the projected UMAP on ``full``.
    refdata             : labels to transfer. A ``dict`` ``{new_col: sketch_col}``
                          writes each transferred label under ``new_col`` (plus
                          ``new_col.score``); a ``str`` column / 1-D array writes
                          the full ``predicted.id`` / ``prediction.score.*`` frame.
                          ``None`` skips transfer.
    project_umap        : also project the sketch's UMAP when a fitted model exists.
    dims                : PCA dimensions to feed the UMAP model (default: all).
    k_weight, sd_weight : anchor-weighting knobs for :func:`shanuz.transfer_data`.
    layer               : layer to draw the full data's expression from.
    seed                : random seed for the anchor step.

    Returns
    -------
    Shanuz
        ``full``, now carrying the projected reduction(s) and any transferred labels.
    """
    from .mapping import _project_into_reference_pca, project_umap as _project_umap

    if reduction not in sketch.reductions:
        raise KeyError(
            f"Sketch reduction {reduction!r} not found; run run_pca(sketch) first."
        )

    sketch_pca = sketch.reductions[reduction]
    full_pca = _project_into_reference_pca(full, sketch, reduction, layer)

    key = getattr(sketch_pca, "_key", "PC_") or "PC_"
    dim_names = [f"{key}{i + 1}" for i in range(full_pca.shape[1])]
    full.reductions[full_reduction] = DimReduc(
        cell_embeddings=full_pca,
        cell_names=full.cell_names(),
        feature_loadings=np.asarray(sketch_pca.feature_loadings),
        feature_names=dim_names if not sketch_pca.features() else list(sketch_pca.features()),
        assay_used=full.active_assay,
        key=key,
        misc={"projected_from": reduction, "sketch_reduction": reduction},
    )

    if (
        project_umap
        and umap_reduction in sketch.reductions
        and sketch.reductions[umap_reduction].misc.get("umap_model") is not None
    ):
        _project_umap(
            full,
            sketch,
            reduction=reduction,
            umap_reduction=umap_reduction,
            dims=dims,
            reduction_name=full_umap_reduction,
            layer=layer,
        )

    if refdata is not None:
        from .transfer import find_transfer_anchors

        anchors = find_transfer_anchors(
            sketch, full, reduction="pcaproject", layer=layer, seed=seed
        )
        _write_refdata(anchors, full, refdata, k_weight, sd_weight)

    return full


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _count_sketch(n: int, nsketch: int, rng: np.random.Generator) -> sp.csr_matrix:
    """A CountSketch subspace embedding — an ``nsketch × n`` sparse sign matrix.

    Each of the ``n`` columns has exactly one non-zero, a random ``±1`` in a random
    row, so ``S @ A`` hashes the cells into ``nsketch`` signed buckets in a single
    pass over the data's non-zeros.
    """
    rows = rng.integers(0, nsketch, size=n)
    signs = rng.choice(np.array([-1.0, 1.0]), size=n)
    cols = np.arange(n)
    return sp.csr_matrix((signs, (rows, cols)), shape=(nsketch, n))


def _leverage_from_matrix(
    A: np.ndarray, nsketch: int, rng: np.random.Generator, eps: float
) -> np.ndarray:
    """Leverage scores of the rows of ``A`` (``n × d``), via a sketched whitening.

    ``B = S A`` is a small subspace embedding (or ``A`` itself when no sketch is
    taken). Whitening ``A`` with ``B``'s right singular vectors, ``Z = A V Σ⁻¹``,
    makes ``Z`` approximately orthonormal, so ``ℓ_i = ‖Z_i‖²`` approximates the
    exact leverage ``a_iᵀ (AᵀA)⁻¹ a_i`` — and equals it when ``B = A``.
    """
    n, d = A.shape
    if 0 < nsketch < n and nsketch >= d:
        B = np.asarray(_count_sketch(n, nsketch, rng) @ A)
    else:
        B = A

    _, s, vt = np.linalg.svd(B, full_matrices=False)
    if s.size == 0 or s[0] == 0:
        return np.zeros(n)

    tol = s[0] * max(B.shape) * eps
    keep = s > tol
    whiten = vt[keep].T / s[keep]          # (d × r)
    z = A @ whiten                          # (n × r)
    return np.einsum("ij,ij->i", z, z)


def _write_refdata(anchors, full, refdata, k_weight: int, sd_weight: float) -> None:
    """Carry ``refdata`` from the sketch onto the full object's metadata."""
    from .transfer import transfer_data

    meta = full.meta_data
    if isinstance(refdata, dict):
        for out_col, src_col in refdata.items():
            preds = transfer_data(
                anchors, src_col, k_weight=k_weight, sd_weight=sd_weight
            )
            meta[out_col] = preds["predicted.id"].reindex(meta.index).to_numpy()
            meta[f"{out_col}.score"] = (
                preds["prediction.score.max"].reindex(meta.index).to_numpy()
            )
    else:
        preds = transfer_data(
            anchors, refdata, k_weight=k_weight, sd_weight=sd_weight
        )
        for col in preds.columns:
            meta[col] = preds[col].reindex(meta.index).to_numpy()
