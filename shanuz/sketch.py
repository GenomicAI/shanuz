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
classical fix (Drineas et al. 2012) is a **subspace embedding**: hit the *n* rows
with a sparse random ``CountSketch`` ``S`` (one signed entry per column, so ``SA``
costs a single pass over the non-zeros) to get a small ``nsketch × d`` matrix with
the same Gram matrix up to distortion, whiten ``A`` against *its* triangular
factor, and read the row norms off a further random projection.

Seurat picks between two regimes, and :func:`leverage_score` mirrors both:

* **Few cells** (``ncells < nsketch * 1.5``) — no sketch is worth taking, so the
  scores come from a **rank-50 truncated SVD** of the *(features × cells)* matrix:
  ``ℓ = rowSums(V²)`` over the leading 50 right singular vectors. The scores then
  sum to 50, not to the rank.
* **Many cells** — ``CountSketch`` → ``QR`` → back-substitution → a
  Johnson–Lindenstrauss projection (``JLEmbed``), with ``ℓ_i = ‖Z_i‖²``.

**The truncation is the whole point, and getting it wrong is silent.** Leverage
scores sum to whatever rank you whiten against. Whiten against the *full* rank of
a typical single-cell matrix — 2000 variable features over a few thousand cells —
and every score is crushed towards ``d/n`` with a ceiling of 1: on PBMC 3k that is
a mean of 0.74 spanning only 0.54–0.99. Sampling weights that flat are
indistinguishable from uniform, which defeats the entire method, and nothing about
the output *looks* broken. shanuz did exactly this until the sketching tutorial
compared it with R (see ``sketch_vignette.md``); truncating to 50 restores R's
spread (max/median 6.5 rather than 1.3).

:func:`leverage_score` is that computation on its own; :func:`sketch_data` draws
the weighted subset and returns it as a standalone object; :func:`project_data`
is the inverse map — it pushes every full-dataset cell through the *sketch's*
PCA loadings (and, if present, its fitted UMAP model, reusing
:func:`shanuz.project_umap`), and optionally carries the sketch's labels out to
every cell by the weighted k-nearest-neighbour vote Seurat uses, **not** through
integration anchors: anchor finding against the full dataset costs precisely what
sketching is for.
"""
from __future__ import annotations

import warnings
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
    ndims: Optional[int] = None,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    layer: str = "data",
    var_name: Optional[str] = "leverage.score",
    eps: float = 0.5,
    seed: int = 123,
) -> np.ndarray:
    """Per-cell statistical leverage (Seurat's ``LeverageScore``).

    A cell's leverage is how much it influences the column space of the data — low
    in a dense, redundant cloud, high in a sparse, distinctive corner — so sampling
    proportional to it keeps the rare states a uniform draw would lose.

    Which of the two regimes runs is decided exactly as Seurat decides it. With
    fewer than ``nsketch * 1.5`` cells the scores come from a rank-``50`` truncated
    SVD and sum to 50; above that a ``CountSketch`` embedding, a ``QR`` and a
    Johnson–Lindenstrauss projection stand in for the SVD, and the scores are on
    the projection's scale rather than summing to 50. Compare scores *within* one
    call, never across the two regimes.

    Parameters
    ----------
    obj      : a :class:`~shanuz.Shanuz` object (normalized).
    nsketch  : rows of the random sketch, and the threshold that picks the regime.
    ndims    : dimension the JL projection targets before ``eps`` shrinks it
               (default: the cell count, as in Seurat). Sketched regime only.
    features : features to score on (default: the assay's variable features).
    assay    : assay to use (default: active assay).
    layer    : layer to draw the data from. Defaults to ``"data"`` — the
               log-normalized values, which is what Seurat scores — *not*
               ``"scale.data"``.
    var_name : if given, the scores are also written to ``obj.meta_data[var_name]``.
    eps      : Johnson–Lindenstrauss distortion, ``0 < eps <= 1`` (Seurat's 0.5).
               Smaller keeps more projected dimensions. Sketched regime only.
    seed     : random seed for the sketch and the projection.

    Returns
    -------
    numpy.ndarray
        One leverage score per cell, in ``obj.cell_names()`` order.
    """
    assay_name = assay or obj.active_assay
    assay_obj = obj.assays[assay_name]
    feats = _default_features(assay_obj, features)

    # (n_features × n_cells), left sparse when it already is: sketching exists to
    # keep large data cheap, and densifying here would give that away up front.
    A = _leverage_matrix(assay_obj, feats, layer)
    n_cells = A.shape[1]

    if n_cells < nsketch * 1.5:
        scores = _leverage_exact(A)
    else:
        rng = np.random.default_rng(seed)
        scores = _leverage_sketched(A, nsketch, ndims or n_cells, eps, rng)

    if var_name is not None:
        obj.meta_data[var_name] = scores
    return scores


def sketch_data(
    obj,
    ncells: int = 5000,
    method: str = "LeverageScore",
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    layer: str = "data",
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
    obj            : a :class:`~shanuz.Shanuz` object (normalized).
    ncells         : cells to keep (capped at the number available).
    method         : ``"LeverageScore"`` (leverage-weighted) or ``"Uniform"``
                     (equal weights), as in Seurat. ``"Uniform"`` is the control
                     that shows what leverage weighting is buying.
    features       : features to score on (default: variable features).
    assay          : assay to use (default: active assay).
    layer          : layer to draw the data from (default ``"data"``).
    nsketch        : sketch size passed to :func:`leverage_score`.
    sketched_assay : name the returned object's active assay is renamed to.
    var_name       : metadata column the scores are written to on ``obj``.
    seed           : random seed for scoring and sampling.

    Returns
    -------
    Shanuz
        The sketched subset (a new object).
    """
    if method not in ("LeverageScore", "Uniform"):
        raise ValueError(
            f"Unknown sketch method {method!r}; expected 'LeverageScore' or 'Uniform'."
        )

    n = len(obj)
    if method == "Uniform":
        scores = np.ones(n)
        if var_name is not None:
            obj.meta_data[var_name] = scores
    else:
        scores = leverage_score(
            obj, nsketch=nsketch, features=features, assay=assay, layer=layer,
            var_name=var_name, seed=seed,
        )

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
    refdata: Optional[Union[str, dict]] = None,
    project_umap: bool = True,
    dims: Optional[Union[list[int], range]] = None,
    k_weight: int = 50,
    sd_weight: float = 1.0,
    layer: str = "scale.data",
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
    3. **Labels** (when ``refdata`` is given): a weighted k-nearest-neighbour vote
       *inside the projected reduction*, where the sketch's own rows are the
       reference — Seurat's ``TransferSketchLabels``. Written onto
       ``full.meta_data``.

    Step 3 is deliberately **not** the :mod:`shanuz.transfer` anchor path, which
    is what an earlier version of this function used. Seurat does not use anchors
    here, and the difference is not academic: finding anchors between the sketch
    and the full dataset costs exactly what sketching exists to avoid, so on the
    million-cell objects this is written for the anchor route is unusable rather
    than merely different. On ifnb the two now agree per-cell **98.1 %** of the
    time, at matching accuracy.

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
    refdata             : sketch metadata to transfer, as in Seurat: a ``dict``
                          ``{new_col: sketch_col}`` writes each label under
                          ``new_col`` plus ``new_col.score``, and a bare ``str``
                          is shorthand for ``{col: col}``. Must name a column on
                          the *sketch* — like R, raw label arrays are not taken.
                          ``None`` skips transfer.
    project_umap        : also project the sketch's UMAP when a fitted model exists.
    dims                : reduction dimensions used for the UMAP model and the
                          label vote (default: all).
    k_weight            : neighbours each cell votes over (Seurat's ``k.weight``).
    sd_weight           : bandwidth of the distance kernel (Seurat fixes this at 1).
    layer               : layer to draw the full data's expression from.

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
        _transfer_from_sketch(
            full, sketch, full_reduction, refdata, k_weight, sd_weight, dims
        )

    return full


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


# Seurat truncates the exact regime's SVD at this many components, and the choice
# is load-bearing: the scores sum to it, so it sets how much spread the sampling
# weights can have. See the module docstring.
_EXACT_NDIMS = 50


def _leverage_matrix(assay_obj, features: list[str], layer: str):
    """The ``(features × cells)`` matrix to score, sparse if the layer is sparse."""
    from .assay5 import Assay5

    if isinstance(assay_obj, Assay5) and layer in assay_obj.layers:
        mat = assay_obj.layers[layer]
        # Only scale.data is stored against a feature subset; every other layer is
        # indexed by the assay's full feature list. Mixing the two silently reads
        # the wrong rows.
        stored = assay_obj._all_feature_names
        if layer == "scale.data":
            stored = getattr(assay_obj, "_scaled_features", None) or stored
        pos = {f: i for i, f in enumerate(stored)}
        idx = [pos[f] for f in features if f in pos]
    elif not isinstance(assay_obj, Assay5) and layer in ("data", "counts"):
        mat = assay_obj.data if layer == "data" else assay_obj.counts
        pos = {f: i for i, f in enumerate(assay_obj._feature_names)}
        idx = [pos[f] for f in features if f in pos]
    else:
        # Anything else (notably "scale.data") goes through the shared accessor,
        # which handles the fallbacks and returns dense.
        return np.asarray(_get_scaled_data(assay_obj, features, layer), dtype=float)

    if not idx:
        raise ValueError(
            f"None of the {len(features)} requested features are present in layer "
            f"{layer!r}; nothing to score."
        )
    if sp.issparse(mat):
        return sp.csr_matrix(mat)[idx, :].astype(float)
    return np.asarray(mat, dtype=float)[idx, :]


def _count_sketch(n: int, nsketch: int, rng: np.random.Generator) -> sp.csr_matrix:
    """A CountSketch subspace embedding — an ``nsketch × n`` sparse sign matrix.

    Each of the ``n`` columns has exactly one non-zero, a random ``±1`` in a random
    row, so ``S @ A`` hashes the cells into ``nsketch`` signed buckets in a single
    pass over the data's non-zeros. Mirrors Seurat's ``CountSketch``.
    """
    rows = rng.integers(0, nsketch, size=n)
    signs = rng.choice(np.array([-1.0, 1.0]), size=n)
    cols = np.arange(n)
    return sp.csr_matrix((signs, (rows, cols)), shape=(nsketch, n))


def _jl_embed(nrow: int, ncol: int, eps: float, rng: np.random.Generator) -> np.ndarray:
    """Seurat's ``JLEmbed`` — Li's sparse Johnson–Lindenstrauss projection.

    ``eps`` first *replaces* the requested width with the JL bound
    ``4 log(ncol) / (eps²/2 - eps³/3)``, then entries are drawn from
    ``{-1, 0, +1}`` with ``P(±1) = 1/2s`` for ``s = ceil(sqrt(ncol))``. As in
    Seurat the result is left unscaled, so the row norms it produces are on the
    projection's own scale rather than the exact regime's.
    """
    if not 0.0 < eps <= 1.0:
        raise ValueError(f"eps must satisfy 0 < eps <= 1; got {eps!r}")
    ncol = int(np.floor(4.0 * np.log(ncol) / ((eps ** 2) / 2.0 - (eps ** 3) / 3.0)))
    s = int(np.ceil(np.sqrt(ncol)))
    probs = [1.0 / (2 * s), 1.0 - 1.0 / s, 1.0 / (2 * s)]
    return rng.choice(np.array([-1.0, 0.0, 1.0]), size=(nrow, ncol), p=probs)


def _leverage_exact(A) -> np.ndarray:
    """Leverage from a rank-50 truncated SVD of ``A`` (``features × cells``).

    Seurat's small-data branch: ``rowSums(irlba(A, nv = 50)$v ^ 2)``. The right
    singular vectors are per-*cell*, so no transpose is needed and the scores sum
    to the number of components kept.

    Seurat hardcodes 50 and would fail on a matrix with fewer rows than that (it
    computes the smaller fallback and then does not use it); the fallback its own
    code names is applied here instead.
    """
    n_features, n_cells = A.shape
    k = min(_EXACT_NDIMS, min(n_features, n_cells) - 1)
    if k < 1:
        return np.zeros(n_cells)

    if sp.issparse(A):
        from scipy.sparse.linalg import svds

        # svds returns ascending singular values; the leading k are the last k.
        _, _, vt = svds(A, k=k)
    else:
        _, _, vt = np.linalg.svd(A, full_matrices=False)
        vt = vt[:k]
    return np.einsum("ij,ij->j", vt, vt)


def _leverage_sketched(A, nsketch: int, ndims: int, eps: float,
                       rng: np.random.Generator) -> np.ndarray:
    """Leverage via CountSketch → QR → JL, Seurat's large-data branch.

    ``S A ᵀ`` is a small subspace embedding of the cells-by-features matrix, so its
    triangular factor ``R`` satisfies ``RᵀR ≈ AᵀA`` and ``Aᵀ R⁻¹`` is approximately
    orthonormal. Its row norms are then read off a JL projection rather than formed
    directly, which is what keeps the cost independent of the cell count.
    """
    n_features, n_cells = A.shape

    # Seurat's two hard stops, ported as-is: the sketch has to be a *reduction*,
    # and the QR needs more rows than columns to give a square R.
    if n_features > 5000:
        raise ValueError(
            f"Scoring {n_features} features in the sketched regime is too slow; "
            f"Seurat refuses above 5000. Pass a smaller `features` list "
            f"(the variable features are the usual choice)."
        )
    if n_features > n_cells / 1.1:
        raise ValueError(
            f"Matrix is too square to sketch: {n_features} features against "
            f"{n_cells} cells. Sketching needs many more cells than features."
        )
    if nsketch < 1.1 * n_features:
        bumped = int(np.ceil(1.1 * n_features))
        warnings.warn(
            f"nsketch ({nsketch}) is too close to the number of features "
            f"({n_features}); raising it to {bumped}.",
            stacklevel=3,
        )
        nsketch = bumped
    nsketch = min(nsketch, ndims)

    At = A.T                                            # (cells × features)
    sa = _count_sketch(At.shape[0], nsketch, rng) @ At  # (nsketch × features)
    sa = np.asarray(sa.todense() if sp.issparse(sa) else sa, dtype=float)

    # mode="r" returns the R factor alone; asarray pins that for the type checker,
    # which types np.linalg.qr by its multi-mode signature.
    r = np.asarray(np.linalg.qr(sa, mode="r"))

    jl = _jl_embed(r.shape[1], ndims, eps, rng)
    # solve_triangular(R, JL) is R⁻¹ @ JL without ever forming the inverse.
    from scipy.linalg import solve_triangular

    proj = solve_triangular(r, jl, lower=False)
    z = np.asarray(At @ proj)
    return np.einsum("ij,ij->i", z, z)


def _sketch_weight_matrix(query_emb, ref_emb, k_weight: int, sd_weight: float):
    """Seurat's ``FindWeightsNN`` — each full cell weighted over its nearest sketch cells.

    Three steps, matching R exactly:

    1. the ``k`` nearest *sketch* cells for every cell, in the projected space;
    2. distances rescaled per row against the furthest of those k,
       ``d' = 1 - d/d_k``, so the nearest neighbour scores 1 and the furthest 0;
    3. the Seurat kernel ``w = 1 - exp(-d' / (2/sd)²)``, normalized to sum to 1.

    Verified against ``Seurat:::FindWeightsC`` term for term: distances
    ``(0, 0.25, 0.5, 1.0)`` give ``(0, 0.151737, 0.2942806, 0.5539824)`` in both.
    """
    from sklearn.neighbors import NearestNeighbors

    k = int(min(k_weight, ref_emb.shape[0]))
    if k < 1:
        raise ValueError("the sketch has no cells to transfer from")

    nn = NearestNeighbors(n_neighbors=k).fit(ref_emb)
    dist, idx = nn.kneighbors(query_emb)

    furthest = dist[:, -1:].copy()
    furthest[furthest == 0] = 1.0            # a cell sitting on all k neighbours
    scaled = 1.0 - dist / furthest
    w = 1.0 - np.exp(-scaled / (2.0 / sd_weight) ** 2)
    total = w.sum(axis=1, keepdims=True)
    total[total == 0] = 1.0
    return w / total, idx


def _transfer_from_sketch(full, sketch, full_reduction: str, refdata,
                          k_weight: int, sd_weight: float, dims) -> None:
    """Carry ``refdata`` from the sketch onto every cell (``TransferSketchLabels``).

    Seurat's ``ProjectData`` does **not** transfer through integration anchors — it
    runs a weighted k-nearest-neighbour vote inside the projected reduction, where
    the sketch's own rows act as the reference. That distinction is not cosmetic:
    anchor finding between the sketch and the full dataset costs what sketching
    exists to avoid, so on the million-cell objects this is written for the anchor
    route is not merely different but unusable.
    """
    emb = np.asarray(full.reductions[full_reduction].cell_embeddings, dtype=float)
    if dims is not None:
        emb = emb[:, list(dims)]

    full_cells = full.cell_names()
    position = {c: i for i, c in enumerate(full_cells)}
    ref_rows = [position[c] for c in sketch.cell_names() if c in position]
    if not ref_rows:
        raise ValueError(
            "none of the sketch's cells are present in the full object — the two "
            "were not drawn from the same data."
        )

    weights, neighbours = _sketch_weight_matrix(emb, emb[ref_rows], k_weight, sd_weight)
    sketch_meta = sketch.meta_data
    meta = full.meta_data

    items = refdata.items() if isinstance(refdata, dict) else [(str(refdata), refdata)]
    for out_col, src_col in items:
        if src_col not in sketch_meta.columns:
            raise KeyError(f"refdata column {src_col!r} not found in the sketch's meta_data.")
        labels = sketch_meta.loc[sketch.cell_names(), src_col].to_numpy().astype(str)

        classes = sorted(set(labels))
        index = {c: i for i, c in enumerate(classes)}
        onehot = np.zeros((len(labels), len(classes)))
        onehot[np.arange(len(labels)), [index[x] for x in labels]] = 1.0

        # Weighted vote: sum each neighbour's weight into its label's column.
        scores = np.zeros((emb.shape[0], len(classes)))
        for j in range(neighbours.shape[1]):
            np.add.at(scores, (np.arange(emb.shape[0]), onehot[neighbours[:, j]].argmax(1)),
                      weights[:, j])

        meta[out_col] = np.asarray(classes, dtype=object)[scores.argmax(axis=1)]
        meta[f"{out_col}.score"] = scores.max(axis=1)
