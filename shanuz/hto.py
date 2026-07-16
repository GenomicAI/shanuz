"""Cell-hashing demultiplexing — Seurat's ``HTODemux``.

Cell Hashing (Stoeckius et al. 2018) tags each *sample* with a distinct
antibody-conjugated oligo — a **hashtag** (HTO) — before the samples are pooled
and run together on one droplet lane. Every droplet then carries, on top of its
mRNA, a little barcode of hashtag counts saying which sample the cell came from.
Demultiplexing turns that hashtag matrix back into a per-cell assignment: this
cell is a **singlet** from sample 3, that droplet caught two cells and is a
**doublet** of samples 1 and 2, this empty-ish droplet is **negative**. Pooling
many samples this way is cheap and cancels batch effects, but only if the
demultiplexing is trustworthy — that is what :func:`hto_demux` provides.

The hard part is picking, per hashtag, the count above which a cell is genuinely
"positive" for that tag. A fixed threshold fails because every antibody has its
own staining background and every experiment its own depth. Seurat's ``HTODemux``
learns the threshold from the data:

1. **Normalize.** Hashtag counts are compositional (what matters is a cell's
   *relative* tag levels, not its total), so they are centered-log-ratio (CLR)
   normalized — the same transform :func:`shanuz.normalize_data` applies with
   ``method="CLR"``. Margin 1 — per hashtag, across cells — is what the hashing
   vignette uses and is this module's default; it is Seurat's default too.
2. **Cluster.** ``k = n_hashtags + 1`` groups split the cells into one high
   cluster per hashtag plus a background cluster, by either k-means or ``clara``
   (k-medoids; Seurat's default — see :mod:`shanuz._clara`). For any given
   hashtag, the cluster with the *lowest* average expression of that tag is its
   **negative** population — real cells that were never stained by this antibody,
   i.e. a clean sample of the tag's background. Only that ranking feeds the next
   step, which is why the choice of clustering rarely changes the calls.
3. **Fit a background model.** A negative binomial is fit (maximum likelihood) to
   the raw counts of that hashtag in its negative cluster, and the
   ``positive.quantile`` (default 0.99) of the fitted distribution becomes the
   cutoff. A cell is positive for the hashtag when its raw count exceeds the
   cutoff. The negative binomial — not Poisson — is used because antibody
   background is overdispersed; when the data happen to be Poisson-like the fit
   degenerates to Poisson on its own.
4. **Classify.** Count how many hashtags each cell is positive for: zero →
   ``Negative``, one → ``Singlet`` (assigned to that hashtag), more than one →
   ``Doublet`` (labelled by its top two tags).

The per-cell results are written to ``obj.meta_data`` under the Seurat column
names (``<assay>_maxID``, ``<assay>_secondID``, ``<assay>_margin``,
``<assay>_classification``, ``<assay>_classification.global``) plus a convenient
``hash.ID`` — the hashtag name for singlets, ``"Doublet"``, or ``"Negative"`` —
which is also set as the active identity so ``subset`` can immediately pull the
singlets of a given sample.

Seurat also ships ``MULTIseqDemux`` (McGinnis et al., a quantile-sweep
alternative to the negative-binomial threshold); it shares this module's
normalize-then-threshold skeleton and is the natural follow-on.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .preprocessing import _clr_normalize


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def hto_demux(
    seurat,
    assay: str = "HTO",
    positive_quantile: float = 0.99,
    init: Optional[int] = None,
    nstarts: int = 10,
    kfunc: str = "kmeans",
    nsamples: int = 100,
    normalize: bool = True,
    margin: int = 1,
    seed: int = 42,
    verbose: bool = False,
):
    """Demultiplex pooled samples from hashtag counts (Seurat's ``HTODemux``).

    Mirrors ``HTODemux(object, assay = "HTO", positive.quantile = 0.99)``. Each
    hashtag's positive/negative cutoff is learned by fitting a negative binomial
    to the tag's background — the k-means cluster in which it is least expressed —
    and thresholding at ``positive_quantile``. Cells positive for zero / one /
    many hashtags are called ``Negative`` / ``Singlet`` / ``Doublet``.

    The object is mutated in place: five ``<assay>_*`` metadata columns plus
    ``hash.ID`` are written (matching Seurat), the active identity is set to
    ``hash.ID``, and the learned cutoffs are stashed in ``obj.misc["hto_demux"]``.

    Parameters
    ----------
    seurat            : a :class:`~shanuz.Shanuz` object carrying a hashtag assay.
    assay             : the hashtag assay to demultiplex (default ``"HTO"``).
    positive_quantile : quantile of each tag's fitted background at which the
                        positive cutoff is set (Seurat default 0.99).
    init              : number of k-means centers; default ``n_hashtags + 1``.
    nstarts           : k-means restarts (``n_init``). Seurat uses 100; 10 is a
                        faster, usually-equivalent default. Ignored by ``clara``.
    kfunc             : ``"kmeans"`` (default) or ``"clara"``, Seurat's k-medoids.
                        **Note shanuz's default differs from Seurat's**, which is
                        ``"clara"``; pass ``kfunc="clara"`` to follow R. The two
                        rarely disagree on which cluster is a tag's background, so
                        the calls usually match either way.
    nsamples          : ``clara`` only — sub-samples to draw (Seurat's default,
                        100). Ignored by ``kmeans``.
    normalize         : CLR-normalize the counts internally for clustering and
                        margins (default). Set False to use the assay's existing
                        ``data`` layer (e.g. a prior ``normalize_data(method="CLR")``).
    margin            : CLR margin when ``normalize`` is True — 1 (per hashtag
                        across cells; Seurat's default, and what the hashing
                        vignette normalizes with) or 2 (per cell across hashtags).
    seed              : random seed for k-means. Has no effect on ``clara``, which
                        draws from its own generator that R cannot seed either —
                        see :mod:`shanuz._clara`.
    verbose           : print each hashtag's learned cutoff.

    Returns
    -------
    Shanuz
        ``seurat``, with the classification metadata and ``hash.ID`` identity.
    """
    if kfunc not in ("kmeans", "clara"):
        raise NotImplementedError(
            f"kfunc={kfunc!r} is not supported; choose from 'kmeans' or 'clara'."
        )

    counts, data, feats, cells = _hto_matrices(seurat, assay, normalize, margin)
    n_htos, n_cells = data.shape
    if n_htos < 2:
        raise ValueError(
            f"HTODemux needs at least 2 hashtags; assay {assay!r} has {n_htos}."
        )

    labels, ncenters = _cluster_cells(data, init, nstarts, seed, kfunc, nsamples)

    # Average (de-logged) expression of each hashtag within each cluster, so the
    # least-expressing cluster can be read off as that hashtag's background.
    expd = np.expm1(data)
    avg = np.full((n_htos, ncenters), np.inf)
    for c in range(ncenters):
        mask = labels == c
        if mask.any():
            avg[:, c] = expd[:, mask].mean(axis=1)

    # Per-hashtag negative-binomial threshold on its background cluster.
    discrete = np.zeros((n_htos, n_cells), dtype=bool)
    cutoffs: dict[str, float] = {}
    for i in range(n_htos):
        neg_cluster = int(np.argmin(avg[i]))
        values_use = counts[i, labels == neg_cluster]
        cutoff = _positive_cutoff(values_use, positive_quantile)
        cutoffs[feats[i]] = cutoff
        discrete[i] = counts[i] > cutoff
        if verbose:
            print(f"Cutoff for {feats[i]}: {cutoff:g} reads")

    _write_classification(
        seurat, assay, data, discrete, feats, cells,
    )
    seurat.misc.setdefault("hto_demux", {})[assay] = {
        "cutoffs": cutoffs,
        "ncenters": ncenters,
        "positive_quantile": positive_quantile,
    }
    return seurat


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _dense(mat) -> np.ndarray:
    """Densify a possibly-sparse matrix to a float ``ndarray`` (hashtag panels
    are only a handful of features, so densifying is cheap)."""
    if sp.issparse(mat):
        return np.asarray(mat.toarray(), dtype=float)
    return np.asarray(mat, dtype=float)


def _hto_matrices(seurat, assay, normalize, margin):
    """Return ``(counts, data, feature_names, cell_names)`` for the hashtag assay.

    ``counts`` are the raw hashtag counts (for the background fit) and ``data`` is
    the CLR-normalized matrix (for clustering and margins) — either recomputed
    from the counts (``normalize``) or taken from the assay's ``data`` layer.
    """
    from .assay5 import Assay5
    from ._sparse import is_matrix_empty

    if assay not in seurat.assays:
        raise KeyError(f"Assay {assay!r} not found on the object.")
    assay_obj = seurat.assays[assay]

    if isinstance(assay_obj, Assay5):
        feats = list(assay_obj._all_feature_names)
        cells = list(assay_obj._all_cell_names)
        counts = assay_obj.layers.get("counts")
        if counts is None:
            counts = assay_obj.layers.get("data")
        data_layer = assay_obj.layers.get("data")
    else:
        feats = list(assay_obj._feature_names)
        cells = list(assay_obj._cell_names)
        counts = assay_obj.counts
        data_layer = None if is_matrix_empty(assay_obj.data) else assay_obj.data

    if counts is None:
        raise ValueError(f"Assay {assay!r} has no counts to demultiplex.")

    counts = _dense(counts)
    if normalize or data_layer is None:
        data = _clr_normalize(counts, margin=margin)
    else:
        data = _dense(data_layer)
    return counts, data, feats, cells


def _cluster_cells(
    data: np.ndarray, init, nstarts: int, seed: int, kfunc: str, nsamples: int,
):
    """Cluster the cells (columns) into ``init or n_hashtags + 1`` groups.

    Returns ``(labels, ncenters)`` with one integer cluster label per cell. Only
    the split matters downstream: each hashtag's least-expressing cluster becomes
    its background, so ``kmeans`` and ``clara`` usually agree on the cutoffs even
    where they draw slightly different boundaries.
    """
    n_htos, n_cells = data.shape
    ncenters = init if init is not None else n_htos + 1
    ncenters = int(min(ncenters, n_cells))
    if ncenters < 2:
        raise ValueError("Too few cells to cluster for HTODemux.")

    if kfunc == "clara":
        from ._clara import clara

        # R: clara(x = t(GetAssayData(object, assay)), k = ncenters,
        #          samples = nsamples), leaving sampsize at its default.
        # clara needs k <= n - 1, which is tighter than k <= n above.
        if ncenters > n_cells - 1:
            raise ValueError(
                f"clara needs at least {ncenters + 1} cells to find {ncenters} "
                f"clusters; got {n_cells}."
            )
        return clara(data.T, k=ncenters, samples=nsamples), ncenters

    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=ncenters, n_init=nstarts, random_state=seed)
    return km.fit_predict(data.T), ncenters


def _fit_nbinom(x: np.ndarray) -> Optional[tuple[float, float]]:
    """Maximum-likelihood negative binomial fit, returning ``(size, mu)``.

    The MLE of the mean is the sample mean; the dispersion ``size`` (``r``) is
    found by a 1-D optimization of the profile log-likelihood. Under-dispersed
    (Poisson-like) data drive ``size`` to the upper bound, where the negative
    binomial converges to Poisson — exactly the desired degenerate behaviour.
    Returns ``None`` when the background mean is zero (no positive threshold to
    learn).
    """
    from scipy.optimize import minimize_scalar
    from scipy.special import gammaln

    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return None
    mu = float(x.mean())
    if mu <= 0:
        return None

    def neg_ll(log_r: float) -> float:
        r = np.exp(log_r)
        p = r / (r + mu)
        return -np.sum(
            gammaln(x + r) - gammaln(r) - gammaln(x + 1)
            + r * np.log(p) + x * np.log1p(-p)
        )

    res = minimize_scalar(
        neg_ll, bounds=(np.log(1e-4), np.log(1e8)), method="bounded"
    )
    return float(np.exp(res.x)), mu


def _positive_cutoff(values_use: np.ndarray, positive_quantile: float) -> float:
    """Positive cutoff for one hashtag: the ``positive_quantile`` of the negative
    binomial fit to its background counts (falls back sanely when degenerate)."""
    from scipy.stats import nbinom

    fit = _fit_nbinom(values_use)
    if fit is None:
        return 0.0
    r, mu = fit
    p = r / (r + mu)
    cutoff = float(nbinom.ppf(positive_quantile, r, p))
    if not np.isfinite(cutoff):
        return float(values_use.max()) if values_use.size else 0.0
    return cutoff


def _write_classification(seurat, assay, data, discrete, feats, cells) -> None:
    """Turn the per-hashtag positive calls into the Seurat metadata columns."""
    n_cells = data.shape[1]
    cols = np.arange(n_cells)

    # Top-two hashtags per cell (by normalized expression) and their margin.
    order = np.argsort(-data, axis=0)
    max_idx = order[0]
    second_idx = order[1]
    hash_max = data[max_idx, cols]
    hash_second = data[second_idx, cols]
    max_id = [feats[i] for i in max_idx]
    second_id = [feats[i] for i in second_idx]
    margin = hash_max - hash_second

    npos = discrete.sum(axis=0)
    global_class = np.where(
        npos == 0, "Negative", np.where(npos == 1, "Singlet", "Doublet")
    )

    classification = []
    hash_id = []
    for j in range(n_cells):
        g = global_class[j]
        if g == "Singlet":
            classification.append(max_id[j])
            hash_id.append(max_id[j])
        elif g == "Doublet":
            classification.append("_".join(sorted((max_id[j], second_id[j]))))
            hash_id.append("Doublet")
        else:
            classification.append("Negative")
            hash_id.append("Negative")

    target = seurat.cell_names()

    def put(col, values):
        seurat.meta_data[col] = (
            pd.Series(list(values), index=cells).reindex(target).values
        )

    put(f"{assay}_maxID", max_id)
    put(f"{assay}_secondID", second_id)
    put(f"{assay}_margin", margin)
    put(f"{assay}_classification", classification)
    put(f"{assay}_classification.global", list(global_class))
    put("hash.ID", hash_id)

    seurat.idents = list(
        pd.Series(hash_id, index=cells).reindex(target).values
    )
