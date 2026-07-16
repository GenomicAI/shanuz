"""MULTI-seq demultiplexing — Seurat's ``MULTIseqDemux``.

MULTI-seq (McGinnis, Patterson et al. 2019) is a second cell-hashing chemistry:
lipid- or cholesterol-anchored barcode oligos label each *sample* before pooling,
exactly as antibody hashtags do in Cell Hashing. Demultiplexing is the same
problem as :func:`shanuz.hto_demux` solves — turn the per-cell barcode count
matrix back into singlet / doublet / negative calls — but MULTI-seq reaches the
answer a different way, and the two methods often disagree at the margins, so
having both is useful.

Where ``HTODemux`` learns each tag's cutoff by fitting a **negative binomial** to
its background cluster, ``MULTIseqDemux`` reads the cutoff straight off the shape
of each barcode's distribution:

1. **Normalize.** Barcode counts are compositional, so they are centered-log-ratio
   (CLR) normalized — the same transform :func:`shanuz.normalize_data` applies with
   ``method="CLR"`` (margin 1, per barcode across cells, as Seurat's default does). As with
   :func:`shanuz.hto_demux`, ``normalize=True`` (default) recomputes this from the
   raw counts so the function works straight after object creation; ``normalize=False``
   uses the assay's existing ``data`` layer, matching Seurat's expectation that you
   normalize first.
2. **Find each barcode's threshold.** For one barcode, the CLR values across all
   cells form a bimodal distribution — a tall low mode (cells the barcode never
   stained: its background) and a shorter high mode (cells it did stain). A smooth
   Gaussian kernel density estimate over a 100-point grid exposes those two modes
   as the two tallest local maxima. The threshold is placed a fraction ``quantile``
   of the way from the low mode to the high mode (``quantile=0.7`` sits 70 % toward
   the positive peak). A cell is positive for the barcode when its CLR value clears
   that threshold. Barcodes with no discernible second mode get no threshold and
   mark no cells (Seurat likewise skips them).
3. **Classify.** Count how many barcodes each cell clears: zero → ``Negative``,
   one → that barcode (``Singlet``), more than one → ``Doublet``.

The ``quantile`` that best separates the data is not known a priori. With
``autothresh=True`` the classifier sweeps ``qrange`` (Seurat's default
``0.1 … 0.9``), picks the ``q`` that yields the most singlets, strips the cells it
calls negative, and repeats on the remainder — up to ``maxiter`` rounds — so that
peeling away empty droplets sharpens the modes for the cells that remain. This is
the iterative "semi-supervised" thresholding from McGinnis's ``deMULTIplex``.

Results are written to ``obj.meta_data`` under Seurat's names: ``MULTI_ID`` — the
barcode name for singlets, ``"Doublet"``, or ``"Negative"`` — which is also set as
the active identity (so ``subset`` can pull a sample's singlets immediately), and a
character copy in ``MULTI_classification``. The learned per-barcode thresholds are
stashed in ``obj.misc["multiseq_demux"]``.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .hto import _hto_matrices


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def multiseq_demux(
    seurat,
    assay: str = "HTO",
    quantile: float = 0.7,
    autothresh: bool = False,
    maxiter: int = 5,
    qrange: Optional[Sequence[float]] = None,
    normalize: bool = True,
    margin: int = 1,
    verbose: bool = False,
):
    """Demultiplex pooled samples from barcode counts (Seurat's ``MULTIseqDemux``).

    Mirrors ``MULTIseqDemux(object, assay = "HTO", quantile = 0.7)``. For each
    barcode a Gaussian-kernel-density threshold is placed a fraction ``quantile``
    of the way between its background and positive modes; cells positive for zero /
    one / many barcodes are called ``Negative`` / ``Singlet`` / ``Doublet``.

    The object is mutated in place: ``MULTI_ID`` and ``MULTI_classification``
    metadata columns are written, the active identity is set to ``MULTI_ID``, and
    the learned thresholds are stashed in ``obj.misc["multiseq_demux"]``.

    Parameters
    ----------
    seurat     : a :class:`~shanuz.Shanuz` object carrying a hashtag/barcode assay.
    assay      : the barcode assay to demultiplex (default ``"HTO"``).
    quantile   : fraction between each barcode's background and positive modes at
                 which its positive cutoff is placed (Seurat default 0.7). Ignored
                 when ``autothresh`` is True.
    autothresh : sweep ``qrange`` for the quantile that maximizes the singlet rate,
                 iteratively removing negatives and re-thresholding the remainder
                 (up to ``maxiter`` rounds). Overrides ``quantile``.
    maxiter    : maximum auto-threshold rounds (default 5).
    qrange     : quantiles swept when ``autothresh`` is True (default
                 ``0.1, 0.15, … 0.9``).
    normalize  : CLR-normalize the counts internally (default). Set False to use the
                 assay's existing ``data`` layer (e.g. a prior
                 ``normalize_data(method="CLR")``).
    margin     : CLR margin when ``normalize`` is True — 1 (per barcode across cells;
                 Seurat's default) or 2 (per cell across barcodes).
    verbose    : print each auto-threshold round's chosen quantile.

    Returns
    -------
    Shanuz
        ``seurat``, with the ``MULTI_ID`` classification and identity.
    """
    _, data, feats, cells = _hto_matrices(seurat, assay, normalize, margin)
    n_bc, n_cells = data.shape
    if n_bc < 2:
        raise ValueError(
            f"MULTIseqDemux needs at least 2 barcodes; assay {assay!r} has {n_bc}."
        )

    if qrange is None:
        qrange = np.round(np.arange(0.1, 0.9 + 1e-9, 0.05), 4)
    else:
        qrange = np.asarray(qrange, dtype=float)

    if autothresh:
        calls, thresholds, q_used = _auto_classify(
            data, feats, qrange, maxiter, verbose
        )
    else:
        calls, thresholds = _classify_cells(data, feats, quantile)
        q_used = quantile

    _write_calls(seurat, calls, cells)
    seurat.misc.setdefault("multiseq_demux", {})[assay] = {
        "thresholds": thresholds,
        "quantile": float(q_used),
        "autothresh": bool(autothresh),
    }
    return seurat


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _local_maxima(y: np.ndarray) -> np.ndarray:
    """Indices of local maxima of ``y`` (mirrors Seurat's ``LocalMaxima``, which
    pads with ``-inf`` so a peak sitting on either boundary is still found)."""
    padded = np.concatenate(([-np.inf], np.asarray(y, dtype=float)))
    rising = np.diff(padded) > 0                # rising[i] True if y rose into i
    n = y.size
    peaks = [
        i for i in range(n)
        if rising[i] and (i == n - 1 or not rising[i + 1])
    ]
    return np.asarray(peaks, dtype=int)


def _barcode_threshold(x: np.ndarray, q: float) -> Optional[float]:
    """Positive cutoff for one barcode: a fraction ``q`` of the way from its low
    (background) mode to its high (positive) mode, read off a Gaussian KDE.

    Returns ``None`` when the barcode has no spread or shows fewer than two modes —
    i.e. no threshold can be learned — in which case it marks no cells positive.
    """
    from scipy.stats import gaussian_kde

    x = np.asarray(x, dtype=float)
    if x.size < 2 or np.ptp(x) <= 0:
        return None
    try:
        kde = gaussian_kde(x)
    except Exception:
        return None

    grid = np.linspace(x.min(), x.max(), 100)
    dens = kde(grid)
    peaks = _local_maxima(dens)
    if peaks.size < 2:
        return None

    # The two tallest modes are the background and positive populations; place the
    # cutoff a fraction q between their positions (quantile of the two, as Seurat).
    top2 = peaks[np.argsort(dens[peaks])[-2:]]
    lo, hi = np.sort(grid[top2])
    return float(lo + q * (hi - lo))


def _classify_cells(data: np.ndarray, feats, q: float):
    """One thresholding pass at quantile ``q``.

    Returns ``(calls, thresholds)`` — a per-cell array of barcode name /
    ``"Doublet"`` / ``"Negative"``, and the per-barcode cutoff dict (``inf`` for a
    barcode that yielded no threshold).
    """
    n_bc, n_cells = data.shape
    thresholds: dict[str, float] = {}
    discrete = np.zeros((n_bc, n_cells), dtype=bool)
    for i in range(n_bc):
        thr = _barcode_threshold(data[i], q)
        if thr is None:
            thresholds[feats[i]] = float("inf")
            continue
        thresholds[feats[i]] = thr
        discrete[i] = data[i] > thr

    npos = discrete.sum(axis=0)
    calls = np.empty(n_cells, dtype=object)
    for j in range(n_cells):
        if npos[j] == 0:
            calls[j] = "Negative"
        elif npos[j] == 1:
            calls[j] = feats[int(np.argmax(discrete[:, j]))]
        else:
            calls[j] = "Doublet"
    return calls, thresholds


def _find_best_q(data: np.ndarray, feats, qrange) -> float:
    """The quantile in ``qrange`` that classifies the most cells as singlets
    (Seurat maximizes ``pSinglet`` over the sweep; ties keep the earliest ``q``)."""
    n_cells = data.shape[1]
    best_q = float(qrange[0])
    best_frac = -1.0
    for q in qrange:
        calls, _ = _classify_cells(data, feats, float(q))
        nsing = int(np.sum((calls != "Negative") & (calls != "Doublet")))
        frac = nsing / n_cells if n_cells else 0.0
        if frac > best_frac:
            best_frac = frac
            best_q = float(q)
    return best_q


def _auto_classify(data: np.ndarray, feats, qrange, maxiter: int, verbose: bool):
    """Iterative auto-thresholding: sweep for the best quantile, peel off the
    cells it calls negative, and re-threshold the remainder until nothing new is
    negative or ``maxiter`` rounds elapse.

    Returns ``(calls, thresholds, q_used)`` where ``thresholds``/``q_used`` are
    from the final round.
    """
    n_bc, n_cells = data.shape
    remaining = np.arange(n_cells)
    final = np.array(["Negative"] * n_cells, dtype=object)
    thresholds: dict[str, float] = {feats[i]: float("inf") for i in range(n_bc)}
    q_used = float(qrange[0])

    for it in range(1, maxiter + 1):
        sub = data[:, remaining]
        q_used = _find_best_q(sub, feats, qrange)
        calls, thresholds = _classify_cells(sub, feats, q_used)
        final[remaining] = calls

        neg_local = np.where(calls == "Negative")[0]
        if verbose:
            print(
                f"[multiseq] round {it}: q={q_used:g}, "
                f"{neg_local.size} new negatives"
            )
        if neg_local.size == 0:
            break
        remaining = np.setdiff1d(remaining, remaining[neg_local])
        if remaining.size == 0:
            break

    return final, thresholds, q_used


def _write_calls(seurat, calls, cells) -> None:
    """Write ``MULTI_ID`` / ``MULTI_classification`` and set the active identity."""
    target = seurat.cell_names()
    aligned = pd.Series(list(calls), index=cells).reindex(target).values
    seurat.meta_data["MULTI_ID"] = list(aligned)
    seurat.meta_data["MULTI_classification"] = list(aligned)
    seurat.idents = list(aligned)
