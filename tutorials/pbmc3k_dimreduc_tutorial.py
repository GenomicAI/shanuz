"""Dimensional-reduction extras — Shanuz vs R Seurat on PBMC 3k.

The three reductions that sit *beside* the standard PCA → UMAP path, and the one
question every guided-clustering run has to answer: **how many PCs are real?**

  * :func:`shanuz.jack_straw` ↔ ``JackStraw`` — a permutation test for PC
    significance. A small fraction of features is scrambled across cells, the
    loadings are recomputed, and each observed loading gets an empirical p-value
    against that null.
  * :func:`shanuz.score_jackstraw` ↔ ``ScoreJackStraw`` — aggregates the
    per-feature p-values into one score per PC, so you can read off a cutoff.
  * :func:`shanuz.run_ica` ↔ ``RunICA`` — independent components instead of
    principal ones.
  * :func:`shanuz.run_tsne` ↔ ``RunTSNE`` — the older 2-D embedding, still the
    reference point for a lot of published figures.

Why this tutorial exists — and what it found
--------------------------------------------
All four landed with synthetic fixtures only and had never met a Seurat
reference. Building this comparison turned up **two real defects in JackStraw**,
both fixed in the same pull request as this tutorial:

1. **The null was too tight.** R's ``JackRandom`` permutes the chosen rows and
   **re-runs a full PCA**, taking the null loadings from that refit basis.
   Shanuz projected the permuted rows onto the **fixed** original embedding — a
   fixed basis cannot rotate to absorb the scrambled signal, so the permuted
   loadings came out too small and ordinary noise looked extreme against them.
   On the pure-noise PCs 14-20, that put **109-203** of 2000 features below
   p ≤ 1e-5 where R finds **0-5**.
2. **The aggregation was the wrong test.** R's ``ScoreJackStraw`` runs
   ``prop.test`` on the count of features below ``score.thresh`` against the
   count expected under a uniform null. Shanuz ran a one-sided KS test against
   Uniform(0, 1), which is enormously more sensitive: its largest score across
   all 20 pbmc3k PCs was **8.1e-112**, so no PC ever failed.

Together they made the function unable to do its one job: shanuz kept **all 20**
PCs where Seurat kept 13. After the fix both tools keep **13**. A third, smaller
gap is closed too — ``JackStrawData.fake_reduction_scores`` was declared but
never populated, where R stores the null.

The comparison that matters is the one an analyst acts on — the **set of PCs
kept** — so that is the headline, with the per-feature p-values underneath it to
show where any residual difference comes from.

What remains is permutation scatter, not a defect. R fixes each replicate's seed
to its loop index, so ``JackStraw`` there returns the same answer on every run;
shanuz seeds from its ``seed`` argument, and across seeds 0/1/7/42/2024 it keeps
13, 14, 15, 13 and 13 PCs. R's deterministic 13 sits at the bottom of that
spread — the scatter is the method's, not the port's.

ICA and t-SNE need a different yardstick. Independent components are defined
only up to sign and order, and t-SNE coordinates are not comparable across
implementations at all (R's ``Rtsne`` is Barnes-Hut, shanuz calls scikit-learn).
So neither is compared coordinate-wise:

  * **ICA** — components are Hungarian-matched between the tools by \\|Pearson r\\|,
    and we report the mean matched \\|r\\|. A faithful port recovers the same
    subspace, whatever order it names it in.
  * **t-SNE** — compared on *structure*: the shared-k-nearest-neighbour fraction
    between the two embeddings, and each embedding's kNN agreement with the PCA
    space it was built from. A t-SNE that preserves its input's neighbourhoods
    as well as R's does is doing its job, even at different coordinates.

Note on the shared basis
------------------------
JackStraw's null is built from the scaled matrix and the PCA basis, so the two
tools must start from the *same* PCA or a divergence downstream says nothing.
The Python run writes the exact cell barcodes and HVG list it used to
``figures_dimreduc/cells.txt`` and ``hvg_features.txt``; the R script subsets to
those cells and scales and runs PCA on those features. Step 0 of the comparison
then checks the two PCA bases actually agree (per-PC \\|correlation\\|) **before**
any JackStraw number is interpreted — same discipline as the Mixscape and
reference-mapping tutorials.

Usage
-----
    python tutorials/pbmc3k_dimreduc_tutorial.py    # downloads ~8 MB, writes hvg_features.txt

Then, for the side-by-side numbers and figures:

    Rscript tutorials/pbmc3k_dimreduc_verify.R      # slow: 100 PCA refits
    python  tutorials/generate_dimreduc_plots.py

References
----------
Chung NC, Storey JD (2015) **Statistical significance of variables driving
systematic variation in high-dimensional data.** Bioinformatics 31, 545-554.
https://doi.org/10.1093/bioinformatics/btu674

Macosko EZ, Basu A, Satija R, et al. (2015) **Highly parallel genome-wide
expression profiling of individual cells using nanoliter droplets.** Cell 161,
1202-1214. https://doi.org/10.1016/j.cell.2015.05.002
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import pbmc3k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import normalize_data, find_variable_features, scale_data
from shanuz.reduction import run_pca, run_ica, run_tsne
from shanuz.jackstraw import jack_straw, score_jackstraw

FIGURES = Path(__file__).parent / "figures_dimreduc"

# Preprocessing — the standard pbmc3k recipe, so the PCA basis is the familiar one.
N_HVG = 2000
N_PCS = 50

# JackStraw — Seurat's defaults, matched exactly on both sides.
JS_DIMS = 20
JS_REPLICATES = 100
JS_PROP_FREQ = 0.01
SCORE_THRESH = 1e-5          # R's ScoreJackStraw(score.thresh=)
ALPHA = 0.05                 # the cutoff an analyst actually applies

# ICA / t-SNE. 20 components rather than Seurat's default 50: FastICA on 2000
# features is the slowest step here, and 20 is plenty to establish the subspace.
N_ICS = 20
TSNE_DIMS = 10               # t-SNE runs on PC 1..10, as in the Seurat vignette
KNN_K = 30                   # neighbourhood size for the structure comparisons


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_dimreduc_tutorial.py)
# ---------------------------------------------------------------------------

def significant_dims(pvals, alpha: float = ALPHA) -> np.ndarray:
    """1-based indices of the dimensions whose score clears ``alpha``.

    Both tools return "small score ⇒ significant PC", so this is the same read
    on either side even though the underlying statistics differ.
    """
    p = np.asarray(pvals, dtype=float)
    if p.ndim != 1 or p.size == 0:
        raise ValueError("pvals must be a non-empty 1-D vector")
    return np.flatnonzero(p <= alpha) + 1


def n_leading_significant(pvals, alpha: float = ALPHA) -> int:
    """How many PCs to keep: the run of significant dims before the first failure.

    The number an analyst actually uses. Seurat's guidance is to cut at the
    drop-off rather than to keep every scattered PC that happens to clear the
    threshold, so a trailing significant PC after a gap does not extend the run.
    """
    p = np.asarray(pvals, dtype=float)
    if p.ndim != 1 or p.size == 0:
        raise ValueError("pvals must be a non-empty 1-D vector")
    failed = np.flatnonzero(p > alpha)
    return int(p.size if failed.size == 0 else failed[0])


def _corr_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Column-wise Pearson correlation between two (n × k) matrices."""
    a = a - a.mean(axis=0, keepdims=True)
    b = b - b.mean(axis=0, keepdims=True)
    sa = a.std(axis=0, ddof=0, keepdims=True)
    sb = b.std(axis=0, ddof=0, keepdims=True)
    sa[sa == 0] = 1.0                       # a constant component correlates with nothing
    sb[sb == 0] = 1.0
    return (a / sa).T @ (b / sb) / a.shape[0]


def matched_component_correlation(a, b) -> dict:
    """Match two sets of components one-to-one and report the matched \\|r\\|.

    ICA returns its components in arbitrary order and arbitrary sign, so a
    coordinate-wise comparison is meaningless — component 3 in R may be
    component 11, negated, in Python. Matching on \\|Pearson r\\| with the
    Hungarian algorithm asks the question that *is* well posed: did the two runs
    find the same components at all? Returns ``mean_abs_r`` / ``min_abs_r``, the
    per-pair correlations, and the pairing itself.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.ndim != 2 or b.ndim != 2:
        raise ValueError("both embeddings must be 2-D (cells × components)")
    if a.shape[0] != b.shape[0] or a.shape[0] == 0:
        raise ValueError("both embeddings must cover the same non-zero cell set")
    from scipy.optimize import linear_sum_assignment

    k = min(a.shape[1], b.shape[1])
    corr = np.abs(_corr_matrix(a[:, :k], b[:, :k]))
    rows, cols = linear_sum_assignment(-corr)          # maximise total |r|
    matched = corr[rows, cols]
    return {
        "mean_abs_r": float(matched.mean()),
        "min_abs_r": float(matched.min()),
        "matched_r": matched,
        "pairs": list(zip((rows + 1).tolist(), (cols + 1).tolist())),
    }


def knn_overlap(a, b, k: int = KNN_K) -> float:
    """Mean fraction of each cell's k nearest neighbours shared by two embeddings.

    The comparison for embeddings whose coordinates are not comparable. Used
    twice: between R's and shanuz's t-SNE (do they place the same cells
    together?), and between a t-SNE and the PCA space it came from (does the
    embedding preserve its input's neighbourhoods?). 1.0 is identical
    neighbourhoods; for k ≪ n, random embeddings score ≈ k/n.
    """
    from sklearn.neighbors import NearestNeighbors

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.ndim != 2 or b.ndim != 2:
        raise ValueError("both embeddings must be 2-D (cells × dims)")
    n = a.shape[0]
    if n != b.shape[0] or n == 0:
        raise ValueError("both embeddings must cover the same non-zero cell set")
    k = int(min(k, n - 1))
    if k < 1:
        raise ValueError("need at least 2 cells to have a neighbourhood")

    def _neighbours(x):
        nn = NearestNeighbors(n_neighbors=k + 1).fit(x)
        idx = nn.kneighbors(x, return_distance=False)
        return [set(row[row != i][:k]) for i, row in enumerate(idx)]

    na, nb = _neighbours(a), _neighbours(b)
    return float(np.mean([len(x & y) / k for x, y in zip(na, nb)]))


def basis_agreement(py_emb, r_emb) -> dict:
    """How far the two PCA bases agree, and where they stop agreeing *in order*.

    A bare "min correlation over 20 PCs" is misleading here. Leading PCs carry
    real structure and match one-to-one; the noise tail does not have a stable
    ordering, so the two tools routinely emit the same subspace with PC 16-19
    permuted. A diagonal comparison reads that permutation as total
    disagreement, when the components are in fact still matched.

    So this reports both: the diagonal \\|r\\| per PC, *and* which R PC each Python
    PC actually matches best. ``aligned_through`` is the last PC before the first
    reordering — the range over which a downstream per-PC comparison (JackStraw)
    is comparing like with like, and therefore the range over which any
    divergence is attributable to the function under test rather than to PCA.
    """
    py_emb = np.asarray(py_emb, dtype=float)
    r_emb = np.asarray(r_emb, dtype=float)
    if py_emb.shape[0] != r_emb.shape[0] or py_emb.shape[0] == 0:
        raise ValueError("both embeddings must cover the same non-zero cell set")

    n = min(py_emb.shape[1], r_emb.shape[1])
    corr = np.abs(_corr_matrix(py_emb[:, :n], r_emb[:, :n]))
    diag = np.diag(corr).copy()
    best = corr.argmax(axis=0) + 1                      # best Python PC per R PC
    in_order = best == np.arange(1, n + 1)
    broke = np.flatnonzero(~in_order)
    aligned = int(n if broke.size == 0 else broke[0])
    return {
        "pca_basis_abs_r": diag,
        "pca_basis_best_match": best,
        "pca_basis_aligned_through": aligned,
        "pca_basis_min_aligned_r": float(diag[:aligned].min()) if aligned else float("nan"),
    }


def pc_table(py_pvals, py_emp, r_pvals=None, r_emp=None,
             score_thresh: float = SCORE_THRESH) -> pd.DataFrame:
    """Per-PC JackStraw comparison: the aggregate score and what drives it.

    ``n_sig_features`` — how many features fall below ``score_thresh`` on a PC —
    is the statistic R's ``prop.test`` consumes, and it is computable from either
    tool's per-feature p-value matrix. That makes it the apples-to-apples column
    even though the two aggregate *scores* come from different tests: if the
    tools agree here but not on the score, the difference is the aggregation; if
    they disagree here too, it is the null.
    """
    py_pvals = np.asarray(py_pvals, dtype=float)
    py_emp = np.asarray(py_emp, dtype=float)
    if py_emp.shape[1] < py_pvals.size:
        raise ValueError("per-feature p-value matrix has fewer dims than scores")
    if r_pvals is not None:
        r_pvals = np.asarray(r_pvals, dtype=float)
    if r_emp is not None:
        r_emp = np.asarray(r_emp, dtype=float)

    rows = []
    for j in range(py_pvals.size):
        row = {
            "PC": j + 1,
            "py_score": py_pvals[j],
            "py_n_sig_features": int((py_emp[:, j] <= score_thresh).sum()),
            "py_median_p": float(np.median(py_emp[:, j])),
        }
        if r_pvals is not None:
            row["r_score"] = r_pvals[j] if j < r_pvals.size else np.nan
        if r_emp is not None:
            row["r_n_sig_features"] = (int((r_emp[:, j] <= score_thresh).sum())
                                       if j < r_emp.shape[1] else -1)
            row["r_median_p"] = (float(np.median(r_emp[:, j]))
                                 if j < r_emp.shape[1] else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def build_scoreboard(rows: list[dict]) -> pd.DataFrame:
    """Assemble the per-method comparison metrics into a tidy, ordered table."""
    cols = ["method", "metric", "python", "r", "agreement"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df.columns]]


def _r_feature_key(name: str) -> str:
    """Normalise a feature symbol to how R's ``Read10X`` would have stored it.

    R rewrites underscores in feature names to dashes, so pbmc3k's ``Y_RNA``
    arrives as ``Y-RNA`` on that side. Comparing the two p-value matrices *by
    name* rather than by row position means this — and any future reordering or
    dropped feature — can never silently misalign the comparison.
    """
    return str(name).replace("_", "-")


def _read_feature_matrix(path, features):
    """Read an R feature × dim CSV and align its rows to ``features`` by name."""
    if not Path(path).exists():
        return None
    df = pd.read_csv(path).set_index("feature")
    df.index = [_r_feature_key(i) for i in df.index]
    aligned = df.reindex([_r_feature_key(f) for f in features])
    if aligned.isna().all(axis=1).any():
        missing = int(aligned.isna().all(axis=1).sum())
        raise ValueError(
            f"{missing} of the {len(aligned)} Python features are absent from "
            f"{Path(path).name} — the two runs are not on the same feature basis")
    return aligned.to_numpy()


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_pbmc3k_object(data_dir=None):
    """The PBMC 3k object, at the plain ``min_cells`` / ``min_features`` gates.

    Deliberately *without* the guided tutorial's mitochondrial-percentage and
    upper-nFeature filters. Nothing here is about QC, and every extra filtering
    step is one more place the two tools could part company before the functions
    under test are even reached; :func:`prep` writes the surviving cell list out
    so the R script subsets to exactly this set.
    """
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    return create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc3k_dimreduc", feature_names=list(genes), cell_names=list(cells),
    )


def prep(obj, n_hvg=N_HVG, n_pcs=N_PCS):
    """Normalize → HVG → scale → PCA, and write the HVG list for the R reference.

    Scales only the variable features (Seurat's ``ScaleData`` default): JackStraw
    subsets ``scale.data`` to the reduction's features anyway, so scaling the
    other 11k genes would cost time and change nothing.
    """
    normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(obj, selection_method="vst", nfeatures=n_hvg)
    hvg = list(obj.assays["RNA"].variable_features)
    scale_data(obj, features=hvg)
    run_pca(obj, n_pcs=n_pcs, features=hvg, reduction_name="pca")

    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    (FIGURES / "cells.txt").write_text("\n".join(obj.cell_names()) + "\n")
    return hvg


def run_jackstraw(obj, dims=JS_DIMS, num_replicate=JS_REPLICATES,
                  prop_freq=JS_PROP_FREQ, seed=42):
    """The permutation test, then the per-PC aggregation. Returns (js, scores)."""
    js = jack_straw(obj, reduction="pca", dims=dims,
                    num_replicate=num_replicate, prop_freq=prop_freq, seed=seed)
    scores = score_jackstraw(obj, reduction="pca", dims=dims,
                             score_thresh=SCORE_THRESH)
    return js, scores


def run_reductions(obj, n_ics=N_ICS, tsne_dims=TSNE_DIMS, seed=42):
    """ICA on the scaled HVGs, t-SNE on the leading PCs."""
    run_ica(obj, nics=n_ics, seed=seed, reduction_name="ica")
    run_tsne(obj, dims=list(range(tsne_dims)), reduction="pca",
             seed=seed, reduction_name="tsne")


def summarize(obj, scores, js, verbose=True) -> dict:
    """Summarise the Python side: PC cutoff, ICA/t-SNE shapes, structure retained."""
    sig = significant_dims(scores)
    keep = n_leading_significant(scores)
    pca = obj.reductions["pca"].cell_embeddings

    out = {
        "n_cells": int(len(obj.cell_names())),
        "n_features": int(js.empirical_p_values.shape[0]),
        "pc_scores": np.asarray(scores, dtype=float),
        "significant_dims": sig,
        "n_leading_significant": keep,
        "pc_table": pc_table(scores, js.empirical_p_values),
    }
    # ICA / t-SNE are optional (run_full(do_reductions=False) skips them, and the
    # JackStraw half of the tutorial stands on its own).
    if "ica" in obj.reductions:
        out["n_ics"] = int(obj.reductions["ica"].cell_embeddings.shape[1])
    if "tsne" in obj.reductions:
        out["tsne_knn_vs_pca"] = knn_overlap(
            obj.reductions["tsne"].cell_embeddings, pca[:, :TSNE_DIMS], k=KNN_K)

    if not verbose:
        return out

    section("JackStraw — which PCs are significant? (shanuz)")
    board = out["pc_table"][["PC", "py_score", "py_n_sig_features", "py_median_p"]]
    print("    " + board.head(JS_DIMS).round(6).to_string(index=False)
          .replace("\n", "\n    "))
    print(f"\n  PCs clearing alpha={ALPHA}: {sig.tolist()}")
    print(f"  PCs to keep (run before the first drop-off): {keep}")
    if "tsne_knn_vs_pca" in out:
        print(f"\n  t-SNE keeps {out['tsne_knn_vs_pca']:.3f} of each cell's "
              f"{KNN_K} PCA neighbours.")
    return out


def report_concordance(obj, summary, verbose=True) -> dict | None:
    """Compare Python vs R, if ``pbmc3k_dimreduc_verify.R`` has run.

    Step 0 checks the shared PCA basis; only then are the JackStraw, ICA and
    t-SNE numbers meaningful. Returns ``None`` (with a hint) when R has not run.
    """
    pca_path = FIGURES / "r_pca.csv"
    if not pca_path.exists():
        if verbose:
            print(f"\n  [concordance] {pca_path.name} not found — run "
                  "`Rscript tutorials/pbmc3k_dimreduc_verify.R` first.")
        return None

    cells = obj.cell_names()

    def _read(name):
        path = FIGURES / name
        if not path.exists():
            return None
        return pd.read_csv(path).set_index("cell").reindex(cells)

    r_pca = _read("r_pca.csv")
    py_pca = obj.reductions["pca"].cell_embeddings

    # --- step 0: the shared basis -----------------------------------------
    out = basis_agreement(py_pca, r_pca.to_numpy())

    # --- JackStraw ---------------------------------------------------------
    pcs_path = FIGURES / "r_jackstraw_pcs.csv"
    if pcs_path.exists():
        r_pcs = pd.read_csv(pcs_path)
        r_emp = _read_feature_matrix(FIGURES / "r_jackstraw_p.csv",
                                     obj.reductions["pca"].features())
        py_scores = summary["pc_scores"]
        r_scores = r_pcs["R_Score"].to_numpy()
        table = pc_table(py_scores, obj.reductions["pca"].jackstraw.empirical_p_values,
                         r_pvals=r_scores, r_emp=r_emp)
        py_sig, r_sig = significant_dims(py_scores), significant_dims(r_scores)
        out.update({
            "pc_table": table,
            "py_significant": py_sig,
            "r_significant": r_sig,
            "significant_agreement": float(
                len(set(py_sig.tolist()) & set(r_sig.tolist()))
                / max(len(set(py_sig.tolist()) | set(r_sig.tolist())), 1)),
            "py_keep": n_leading_significant(py_scores),
            "r_keep": n_leading_significant(r_scores),
        })
        if r_emp is not None:
            from scipy.stats import spearmanr
            n = min(r_emp.shape[1], obj.reductions["pca"].jackstraw
                    .empirical_p_values.shape[1])
            py_emp = obj.reductions["pca"].jackstraw.empirical_p_values
            out["feature_p_spearman"] = np.array([
                float(spearmanr(py_emp[:, j], r_emp[:, j]).statistic)
                for j in range(n)
            ])

    # --- ICA ---------------------------------------------------------------
    r_ica = _read("r_ica.csv")
    if r_ica is not None and "ica" in obj.reductions:
        out["ica"] = matched_component_correlation(
            obj.reductions["ica"].cell_embeddings, r_ica.to_numpy())

    # --- t-SNE -------------------------------------------------------------
    r_tsne = _read("r_tsne.csv")
    if r_tsne is not None and "tsne" in obj.reductions:
        py_tsne = obj.reductions["tsne"].cell_embeddings
        out["tsne_cross_knn"] = knn_overlap(py_tsne, r_tsne.to_numpy(), k=KNN_K)
        out["tsne_r_knn_vs_pca"] = knn_overlap(
            r_tsne.to_numpy(), r_pca.to_numpy()[:, :TSNE_DIMS], k=KNN_K)

    if verbose:
        _print_concordance(out, summary)
    return out


def _print_concordance(out, summary) -> None:
    """Render the R-vs-Python report — basis check first, then each method."""
    section("Step 0 — do the two PCA bases agree?")
    basis = out["pca_basis_abs_r"]
    aligned = out["pca_basis_aligned_through"]
    print(f"  Per-PC |correlation| over {len(basis)} PCs: "
          f"median {np.median(basis):.4f}")
    print(f"  Matched one-to-one and in the same order through PC {aligned} "
          f"(min |r| there: {out['pca_basis_min_aligned_r']:.4f}).")
    if aligned < len(basis):
        swapped = [f"py{i + 1}~R{out['pca_basis_best_match'][i]}"
                   for i in range(aligned, len(basis))
                   if out["pca_basis_best_match"][i] != i + 1]
        print(f"  Beyond PC {aligned} the noise-tail PCs reorder "
              f"({', '.join(swapped[:4])}, ...):")
        print("  the same subspace, differently numbered. Per-PC rows past that")
        print("  point are not like-for-like — the comparison below is decisive")
        print(f"  through PC {aligned}.")

    if "pc_table" in out:
        section("JackStraw — R vs shanuz")
        cols = ["PC", "py_score", "r_score", "py_n_sig_features",
                "r_n_sig_features", "py_median_p", "r_median_p"]
        table = out["pc_table"]
        show = [c for c in cols if c in table.columns]
        print("    " + table[show].round(6).to_string(index=False)
              .replace("\n", "\n    "))
        print(f"\n  PCs significant at alpha={ALPHA}:")
        print(f"    shanuz : {out['py_significant'].tolist()}")
        print(f"    R      : {out['r_significant'].tolist()}")
        print(f"    Jaccard agreement: {out['significant_agreement']:.4f}")
        print(f"  PCs kept (run before the first drop-off): "
              f"shanuz {out['py_keep']}, R {out['r_keep']}")
        if "feature_p_spearman" in out:
            fs = out["feature_p_spearman"]
            print(f"\n  Per-feature empirical-p Spearman, per PC: "
                  f"min {fs.min():.4f}, median {np.median(fs):.4f}")
            print("  py_n_sig_features vs r_n_sig_features is the like-for-like")
            print("  column — it is computed the same way from either tool's")
            print("  p-value matrix, so it separates a null difference (these")
            print("  differ) from an aggregation difference (only the scores do).")

    rows = []
    if "ica" in out:
        rows.append({"method": "ICA", "metric": f"matched |r| over {N_ICS} comps",
                     "agreement": round(out["ica"]["mean_abs_r"], 4)})
    if "tsne_cross_knn" in out:
        rows.append({"method": "t-SNE", "metric": f"shared {KNN_K}-NN, R vs py",
                     "agreement": round(out["tsne_cross_knn"], 4)})
        if "tsne_knn_vs_pca" in summary:
            rows.append({"method": "t-SNE", "metric": f"{KNN_K}-NN kept from PCA",
                         "python": round(summary["tsne_knn_vs_pca"], 4),
                         "r": round(out["tsne_r_knn_vs_pca"], 4)})
    if rows:
        section("ICA & t-SNE — structure, not coordinates")
        # fillna: not every row has both a python and an r column (the shared
        # metrics are a single agreement number), and "NaN" reads as a failure.
        print("    " + build_scoreboard(rows).fillna("").to_string(index=False)
              .replace("\n", "\n    "))
        if "ica" in out:
            print(f"\n  ICA components are matched one-to-one by |Pearson r| "
                  f"(sign- and\n  order-free); worst matched pair "
                  f"{out['ica']['min_abs_r']:.4f}.")
        print("  t-SNE coordinates are not comparable across implementations —")
        print("  'kept from PCA' is each tool judged against its own input, so")
        print("  the two columns are directly comparable to each other.")


def run_full(data_dir=None, verbose=True, seed=42,
             num_replicate=JS_REPLICATES, do_reductions=True):
    """Load, preprocess, JackStraw, ICA/t-SNE, summarize, compare to R."""
    t0 = time.time()
    if verbose:
        section("Loading PBMC 3k (dimensional-reduction extras)")
    obj = load_pbmc3k_object(data_dir=data_dir)

    hvg = prep(obj)
    if verbose:
        print(f"  {len(obj.assays['RNA'].features())} genes x "
              f"{len(obj.cell_names())} cells; {len(hvg)} HVGs -> {N_PCS} PCs")
        print("  HVG list written to figures_dimreduc/ for the R reference")

    t1 = time.time()
    js, scores = run_jackstraw(obj, num_replicate=num_replicate, seed=seed)
    if verbose:
        print(f"  JackStraw: {num_replicate} replicates x {JS_DIMS} dims "
              f"({time.time() - t1:.1f}s)")

    if do_reductions:
        t1 = time.time()
        run_reductions(obj, seed=seed)
        if verbose:
            print(f"  ICA ({N_ICS} comps) + t-SNE (PC 1-{TSNE_DIMS}) "
                  f"({time.time() - t1:.1f}s)")

    summary = summarize(obj, scores, js, verbose=verbose)
    report_concordance(obj, summary, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return obj, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PBMC 3k dimensional-reduction extras tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--replicates", type=int, default=JS_REPLICATES,
                        help="JackStraw permutation replicates (Seurat default 100)")
    args = parser.parse_args()
    run_full(data_dir=args.data_dir, num_replicate=args.replicates)
