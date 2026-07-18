"""Leverage-score sketching — Shanuz vs R Seurat on ifnb.

The Seurat v5 answer to "this atlas has two million cells and I have a laptop":
don't subsample uniformly, subsample by **statistical leverage**, analyse the
small subset, then push the answers back onto every cell.

  * :func:`shanuz.leverage_score` ↔ ``LeverageScore`` — per-cell influence on the
    column space of the data. High for a cell in a sparse, distinctive corner;
    low for one of ten thousand near-identical cells in a common state.
  * :func:`shanuz.sketch_data` ↔ ``SketchData`` — draws the leverage-weighted
    subset (or a uniform one, as the control).
  * :func:`shanuz.project_data` ↔ ``ProjectData`` — the inverse map: every
    full-dataset cell through the sketch's PCA, its fitted UMAP, and its labels.

Why this tutorial exists — and what it found
--------------------------------------------
All three landed with synthetic fixtures only and had never met a Seurat
reference. Building this comparison turned up **two real defects**, both fixed in
the same pull request as this tutorial.

**Defect 2, in one paragraph, because it is the more instructive one.**
``project_data`` transferred the sketch's labels through the
:mod:`shanuz.transfer` integration anchors. Seurat's ``ProjectData`` does not:
it calls ``TransferSketchLabels``, a weighted k-nearest-neighbour vote *inside
the projected reduction*, with the sketch's own rows as the reference. The
anchor route scored **better** on ifnb — 0.936 against Seurat's 0.905 — which is
exactly why it survived review. It is still wrong, and not only on fidelity:
finding anchors between the sketch and the full dataset costs precisely what
sketching exists to remove, so on the million-cell objects this API is written
for the anchor route is unusable rather than merely different.

Matching Seurat moved accuracy **down**, to Seurat's. Measured two ways, because
they answer different questions:

* handed **R's own sketch cells**, so the transfer is the only variable, the two
  agree on **98.1 %** of cells (accuracy 0.9031 against R's 0.9050);
* end to end, each tool drawing its own sketch, both land on **exactly 0.9050**
  at 94.9 % per-cell agreement — the extra disagreement is the differing draw,
  not the transfer.

A regression here would look like an improvement, so the test that guards it
checks the mechanism rather than the score.

Defect 1 is ``leverage_score`` itself.

Seurat computes leverage from a **rank-50 truncated SVD** — ``ℓ = rowSums(V²)``
over the leading 50 right singular vectors — so the scores sum to 50. Shanuz
whitened against the **full rank** instead. Both are defensible as "leverage" in
the textbook sense, and that is exactly what made it hard to see: the full-rank
version is the classical hat-matrix definition, and shanuz's own test asserted it
to six decimals.

But full-rank leverage is useless on data shaped like this. Scores sum to the
rank and are capped at 1, so with 2000 variable features over a few thousand
cells every score is crushed towards ``d/n``:

===========================  ==========  ==========
PBMC 3k, 2000 HVGs           shanuz      R Seurat
===========================  ==========  ==========
sum of scores                2000.0      50.0
coefficient of variation     0.085       0.384
max / median                 1.34        6.48
Spearman vs R                0.238       --
===========================  ==========  ==========

A max/median of 1.34 means the single most influential cell in the dataset was
only 34 % more likely to be drawn than a median one — against 6.48× in Seurat,
and 1.00 for uniform sampling. **Sketching was therefore doing almost nothing**:
``sketch_data`` was an expensive way to sample uniformly, which is precisely the
thing leverage sampling exists to avoid. Nothing about the output looked wrong.

After the fix the exact regime agrees with Seurat to **2.4e-6** — below R's own
run-to-run noise, since ``irlba`` starts from a random vector and R never seeds
it.

Why ifnb, and why two regimes
-----------------------------
Seurat picks between two algorithms on cell count, and the tutorial exercises
both on the same data by moving the threshold rather than the dataset:

  * ``nsketch = 10000`` → **exact**: 13,999 cells < 1.5 × 10000, so R takes the
    truncated SVD. Deterministic up to ``irlba``, and shanuz matches it exactly.
  * ``nsketch = 5000``  → **sketched**: ``CountSketch`` → ``QR`` →
    Johnson–Lindenstrauss. Both tools draw their random matrices from their own
    RNG, so this one is compared statistically, never bit-for-bit.

ifnb also supplies the thing the method is *for*: 13 annotated cell types
spanning 4,362 cells (CD14 Mono) down to 55 (Eryth). That makes the headline
metric an honest one — not "do the numbers match" but **does the method keep the
rare cells**, measured against held-out annotations both tools can see.

The residuals, and why each is RNG rather than a gap
----------------------------------------------------
Three differences survive, and all three were checked distribution-against-
distribution over matched seeds rather than argued from a single pair of runs:

* **The sketched regime** differs from R at Spearman 0.946. Over seeds
  123/1/7/42/2024 shanuz's score sums span 69,611-70,968 and R's 70,278-70,867;
  the CVs span 0.5600-0.5770 and 0.5657-0.5776. The distributions overlap, so
  the difference is the two RNGs, not the algorithm.
* **Sketch composition.** Over the same seeds the rarest types land at
  Eryth 1.44x (1.15-1.91) / pDC 2.26x / Mk 1.58x in shanuz against
  1.29x (0.64-2.04) / 2.23x / 1.47x in R — overlapping, and R's Eryth spread is
  the wider of the two. A single-seed comparison made this look like a real
  divergence (0.89x against 1.65x); it was two draws from overlapping
  distributions landing at opposite ends. With ~11 erythrocytes in a 2,000-cell
  sketch, that spread is Poisson noise on a small count.
* **The two regimes against each other.** Seurat's sketched approximation tracks
  its own exact scores at Spearman only 0.31 — the JL projection is a large
  approximation, and Seurat leaves it unscaled, so the sketched scores are not
  even on the same scale (they sum to ~70,000, not 50). shanuz reproduces that
  gap at 0.31 too (R 0.309, shanuz 0.307), which is the useful evidence that the
  sketched path is faithful. **Compare scores within one regime, never across
  the two**, and read the sketched one as "cheap and roughly right".

Lazy matrices have no R side here
---------------------------------
``LazyMatrix`` / ``open_lazy_matrix`` are shanuz's on-disk backing, R's
equivalent being BPCells, which is not installed. So that section is **not** a
side-by-side: it checks the one property that matters, that going through disk
does not change the answer, and is reported separately so it is not mistaken for
an R comparison.

Note on the shared basis
------------------------
Leverage is a property of a *matrix*, so the two tools must score the same one.
The Python run writes the exact cell barcodes and HVG list it used to
``figures_sketch/cells.txt`` and ``hvg_features.txt``; the R script subsets to
those cells and passes that feature list to ``LeverageScore`` explicitly.

**That last word is load-bearing.** ``VariableFeatures(obj) <- hvg`` does *not*
register for ``layer = "data"``, so ``LeverageScore`` silently finds no variable
features, falls back, and scores **all 13,714 genes** — a different matrix, and
no error. Pass ``features =`` explicitly. This cost an hour of chasing a
divergence that was never there.

Usage
-----
    Rscript tutorials/export_seuratdata.R ifnb      # one-time counts export
    python  tutorials/ifnb_sketch_tutorial.py       # writes the shared lists

Then, for the side-by-side numbers and figures:

    Rscript tutorials/ifnb_sketch_verify.R
    python  tutorials/generate_sketch_plots.py

References
----------
Hao Y, Stuart T, Kowalski MH, et al. (2024) **Dictionary learning for integrative,
multimodal and scalable single-cell analysis.** Nature Biotechnology 42, 293-304.
https://doi.org/10.1038/s41587-023-01767-y

Drineas P, Magdon-Ismail M, Mahoney MW, Woodruff DP (2012) **Fast approximation
of matrix coherence and statistical leverage.** JMLR 13, 3475-3506.

Kang HM, Subramaniam M, Targ S, et al. (2018) **Multiplexed droplet single-cell
RNA-sequencing using natural genetic variation.** Nature Biotechnology 36, 89-94.
https://doi.org/10.1038/nbt.4042
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

from shanuz.datasets import ifnb
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import normalize_data, find_variable_features, scale_data
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.sketch import leverage_score, sketch_data, project_data

FIGURES = Path(__file__).parent / "figures_sketch"

CELLTYPE = "seurat_annotations"
BATCH = "stim"

N_HVG = 2000

# The two regimes, chosen so one dataset exercises both of Seurat's branches.
# R switches at ncells < nsketch * 1.5, and ifnb has 13,999 cells.
NSKETCH_EXACT = 10_000      # 13999 < 15000  → truncated SVD
NSKETCH_SKETCHED = 5_000    # 13999 > 7500   → CountSketch + QR + JL

SKETCH_CELLS = 2_000        # cells kept in the sketch (~14 % of the data)
N_PCS = 30
CLUSTER_RESOLUTION = 0.6
SEED = 123

# Ranks compared when asking whether the two tools agree about *which* cells
# matter, rather than about the score values.
TOP_K = (200, 1000)


# ----------------------------------------------------------------------
# Pure helpers — no data, no I/O, unit-testable
# ----------------------------------------------------------------------


def group_means(values, groups):
    """Mean of ``values`` within each group, as a Series indexed by group."""
    return pd.Series(np.asarray(values, dtype=float)).groupby(np.asarray(groups)).mean()


def group_enrichment(scores, groups):
    """Per-group mean leverage relative to the overall mean.

    The readable form of "does this group punch above its weight": 1.0 is
    exactly average, 2.9 means cells of that type carry nearly three times the
    average influence on the data's column space.
    """
    scores = np.asarray(scores, dtype=float)
    groups = np.asarray(groups)
    overall = scores.mean()
    sizes = pd.Series(groups).value_counts()
    means = group_means(scores, groups)
    out = pd.DataFrame({
        "n": sizes.reindex(means.index).to_numpy(),
        "mean_leverage": means.to_numpy(),
        "vs_overall": (means / overall).to_numpy() if overall else np.zeros(means.size),
    }, index=means.index)
    return out.sort_values("n", ascending=False)


def rarity_correlation(scores, groups):
    """Spearman between a group's mean leverage and how big the group is.

    The single number that says whether leverage is doing its job. Strongly
    **negative** is the intended behaviour: the rarer the population, the higher
    its cells score. Returns ``nan`` when there are fewer than three groups,
    which is not enough to rank.
    """
    from scipy.stats import spearmanr

    table = group_enrichment(scores, groups)
    if len(table) < 3:
        return float("nan")
    return float(spearmanr(table["mean_leverage"].to_numpy(),
                           table["n"].to_numpy()).statistic)


def composition(groups, index=None):
    """Fraction of cells in each group, optionally reindexed to a fixed order."""
    frac = pd.Series(np.asarray(groups)).value_counts(normalize=True)
    return frac if index is None else frac.reindex(index).fillna(0.0)


def sketch_fold_change(full_groups, sketch_groups):
    """How each group's share changes between the full data and the sketch.

    ``fold > 1`` means the sketch over-represents that group. This is the metric
    the method is actually sold on, and it is computed identically from either
    tool's cell list, so it compares like for like even though the two draw
    different cells.
    """
    full = composition(full_groups)
    sk = composition(sketch_groups, index=full.index)
    out = pd.DataFrame({
        "n_full": pd.Series(np.asarray(full_groups)).value_counts().reindex(full.index),
        "frac_full": full.to_numpy(),
        "frac_sketch": sk.to_numpy(),
    }, index=full.index)
    out["fold"] = np.where(out["frac_full"] > 0,
                           out["frac_sketch"] / out["frac_full"], np.nan)
    return out.sort_values("n_full", ascending=False)


def agreement(a, b, top_k=TOP_K):
    """Rank/linear agreement between two score vectors, plus top-k set overlap.

    Correlation alone is a weak summary for a sampling weight: what matters is
    whether the two tools would *draw the same cells*, so the overlap of the
    highest-scoring k is reported alongside.
    """
    from scipy.stats import spearmanr

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"score vectors differ in shape: {a.shape} vs {b.shape}")

    out = {
        "spearman": float(spearmanr(a, b).statistic),
        "pearson": float(np.corrcoef(a, b)[0, 1]),
        "max_abs_diff": float(np.abs(a - b).max()),
    }
    for k in top_k:
        k = min(k, a.size)
        top_a = set(np.argsort(a)[-k:])
        top_b = set(np.argsort(b)[-k:])
        out[f"top{k}_overlap"] = len(top_a & top_b) / k
    return out


def label_accuracy(predicted, truth):
    """Fraction of cells whose projected label matches the held-out annotation."""
    predicted = np.asarray(predicted).astype(str)
    truth = np.asarray(truth).astype(str)
    if predicted.shape != truth.shape:
        raise ValueError("predicted and truth differ in length")
    keep = truth != "nan"
    return float((predicted[keep] == truth[keep]).mean()) if keep.any() else float("nan")


def _r_feature_key(name):
    """R's ``Read10X`` rewrites underscores in feature names to dashes."""
    return name.replace("_", "-")


def _read_r_scores(path, column, cells):
    """Read an R per-cell CSV back, aligned to the Python cell order **by name**.

    Aligning by name rather than position: the R object's column order comes out
    of ``subset``, and trusting it to match would make every downstream number
    quietly wrong rather than loudly absent.
    """
    frame = pd.read_csv(path)
    series = frame.set_index(frame.columns[0])[column]
    missing = [c for c in cells if c not in series.index]
    if missing:
        raise KeyError(
            f"{len(missing)} of the Python cells are absent from {path.name} "
            f"(e.g. {missing[:3]}) — the two runs used different cell sets."
        )
    return series.reindex(cells).to_numpy(dtype=float)


# ----------------------------------------------------------------------
# Pipeline
# ----------------------------------------------------------------------


def load_ifnb_sketch_object(data_dir=None, min_cells=3):
    """The ifnb object carrying the cell-type annotations the metrics need."""
    counts, genes, cells, meta = ifnb(data_dir=data_dir)
    obj = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=min_cells,
        project="ifnb_sketch", feature_names=list(genes), cell_names=list(cells),
        meta_data=meta,
    )
    if CELLTYPE not in obj.meta_data.columns:
        raise KeyError(
            f"ifnb metadata is missing {CELLTYPE!r}; re-export with "
            f"`Rscript tutorials/export_seuratdata.R ifnb`."
        )
    return obj


def prep(obj, n_hvg=N_HVG):
    """Normalize → HVG, and write the shared cell + feature lists for R.

    Deliberately no ``scale_data`` before scoring: Seurat's ``LeverageScore``
    reads the **data** layer, not ``scale.data``, and shanuz now defaults the
    same way. Scaling happens later, for the sketch's PCA.
    """
    normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(obj, selection_method="vst", nfeatures=n_hvg)
    hvg = list(obj.assays[obj.active_assay].variable_features)

    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    (FIGURES / "cells.txt").write_text("\n".join(obj.cell_names()) + "\n")
    return hvg


def run_leverage(obj, hvg, seed=SEED):
    """Leverage in both of Seurat's regimes. Returns ``{regime: scores}``."""
    out = {}
    for name, nsketch in (("exact", NSKETCH_EXACT), ("sketched", NSKETCH_SKETCHED)):
        t0 = time.time()
        out[name] = leverage_score(
            obj, nsketch=nsketch, features=hvg, seed=seed,
            var_name=f"leverage.{name}",
        )
        out[f"{name}_seconds"] = time.time() - t0
    return out


def run_sketch(obj, hvg, ncells=SKETCH_CELLS, seed=SEED):
    """Draw the leverage-weighted sketch and a uniform one as the control.

    The uniform draw is not decoration: "the sketch contains rare cells" means
    nothing without it, because a large enough sketch contains rare cells by
    accident. The comparison is leverage *against* uniform, at the same size.
    """
    sketches = {}
    for method in ("LeverageScore", "Uniform"):
        sketches[method] = sketch_data(
            obj, ncells=ncells, method=method, features=hvg,
            nsketch=NSKETCH_SKETCHED, var_name=None, seed=seed,
        )
    return sketches


def analyse_sketch(sketch, hvg, n_pcs=N_PCS, resolution=CLUSTER_RESOLUTION, seed=SEED):
    """The expensive analysis, run on the sketch only: scale → PCA → clusters → UMAP."""
    present = [f for f in hvg if f in set(sketch.assays[sketch.active_assay].features())]
    scale_data(sketch, features=present)
    run_pca(sketch, n_pcs=n_pcs, features=present, reduction_name="pca", seed=seed)
    find_neighbors(sketch, reduction="pca", dims=list(range(n_pcs)))
    find_clusters(sketch, resolution=resolution, random_seed=seed)
    run_umap(sketch, reduction="pca", dims=list(range(n_pcs)), seed=seed)
    return sketch


def run_projection(obj, sketch, hvg):
    """Push the sketch's analysis back onto every cell (Seurat's ``ProjectData``).

    Labels transferred are the sketch's *own* held-out annotations, so the
    accuracy reported afterwards is against ground truth the projection never
    saw for the 86 % of cells outside the sketch.

    No ``seed``: the label vote is a deterministic weighted k-NN in the projected
    space, so unlike the sketch draw there is nothing here to seed.
    """
    present = [f for f in hvg if f in set(sketch.assays[sketch.active_assay].features())]
    scale_data(obj, features=present)
    project_data(
        obj, sketch, reduction="pca", full_reduction="pca.full",
        umap_reduction="umap", full_umap_reduction="ref.umap",
        refdata={"projected.celltype": CELLTYPE},
    )
    return obj


def lazy_roundtrip(obj, hvg, tmpdir=None, seed=SEED):
    """Does backing the matrix with an on-disk ``LazyMatrix`` change the answer?

    No R counterpart — BPCells is not installed — so this is a shanuz-internal
    invariant, reported separately from every R comparison. The only acceptable
    answer is "not at all".
    """
    from shanuz.lazy import write_lazy_matrix, open_lazy_matrix, is_lazy

    assay = obj.assays[obj.active_assay]
    target = Path(tmpdir or FIGURES) / "ifnb_data.lazy"
    try:
        write_lazy_matrix(assay.layer_data("data"), target, overwrite=True)
        lazy = open_lazy_matrix(target)
    except Exception as exc:                     # pragma: no cover - environment
        return {"available": False, "reason": f"{type(exc).__name__}: {exc}"}

    dense_scores = leverage_score(obj, nsketch=NSKETCH_EXACT, features=hvg,
                                  seed=seed, var_name=None)
    swapped = assay.layers["data"]
    try:
        assay.layers["data"] = lazy
        lazy_scores = leverage_score(obj, nsketch=NSKETCH_EXACT, features=hvg,
                                     seed=seed, var_name=None)
    except Exception as exc:                     # pragma: no cover - environment
        return {"available": False, "reason": f"{type(exc).__name__}: {exc}"}
    finally:
        assay.layers["data"] = swapped

    return {
        "available": True,
        "is_lazy": bool(is_lazy(lazy)),
        "max_abs_diff": float(np.abs(dense_scores - lazy_scores).max()),
        "identical": bool(np.allclose(dense_scores, lazy_scores, rtol=0, atol=1e-10)),
    }


def summarize(obj, scores, sketches, hvg):
    """Every Python-side number the vignette and the figures quote."""
    groups = obj.meta_data[CELLTYPE].astype(str).to_numpy()
    cells = obj.cell_names()

    summary = {
        "n_cells": len(obj),
        "n_hvg": len(hvg),
        "cells": cells,
        "celltype": groups,
        "leverage": {k: v for k, v in scores.items() if not k.endswith("_seconds")},
        "seconds": {k.replace("_seconds", ""): v
                    for k, v in scores.items() if k.endswith("_seconds")},
    }

    summary["enrichment"] = group_enrichment(scores["exact"], groups)
    summary["rarity_spearman"] = {
        regime: rarity_correlation(scores[regime], groups)
        for regime in ("exact", "sketched")
    }
    # The two regimes against each other, within shanuz. R's own gap here is the
    # yardstick this is read against — see the module docstring.
    summary["exact_vs_sketched"] = agreement(scores["exact"], scores["sketched"])

    summary["sketch_composition"] = {
        method: sketch_fold_change(groups,
                                   obj.meta_data.loc[sk.cell_names(), CELLTYPE]
                                   .astype(str).to_numpy())
        for method, sk in sketches.items()
    }
    if "projected.celltype" in obj.meta_data.columns:
        predicted = obj.meta_data["projected.celltype"].to_numpy()
        summary["projected"] = predicted
        summary["projection_accuracy"] = label_accuracy(predicted, groups)
        in_sketch = set(sketches["LeverageScore"].cell_names())
        held_out = np.array([c not in in_sketch for c in cells])
        summary["projection_accuracy_heldout"] = label_accuracy(
            predicted[held_out], groups[held_out]
        )
    return summary


def report_concordance(summary, figures=FIGURES):
    """Compare against the R references, if ``ifnb_sketch_verify.R`` has run."""
    lev_path = figures / "r_leverage.csv"
    if not lev_path.exists():
        print("\n  R reference not found — run `Rscript tutorials/ifnb_sketch_verify.R`")
        print("  for the side-by-side numbers. Python-side results above stand alone.")
        return None

    cells = summary["cells"]
    out = {"leverage": {}}
    for regime in ("exact", "sketched"):
        r_scores = _read_r_scores(lev_path, f"r_leverage_{regime}", cells)
        out["leverage"][regime] = agreement(summary["leverage"][regime], r_scores)
        out["leverage"][regime]["r_sum"] = float(r_scores.sum())
        out["leverage"][regime]["py_sum"] = float(summary["leverage"][regime].sum())
        out["leverage"][regime]["r_rarity_spearman"] = rarity_correlation(
            r_scores, summary["celltype"]
        )

    comp_path = figures / "r_sketch_composition.csv"
    if comp_path.exists():
        out["sketch_composition"] = pd.read_csv(comp_path)
    proj_path = figures / "r_projection.csv"
    if proj_path.exists() and "projected" in summary:
        frame = pd.read_csv(proj_path).set_index("cell")
        shared = [c for c in cells if c in frame.index]
        position = {c: i for i, c in enumerate(cells)}
        py = np.array([str(summary["projected"][position[c]]) for c in shared])
        r_pred = frame.loc[shared, "r_predicted"].astype(str).to_numpy()
        r_truth = frame.loc[shared, "truth"].astype(str).to_numpy()
        # Per-cell agreement, not two accuracy figures: two methods can be equally
        # accurate while labelling different cells, and it is the calls that are
        # the port's output.
        out["projection"] = {
            "n_compared": len(shared),
            "per_cell_concordance": float((py == r_pred).mean()),
            "py_accuracy": float((py == r_truth).mean()),
            "r_accuracy": float((r_pred == r_truth).mean()),
        }
    return out


def _print_report(summary, concordance):
    def section(title):
        print(f"\n{title}\n{'-' * len(title)}")

    section("Leverage scores — the two regimes")
    for regime in ("exact", "sketched"):
        s = summary["leverage"][regime]
        print(f"  {regime:9s} sum={s.sum():10.1f}  CV={s.std() / s.mean():.4f}  "
              f"max/median={s.max() / np.median(s):5.2f}  "
              f"({summary['seconds'][regime]:.1f}s)")
    print(f"\n  the two regimes agree with each other at Spearman "
          f"{summary['exact_vs_sketched']['spearman']:.3f} — low by design, and R's")
    print("  own gap is the same size (see the vignette).")

    section("Does leverage track rarity? (the point of the method)")
    print(f"  Spearman(mean leverage, cell-type size), exact regime: "
          f"{summary['rarity_spearman']['exact']:+.3f}   (want strongly negative)")
    table = summary["enrichment"]
    print(f"\n  {'cell type':22s} {'n':>6s} {'mean lev / overall':>20s}")
    for name, row in table.iterrows():
        print(f"  {str(name):22s} {int(row['n']):6d} {row['vs_overall']:19.2f}x")

    section("What the sketch keeps")
    lev = summary["sketch_composition"]["LeverageScore"]
    uni = summary["sketch_composition"]["Uniform"]
    print(f"  {'cell type':22s} {'n':>6s} {'full':>8s} {'leverage':>10s} {'uniform':>9s}")
    for name in lev.index:
        print(f"  {str(name):22s} {int(lev.loc[name, 'n_full']):6d} "
              f"{lev.loc[name, 'frac_full']:8.4f} {lev.loc[name, 'fold']:9.2f}x "
              f"{uni.loc[name, 'fold']:8.2f}x")

    if "projection_accuracy" in summary:
        section("project_data — the sketch's labels on every cell")
        print(f"  accuracy, all cells        {summary['projection_accuracy']:.4f}")
        print(f"  accuracy, cells NOT in the sketch  "
              f"{summary['projection_accuracy_heldout']:.4f}")

    if concordance is None:
        return
    section("Shanuz vs R Seurat")
    for regime in ("exact", "sketched"):
        c = concordance["leverage"][regime]
        print(f"  {regime:9s} spearman={c['spearman']:.6f}  pearson={c['pearson']:.6f}  "
              f"sum py={c['py_sum']:.1f} R={c['r_sum']:.1f}")
        overlaps = "  ".join(f"top{k}={c[f'top{k}_overlap']:.1%}" for k in TOP_K)
        print(f"            {overlaps}   max|diff|={c['max_abs_diff']:.3e}")
        print(f"            rarity Spearman: py={summary['rarity_spearman'][regime]:+.3f}  "
              f"R={c['r_rarity_spearman']:+.3f}")
    if "projection" in concordance:
        p = concordance["projection"]
        print(f"\n  project_data, {p['n_compared']} cells:")
        print(f"    per-cell concordance  {p['per_cell_concordance']:.4f}")
        print(f"    accuracy vs truth     py {p['py_accuracy']:.4f}   R {p['r_accuracy']:.4f}")


def run_full(data_dir=None, verbose=True, do_projection=True, do_lazy=True):
    """The whole tutorial. Returns ``(obj, summary)``."""
    obj = load_ifnb_sketch_object(data_dir=data_dir)
    hvg = prep(obj)
    if verbose:
        print(f"ifnb: {len(obj)} cells x {len(hvg)} variable features")

    scores = run_leverage(obj, hvg)
    sketches = run_sketch(obj, hvg)
    analyse_sketch(sketches["LeverageScore"], hvg)
    if do_projection:
        run_projection(obj, sketches["LeverageScore"], hvg)

    summary = summarize(obj, scores, sketches, hvg)
    summary["sketches"] = sketches
    if do_lazy:
        summary["lazy"] = lazy_roundtrip(obj, hvg)

    if verbose:
        concordance = report_concordance(summary)
        _print_report(summary, concordance)
        summary["concordance"] = concordance
        lazy = summary.get("lazy", {})
        if lazy.get("available"):
            print(f"\n  lazy matrix round-trip (no R counterpart): "
                  f"identical={lazy['identical']}, max|diff|={lazy['max_abs_diff']:.2e}")
        elif lazy:
            print(f"\n  lazy matrix round-trip skipped: {lazy['reason']}")
    return obj, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ifnb leverage-score sketching tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--no-projection", action="store_true",
                        help="skip project_data (the slow anchor step)")
    args = parser.parse_args()
    run_full(data_dir=args.data_dir, do_projection=not args.no_projection)
