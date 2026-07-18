"""Batch integration tutorial — Harmony / CCA / RPCA with Shanuz on ifnb.

A Python port of Seurat's integration vignettes
(https://satijalab.org/seurat/articles/integration_introduction) on the
Kang et al. 2018 PBMC dataset (``ifnb``): ~14,000 human peripheral blood
mononuclear cells, half left resting (``CTRL``) and half stimulated for six
hours with interferon-beta (``STIM``). Interferon drives a strong, near-global
transcriptional response, so without correction the cells split first by
*condition* and only then by *cell type* — the textbook batch effect.

The task integration solves: **make the same cell type from the two conditions
overlap, without erasing the biology that separates cell types.** A good
integration collapses the CTRL/STIM shift (so a CD14 monocyte looks like a
CD14 monocyte in either condition) while keeping monocytes, T cells, B cells,
NK cells and the rest apart.

It demonstrates Shanuz's three integration paths side by side with Seurat's:
  * :func:`shanuz.run_harmony` — Seurat's ``RunHarmony(group.by.vars=...)``:
    iteratively nudges the PCA embedding so batches mix while clusters hold,
  * :func:`shanuz.integrate_layers` ``method="cca"`` — Seurat's
    ``FindIntegrationAnchors(reduction="cca")`` + ``IntegrateData``: anchors are
    mutual nearest neighbours in a shared canonical-correlation space,
  * :func:`shanuz.integrate_layers` ``method="rpca"`` — the reciprocal-PCA
    variant: each dataset is projected into the other's PCA space before the
    mutual-nearest-neighbour search — faster and more conservative than CCA.

Why this tutorial exists
------------------------
Every integration function landed in v0.2.0 and had only ever been checked
against synthetic ``default_rng`` fixtures with *balanced* batches. This is the
first time ``run_harmony`` / ``integrate_layers`` meet a real dataset with a
Seurat reference and unequal batch sizes (CTRL 6548 vs STIM 7451). Integration
embeddings are *not* expected to be byte-identical — harmonypy and R's harmony
are separate implementations, and CCA/RPCA differ in SVD sign and neighbour
ties — so the comparison target is **not** matching coordinates but *do the two
tools recover the same structure*: the same drop in batch separation, the same
recovery of cell type, and cluster partitions that agree (adjusted Rand index),
all reported by :func:`report_concordance` once ``ifnb_integration_verify.R``
has run.

Note on the data and the R comparison
-------------------------------------
``ifnb`` is a curated SeuratData object with no clean raw source, so both tools
read the *same* counts exported once from R by ``tutorials/export_seuratdata.R``
(a 10x-style matrix folder), guaranteeing byte-identical input and cell order.
To keep the reduction on one shared gene basis, the Python run writes the
variable features it selected to ``figures_integration/hvg_features.txt`` and the
R script reads them back — so the only divergences left are the genuinely
method-level ones (PCA numerics, the Harmony/anchor algorithms, Louvain ties).

Usage
-----
    Rscript tutorials/export_seuratdata.R ifnb        # one-time, writes the counts
    python  tutorials/ifnb_integration_tutorial.py    # writes the shared HVGs

Then, for the side-by-side numbers and figures:

    Rscript tutorials/ifnb_integration_verify.R
    python  tutorials/generate_integration_plots.py

References
----------
Kang HM, Subramaniam M, Targ S, et al. (2018)
**Multiplexed droplet single-cell RNA-sequencing using natural genetic
variation.** Nature Biotechnology 36, 89-94. https://doi.org/10.1038/nbt.4042

Korsunsky I, Millard N, Fan J, et al. (2019)
**Fast, sensitive and accurate integration of single-cell data with Harmony.**
Nature Methods 16, 1289-1296. https://doi.org/10.1038/s41592-019-0619-0
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
from shanuz.integration import run_harmony, integrate_layers
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap

FIGURES = Path(__file__).parent / "figures_integration"

# The two labels the whole tutorial turns on: the batch to remove and the
# cell-type annotation to preserve (both ship with the dataset).
BATCH = "stim"
CELLTYPE = "seurat_annotations"

# Reduction params mirror the Seurat vignette: 2000 variable features, 30 PCs,
# batch integration on those 30 dims, Louvain at resolution 0.5 on all methods
# so the R-vs-Python cluster comparison is like-for-like.
N_HVG = 2000
N_PCS = 30
RESOLUTION = 0.5

# The three integration reductions, each keyed for its clusters/umap/metrics.
METHODS = ("harmony", "cca", "rpca")


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_integration_tutorial.py)
# ---------------------------------------------------------------------------

def mixing_metrics(embedding, batch, celltype) -> dict:
    """Silhouette summaries of an embedding: batch separation and cell-type separation.

    ``sil_batch`` is the average silhouette width using the batch labels — how
    cleanly the two conditions sit apart, so **lower is better** (a good
    integration mixes them). ``sil_celltype`` uses the cell-type labels — how
    cleanly cell types sit apart, so **higher is better** (integration must not
    blur them). Both are rotation-invariant, so they compare fairly across two
    tools whose embeddings are not coordinate-aligned.
    """
    from sklearn.metrics import silhouette_score

    emb = np.asarray(embedding)
    return {
        "sil_batch": float(silhouette_score(emb, np.asarray(batch))),
        "sil_celltype": float(silhouette_score(emb, np.asarray(celltype))),
    }


def cluster_ari(a, b) -> float:
    """Adjusted Rand index between two per-cell label vectors (order-invariant)."""
    from sklearn.metrics import adjusted_rand_score

    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("label vectors must be the same non-zero length")
    return float(adjusted_rand_score(a, b))


def batch_entropy(clusters, batch) -> float:
    """Size-weighted mean per-cluster batch entropy, normalised to [0, 1].

    For each cluster, the Shannon entropy of its batch composition divided by
    ``log(n_batches)``: **1.0** when every cluster holds the batches in their
    global proportions (perfectly mixed), **0.0** when each cluster is a single
    batch (no integration). A partition-only metric — needs just the cluster and
    batch labels, no embedding — so R and Python are compared on identical terms.
    """
    clusters = np.asarray(clusters)
    batch = np.asarray(batch)
    levels = np.unique(batch)
    if levels.size < 2:
        return float("nan")
    norm = np.log(levels.size)
    total = 0.0
    for c in np.unique(clusters):
        mask = clusters == c
        n = int(mask.sum())
        p = np.array([(batch[mask] == lv).mean() for lv in levels])
        p = p[p > 0]
        ent = -(p * np.log(p)).sum() / norm
        total += n * ent
    return float(total / len(clusters))


def build_scoreboard(rows: list[dict]) -> pd.DataFrame:
    """Assemble per-method integration metrics into a tidy, ordered table."""
    cols = ["method", "sil_batch", "sil_celltype", "n_clusters",
            "ari_celltype", "batch_entropy"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df.columns]]


def variable_feature_list(obj, assay: str = "RNA") -> list[str]:
    """The assay's selected variable features, as a plain list of names."""
    a = obj.assays[assay]
    vf = getattr(a, "variable_features", None)
    if vf is None:
        vf = getattr(a, "var_features", None)
    return list(vf) if vf is not None else []


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_ifnb_object(data_dir=None, min_cells=3):
    """Load ``ifnb`` and build the RNA object carrying the stim + cell-type labels.

    Reads the counts exported by ``tutorials/export_seuratdata.R ifnb`` (raises a
    helpful error if that one-time export has not run) and keeps only the two
    metadata columns the tutorial needs: :data:`BATCH` and :data:`CELLTYPE`.
    """
    counts, genes, cells, meta = ifnb(data_dir=data_dir)
    keep = meta[[c for c in (BATCH, CELLTYPE) if c in meta.columns]].copy()
    obj = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=min_cells, min_features=0,
        project="ifnb", feature_names=list(genes), cell_names=list(cells),
        meta_data=keep,
    )
    return obj


def prep_reduction(obj, n_hvg=N_HVG, n_pcs=N_PCS, verbose=False):
    """Standard prep up to PCA — the uncorrected starting point every method shares.

    Normalize (LogNormalize) → variable features → scale → PCA, the vignette's
    ``NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA``. Returns
    the selected variable-feature list (the shared basis written for the R run).
    """
    normalize_data(obj, assay="RNA")
    find_variable_features(obj, assay="RNA", nfeatures=n_hvg)
    scale_data(obj, assay="RNA")
    run_pca(obj, assay="RNA", n_pcs=n_pcs)
    return variable_feature_list(obj, assay="RNA")


def cluster_and_embed(obj, reduction, key, resolution=RESOLUTION, n_pcs=N_PCS,
                      do_umap=True, seed=42):
    """Neighbours + Louvain + (optional) UMAP on one reduction, stored under ``key``.

    ``find_clusters`` overwrites ``seurat_clusters`` each call, so the labels are
    copied out to ``clusters_<key>`` and the UMAP is stored as ``umap_<key>`` —
    letting every method's partition and map coexist on the one object.
    """
    dims = range(min(n_pcs, obj.reductions[reduction].cell_embeddings.shape[1]))
    find_neighbors(obj, reduction=reduction, dims=dims)
    find_clusters(obj, resolution=resolution)
    obj.meta_data[f"clusters_{key}"] = np.asarray(obj.meta_data["seurat_clusters"])
    if do_umap:
        run_umap(obj, reduction=reduction, dims=dims,
                 reduction_name=f"umap_{key}", seed=seed)
    return obj


def run_uncorrected(obj, do_umap=True):
    """Baseline: cluster the raw PCA embedding (no batch correction)."""
    return cluster_and_embed(obj, reduction="pca", key="pca", do_umap=do_umap)


def run_integration(obj, method, group_by=BATCH, do_umap=True, seed=0):
    """One integration method → a corrected reduction, then cluster + embed it.

    ``harmony`` corrects the existing PCA in place; ``cca`` / ``rpca`` split by
    ``group_by``, anchor the batches and rebuild an integrated reduction. Each is
    stored under its method name and clustered identically to the baseline.
    """
    if method == "harmony":
        run_harmony(obj, group_by=group_by, reduction="pca",
                    reduction_name="harmony", seed=seed)
    else:
        integrate_layers(obj, method=method, group_by=group_by,
                         new_reduction=method, seed=42)
    return cluster_and_embed(obj, reduction=method, key=method, do_umap=do_umap)


def _scoreboard_row(obj, name, reduction, key):
    """Metrics for one reduction: batch/celltype silhouette, clusters, ARIs."""
    emb = obj.reductions[reduction].cell_embeddings
    batch = np.asarray(obj.meta_data[BATCH])
    celltype = np.asarray(obj.meta_data[CELLTYPE])
    clusters = np.asarray(obj.meta_data[f"clusters_{key}"])
    metrics = mixing_metrics(emb, batch, celltype)
    return {
        "method": name,
        "sil_batch": round(metrics["sil_batch"], 4),
        "sil_celltype": round(metrics["sil_celltype"], 4),
        "n_clusters": int(len(set(clusters))),
        "ari_celltype": round(cluster_ari(clusters, celltype), 4),
        "batch_entropy": round(batch_entropy(clusters, batch), 4),
    }


def summarize(obj, methods=METHODS, verbose=True) -> dict:
    """Build the integration scoreboard (uncorrected + each method) and return it."""
    rows = [_scoreboard_row(obj, "uncorrected (PCA)", "pca", "pca")]
    for m in methods:
        if m in obj.reductions:
            rows.append(_scoreboard_row(obj, m, m, m))
    board = build_scoreboard(rows)

    out = {
        "n_cells": int(len(obj.cell_names())),
        "n_batches": int(obj.meta_data[BATCH].nunique()),
        "n_celltypes": int(obj.meta_data[CELLTYPE].nunique()),
        "scoreboard": board,
    }
    if not verbose:
        return out

    section("Integration scoreboard (lower sil_batch = better mixed; "
            "higher sil_celltype / ari_celltype = better preserved)")
    print("    " + board.to_string(index=False).replace("\n", "\n    "))
    base = board.iloc[0]
    print(f"\n  Uncorrected: cells separate by {BATCH} "
          f"(sil_batch={base['sil_batch']}). Each method below lowers that while "
          f"holding or raising sil_celltype / ari_celltype.")
    return out


def report_concordance(obj, r_calls_path=None, verbose=True) -> dict | None:
    """Compare Python vs R integration, if ``ifnb_integration_verify.R`` has run.

    Reads ``figures_integration/r_calls.csv`` (per-cell R cluster labels for each
    method, plus stim and cell type) and reports, per method: the adjusted Rand
    index between the two tools' cluster partitions, and each tool's own
    cell-type recovery ARI and batch-mixing entropy — so a difference shows up as
    *both tools integrate about equally well* rather than *identical clusters*.
    Returns ``None`` (with a hint) when the R calls are absent.
    """
    path = Path(r_calls_path) if r_calls_path else FIGURES / "r_calls.csv"
    if not path.exists():
        if verbose:
            print(f"\n  [concordance] {path.name} not found — run "
                  "`Rscript tutorials/ifnb_integration_verify.R` first.")
        return None

    r = pd.read_csv(path).set_index("cell")
    cells = obj.cell_names()
    r = r.reindex(cells)
    meta = obj.meta_data
    celltype = np.asarray(meta[CELLTYPE])
    batch = np.asarray(meta[BATCH])

    out = {}
    rows = []
    for m in ("pca", *METHODS):
        rcol = f"R_{m}"
        pycol = f"clusters_{m}"
        if rcol not in r.columns or pycol not in meta.columns:
            continue
        py = np.asarray(meta[pycol])
        rv = r[rcol].to_numpy()
        out[m] = {
            "ari_partitions": cluster_ari(py, rv),
            "py_ari_celltype": cluster_ari(py, celltype),
            "r_ari_celltype": cluster_ari(rv, celltype),
            "py_batch_entropy": batch_entropy(py, batch),
            "r_batch_entropy": batch_entropy(rv, batch),
        }
        rows.append((m, out[m]["ari_partitions"], out[m]["py_ari_celltype"],
                     out[m]["r_ari_celltype"], out[m]["py_batch_entropy"],
                     out[m]["r_batch_entropy"]))

    if verbose:
        section("R-vs-Python concordance (shared counts & variable-feature basis)")
        tbl = pd.DataFrame(rows, columns=[
            "method", "ARI(py,R)", "py ARI→type", "R ARI→type",
            "py mix", "R mix"]).round(4)
        print("    " + tbl.to_string(index=False).replace("\n", "\n    "))
        print("\n  ARI(py,R): agreement of the two cluster partitions (1 = identical).")
        print("  ARI→type : each tool's clusters vs the known cell types — the")
        print("             biological check; the two columns should track closely.")
        print("  mix      : batch-mixing entropy (1 = fully mixed, 0 = unintegrated).")
    return out


def run_full(data_dir=None, verbose=True, methods=METHODS, do_umap=False):
    """Load, prep, run uncorrected + each integration method, score, compare to R."""
    t0 = time.time()
    if verbose:
        section("Loading ifnb — IFN-beta stimulated vs control PBMCs")
    obj = load_ifnb_object(data_dir=data_dir)
    if verbose:
        n_ctrl = int((obj.meta_data[BATCH] == "CTRL").sum())
        n_stim = int((obj.meta_data[BATCH] == "STIM").sum())
        print(f"  {len(obj.assays['RNA'].features())} genes x "
              f"{len(obj.cell_names())} cells | CTRL {n_ctrl} vs STIM {n_stim} | "
              f"{obj.meta_data[CELLTYPE].nunique()} cell types")

    hvg = prep_reduction(obj, verbose=verbose)
    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    if verbose:
        print(f"  variable features: {len(hvg)} (written to "
              f"figures_integration/hvg_features.txt for the R reference)")

    run_uncorrected(obj, do_umap=do_umap)
    for m in methods:
        t1 = time.time()
        run_integration(obj, m, do_umap=do_umap)
        if verbose:
            print(f"  {m}: integrated + clustered ({time.time() - t1:.1f}s)")

    summary = summarize(obj, methods=methods, verbose=verbose)
    report_concordance(obj, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return obj, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ifnb Harmony/CCA/RPCA integration tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--methods", nargs="+", default=list(METHODS),
                        choices=list(METHODS), help="which integration methods to run")
    parser.add_argument("--umap", action="store_true",
                        help="also compute UMAP embeddings (needed only for figures)")
    args = parser.parse_args()
    run_full(data_dir=args.data_dir, methods=tuple(args.methods), do_umap=args.umap)
