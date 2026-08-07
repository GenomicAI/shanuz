"""PBMC 3k Tutorial — Truecell Python implementation.

Mirrors the official Seurat PBMC 3k guided clustering tutorial step-by-step:
  https://satijalab.org/seurat/articles/pbmc3k_tutorial

Each step prints summary statistics comparable to the R tutorial so you can
validate the Python results directly.

Usage
-----
    python tutorials/pbmc3k_tutorial.py [--data-dir PATH]

If --data-dir is not supplied the PBMC 3k dataset is downloaded automatically
to ~/.truecell_data/pbmc3k (~24 MB, 10X Genomics).

References
----------
Hao et al. (2024) Nature Biotechnology — https://doi.org/10.1038/s41587-023-01767-y
Stuart et al. (2019) Cell — https://doi.org/10.1016/j.cell.2019.05.031
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Ensure the package root is on the path when running the script directly
# ---------------------------------------------------------------------------
_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from truecell.datasets import pbmc3k
from truecell.truecell import create_truecell_object
from truecell.preprocessing import (
    normalize_data,
    find_variable_features,
    scale_data,
    percentage_feature_set,
)
from truecell.reduction import run_pca
from truecell.neighbors import find_neighbors
from truecell.clustering import find_clusters
from truecell.umap import run_umap
from truecell.markers import find_markers, find_all_markers


# ---------------------------------------------------------------------------
# Known expected results from the R tutorial (used for validation)
# ---------------------------------------------------------------------------
EXPECTED = {
    # After CreateSeuratObject(min.cells=3, min.features=200)
    "n_features_raw": 13714,
    "n_cells_raw": 2700,
    # After QC filter (nFeature_RNA 200-2500, percent.mt < 5)
    "n_cells_filtered": 2638,
    # Top 10 HVGs from the R tutorial
    "top10_hvg": {
        "PPBP", "LYZ", "S100A9", "IGLL5", "GNLY",
        "FTL", "PF4", "FTH1", "GNG11", "S100A8",
    },
    # Canonical marker genes expected per cell type
    "canonical_markers": {
        "CD14+ Mono":    ["LYZ", "CD14", "S100A9"],
        "NK":            ["NKG7", "GNLY"],
        "B":             ["MS4A1", "CD79A"],
        "CD8 T":         ["CD8A"],
        "DC":            ["FCER1A"],
        "Platelet":      ["PPBP"],
    },
    # Number of clusters expected (resolution=0.5 → ~9)
    "n_clusters_expected": 9,
}


def validate(label: str, value, expected=None, atol: float = 0.05) -> None:
    """Print a validation line. Green check if matches, red ✗ otherwise."""
    if expected is None:
        print(f"  [INFO] {label}: {value}")
        return
    if isinstance(expected, set) and isinstance(value, (set, list)):
        overlap = set(value) & expected
        pct = len(overlap) / len(expected) * 100
        ok = pct >= 50
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {label}: overlap {len(overlap)}/{len(expected)} "
              f"({pct:.0f}%)  got={sorted(set(value))[:10]}")
    elif isinstance(expected, int):
        ok = abs(value - expected) <= max(1, int(expected * atol))
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {label}: {value}  (expected ~{expected})")
    else:
        print(f"  [INFO] {label}: {value}  (expected {expected})")


def section(title: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}")


def top_markers_per_cluster(all_markers, n: int = 3) -> dict:
    """Map each cluster to its `n` most significant marker genes.

    Filters per cluster rather than grouping: `groupby(...).apply(...)` stopped
    passing the grouping column into the callable in pandas 3, which silently
    drops `cluster` from the result instead of raising. Mirrors
    `pbmc8k_subclustering_tutorial.top_markers_table`.
    """
    return {
        cluster: all_markers[all_markers["cluster"] == cluster]
        .nsmallest(n, "p_val")["gene"]
        .tolist()
        for cluster in sorted(all_markers["cluster"].unique(), key=int)
    }


# ---------------------------------------------------------------------------
# Main tutorial
# ---------------------------------------------------------------------------

def run_tutorial(data_dir: str | None = None) -> None:
    t0_total = time.time()

    # -----------------------------------------------------------------------
    section("1. Load Data")
    # -----------------------------------------------------------------------
    t0 = time.time()
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    print(f"  Raw matrix: {counts.shape[0]} genes × {counts.shape[1]} cells  "
          f"({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("2. Create Truecell Object (min.cells=3, min.features=200)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    pbmc = create_truecell_object(
        counts=counts,
        assay="RNA",
        min_cells=3,
        min_features=200,
        project="pbmc3k",
        feature_names=genes,
        cell_names=cells,
    )
    n_feat = len(pbmc.assays["RNA"]._all_feature_names)
    n_cells = len(pbmc.cell_names())
    print(f"  {n_feat} features × {n_cells} cells  ({time.time() - t0:.1f}s)")
    # Stash the pre-QC dimensions in `misc` (Seurat's own slot for this, and it
    # survives subset()); both are handoff anchors, and once the QC filter runs
    # they cannot be recovered from the object.
    pbmc.misc["n_genes_raw"] = n_feat
    pbmc.misc["n_cells_raw"] = n_cells
    validate("n_features after min.cells=3", n_feat, EXPECTED["n_features_raw"])
    validate("n_cells after min.features=200", n_cells, EXPECTED["n_cells_raw"])

    # -----------------------------------------------------------------------
    section("3. QC Metrics")
    # -----------------------------------------------------------------------
    t0 = time.time()
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    md = pbmc.meta_data
    print(f"  nFeature_RNA: mean={md['nFeature_RNA'].mean():.0f}  "
          f"min={md['nFeature_RNA'].min():.0f}  "
          f"max={md['nFeature_RNA'].max():.0f}")
    print(f"  percent.mt:   mean={md['percent.mt'].mean():.2f}%  "
          f"max={md['percent.mt'].max():.2f}%  ({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("4. Filter Cells (nFeature 200-2500, percent.mt < 5)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    keep = (
        (md["nFeature_RNA"] > 200) &
        (md["nFeature_RNA"] < 2500) &
        (md["percent.mt"] < 5)
    )
    keep_cells = list(md.index[keep])
    pbmc = pbmc.subset(cells=keep_cells)
    n_cells_filt = len(pbmc.cell_names())
    print(f"  {n_cells_filt} cells retained  ({time.time() - t0:.1f}s)")
    validate("n_cells after QC filter", n_cells_filt, EXPECTED["n_cells_filtered"])

    # -----------------------------------------------------------------------
    section("5. Normalize Data (LogNormalize, scale.factor=10000)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    print(f"  Log-normalization complete  ({time.time() - t0:.1f}s)")

    # Quick sanity check: mean of normalized data should be ~1-3 log units
    rna = pbmc.assays["RNA"]
    norm = rna.layers["data"]
    mean_expr = float(np.array(norm.mean()))
    print(f"  Mean log-normalized expression: {mean_expr:.4f}")

    # -----------------------------------------------------------------------
    section("6. Find Variable Features (VST, nfeatures=2000)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    hvg = pbmc.assays["RNA"].variable_features
    top10 = hvg[:10]
    print(f"  {len(hvg)} variable features selected  ({time.time() - t0:.1f}s)")
    print(f"  Top 10 HVGs: {top10}")
    validate("Top 10 HVG overlap with R tutorial", set(top10), EXPECTED["top10_hvg"])

    # -----------------------------------------------------------------------
    section("7. Scale Data")
    # -----------------------------------------------------------------------
    t0 = time.time()
    all_genes = pbmc.assays["RNA"]._all_feature_names
    scale_data(pbmc, features=all_genes)
    print(f"  Scaled {len(all_genes)} genes  ({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("8. Run PCA (npc=50)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    run_pca(pbmc, n_pcs=50, features=hvg, reduction_name="pca")
    pca_emb = pbmc.reductions["pca"].cell_embeddings
    print(f"  PCA: {pca_emb.shape[0]} cells × {pca_emb.shape[1]} PCs  "
          f"({time.time() - t0:.1f}s)")
    stdev = pbmc.reductions["pca"].stdev
    print(f"  PC1 stdev={stdev[0]:.3f}  PC2 stdev={stdev[1]:.3f}  "
          f"PC10 stdev={stdev[9]:.3f}")
    # Top loadings of PC1 (should be ribosomal/mitochondrial or strong cell-type genes)
    loadings = pbmc.reductions["pca"].feature_loadings
    feat_names = pbmc.reductions["pca"]._feature_names
    top_pc1 = [feat_names[i] for i in np.argsort(np.abs(loadings[:, 0]))[::-1][:5]]
    print(f"  Top PC1 loading genes: {top_pc1}")

    # -----------------------------------------------------------------------
    section("9. Find Neighbors (dims=1:10, k=20)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_neighbors(pbmc, dims=range(10), k_param=20)
    print(f"  KNN+SNN graphs built  ({time.time() - t0:.1f}s)")
    print(f"  Graphs: {list(pbmc.graphs)}")

    # -----------------------------------------------------------------------
    section("10. Find Clusters (resolution=0.5)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)
    cluster_counts = pbmc.meta_data["seurat_clusters"].value_counts().sort_index()
    n_clusters = len(cluster_counts)
    print(f"  {n_clusters} clusters found  ({time.time() - t0:.1f}s)")
    print("  Cells per cluster:")
    for c, n in cluster_counts.items():
        print(f"    Cluster {c}: {n} cells")
    validate("Number of clusters", n_clusters, EXPECTED["n_clusters_expected"])

    # -----------------------------------------------------------------------
    section("11. Run UMAP (dims=1:10)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    run_umap(pbmc, dims=range(10), reduction_name="umap", seed=42)
    umap_emb = pbmc.reductions["umap"].cell_embeddings
    print(f"  UMAP: {umap_emb.shape}  ({time.time() - t0:.1f}s)")
    print(f"  UMAP range: x=[{umap_emb[:,0].min():.2f},{umap_emb[:,0].max():.2f}]  "
          f"y=[{umap_emb[:,1].min():.2f},{umap_emb[:,1].max():.2f}]")

    # -----------------------------------------------------------------------
    section("12. Find Cluster 2 Markers")
    # -----------------------------------------------------------------------
    t0 = time.time()
    c2_markers = find_markers(pbmc, ident_1="2", only_pos=True)
    print(f"  Cluster 2 markers: {len(c2_markers)} genes  ({time.time() - t0:.1f}s)")
    if len(c2_markers) > 0:
        top5 = c2_markers.head(5)
        print("  Top 5 cluster 2 markers:")
        print(top5[["avg_log2FC", "pct.1", "pct.2", "p_val_adj"]].to_string())

    # -----------------------------------------------------------------------
    section("13. Find All Markers (only.pos=True, logfc.threshold=0.25)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    all_markers = find_all_markers(
        pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25
    )
    print(f"  Total marker genes: {len(all_markers)}  ({time.time() - t0:.1f}s)")
    print("\n  Top 3 markers per cluster:")
    for cluster, genes in top_markers_per_cluster(all_markers, n=3).items():
        print(f"    Cluster {cluster}: {', '.join(genes)}")

    # -----------------------------------------------------------------------
    section("14. Validate Canonical Marker Expression")
    # -----------------------------------------------------------------------
    found_in_markers = set(all_markers["gene"].tolist())
    for cell_type, canon_genes in EXPECTED["canonical_markers"].items():
        hits = [g for g in canon_genes if g in found_in_markers]
        pct = len(hits) / len(canon_genes) * 100
        mark = "OK" if pct >= 50 else "WARN"
        print(f"  [{mark}] {cell_type}: {hits} ({pct:.0f}% canonical markers found)")

    # -----------------------------------------------------------------------
    section("15. Cell Type Annotation (R tutorial mapping)")
    # -----------------------------------------------------------------------
    # Map cluster labels to cell types using known marker patterns
    # This is R's RenameIdents() step
    cluster_to_celltype = _assign_cell_types(all_markers, pbmc)
    names = {str(k): v for k, v in cluster_to_celltype.items()}
    pbmc.rename_idents(names)
    pbmc.meta_data["celltype"] = [str(i) for i in pbmc.idents]

    celltype_counts = pd.Series(list(pbmc.idents)).value_counts()
    print("\n  Cell type distribution:")
    for ct, n in celltype_counts.sort_values(ascending=False).items():
        print(f"    {ct}: {n} cells")

    # -----------------------------------------------------------------------
    section("Summary")
    # -----------------------------------------------------------------------
    total = time.time() - t0_total
    print(f"\n  Total runtime: {total:.1f}s")
    print(f"\n{pbmc}")

    return pbmc, all_markers


def _assign_cell_types(
    all_markers: pd.DataFrame,
    pbmc,
) -> dict[int, str]:
    """Heuristically assign cell types to clusters based on top markers.

    Mirrors the manual annotation step in the R tutorial.
    """
    # Canonical markers for each cell type (ordered by specificity)
    markers_ref = {
        "Naive CD4 T":   ["IL7R", "CCR7"],
        "CD14+ Mono":    ["CD14", "LYZ"],
        "Memory CD4 T":  ["IL7R", "S100A4"],
        "B":             ["MS4A1"],
        "CD8 T":         ["CD8A"],
        "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
        "NK":            ["GNLY", "NKG7"],
        "DC":            ["FCER1A", "CST3"],
        "Platelet":      ["PPBP"],
    }

    # Sort numerically, not lexicographically: a type is consumed once assigned
    # (`used_types` below), so the loop order is part of the definition, and
    # string order would visit cluster 10 before cluster 2 and hand out the
    # panels in a different sequence. `pbmc3k_verify.R` iterates the same way.
    clusters = sorted(set(str(i) for i in pbmc.idents), key=int)
    cluster_top_genes: dict[str, set] = {}
    for cluster in clusters:
        sub = all_markers[all_markers["cluster"] == cluster].head(50)
        cluster_top_genes[cluster] = set(sub["gene"].tolist())

    assignment: dict[str, str] = {}
    used_types: set[str] = set()

    for cluster in clusters:
        top_genes = cluster_top_genes.get(cluster, set())
        best_type = "Unknown"
        best_score = 0
        for cell_type, canon in markers_ref.items():
            if cell_type in used_types:
                continue
            score = sum(1 for g in canon if g in top_genes)
            if score > best_score:
                best_score = score
                best_type = cell_type
        if best_score > 0:
            used_types.add(best_type)
        assignment[cluster] = best_type

    return assignment


# ---------------------------------------------------------------------------
# Numeric handoff against R
# ---------------------------------------------------------------------------
#
# Every figure in this tutorial links the canonical satijalab.org image, which
# means the whole guided workflow was compared by eye. These files are not.
#
# Nothing is pinned across the two sides: R runs its own pipeline from the same
# 10x bytes. That is deliberate. `pbmc3k_dimreduc_verify.R` pins Python's cells
# and features to isolate the post-PCA machinery; here the question is the one
# the tutorial itself makes — do two independent runs keep the same cells,
# choose the same variable genes, and land on the same clusters and markers?
#
# Per-cell rows are keyed by barcode and per-gene rows by symbol so the report
# can compare on a shared key. Cluster *numbering* is arbitrary in both tools,
# so cluster comparisons match the two partitions one-to-one on overlap before
# scoring; without that step a pure relabelling reads as total disagreement.

FIGURES = Path(__file__).parent / "figures"

# R's Read10X() rewrites underscores in gene symbols to dashes ("RP11-34P13_3"
# -> "RP11-34P13-3"); truecell's loader keeps the file's spelling. Map Python's
# symbols through the same rule before joining, or ~30 genes silently drop out
# of every per-gene comparison.
def _r_symbols(genes) -> list:
    return [str(g).replace("_", "-") for g in genes]


def match_partitions(a, b) -> dict:
    """Match two cluster labellings of the same cells one-to-one on overlap.

    Cluster ids are arbitrary in both tools — truecell's cluster 3 may be
    Seurat's cluster 5 — so a coordinate-wise comparison of the labels answers
    the wrong question. The Hungarian algorithm on the contingency table finds
    the pairing that maximises the number of cells the two agree on, and
    ``concordance`` is the fraction of cells that fall inside a matched pair.
    ``ari`` is label-invariant already and needs no matching; it is reported
    alongside because it also penalises splits and merges, which a
    best-matching concordance quietly forgives.
    """
    from scipy.optimize import linear_sum_assignment
    from sklearn.metrics import adjusted_rand_score

    a = np.asarray([str(x) for x in a])
    b = np.asarray([str(x) for x in b])
    if a.shape != b.shape or a.size == 0:
        raise ValueError("both labellings must cover the same non-empty cells")
    la, lb = sorted(set(a)), sorted(set(b))
    table = np.zeros((len(la), len(lb)), dtype=int)
    for i, x in enumerate(la):
        mask = a == x
        for j, y in enumerate(lb):
            table[i, j] = int(np.count_nonzero(mask & (b == y)))
    rows, cols = linear_sum_assignment(-table)
    return {
        "ari": float(adjusted_rand_score(a, b)),
        "concordance": float(table[rows, cols].sum() / a.size),
        "mapping": {la[i]: lb[j] for i, j in zip(rows, cols)},
        "table": pd.DataFrame(table, index=la, columns=lb),
        "n_a": len(la),
        "n_b": len(lb),
    }


def write_anchors(pbmc, all_markers) -> None:
    """Dump the per-cell, per-gene and per-marker tables plus the scalars."""
    import json

    FIGURES.mkdir(exist_ok=True)
    md = pbmc.meta_data
    emb = pbmc.reductions["pca"].cell_embeddings[:, :10]
    cells = pd.DataFrame({
        "cell": list(pbmc.cell_names()),
        "nCount_RNA": md["nCount_RNA"].to_numpy(),
        "nFeature_RNA": md["nFeature_RNA"].to_numpy(),
        "percent.mt": md["percent.mt"].to_numpy(),
        "cluster": md["seurat_clusters"].astype(str).to_numpy(),
        "celltype": md["celltype"].astype(str).to_numpy(),
    })
    for k in range(10):
        cells[f"PC_{k + 1}"] = emb[:, k]
    cells.to_csv(FIGURES / "py_cell_meta.csv", index=False)

    rna = pbmc.assays["RNA"]
    # truecell stores the VST statistics on the assay's meta_data under the names
    # HVFInfo() uses, so these columns carry straight through to the R side of
    # the handoff without being renamed on the way out.
    hvf = rna.meta_data
    selected = list(rna.variable_features)
    rank = {g: i + 1 for i, g in enumerate(selected)}
    pd.DataFrame({
        "gene": _r_symbols(hvf.index),
        "mean": hvf["mean"].to_numpy(),
        "variance": hvf["variance"].to_numpy(),
        "var.expected": hvf["variance.expected"].to_numpy(),
        "var.std": hvf["variance.standardized"].to_numpy(),
        "hvg_rank": [rank.get(g, np.nan) for g in hvf.index],
    }).to_csv(FIGURES / "py_hvg.csv", index=False)

    markers = all_markers.copy()
    markers["gene"] = _r_symbols(markers["gene"])
    markers[["cluster", "gene", "avg_log2FC", "pct.1", "pct.2",
             "p_val", "p_val_adj"]].to_csv(FIGURES / "py_markers.csv", index=False)

    knn = pbmc.graphs["RNA_nn"]
    snn = pbmc.graphs["RNA_snn"]
    sizes = md["seurat_clusters"].astype(str).value_counts()
    anchors = {
        "n_genes_raw": int(pbmc.misc["n_genes_raw"]),
        "n_cells_raw": int(pbmc.misc["n_cells_raw"]),
        "n_cells_qc": len(pbmc.cell_names()),
        "n_hvg": len(selected),
        "n_clusters": int(md["seurat_clusters"].nunique()),
        "n_markers": int(len(all_markers)),
        "knn_nnz": int(knn.nnz),
        "snn_nnz": int(snn.nnz),
        "snn_weight_sum": float(snn.sum()),
        "data_sum": float(rna.layers["data"].sum()),
        "pc_stdev": [float(s) for s in pbmc.reductions["pca"].stdev[:10]],
        "cluster_sizes": {c: int(sizes[c]) for c in sorted(sizes.index, key=int)},
    }
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    print(f"\n  Wrote py_cell_meta.csv, py_hvg.csv, py_markers.csv and "
          f"py_anchors.json to {FIGURES}")
    print("  Next: Rscript tutorials/pbmc3k_verify.R"
          "  then  python tutorials/pbmc3k_tutorial.py --report")


def report() -> None:
    """Compare the whole guided workflow against the R reference."""
    import json

    from scipy.stats import pearsonr, spearmanr

    need = ["py_cell_meta.csv", "r_cell_meta.csv", "py_hvg.csv", "r_hvg.csv",
            "py_markers.csv", "r_markers.csv", "py_anchors.json", "r_anchors.json"]
    missing = [f for f in need if not (FIGURES / f).exists()]
    if missing:
        print(f"  missing {missing} — run the tutorial and then "
              f"`Rscript tutorials/pbmc3k_verify.R`")
        return

    print("=" * 78)
    print("truecell vs Seurat 5.5.1 — PBMC 3k guided clustering, end to end")
    print("=" * 78)

    pc = pd.read_csv(FIGURES / "py_cell_meta.csv", float_precision="round_trip").set_index("cell")
    rc = pd.read_csv(FIGURES / "r_cell_meta.csv", float_precision="round_trip").set_index("cell")

    # ---- 1. did QC keep the same cells, with the same metrics? -------------
    only_py = pc.index.difference(rc.index)
    only_r = rc.index.difference(pc.index)
    cells = pc.index.intersection(rc.index)
    print(f"\n  QC — cells retained: truecell {len(pc)}, R {len(rc)}, "
          f"shared {len(cells)}")
    print(f"       only truecell {len(only_py)}   only R {len(only_r)}   "
          f"{'IDENTICAL CELL SET' if not len(only_py) and not len(only_r) else 'SETS DIFFER'}")
    p, r = pc.loc[cells], rc.loc[cells]
    print(f"\n  {'per-cell metric':<18}{'max|diff|':>14}")
    print(f"  {'-' * 32}")
    for col in ("nCount_RNA", "nFeature_RNA", "percent.mt"):
        d = np.abs(p[col].to_numpy() - r[col].to_numpy()).max()
        print(f"  {col:<18}{d:>14.3e}")

    # ---- 2. the VST statistics and the 2,000 it selected -------------------
    ph = pd.read_csv(FIGURES / "py_hvg.csv", float_precision="round_trip").set_index("gene")
    rh = pd.read_csv(FIGURES / "r_hvg.csv", float_precision="round_trip").set_index("gene")
    genes = ph.index.intersection(rh.index)
    print(f"\n  Variable features (VST) — {len(genes)} shared genes "
          f"of {len(ph)} / {len(rh)}")
    print(f"  {'per-gene statistic':<20}{'max|diff|':>14}{'max rel':>14}")
    print(f"  {'-' * 48}")
    # `var.expected` is the LOESS fit, and `var.std` is proportional to its
    # reciprocal — so the two columns carry the same disagreement, and having
    # both is what separates "the standardization drifted" from "the fit did".
    #
    # Relative as well as absolute, because these four live on wildly different
    # scales: an expected variance runs into the hundreds for an abundant gene,
    # so its absolute gap reads as enormous next to a mean's and says nothing.
    for col in ("mean", "variance", "var.expected", "var.std"):
        a = ph.loc[genes, col].to_numpy()
        b = rh.loc[genes, col].to_numpy()
        d = np.abs(a - b)
        rel = d / np.maximum(np.abs(b), np.finfo(float).tiny)
        print(f"  {col:<20}{d.max():>14.3e}{rel.max():>14.3e}")
    py_sel = set(ph.index[ph["hvg_rank"].notna()])
    r_sel = set(rh.index[rh["hvg_rank"].notna()])
    both = py_sel & r_sel
    print(f"    selected set  truecell {len(py_sel)}  R {len(r_sel)}  "
          f"shared {len(both)} ({len(both) / max(len(r_sel), 1):.4f})")
    if both:
        rr = spearmanr(ph.loc[sorted(both), "hvg_rank"],
                       rh.loc[sorted(both), "hvg_rank"]).statistic
        print(f"    rank agreement on the shared selection: spearman {rr:.4f}")

    # ---- 3. PCA, matched one-to-one because sign and order are arbitrary ---
    pe = p[[f"PC_{k}" for k in range(1, 11)]].to_numpy()
    re_ = r[[f"PC_{k}" for k in range(1, 11)]].to_numpy()
    from scipy.optimize import linear_sum_assignment
    corr = np.abs(np.corrcoef(pe.T, re_.T)[:10, 10:])
    rows, cols = linear_sum_assignment(-corr)
    matched = corr[rows, cols]
    print("\n  PCA — the 10 dims the clustering runs on, matched on |r|")
    print(f"    mean |r| {matched.mean():.6f}   min |r| {matched.min():.6f}   "
          f"in order: {'yes' if list(rows) == list(cols) else 'NO — reordered'}")

    # ---- 4. clusters: match the partitions, then score ---------------------
    m = match_partitions(p["cluster"], r["cluster"])
    print(f"\n  Clusters — truecell {m['n_a']}, R {m['n_b']}")
    print(f"    adjusted Rand index {m['ari']:.4f}   "
          f"best-match concordance {m['concordance']:.4f} "
          f"({int(round(m['concordance'] * len(cells)))}/{len(cells)} cells)")
    print(f"\n    {'truecell':>8}{'R':>6}{'n py':>8}{'n R':>8}{'shared':>9}"
          f"   cell type py / R")
    print(f"    {'-' * 62}")
    for a, b in m["mapping"].items():
        shared = int(m["table"].loc[a, b])
        ma, mb = p["cluster"].astype(str) == a, r["cluster"].astype(str) == b
        la, lb = p.loc[ma, "celltype"].mode().iat[0], r.loc[mb, "celltype"].mode().iat[0]
        flag = "" if la == lb else "   <-- differs"
        print(f"    {a:>8}{b:>6}{int(ma.sum()):>8}{int(mb.sum()):>8}{shared:>9}"
              f"   {la} / {lb}{flag}")
    # The counts above differ, so one side has clusters with no partner. Those
    # are the whole difference between the two runs and must not be left out of
    # the table just because the matching had nowhere to put them — so name
    # them, and say where the other side put their cells instead.
    for label, side, other in (("truecell", p, r), ("R", r, p)):
        paired = set(m["mapping"]) if label == "truecell" else set(m["mapping"].values())
        for c in sorted(set(side["cluster"].astype(str)) - paired, key=int):
            mask = side["cluster"].astype(str) == c
            lands = other.loc[mask, "cluster"].astype(str).value_counts()
            where = ", ".join(f"{k} ({v})" for k, v in lands.head(3).items())
            print(f"    unmatched {label} cluster {c}: {int(mask.sum())} cells, "
                  f"labelled {side.loc[mask, 'celltype'].mode().iat[0]}")
            print(f"      -> they sit in the other run's cluster {where}")

    # ---- 5. the cell-type labels, which share a vocabulary -----------------
    # Read this against the table above, not on its own. `_assign_cell_types`
    # is a greedy heuristic that consumes each cell type once, scoring clusters
    # on whether canonical genes reach their top 50 markers — so a cluster the
    # two tools agree about, cell for cell, can still come out with different
    # labels if one marker slipped in or out of a top-50 list. Where that
    # happens the "<-- differs" flag above is the honest reading, not this
    # number.
    same = (p["celltype"].astype(str).to_numpy() == r["celltype"].astype(str).to_numpy())
    agree = sum(1 for a, b in m["mapping"].items()
                if p.loc[p["cluster"].astype(str) == a, "celltype"].mode().iat[0]
                == r.loc[r["cluster"].astype(str) == b, "celltype"].mode().iat[0])
    print(f"\n  Cell types — per-cell label concordance {same.mean():.4f} "
          f"({same.sum()}/{len(cells)}); "
          f"{agree}/{len(m['mapping'])} matched clusters carry the same label")
    py_ct = p["celltype"].value_counts()
    r_ct = r["celltype"].value_counts()
    print(f"\n    {'label':<16}{'truecell':>9}{'R':>9}")
    print(f"    {'-' * 34}")
    for ct in sorted(set(py_ct.index) | set(r_ct.index)):
        print(f"    {ct:<16}{py_ct.get(ct, 0):>9}{r_ct.get(ct, 0):>9}")

    # ---- 6. markers, on the matched clusters -------------------------------
    pm = pd.read_csv(FIGURES / "py_markers.csv", float_precision="round_trip")
    rm = pd.read_csv(FIGURES / "r_markers.csv", float_precision="round_trip")
    pm["cluster"] = pm["cluster"].astype(str)
    rm["cluster"] = rm["cluster"].astype(str)
    print(f"\n  Markers — {len(pm)} rows truecell, {len(rm)} R "
          f"(only.pos, min.pct 0.25, logfc 0.25, return.thresh 0.01)")
    print("  Rows are per matched cluster pair. `top10 by p` is scored twice"
          " because R's\n  Wilcoxon p-values underflow to exactly 0 for the"
          " strongest markers, so its\n  own ordering among them is decided by"
          " the log2FC tie-break, not by p.")
    print(f"\n    {'py→R':>8}{'R ties@p=0':>12}{'top10 by p':>12}"
          f"{'top10 by FC':>13}{'genes shared':>14}{'log2FC r':>10}{'max|Δ|':>10}")
    print(f"    {'-' * 79}")
    for a, b in m["mapping"].items():
        pa = pm[pm["cluster"] == a].set_index("gene")
        rb = rm[rm["cluster"] == b].set_index("gene")
        ties = int((rb["p_val"] == 0).sum())
        top_p = set(pa.nsmallest(10, "p_val").index) & set(rb.nsmallest(10, "p_val").index)
        top_f = set(pa.nlargest(10, "avg_log2FC").index) & set(rb.nlargest(10, "avg_log2FC").index)
        shared = pa.index.intersection(rb.index)
        if len(shared) >= 2:
            x = pa.loc[shared, "avg_log2FC"].to_numpy()
            y = rb.loc[shared, "avg_log2FC"].to_numpy()
            rr = f"{pearsonr(x, y).statistic:.4f}"
            mx = f"{np.abs(x - y).max():.2e}"
        else:
            rr, mx = "—", "—"
        denom = max(len(pa.index.union(rb.index)), 1)
        print(f"    {a + '→' + b:>8}{ties:>12}{len(top_p):>9}/10{len(top_f):>10}/10"
              f"{len(shared):>9}/{denom:<4}{rr:>10}{mx:>10}")

    # ---- 7. scalars ---------------------------------------------------------
    pa_ = json.loads((FIGURES / "py_anchors.json").read_text())
    ra_ = json.loads((FIGURES / "r_anchors.json").read_text())
    print(f"\n  {'anchor':<16}{'truecell':>22}{'R Seurat':>22}   verdict")
    print(f"  {'-' * 66}")
    for k in ("n_genes_raw", "n_cells_raw", "n_cells_qc", "n_hvg", "n_clusters",
              "n_markers", "knn_nnz", "snn_nnz", "snn_weight_sum", "data_sum"):
        if k not in pa_ or k not in ra_:
            continue
        x, y = pa_[k], ra_[k]
        if isinstance(x, int) and isinstance(y, int):
            verdict = "MATCH" if x == y else f"differ by {x - y:+d}"
            xs, ys = str(x), str(y)
        else:
            verdict = f"rel {abs(x - y) / max(abs(y), 1e-12):.2e}"
            xs, ys = f"{x:.6f}", f"{y:.6f}"
        print(f"  {k:<16}{xs:>22}{ys:>22}   {verdict}")
    sd = np.abs(np.array(pa_["pc_stdev"]) - np.array(ra_["pc_stdev"]))
    print(f"  {'pc_stdev[1:10]':<16}{'':>44}   max|Δ| {sd.max():.3e}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 3k tutorial")
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Directory for PBMC3k data (default: ~/.truecell_data/pbmc3k)",
    )
    parser.add_argument("--report", action="store_true",
                        help="compare against the R reference and exit")
    args = parser.parse_args()
    if args.report:
        report()
    else:
        pbmc, all_markers = run_tutorial(data_dir=args.data_dir)
        write_anchors(pbmc, all_markers)
