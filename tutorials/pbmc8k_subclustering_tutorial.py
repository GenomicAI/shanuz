"""Advanced PBMC 8k Tutorial — Clustering + Subclustering with Shanuz.

A more complex companion to the PBMC 3k tutorial. It reproduces the standard
Seurat guided-clustering workflow (Satija et al. 2015; Butler et al. 2018) on a
larger 10x Genomics dataset (~8,400 PBMCs, GRCh38) and then adds the advanced
*subclustering* step used throughout the Seurat reference papers to resolve
fine-grained immune subsets: the T/NK lymphoid compartment is isolated and
re-analysed from scratch (HVG -> PCA -> neighbours -> clusters -> UMAP) to
separate naive CD4, memory CD4, CD8, and NK populations that the global
clustering lumps together.

Usage
-----
    python tutorials/pbmc8k_subclustering_tutorial.py [--data-dir PATH]

The PBMC 8k dataset (~38 MB) downloads automatically to ~/.shanuz_data/pbmc8k.

References
----------
Satija R, Farrell JA, Gennert D, et al. (2015) Nature Biotechnology 33, 495-502.
Butler A, Hoffman P, Smibert P, et al. (2018) Nature Biotechnology 36, 411-420.
Hao Y, Hao S, Andersen-Nissen E, et al. (2021) Cell 184, 3573-3587.
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

from shanuz.datasets import pbmc8k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data, percentage_feature_set,
)
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_all_markers
from shanuz.plotting import _get_expression


# ---------------------------------------------------------------------------
# Canonical marker panels
# ---------------------------------------------------------------------------

# Broad lineages — each cluster is assigned to the lineage whose markers it
# expresses most strongly (reuse allowed: several clusters can be "CD4 T").
BROAD_MARKERS = {
    "CD4 T":         ["IL7R", "CD3D", "CCR7"],
    "CD8 T":         ["CD8A", "CD8B", "GZMK"],
    "NK":            ["GNLY", "NKG7", "KLRD1"],
    "B":             ["MS4A1", "CD79A", "CD79B"],
    "CD14+ Mono":    ["CD14", "LYZ", "S100A8"],
    "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
    "DC":            ["FCER1A", "CST3"],
    "Platelet":      ["PPBP", "PF4"],
}

# The lymphoid lineages we re-analyse (subcluster) together.
LYMPHOID_LINEAGES = {"CD4 T", "CD8 T", "NK"}

# Genes used by the hierarchical T/NK subset annotator.
TNK_PANEL = ["CD3D", "CD3E", "CD8A", "CD8B", "GNLY", "NKG7", "KLRD1",
             "CCR7", "SELL", "LEF1", "IL7R", "S100A4", "GZMK"]


# ---------------------------------------------------------------------------
# Core pipeline (shared with generate_advanced_plots.py)
# ---------------------------------------------------------------------------

def run_pipeline(pbmc, dims=range(10), resolution=0.5, n_pcs=50, k_param=20,
                 nfeatures=2000, normalize=True):
    """Run the standard workflow on a (possibly already-subset) Shanuz object."""
    if normalize:
        normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(pbmc, selection_method="vst", nfeatures=nfeatures)
    scale_data(pbmc, features=pbmc.assays["RNA"]._all_feature_names)
    run_pca(pbmc, n_pcs=n_pcs, features=pbmc.assays["RNA"].variable_features,
            reduction_name="pca")
    find_neighbors(pbmc, dims=dims, k_param=k_param)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=0)
    run_umap(pbmc, dims=dims, reduction_name="umap", seed=42)
    return pbmc


def load_object(data_dir=None):
    counts, genes, cells = pbmc8k(data_dir=data_dir)
    pbmc = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc8k", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    return pbmc


def qc_filter(pbmc, max_features=2500, max_mt=5.0):
    md = pbmc.meta_data
    keep = (
        (md["nFeature_RNA"] > 200) &
        (md["nFeature_RNA"] < max_features) &
        (md["percent.mt"] < max_mt)
    )
    return pbmc.subset(cells=list(md.index[keep]))


def annotate_clusters(pbmc, marker_sets):
    """Assign each cluster to the lineage whose markers are most *enriched*.

    Each marker's per-cluster mean expression is z-scored across clusters, so a
    cluster scores on a lineage by how relatively enriched its markers are
    (CD8A is lower-magnitude than IL7R but still flags the CD8 cluster). Each
    lineage's score is the mean z-score of its present markers; argmax wins.

    Returns {cluster_label: lineage_name}. Reuse is allowed so several clusters
    can map to the same lineage (e.g. multiple T-cell clusters).
    """
    idents = np.array([str(i) for i in pbmc.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(pbmc.assays["RNA"]._all_feature_names)

    needed = {g for gs in marker_sets.values() for g in gs if g in feats}
    # z-score of each marker's per-cluster mean across clusters
    zmean = {}
    for g in needed:
        expr = _get_expression(pbmc, g)
        per_cluster = np.array([expr[idents == c].mean() for c in clusters])
        sd = per_cluster.std()
        zmean[g] = (per_cluster - per_cluster.mean()) / sd if sd > 1e-9 \
            else np.zeros(len(clusters))

    assignment = {}
    for ci, c in enumerate(clusters):
        best, best_score = "Unknown", -np.inf
        for lineage, genes in marker_sets.items():
            present = [g for g in genes if g in zmean]
            if not present:
                continue
            score = float(np.mean([zmean[g][ci] for g in present]))
            if score > best_score:
                best_score, best = score, lineage
        assignment[c] = best
    return assignment


def annotate_tnk_subsets(sub):
    """Annotate T/NK subclusters with biologically-ordered gating.

    Flat argmax over marker panels fails here because the high-magnitude naive
    markers (CCR7/SELL) outscore the lower-magnitude but definitive CD8 markers.
    We instead gate in lineage-priority order on mean expression:

      1. NK            — CD3 low and NKG7/GNLY high (NK cells are CD3-negative)
      2. CD8 T         — CD8A/B detectable, or a CD3+ cytotoxic (NKG7/GZMK)
                         program (captures CD8 effector plus MAIT/gamma-delta T,
                         which are CD8-lineage cytotoxic cells, never CD4)
      3. CD4 Naive     — naive markers (CCR7/SELL/LEF1) high
      4. CD4 Memory    — otherwise (IL7R / S100A4)
    """
    idents = np.array([str(i) for i in sub.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(sub.assays["RNA"]._all_feature_names)
    expr = {g: _get_expression(sub, g) for g in TNK_PANEL if g in feats}

    def m(genes, mask):
        vals = [expr[g][mask].mean() for g in genes if g in expr]
        return float(np.mean(vals)) if vals else 0.0

    assignment = {}
    for c in clusters:
        mask = idents == c
        cd3 = m(["CD3D", "CD3E"], mask)
        cd8 = m(["CD8A", "CD8B"], mask)
        nk = m(["GNLY", "NKG7"], mask)
        cyto = m(["NKG7", "GZMK"], mask)
        naive = m(["CCR7", "SELL", "LEF1"], mask)
        if cd3 < 0.75 and nk > 1.5:
            assignment[c] = "NK"
        elif cd8 > 0.6 or (cd3 >= 0.75 and cyto > 1.5):
            assignment[c] = "CD8 T"
        elif naive > 0.9:
            assignment[c] = "CD4 Naive"
        else:
            assignment[c] = "CD4 Memory"
    return assignment


# ---------------------------------------------------------------------------
# Pretty-printing helpers
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 64}\n  {title}\n{'=' * 64}")


def top_markers_table(all_markers, n=3):
    out = []
    for cl in sorted(all_markers["cluster"].unique(), key=lambda x: int(x)):
        sub = all_markers[all_markers["cluster"] == cl].nsmallest(n, "p_val")
        out.append(f"    Cluster {cl}: " + ", ".join(sub["gene"].tolist()))
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------

def run_full(data_dir=None, verbose=True):
    """Execute the entire advanced workflow and return all artefacts."""
    t0 = time.time()

    if verbose:
        section("1. Load PBMC 8k + QC")
    pbmc = load_object(data_dir)
    n_raw = len(pbmc)
    pbmc.misc["n_cells_raw"] = n_raw   # not recoverable after qc_filter
    pbmc = qc_filter(pbmc)
    if verbose:
        print(f"  {n_raw} cells -> {len(pbmc)} after QC "
              f"(nFeature 200-2500, percent.mt < 5)")

    if verbose:
        section("2. Standard workflow (normalize -> HVG -> PCA -> cluster -> UMAP)")
    run_pipeline(pbmc, dims=range(10), resolution=0.5)
    n_clusters = pbmc.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_clusters} global clusters at resolution 0.5")
        print(f"  Top HVGs: {pbmc.assays['RNA'].variable_features[:8]}")

    if verbose:
        section("3. Marker genes per global cluster")
    all_markers = find_all_markers(pbmc, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    if verbose:
        print(top_markers_table(all_markers, n=3))

    if verbose:
        section("4. Broad lineage annotation")
    broad = annotate_clusters(pbmc, BROAD_MARKERS)
    pbmc.meta_data["broad_cluster"] = [str(i) for i in pbmc.idents]
    pbmc.stash_ident("global_clusters")
    pbmc.rename_idents(broad)
    pbmc.meta_data["broad_celltype"] = [str(i) for i in pbmc.idents]
    if verbose:
        for c, lin in broad.items():
            print(f"    cluster {c:>2} -> {lin}")
        dist = pd.Series(list(pbmc.idents)).value_counts()
        print("\n  Lineage sizes:")
        for ct, k in dist.items():
            print(f"    {ct}: {k}")

    # ----- Subclustering the T/NK lymphoid compartment -----
    if verbose:
        section("5. Subcluster the T/NK lymphoid compartment")
    lymphoid_clusters = [c for c, lin in broad.items() if lin in LYMPHOID_LINEAGES]
    pbmc.meta_data["global_clusters"]  # ensure present
    global_idents = pbmc.meta_data["global_clusters"].astype(str).values
    cells = pbmc.cell_names()
    lymphoid_cells = [c for c, g in zip(cells, global_idents) if g in set(lymphoid_clusters)]
    pbmc.misc["n_lymphoid_clusters"] = len(lymphoid_clusters)
    if verbose:
        # map(str, ...) because these labels come back as numpy str_, whose repr
        # would render the list as [np.str_('1'), np.str_('2'), ...].
        labels = ", ".join(sorted(map(str, lymphoid_clusters), key=int))
        print(f"  Global clusters {labels} -> {len(lymphoid_cells)} T/NK cells")

    sub = pbmc.subset(cells=lymphoid_cells)
    # Re-analyse from counts; data layer is already normalised, so skip renorm.
    run_pipeline(sub, dims=range(10), resolution=0.6, n_pcs=30, normalize=False)
    n_sub = sub.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_sub} subclusters at resolution 0.6")

    if verbose:
        section("6. Annotate T/NK subclusters")
    sub_markers = find_all_markers(sub, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    sub_anno = annotate_tnk_subsets(sub)
    sub.stash_ident("sub_clusters")
    sub.rename_idents(sub_anno)
    sub.meta_data["tnk_subset"] = [str(i) for i in sub.idents]
    if verbose:
        for c, lin in sub_anno.items():
            print(f"    subcluster {c:>2} -> {lin}")
        dist = pd.Series(list(sub.idents)).value_counts()
        print("\n  Subset sizes:")
        for ct, k in dist.items():
            print(f"    {ct}: {k}")
        print("\n  Subcluster top markers:")
        print(top_markers_table(sub_markers, n=4))

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n  Global:  {pbmc}")
        print(f"\n  T/NK:    {sub}")

    return pbmc, sub, all_markers, sub_markers, broad, sub_anno


# ---------------------------------------------------------------------------
# Numeric handoff against R
# ---------------------------------------------------------------------------
#
# This tutorial's whole point is the second stage: isolate the T/NK compartment
# and re-analyse it. So the handoff is keyed by barcode rather than summarised,
# because the interesting failure mode here is silent — a subclustering fed
# from the wrong global clusters still yields a compartment, still yields
# subclusters, and a count check calls it a match. Only a barcode-level set
# comparison can tell "the same 4,000 cells" from "4,000 cells".
#
# `match_partitions` is shared with the PBMC 3k tutorial (this one is its
# advanced companion), so both tutorials score cluster agreement the same way:
# Hungarian matching on the contingency table, then ARI alongside it.

FIGURES = Path(__file__).parent / "figures_advanced"


def write_anchors(pbmc, sub, all_markers, sub_markers) -> None:
    """Dump the global and T/NK cell tables plus the scalars."""
    import json

    FIGURES.mkdir(exist_ok=True)
    md = pbmc.meta_data
    pd.DataFrame({
        "cell": list(pbmc.cell_names()),
        "nCount_RNA": md["nCount_RNA"].to_numpy(),
        "nFeature_RNA": md["nFeature_RNA"].to_numpy(),
        "percent.mt": md["percent.mt"].to_numpy(),
        "global_cluster": md["global_clusters"].astype(str).to_numpy(),
        "broad_celltype": md["broad_celltype"].astype(str).to_numpy(),
    }).to_csv(FIGURES / "py_cell_meta.csv", index=False)

    smd = sub.meta_data
    pd.DataFrame({
        "cell": list(sub.cell_names()),
        "sub_cluster": smd["sub_clusters"].astype(str).to_numpy(),
        "tnk_subset": smd["tnk_subset"].astype(str).to_numpy(),
    }).to_csv(FIGURES / "py_tnk_cells.csv", index=False)

    rna = pbmc.assays["RNA"]
    anchors = {
        "n_cells_raw": int(pbmc.misc["n_cells_raw"]),
        "n_cells_qc": len(pbmc.cell_names()),
        "n_genes": len(rna._all_feature_names),
        "n_hvg": len(rna.variable_features),
        "n_global_clusters": int(md["global_clusters"].astype(str).nunique()),
        "n_markers": int(len(all_markers)),
        "knn_nnz": int(pbmc.graphs["RNA_nn"].nnz),
        "snn_nnz": int(pbmc.graphs["RNA_snn"].nnz),
        "snn_weight_sum": float(pbmc.graphs["RNA_snn"].sum()),
        "data_sum": float(rna.layers["data"].sum()),
        "n_lymphoid_clusters": int(pbmc.misc["n_lymphoid_clusters"]),
        "n_tnk_cells": len(sub.cell_names()),
        "n_tnk_subclusters": int(smd["sub_clusters"].astype(str).nunique()),
        "n_tnk_markers": int(len(sub_markers)),
        "pc_stdev": [float(s) for s in pbmc.reductions["pca"].stdev[:10]],
    }
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    print(f"\n  Wrote py_cell_meta.csv, py_tnk_cells.csv and py_anchors.json "
          f"to {FIGURES}")
    print("  Next: Rscript tutorials/pbmc8k_subclustering_verify.R"
          "  then  python tutorials/pbmc8k_subclustering_tutorial.py --report")


def report() -> None:
    """Compare the two-stage workflow against the R reference."""
    import json

    from tutorials.pbmc3k_tutorial import match_partitions

    need = ["py_cell_meta.csv", "r_cell_meta.csv", "py_tnk_cells.csv",
            "r_tnk_cells.csv", "py_anchors.json", "r_anchors.json"]
    missing = [f for f in need if not (FIGURES / f).exists()]
    if missing:
        print(f"  missing {missing} — run the tutorial and then "
              f"`Rscript tutorials/pbmc8k_subclustering_verify.R`")
        return

    print("=" * 78)
    print("shanuz vs Seurat 5.5.1 — PBMC 8k, global clustering then T/NK subclustering")
    print("=" * 78)

    pc = pd.read_csv(FIGURES / "py_cell_meta.csv").set_index("cell")
    rc = pd.read_csv(FIGURES / "r_cell_meta.csv").set_index("cell")

    # ---- 1. did QC keep the same cells? ------------------------------------
    only_py, only_r = pc.index.difference(rc.index), rc.index.difference(pc.index)
    cells = pc.index.intersection(rc.index)
    print(f"\n  QC — cells retained: shanuz {len(pc)}, R {len(rc)}, shared {len(cells)}")
    print(f"       only shanuz {len(only_py)}   only R {len(only_r)}   "
          f"{'IDENTICAL CELL SET' if not len(only_py) and not len(only_r) else 'SETS DIFFER'}")
    p, r = pc.loc[cells], rc.loc[cells]
    print(f"\n  {'per-cell metric':<18}{'max|diff|':>14}")
    print(f"  {'-' * 32}")
    for col in ("nCount_RNA", "nFeature_RNA", "percent.mt"):
        d = np.abs(p[col].to_numpy() - r[col].to_numpy()).max()
        print(f"  {col:<18}{d:>14.3e}")

    # ---- 2. stage one: the global clusters ---------------------------------
    m = match_partitions(p["global_cluster"], r["global_cluster"])
    print(f"\n  Stage 1 — global clusters: shanuz {m['n_a']}, R {m['n_b']}")
    print(f"    adjusted Rand index {m['ari']:.4f}   "
          f"best-match concordance {m['concordance']:.4f} "
          f"({int(round(m['concordance'] * len(cells)))}/{len(cells)} cells)")
    print(f"\n    {'shanuz':>8}{'R':>6}{'n py':>8}{'n R':>8}{'shared':>9}"
          f"   lineage py / R")
    print(f"    {'-' * 62}")
    for a, b in m["mapping"].items():
        ma, mb = p["global_cluster"].astype(str) == a, r["global_cluster"].astype(str) == b
        la = p.loc[ma, "broad_celltype"].mode().iat[0]
        lb = r.loc[mb, "broad_celltype"].mode().iat[0]
        flag = "" if la == lb else "   <-- differs"
        print(f"    {a:>8}{b:>6}{int(ma.sum()):>8}{int(mb.sum()):>8}"
              f"{int(m['table'].loc[a, b]):>9}   {la} / {lb}{flag}")
    # Whichever side has more clusters leaves some unpaired. Those unpaired
    # clusters *are* the difference between the two runs, so print them rather
    # than letting the matching quietly drop them off the table.
    for label, side, other in (("shanuz", p, r), ("R", r, p)):
        paired = set(m["mapping"]) if label == "shanuz" else set(m["mapping"].values())
        for c in sorted(set(side["global_cluster"].astype(str)) - paired, key=int):
            mask = side["global_cluster"].astype(str) == c
            lands = other.loc[mask, "global_cluster"].astype(str).value_counts()
            where = ", ".join(f"{k} ({v})" for k, v in lands.head(3).items())
            print(f"    unmatched {label} cluster {c}: {int(mask.sum())} cells, "
                  f"lineage {side.loc[mask, 'broad_celltype'].mode().iat[0]}")
            print(f"      -> they sit in the other run's cluster {where}")

    same = (p["broad_celltype"].astype(str).to_numpy()
            == r["broad_celltype"].astype(str).to_numpy())
    print(f"\n    Broad lineage label, per cell: {same.mean():.4f} "
          f"({same.sum()}/{len(cells)})")

    # ---- 3. the compartment handed to stage two ----------------------------
    # The load-bearing check. Everything downstream is conditioned on this set,
    # so two compartments of the same size drawn from different clusters would
    # make every later number incomparable while looking fine.
    pt = pd.read_csv(FIGURES / "py_tnk_cells.csv").set_index("cell")
    rt = pd.read_csv(FIGURES / "r_tnk_cells.csv").set_index("cell")
    sp_, sr = set(pt.index), set(rt.index)
    inter, union = sp_ & sr, sp_ | sr
    print(f"\n  The T/NK compartment — shanuz {len(sp_)} cells, R {len(sr)}")
    print(f"    shared {len(inter)}   only shanuz {len(sp_ - sr)}   "
          f"only R {len(sr - sp_)}   Jaccard {len(inter) / max(len(union), 1):.4f}")

    # ---- 4. stage two: the subclusters, on the shared cells ----------------
    shared = sorted(inter)
    if len(shared) >= 2:
        ms = match_partitions(pt.loc[shared, "sub_cluster"], rt.loc[shared, "sub_cluster"])
        print(f"\n  Stage 2 — T/NK subclusters on the {len(shared)} shared cells: "
              f"shanuz {ms['n_a']}, R {ms['n_b']}")
        print(f"    adjusted Rand index {ms['ari']:.4f}   "
              f"best-match concordance {ms['concordance']:.4f}")
        sub_same = (pt.loc[shared, "tnk_subset"].astype(str).to_numpy()
                    == rt.loc[shared, "tnk_subset"].astype(str).to_numpy())
        print(f"    T/NK subset label, per cell: {sub_same.mean():.4f} "
              f"({sub_same.sum()}/{len(shared)})")
        py_n = pt.loc[shared, "tnk_subset"].value_counts()
        r_n = rt.loc[shared, "tnk_subset"].value_counts()
        print(f"\n    {'subset':<14}{'shanuz':>9}{'R':>9}")
        print(f"    {'-' * 32}")
        for s in sorted(set(py_n.index) | set(r_n.index)):
            print(f"    {s:<14}{py_n.get(s, 0):>9}{r_n.get(s, 0):>9}")

    # ---- 5. scalars ---------------------------------------------------------
    pa, ra = (json.loads((FIGURES / f"{s}_anchors.json").read_text())
              for s in ("py", "r"))
    print(f"\n  {'anchor':<22}{'shanuz':>20}{'R Seurat':>20}   verdict")
    print(f"  {'-' * 68}")
    for k in ("n_cells_raw", "n_cells_qc", "n_genes", "n_hvg", "n_global_clusters",
              "n_markers", "knn_nnz", "snn_nnz", "snn_weight_sum", "data_sum",
              "n_lymphoid_clusters", "n_tnk_cells", "n_tnk_subclusters",
              "n_tnk_markers"):
        if k not in pa or k not in ra:
            continue
        x, y = pa[k], ra[k]
        if isinstance(x, int) and isinstance(y, int):
            verdict = "MATCH" if x == y else f"differ by {x - y:+d}"
            xs, ys = str(x), str(y)
        else:
            verdict = f"rel {abs(x - y) / max(abs(y), 1e-12):.2e}"
            xs, ys = f"{x:.6f}", f"{y:.6f}"
        print(f"  {k:<22}{xs:>20}{ys:>20}   {verdict}")
    sd = np.abs(np.array(pa["pc_stdev"]) - np.array(ra["pc_stdev"]))
    print(f"  {'pc_stdev[1:10]':<22}{'':>40}   max|Δ| {sd.max():.3e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 8k subclustering tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--report", action="store_true",
                        help="compare against the R reference and exit")
    args = parser.parse_args()
    if args.report:
        report()
    else:
        pbmc, sub, all_markers, sub_markers, _, _ = run_full(data_dir=args.data_dir)
        write_anchors(pbmc, sub, all_markers, sub_markers)
