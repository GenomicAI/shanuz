"""Multimodal CITE-seq Tutorial — RNA + surface protein (ADT) with Truecell.

A Python port of Seurat's multimodal vignette
(https://satijalab.org/seurat/articles/multimodal_vignette) using the CBMC
CITE-seq dataset (GSE100866): ~8,600 cord-blood mononuclear cells measured for
both the transcriptome and 13 surface proteins.

It demonstrates Truecell's multi-assay support:
  * build the object from RNA and run the standard clustering workflow,
  * attach the antibody-capture counts as a second ("ADT") assay,
  * CLR-normalise the proteins (margin=2, per-cell across the panel),
  * visualise protein levels on the RNA-derived UMAP, comparing each protein
    to its encoding gene, and
  * jointly cluster both modalities with Weighted Nearest Neighbor (WNN)
    analysis (learns per-cell RNA-vs-protein weights).

Usage
-----
    python tutorials/cbmc_citeseq_tutorial.py [--data-dir PATH]

The CBMC dataset (~15 MB) downloads automatically to ~/.truecell_data/cbmc.

References
----------
Stoeckius M, Hafemeister C, Stephenson W, et al. (2017)
**Simultaneous epitope and transcriptome measurement in single cells.**
Nature Methods 14, 865-868. https://doi.org/10.1038/nmeth.4380
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from truecell.datasets import cbmc_citeseq
from truecell.truecell import create_truecell_object
from truecell.assay5 import create_assay5_object
from truecell.preprocessing import (
    normalize_data, find_variable_features, scale_data,
)
from truecell.reduction import run_pca
from truecell.neighbors import find_neighbors
from truecell.multimodal import find_multi_modal_neighbors
from truecell.clustering import find_clusters
from truecell.umap import run_umap
from truecell.markers import find_all_markers
from truecell.plotting import _get_expression


# Surface proteins in the CBMC panel, mapped to their encoding gene(s) for the
# protein-vs-RNA comparison plots.
PROTEIN_TO_GENE = {
    "CD3": "CD3E", "CD4": "CD4", "CD8": "CD8A", "CD19": "CD19",
    "CD14": "CD14", "CD16": "FCGR3A", "CD56": "NCAM1", "CD11c": "ITGAX",
    "CD34": "CD34", "CD45RA": "PTPRC", "CD10": "MME", "CCR7": "CCR7",
}


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_object(data_dir=None, clr_margin=2):
    """Load CBMC, build the RNA object, and attach the CLR-normalised ADT assay."""
    rna, genes, adt, proteins, cells = cbmc_citeseq(data_dir=data_dir)

    obj = create_truecell_object(
        counts=rna, assay="RNA", min_cells=3, min_features=0,
        project="cbmc", feature_names=genes, cell_names=cells,
    )
    kept = obj.cell_names()

    # Attach the antibody-capture counts as a second assay, aligned to `kept`.
    cpos = {c: i for i, c in enumerate(cells)}
    adt_aligned = adt[:, [cpos[c] for c in kept]].tocsc()
    obj.assays["ADT"] = create_assay5_object(
        counts=adt_aligned, feature_names=proteins, cell_names=kept, key="adt_",
    )
    # CLR per cell across the protein panel (Seurat's recommended ADT margin).
    normalize_data(obj, assay="ADT", normalization_method="CLR", margin=clr_margin)
    return obj


def run_rna_workflow(obj, dims=range(15), resolution=0.6):
    normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(obj, selection_method="vst", nfeatures=2000)
    scale_data(obj, features=obj.assays["RNA"]._all_feature_names)
    run_pca(obj, n_pcs=30, features=obj.assays["RNA"].variable_features)
    find_neighbors(obj, dims=dims, k_param=20)
    find_clusters(obj, resolution=resolution, algorithm=1, random_seed=0)
    run_umap(obj, dims=dims, seed=42)
    return obj


def run_wnn(obj, rna_dims=range(15), resolution=0.6):
    """Weighted Nearest Neighbor multimodal clustering (Hao et al., Cell 2021).

    Where the RNA workflow clusters on transcriptome alone, WNN learns a
    per-cell weight for each modality — how much to trust RNA vs protein for
    that particular cell — and clusters/embeds on a *joint* graph. Cells whose
    lineage is sharper in protein space (e.g. the CD4/CD8 T split) lean on ADT;
    cells the 13-protein panel can't resolve lean on RNA.

    Mirrors Seurat's ADT `RunPCA` (reduction "apca") →
    `FindMultiModalNeighbors` → `FindClusters(graph = "wsnn")` →
    `RunUMAP(nn.name = "weighted.nn")`. Assumes `run_rna_workflow` has already
    populated the RNA `"pca"` reduction and the CLR-normalised ADT assay.
    """
    # Protein modality needs its own reduction. The ADT panel is small (13
    # proteins), so PCA keeps every informative component.
    adt_features = obj.assays["ADT"]._all_feature_names
    scale_data(obj, assay="ADT", features=adt_features)
    n_apca = min(18, len(adt_features) - 1)
    run_pca(obj, assay="ADT", reduction_name="apca", reduction_key="apca_",
            n_pcs=n_apca, features=adt_features)

    # Learn per-cell modality weights and build the joint wknn/wsnn graphs.
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[rna_dims, range(n_apca)], k_nn=20,
    )

    # Cluster and embed on the joint graph (not the RNA reduction).
    find_clusters(obj, resolution=resolution, graph_name="wsnn", random_seed=0)
    obj.meta_data["wnn_clusters"] = obj.meta_data["seurat_clusters"].astype(str).tolist()
    run_umap(obj, graph="wsnn", reduction_name="wnn_umap", seed=42)
    return obj


# RNA markers for the populations the 13-protein ADT panel can't resolve.
_RNA_FALLBACK = {
    "Platelet":  ["PPBP", "PF4"],
    "Erythroid": ["HBB", "HBA1"],
    "pDC":       ["IGJ", "PLD4", "SERPINF1"],
    "Cycling":   ["STMN1", "MKI67", "TUBB"],
}


def annotate_cells(obj):
    """Annotate RNA clusters using surface protein first, RNA as a fallback.

    CITE-seq proteins are cleaner lineage markers, so immune lineages are gated
    on the ADT assay in priority order — T (CD3) split into CD4/CD8, then NK
    (CD16 & CD56 high), B (CD19), progenitors (CD34), monocytes (CD14), DC
    (CD11c). Populations outside the 13-protein panel (platelets, erythroid,
    pDC, cycling) carry no ADT signal, so they are resolved from RNA markers —
    the same protein+RNA reasoning the Seurat vignette uses by eye.

    The CLR cut-offs match cbmc_citeseq_verify.R exactly: both sides run the
    same transform, so the same thresholds read the same lineages.
    """
    idents = np.array([str(i) for i in obj.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    prot = {p: _get_expression(obj, p, assay="ADT") for p in obj.assays["ADT"]._all_feature_names}
    rna_feats = set(obj.assays["RNA"]._all_feature_names)
    rna = {g: _get_expression(obj, g, assay="RNA")
           for gs in _RNA_FALLBACK.values() for g in gs if g in rna_feats}

    def pm(p, mask):
        return float(prot[p][mask].mean()) if p in prot else -np.inf

    def rna_fallback(mask):
        best, best_score = "Other", 0.30
        for label, genes in _RNA_FALLBACK.items():
            present = [g for g in genes if g in rna]
            if not present:
                continue
            score = float(np.mean([rna[g][mask].mean() for g in present]))
            if score > best_score:
                best_score, best = score, label
        return best

    def rmean(genes, mask):
        present = [g for g in genes if g in rna]
        return float(np.mean([rna[g][mask].mean() for g in present])) if present else 0.0

    assignment = {}
    for c in clusters:
        mask = idents == c
        cd3, cd8 = pm("CD3", mask), pm("CD8", mask)
        cd19, cd14 = pm("CD19", mask), pm("CD14", mask)
        cd16, cd56, cd11c, cd34 = pm("CD16", mask), pm("CD56", mask), pm("CD11c", mask), pm("CD34", mask)
        # Unambiguous RNA-only lineages (no protein in the panel) take priority.
        if rmean(["PPBP", "PF4"], mask) > 2.0:
            assignment[c] = "Platelet"
        elif rmean(["HBB", "HBA1"], mask) > 2.5:
            assignment[c] = "Erythroid"
        elif cd3 > 1.0:
            assignment[c] = "CD8 T" if cd8 > 1.0 else "CD4 T"
        elif cd16 > 0.8 and cd56 > 0.8:
            assignment[c] = "NK"
        elif cd19 > 1.5:
            assignment[c] = "B"
        elif cd34 > 1.0:
            assignment[c] = "Progenitor"
        elif cd14 > 0.9:
            assignment[c] = "CD14+ Mono"
        elif cd11c > 1.0:
            assignment[c] = "DC / Mono"
        else:
            assignment[c] = rna_fallback(mask)
    return assignment


# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 64}\n  {title}\n{'=' * 64}")


def run_full(data_dir=None, verbose=True):
    t0 = time.time()

    if verbose:
        section("1. Load CBMC CITE-seq (RNA + ADT)")
    obj = load_object(data_dir)
    if verbose:
        print(f"  {len(obj)} cells | RNA {len(obj.assays['RNA']._all_feature_names)} genes "
              f"| ADT {len(obj.assays['ADT']._all_feature_names)} proteins")
        print(f"  proteins: {obj.assays['ADT']._all_feature_names}")

    if verbose:
        section("2. RNA workflow (normalize -> HVG -> PCA -> cluster -> UMAP)")
    run_rna_workflow(obj)
    n_clusters = obj.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_clusters} RNA clusters at resolution 0.6")

    if verbose:
        section("3. Surface-protein levels per cluster (ADT, CLR)")
    idents = np.array([str(i) for i in obj.idents])
    cl = sorted(set(idents), key=int)
    rows = []
    panel = ["CD3", "CD4", "CD8", "CD19", "CD14", "CD16", "CD56", "CD11c", "CD34"]
    for p in panel:
        e = _get_expression(obj, p, assay="ADT")
        rows.append([f"{e[idents == c].mean():+.1f}" for c in cl])
    if verbose:
        print(pd.DataFrame(rows, index=panel, columns=[f"c{c}" for c in cl]).to_string())

    if verbose:
        section("4. Annotate clusters by surface protein")
    anno = annotate_cells(obj)
    obj.stash_ident("rna_clusters")
    obj.rename_idents(anno)
    obj.meta_data["protein_celltype"] = [str(i) for i in obj.idents]
    if verbose:
        for c, lab in anno.items():
            print(f"    cluster {c:>2} -> {lab}")
        dist = pd.Series(list(obj.idents)).value_counts()
        print("\n  Cell-type sizes:")
        for ct, k in dist.items():
            print(f"    {ct}: {k}")

    if verbose:
        section("5. RNA markers per cluster (sanity check)")
    obj.idents = obj.meta_data["rna_clusters"].astype(str).tolist()
    all_markers = find_all_markers(obj, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
    obj.rename_idents(anno)  # restore cell-type labels
    if verbose:
        for clid in sorted(all_markers["cluster"].unique(), key=int):
            top = all_markers[all_markers["cluster"] == clid].nsmallest(3, "p_val")
            print(f"    cluster {clid}: " + ", ".join(top["gene"].tolist()))

    if verbose:
        section("6. Weighted Nearest Neighbor (WNN) multimodal clustering")
    run_wnn(obj)
    if verbose:
        n_wnn = obj.meta_data["wnn_clusters"].nunique()
        rna_w = obj.meta_data["RNA.weight"]
        adt_w = obj.meta_data["ADT.weight"]
        print(f"  {n_wnn} WNN clusters (RNA-only gave {n_clusters}) on the joint graph")
        print(f"  mean per-cell modality weight — RNA {rna_w.mean():.2f} | ADT {adt_w.mean():.2f}")
        # Which annotated cell types lean on protein vs RNA? Group by the
        # protein cell-type labels stashed in step 4 (find_clusters has since
        # set the active ident to the new WNN cluster ids).
        print("\n  Mean ADT weight by protein cell type (higher = protein-driven):")
        by_ct = adt_w.groupby(obj.meta_data["protein_celltype"]).mean().sort_values(ascending=False)
        for cell_type, w in by_ct.items():
            print(f"    {cell_type:<14} ADT {w:.2f}")

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n{obj}")

    return obj, all_markers, anno


# ---------------------------------------------------------------------------
# Numeric handoff against R
# ---------------------------------------------------------------------------
#
# The per-cell weights are written keyed by *barcode*, not summarised by cell
# type. Mean ADT weight per cell type — the table printed above — asks two
# different cluster partitions the same question, so it cannot separate "the
# weights differ" from "the labels differ". Same-barcode weights can, and this
# distinction is not hypothetical: the per-cell-type table showed Progenitor at
# 0.35 against Seurat's 0.285 while every other type agreed to 0.02.

FIGURES = Path(__file__).parent / "figures_multimodal"


def write_anchors(obj):
    """Dump the CLR summary, the per-cell weights and the scalars."""
    import json

    FIGURES.mkdir(exist_ok=True)
    adt = obj.assays["ADT"]
    data = adt.layers["data"]
    dense = data.toarray() if sp.issparse(data) else np.asarray(data)
    pd.DataFrame({
        "protein": list(adt._all_feature_names),
        "mean": dense.mean(axis=1),
        "sd": dense.std(axis=1, ddof=1),
        "min": dense.min(axis=1),
        "max": dense.max(axis=1),
    }).to_csv(FIGURES / "py_adt_clr.csv", index=False)

    md = obj.meta_data
    pd.DataFrame({
        "cell": list(obj.cell_names()),
        "RNA.weight": md["RNA.weight"].to_numpy(),
        "ADT.weight": md["ADT.weight"].to_numpy(),
        "rna_cluster": md["rna_clusters"].astype(str).to_numpy(),
        "protein_celltype": md["protein_celltype"].astype(str).to_numpy(),
    }).to_csv(FIGURES / "py_cell_weights.csv", index=False)

    anchors = {
        "n_cells": len(obj),
        "n_genes": len(obj.assays["RNA"]._all_feature_names),
        "n_proteins": len(adt._all_feature_names),
        "n_rna_clusters": int(md["rna_clusters"].nunique()),
        "n_wnn_clusters": int(md["wnn_clusters"].nunique()),
        "mean_rna_weight": float(md["RNA.weight"].mean()),
        "mean_adt_weight": float(md["ADT.weight"].mean()),
        "adt_weight_sum": float(md["ADT.weight"].sum()),
    }
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    print(f"\n  Wrote py_adt_clr.csv, py_cell_weights.csv and py_anchors.json "
          f"to {FIGURES}")
    print("  Next: Rscript tutorials/cbmc_citeseq_verify.R"
          "  then  python tutorials/cbmc_citeseq_tutorial.py --report")


def report():
    """Compare the CLR transform, the WNN weights and the annotation."""
    import json

    from scipy.stats import pearsonr, spearmanr

    need = ["py_adt_clr.csv", "r_adt_clr.csv", "py_cell_weights.csv",
            "r_cell_weights.csv", "py_anchors.json", "r_anchors.json"]
    missing = [f for f in need if not (FIGURES / f).exists()]
    if missing:
        print(f"  missing {missing} — run the tutorial and then "
              f"`Rscript tutorials/cbmc_citeseq_verify.R`")
        return

    print("=" * 74)
    print("truecell vs Seurat 5.5.1 — CITE-seq: CLR, WNN weights, annotation")
    print("=" * 74)

    # ---- 1. the CLR transform, which depends on nothing stochastic ----------
    py = pd.read_csv(FIGURES / "py_adt_clr.csv").set_index("protein")
    r = pd.read_csv(FIGURES / "r_adt_clr.csv").set_index("protein")
    prot = py.index.intersection(r.index)
    print(f"\n  ADT CLR normalisation — {len(prot)} proteins, no RNG anywhere")
    print(f"  {'statistic':<12}{'max|diff|':>14}")
    print(f"  {'-' * 26}")
    for col in ("mean", "sd", "min", "max"):
        d = np.abs(py.loc[prot, col].to_numpy() - r.loc[prot, col].to_numpy()).max()
        print(f"  {col:<12}{d:>14.3e}")

    # ---- 2. the WNN weights, per cell, on shared barcodes -------------------
    pw = pd.read_csv(FIGURES / "py_cell_weights.csv").set_index("cell")
    rw = pd.read_csv(FIGURES / "r_cell_weights.csv").set_index("cell")
    cells = pw.index.intersection(rw.index)
    pw, rw = pw.loc[cells], rw.loc[cells]
    a, b = pw["ADT.weight"].to_numpy(), rw["ADT.weight"].to_numpy()
    print(f"\n  WNN modality weights — {len(cells)} shared barcodes")
    print(f"    ADT weight  pearson {pearsonr(a, b).statistic:.4f}   "
          f"spearman {spearmanr(a, b).statistic:.4f}")
    print(f"    mean        truecell {a.mean():.4f}   R {b.mean():.4f}   "
          f"diff {abs(a.mean() - b.mean()):.4f}")
    print(f"    max|diff| per cell {np.abs(a - b).max():.4f}   "
          f"median|diff| {np.median(np.abs(a - b)):.4f}")

    # ---- 3. is a per-cell-type gap the weights or the labels? ---------------
    # Score each tool's weights under *both* labellings. If the gap follows the
    # labels rather than the weights, it was never a WNN difference.
    print("\n  Mean ADT weight by cell type — under each tool's own labels,")
    print("  then truecell's weights re-grouped by R's labels. A row that moves")
    print("  in the third column was a labelling difference, not a weight one.")
    print(f"\n  {'cell type':<14}{'truecell':>9}{'R':>9}{'truecell w/ R labels':>21}"
          f"{'n py':>8}{'n R':>8}")
    print(f"  {'-' * 69}")
    py_by = pw.groupby("protein_celltype")["ADT.weight"].mean()
    r_by = rw.groupby("protein_celltype")["ADT.weight"].mean()
    cross = pw["ADT.weight"].groupby(rw["protein_celltype"]).mean()
    py_n = pw["protein_celltype"].value_counts()
    r_n = rw["protein_celltype"].value_counts()
    for ct in sorted(set(py_by.index) | set(r_by.index)):
        def f(s, w=9, fmt=".3f"):
            return f"{s[ct]:{fmt}}".rjust(w) if ct in s.index else "—".rjust(w)
        print(f"  {ct:<14}{f(py_by)}{f(r_by)}{f(cross, 21)}"
              f"{f(py_n, 8, 'd')}{f(r_n, 8, 'd')}")

    same = (pw["protein_celltype"].to_numpy() == rw["protein_celltype"].to_numpy())
    print(f"\n  Per-cell cell-type concordance: {same.mean():.4f} "
          f"({same.sum()}/{len(cells)})")

    # ---- 4. scalars ---------------------------------------------------------
    pa = json.loads((FIGURES / "py_anchors.json").read_text())
    ra = json.loads((FIGURES / "r_anchors.json").read_text())
    # Counts must be identical; the aggregate weights are means over a WNN fit
    # that neither tool reproduces bit-for-bit, so they get a relative
    # difference rather than a MATCH/differ verdict that would always read
    # "differ" and tell you nothing.
    print(f"\n  {'anchor':<20}{'truecell':>20}{'R Seurat':>20}   verdict")
    print(f"  {'-' * 68}")
    for k in sorted(set(pa) & set(ra)):
        x, y = pa[k], ra[k]
        if isinstance(x, int) and isinstance(y, int):
            verdict = "MATCH" if x == y else "differ"
        else:
            verdict = f"rel {abs(x - y) / max(abs(y), 1e-12):.2e}"
        print(f"  {k:<20}{x!s:>20}{y!s:>20}   {verdict}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CBMC CITE-seq multimodal tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--report", action="store_true",
                        help="compare against the R reference and exit")
    args = parser.parse_args()
    if args.report:
        report()
    else:
        obj, _, _ = run_full(data_dir=args.data_dir)
        write_anchors(obj)
