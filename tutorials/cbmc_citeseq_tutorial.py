"""Multimodal CITE-seq Tutorial — RNA + surface protein (ADT) with Shanuz.

A Python port of Seurat's multimodal vignette
(https://satijalab.org/seurat/articles/multimodal_vignette) using the CBMC
CITE-seq dataset (GSE100866): ~8,600 cord-blood mononuclear cells measured for
both the transcriptome and 13 surface proteins.

It demonstrates Shanuz's multi-assay support:
  * build the object from RNA and run the standard clustering workflow,
  * attach the antibody-capture counts as a second ("ADT") assay,
  * CLR-normalise the proteins (margin=2, per-protein across cells), and
  * visualise protein levels on the RNA-derived UMAP, comparing each protein
    to its encoding gene.

Usage
-----
    python tutorials/cbmc_citeseq_tutorial.py [--data-dir PATH]

The CBMC dataset (~15 MB) downloads automatically to ~/.shanuz_data/cbmc.

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

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import cbmc_citeseq
from shanuz.shanuz import create_shanuz_object
from shanuz.assay5 import create_assay5_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data,
)
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_all_markers
from shanuz.plotting import _get_expression


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

    obj = create_shanuz_object(
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
    # CLR across cells per protein (Seurat's recommended ADT margin).
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
    (CD16 & CD56 high, CD3-), B (CD19), monocytes (CD14), DC (CD11c),
    progenitors (CD34). Populations outside the 13-protein panel (platelets,
    erythroid, pDC, cycling) carry no ADT signal, so they are resolved from RNA
    markers — the same protein+RNA reasoning the Seurat vignette uses by eye.
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
        cd3, cd4, cd8 = pm("CD3", mask), pm("CD4", mask), pm("CD8", mask)
        cd19, cd14 = pm("CD19", mask), pm("CD14", mask)
        cd16, cd56, cd11c, cd34 = pm("CD16", mask), pm("CD56", mask), pm("CD11c", mask), pm("CD34", mask)
        # Unambiguous RNA-only lineages (no protein in the panel) take priority.
        if rmean(["PPBP", "PF4"], mask) > 2.0:
            assignment[c] = "Platelet"
        elif rmean(["HBB", "HBA1"], mask) > 2.5:
            assignment[c] = "Erythroid"
        elif cd3 > 0.5:
            assignment[c] = "CD8 T" if cd8 > cd4 else "CD4 T"
        elif cd16 > 0.8 and cd56 > 0.8 and cd3 < 0.5:
            assignment[c] = "NK"
        elif cd19 > 1.0 and cd3 < 0.5:
            assignment[c] = "B"
        elif cd14 > 0.5:
            assignment[c] = "CD14+ Mono"
        elif cd11c > 0.8 and cd3 < 0.5:
            assignment[c] = "DC / Mono"
        elif cd34 > 0.8:
            assignment[c] = "Progenitor"
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
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n{obj}")

    return obj, all_markers, anno


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CBMC CITE-seq multimodal tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
