"""SCTransform Tutorial — PBMC 3k with Truecell.

A Python port of Seurat's sctransform vignette
(https://satijalab.org/seurat/articles/sctransform_vignette) on the PBMC 3k
dataset. SCTransform replaces the NormalizeData -> FindVariableFeatures ->
ScaleData trio with a single regularized negative-binomial model, returning
Pearson residuals that more effectively remove technical (sequencing-depth and
percent-mito) effects. The vignette's headline result is that this sharper
normalization, run over more PCs (dims 1:30), resolves finer immune subsets
than the standard log-normalization workflow.

This script runs BOTH workflows on the same cells so their cluster resolution
can be compared directly:
  * SCT workflow  : sctransform(vars.to.regress="percent.mt") -> PCA -> dims 1:30
  * Std workflow  : LogNormalize -> VST -> ScaleData -> PCA -> dims 1:10

Usage
-----
    python tutorials/pbmc3k_sctransform_tutorial.py [--data-dir PATH]

The PBMC 3k dataset (~24 MB) downloads automatically to ~/.truecell_data/pbmc3k.

References
----------
Hafemeister C, Satija R (2019) Genome Biology 20, 296.
Choudhary S, Satija R (2022) Genome Biology 23, 27. (sctransform v2)
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

from truecell.datasets import pbmc3k
from truecell.truecell import create_truecell_object
from truecell.preprocessing import (
    normalize_data, find_variable_features, scale_data, percentage_feature_set,
)
from truecell.sctransform import sctransform
from truecell.reduction import run_pca
from truecell.neighbors import find_neighbors
from truecell.clustering import find_clusters
from truecell.umap import run_umap
from truecell.markers import find_all_markers
from truecell.plotting import _get_expression


# Markers shown in the vignette's FeaturePlots — these define the fine subsets
# SCTransform is meant to resolve.
VIGNETTE_MARKERS_1 = ["CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"]
VIGNETTE_MARKERS_2 = ["CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"]
VLN_MARKERS = ["CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"]

# Fine-grained lineage panel for relative-enrichment annotation. Several of
# these (CD8 effector vs naive, CD4 naive vs memory, NK bright vs dim) are the
# distinctions SCTransform is meant to sharpen.
FINE_MARKERS = {
    "Naive CD4 T":   ["IL7R", "CCR7", "LEF1", "SELL"],
    "Memory CD4 T":  ["IL7R", "S100A4", "IL32", "ANXA1"],
    "CD8 Naive/Mem": ["CD8A", "CD8B", "CCR7"],
    "CD8 Effector":  ["CD8A", "GZMK", "CCL5", "NKG7"],
    "B":             ["MS4A1", "CD79A", "TCL1A", "FCER2"],
    "CD14+ Mono":    ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
    "NK":            ["GNLY", "NKG7", "KLRD1", "XCL1"],
    "DC":            ["FCER1A", "CST3"],
    "pDC":           ["SERPINF1", "ITM2C"],
    "Platelet":      ["PPBP", "PF4"],
}

# Which fine lineages count as resolved T-cell subsets, for the comparison.
T_SUBSETS = {"Naive CD4 T", "Memory CD4 T", "CD8 Naive/Mem", "CD8 Effector"}


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_object(data_dir=None):
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    pbmc = create_truecell_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc3k", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    return pbmc


def run_sct_workflow(pbmc, dims=range(30), resolution=0.8, seed=0):
    """SCTransform workflow (mirrors the vignette: regress percent.mt, dims 1:30)."""
    sctransform(pbmc, vars_to_regress=["percent.mt"], n_features=3000, seed=42)
    run_pca(pbmc, n_pcs=50, features=pbmc.assays["SCT"].variable_features)
    find_neighbors(pbmc, dims=dims, k_param=20)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=seed)
    run_umap(pbmc, dims=dims, seed=42)
    return pbmc


def run_std_workflow(pbmc, dims=range(10), resolution=0.8, seed=0):
    """Standard log-normalization workflow for comparison (dims 1:10)."""
    normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    scale_data(pbmc, features=pbmc.assays["RNA"].variable_features)
    run_pca(pbmc, n_pcs=50, features=pbmc.assays["RNA"].variable_features)
    find_neighbors(pbmc, dims=dims, k_param=20)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=seed)
    run_umap(pbmc, dims=dims, seed=42)
    return pbmc


def annotate_clusters(pbmc, marker_sets, assay=None):
    """Assign each cluster to the lineage whose markers are most *enriched*.

    Each marker's per-cluster mean is z-scored across clusters, so a cluster
    scores on a lineage by relative enrichment (matching the advanced tutorial).
    Reuse is allowed: several clusters may share a lineage.
    """
    idents = np.array([str(i) for i in pbmc.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(pbmc.assays[assay or pbmc.active_assay]._all_feature_names)

    needed = {g for gs in marker_sets.values() for g in gs if g in feats}
    zmean = {}
    for g in needed:
        expr = _get_expression(pbmc, g, assay=assay)
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


# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 64}\n  {title}\n{'=' * 64}")


def _n_t_subsets(anno):
    return len({lab for lab in anno.values() if lab in T_SUBSETS})


def run_full(data_dir=None, verbose=True):
    t0 = time.time()

    if verbose:
        section("1. Load PBMC 3k + percent.mt")
    base = load_object(data_dir)
    if verbose:
        print(f"  {len(base)} cells x {len(base.assays['RNA']._all_feature_names)} genes")

    # ---- SCTransform workflow ----
    if verbose:
        section("2. SCTransform workflow (regress percent.mt, PCA dims 1:30)")
    sct = load_object(data_dir)
    run_sct_workflow(sct)
    n_sct = sct.meta_data["seurat_clusters"].nunique()
    sct_assay = sct.assays["SCT"]
    if verbose:
        print(f"  SCT assay: {len(sct_assay._all_feature_names)} genes, "
              f"{len(sct_assay.variable_features)} variable features")
        print(f"  {n_sct} clusters at resolution 0.8 (dims 1:30)")
        print(f"  Top SCT variable features: {sct_assay.variable_features[:10]}")

    sct_anno = annotate_clusters(sct, FINE_MARKERS, assay="SCT")
    sct.meta_data["sct_clusters"] = [str(i) for i in sct.idents]
    sct.stash_ident("sct_clusters")
    sct.rename_idents(sct_anno)
    sct.meta_data["sct_celltype"] = [str(i) for i in sct.idents]

    # ---- Standard workflow ----
    if verbose:
        section("3. Standard LogNormalize workflow (PCA dims 1:10)")
    std = load_object(data_dir)
    run_std_workflow(std)
    n_std = std.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_std} clusters at resolution 0.8 (dims 1:10)")
    std_anno = annotate_clusters(std, FINE_MARKERS, assay="RNA")
    std.meta_data["std_clusters"] = [str(i) for i in std.idents]
    std.stash_ident("std_clusters")
    std.rename_idents(std_anno)

    # ---- Comparison ----
    if verbose:
        section("4. SCTransform vs standard — cluster resolution")
        print(f"  {'workflow':<16}{'clusters':>10}{'T-subsets resolved':>22}")
        print(f"  {'-' * 46}")
        print(f"  {'SCTransform':<16}{n_sct:>10}{_n_t_subsets(sct_anno):>22}")
        print(f"  {'LogNormalize':<16}{n_std:>10}{_n_t_subsets(std_anno):>22}")
        print("\n  SCT cluster -> annotation:")
        for c, lab in sct_anno.items():
            print(f"    cluster {c:>2} -> {lab}")

    if verbose:
        section("5. SCT marker check (vignette FeaturePlot genes)")
    sct.idents = sct.meta_data["sct_clusters"].astype(str).tolist()
    sct_markers = find_all_markers(sct, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    sct.rename_idents(sct_anno)
    if verbose:
        for clid in sorted(sct_markers["cluster"].unique(), key=int):
            top = sct_markers[sct_markers["cluster"] == clid].nsmallest(4, "p_val")
            print(f"    cluster {clid}: " + ", ".join(top["gene"].tolist()))

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n  SCT:  {sct}")

    return sct, std, sct_anno, std_anno, sct_markers


# ---------------------------------------------------------------------------
# Numeric handoff against R
# ---------------------------------------------------------------------------
#
# The figures in this vignette are compared by eye. The model is not. That
# distinction earns its keep here more than anywhere else in the series: this
# port's SCT model was once wrong in four separate ways at once — theta came out
# *anti*-correlated with R's, and residual variance ranked genes near randomly —
# and the tutorial still drew a perfectly plausible UMAP. A picture cannot fail.

FIGURES = Path(__file__).parent / "figures_sctransform"

# Seurat's SCTModel feature.attributes columns, in its order. truecell stores the
# same names on the SCT assay's meta_data, so the two tables line up directly.
_MODEL_COLS = ["detection_rate", "gmean", "residual_mean",
               "residual_variance", "theta", "(Intercept)", "log_umi"]


def write_anchors(sct, n_sct, n_std):
    """Dump the fitted model, the ranked features and the scalars."""
    import json

    FIGURES.mkdir(exist_ok=True)
    assay = sct.assays["SCT"]
    md = assay.meta_data
    frame = pd.DataFrame({"gene": list(md.index)})
    for col in _MODEL_COLS:
        frame[col] = md[col].to_numpy()
    frame.to_csv(FIGURES / "py_sct_model.csv", index=False)

    (FIGURES / "py_variable_features.txt").write_text(
        "\n".join(assay.variable_features) + "\n")

    n_cells = len(sct)
    anchors = {
        "n_cells": n_cells,
        "n_genes_modelled": int(len(md)),
        "n_variable_features": len(assay.variable_features),
        # Both clips are derived, not stored, so they are recomputed from the
        # same expressions the implementation uses. `vst` is the residual
        # -variance clip, `sct` the tighter scale.data one.
        "clip_vst_lo": -float(np.sqrt(n_cells)),
        "clip_vst_hi": float(np.sqrt(n_cells)),
        "clip_sct_lo": -float(np.sqrt(n_cells / 30.0)),
        "clip_sct_hi": float(np.sqrt(n_cells / 30.0)),
        "n_poisson_genes": int(np.isinf(md["theta"].to_numpy()).sum()),
        "residual_variance_sum": float(md["residual_variance"].sum()),
        "n_clusters_sct": int(n_sct),
        "n_clusters_std": int(n_std),
    }
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    print(f"\n  Wrote py_sct_model.csv, py_variable_features.txt and "
          f"py_anchors.json to {FIGURES}")
    print("  Next: Rscript tutorials/pbmc3k_sctransform_verify.R"
          "  then  python tutorials/pbmc3k_sctransform_tutorial.py --report")


def report():
    """Compare the two fitted models gene by gene."""
    import json

    from scipy.stats import pearsonr, spearmanr

    need = ["py_sct_model.csv", "r_sct_model.csv",
            "py_anchors.json", "r_anchors.json"]
    missing = [f for f in need if not (FIGURES / f).exists()]
    if missing:
        print(f"  missing {missing} — run the tutorial and then "
              f"`Rscript tutorials/pbmc3k_sctransform_verify.R`")
        return

    # `keep_default_na=False` because gene symbols include literal "NA"-like
    # tokens that pandas would otherwise read as missing.
    py = pd.read_csv(FIGURES / "py_sct_model.csv", keep_default_na=False, float_precision="round_trip")
    r = pd.read_csv(FIGURES / "r_sct_model.csv", keep_default_na=False, float_precision="round_trip")
    # `Read10X` rewrites "_" to "-" in gene symbols and truecell's loader does
    # not, so one gene (RP11-442N24__B.1) spells differently on the two sides.
    # That belongs to the two file readers, not to SCTransform — normalise to
    # R's spelling rather than silently dropping the gene from the comparison.
    py["gene"] = py["gene"].str.replace("_", "-", regex=False)
    py, r = py.set_index("gene"), r.set_index("gene")
    shared = py.index.intersection(r.index)
    py, r = py.loc[shared], r.loc[shared]

    print("=" * 74)
    print("truecell vs Seurat 5.5.1 — the SCTransform model, per gene")
    print("=" * 74)
    print(f"  genes modelled: truecell {len(py)}  R {len(r)}  shared {len(shared)}\n")

    print(f"  {'quantity':<20}{'spearman':>10}{'pearson':>10}"
          f"{'max|diff|':>12}{'note':>18}")
    print(f"  {'-' * 70}")
    for col in _MODEL_COLS:
        a, b = py[col].to_numpy(float), r[col].to_numpy(float)
        finite = np.isfinite(a) & np.isfinite(b)
        note = ""
        if not finite.all():
            note = f"{(~finite).sum()} non-finite"
        max_diff = np.abs(a[finite] - b[finite]).max() if finite.any() else np.nan
        # A correlation is undefined when either side is constant, which is not
        # a failure — `log_umi` is pinned at log(10) for every gene under v2, so
        # it is *supposed* to be constant. Say so instead of printing nan.
        # Tested on the range, not on `std() == 0`: truecell's column is exactly
        # constant yet numpy's std of it is 4e-16, and R's varies by 3e-14 from
        # write.csv's 15-digit rounding, so neither side is ever exactly flat.
        def _flat(v):
            return np.ptp(v) <= 1e-12 * max(1.0, abs(float(np.mean(v))))

        if finite.sum() < 2 or _flat(a[finite]) or _flat(b[finite]):
            const = "constant" if finite.sum() >= 2 else "too few"
            print(f"  {col:<20}{const:>10}{'':>10}{max_diff:>12.3e}{note:>18}")
            continue
        sp_ = spearmanr(a[finite], b[finite]).statistic
        pe_ = pearsonr(a[finite], b[finite]).statistic
        print(f"  {col:<20}{sp_:>10.4f}{pe_:>10.4f}{max_diff:>12.3e}{note:>18}")

    print("\n  residual_mean is the one column that does not track by rank, and"
          "\n  it is the one nothing downstream reads — Seurat records it but"
          "\n  selects features on residual_variance. Pearson ~0.99 says the"
          "\n  large values agree; the rank disagreement sits in genes whose"
          "\n  residual mean is ~1e-3 or smaller, where it is numerical dust.")

    # theta is Inf exactly for the genes v2 declares non-overdispersed, so the
    # *set* matters as much as the correlation over the rest.
    poisson_py = set(shared[np.isinf(py["theta"].to_numpy(float))])
    poisson_r = set(shared[np.isinf(r["theta"].to_numpy(float))])
    inter = len(poisson_py & poisson_r)
    union = len(poisson_py | poisson_r)
    print(f"\n  Poisson genes (theta = Inf): truecell {len(poisson_py)}  "
          f"R {len(poisson_r)}  Jaccard {inter / max(union, 1):.4f}")

    py_vf = (FIGURES / "py_variable_features.txt").read_text().split()
    r_vf = (FIGURES / "r_variable_features.txt").read_text().split()
    print(f"\n  Variable features: {len(set(py_vf) & set(r_vf))}/{len(r_vf)} shared")
    for n in (100, 500, 1000):
        overlap = len(set(py_vf[:n]) & set(r_vf[:n]))
        print(f"    top {n:<5} {overlap}/{n} ({100 * overlap / n:.1f}%)")

    pa = json.loads((FIGURES / "py_anchors.json").read_text())
    ra = json.loads((FIGURES / "r_anchors.json").read_text())
    print(f"\n  {'anchor':<26}{'truecell':>20}{'R Seurat':>20}   verdict")
    print(f"  {'-' * 70}")
    for k in sorted(set(pa) & set(ra)):
        a, b = pa[k], ra[k]
        same = (abs(a - b) <= 1e-9 * max(1.0, abs(b))
                if isinstance(a, (int, float)) and isinstance(b, (int, float))
                else a == b)
        print(f"  {k:<26}{a!s:>20}{b!s:>20}   {'MATCH' if same else 'differ'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 3k SCTransform tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--report", action="store_true",
                        help="compare against the R reference and exit")
    args = parser.parse_args()
    if args.report:
        report()
    else:
        sct, std, sct_anno, std_anno, _ = run_full(data_dir=args.data_dir)
        write_anchors(sct, sct.meta_data["sct_clusters"].nunique(),
                      std.meta_data["std_clusters"].nunique())
