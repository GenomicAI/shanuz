"""Xenium spatial tutorial — Shanuz Python implementation & figure generator.

Reproduces, in pure Python with **Shanuz**, the same style of Xenium mast/spatial
analysis used on internal celiac data — here on a public 10x Genomics dataset
(mouse-brain coronal CTX+HP subset, the section featured in Seurat's spatial
vignette). Every step mirrors an R Seurat call so the two can be shown side by
side (see ``xenium_spatial_tutorial.md``).

It exercises the spatial Seurat-parity layer added in this branch:
    load_xenium · get_tissue_coordinates · image_dim_plot · image_feature_plot
    nearest_neighbor_distance · local_neighborhood · build_niche_assay
    composition_test

Usage
-----
    python tutorials/generate_spatial_plots.py [--data-dir PATH]

With no --data-dir the dataset auto-downloads (~20 MB) to
~/.shanuz_data/xenium_mouse_brain via shanuz.datasets.xenium_mouse_brain().

Deterministic anchors (cell counts, marker-defined cell-type counts, spatial
nearest-neighbour distances, region composition test) are written to
figures_spatial/anchors.json so they can be checked against the R reference.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import xenium_mouse_brain
from shanuz.spatial import load_xenium
from shanuz.preprocessing import normalize_data, find_variable_features, scale_data
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz import (
    get_tissue_coordinates, nearest_neighbor_distance, local_neighborhood,
    build_niche_assay, composition_test,
)
from shanuz.plotting import (
    vln_plot, dim_plot, image_dim_plot, image_feature_plot,
)

FIG = _ROOT / "tutorials" / "figures_spatial"
FIG.mkdir(parents=True, exist_ok=True)

# --- deterministic marker panels (all present in the 248-gene panel) ---------
# Ordered: ties in the argmax break toward the earlier type (matches R which.max).
CELLTYPE_MARKERS: dict[str, list[str]] = {
    "Excitatory":      ["Slc17a7", "Slc17a6", "Satb2", "Fezf2", "Neurod6"],
    "Inhibitory":      ["Gad1", "Gad2", "Pvalb", "Sst", "Vip", "Lamp5"],
    "Astrocyte":       ["Aqp4", "Gfap", "Ntsr2"],
    "Oligodendrocyte": ["Sox10", "Opalin", "Gjc3"],
    "OPC":             ["Pdgfra", "Cspg4"],
    "Vascular":        ["Cldn5", "Pecam1", "Kdr", "Emcn", "Adgrl4"],
    "Immune":          ["Cd68", "Trem2", "Siglech", "Laptm5", "Cd53"],
}
FOCAL = "Vascular"          # the sparse, spatially-patterned focal type
NICHE_K = 20
NICHES = 6


def assign_cell_types(obj, assay="Xenium"):
    """Deterministic cell type = argmax of summed raw marker counts per panel.

    Purely a function of the (identical) raw count matrix, so R and Python
    produce byte-identical labels — the anchor for the spatial comparisons.
    """
    a = obj.assays[assay]
    feats = list(a.features())
    idx = {g: i for i, g in enumerate(feats)}
    counts = a.layer_data("counts").tocsr()
    types = list(CELLTYPE_MARKERS)
    score = np.zeros((counts.shape[1], len(types)))
    for j, t in enumerate(types):
        rows = [idx[g] for g in CELLTYPE_MARKERS[t] if g in idx]
        if rows:
            score[:, j] = np.asarray(counts[rows, :].sum(axis=0)).ravel()
    best = score.argmax(axis=1)
    labels = np.where(score.max(axis=1) > 0,
                      np.array(types)[best], "Other").astype(object)
    return pd.Series(labels, index=obj.cell_names(), name="cell_type")


def main(data_dir: str | None):
    t0 = time.time()
    path = data_dir or xenium_mouse_brain()
    print(f"Data: {path}")

    # --- 1. Load Xenium (assay + per-FOV centroids) --------------------------
    obj = load_xenium(path, assay="Xenium")
    n_cells_raw = len(obj.cell_names())
    n_genes = len(obj.assays["Xenium"].features())
    print(f"Loaded {n_cells_raw} cells x {n_genes} genes; images={obj.image_names()}")

    # --- 2. QC filter (transcript-poor cells) --------------------------------
    md = obj.meta_data
    keep = list(md.index[md["nCount_Xenium"] >= 10])
    obj = obj.subset(cells=keep)
    n_cells_qc = len(obj.cell_names())
    print(f"QC nCount_Xenium>=10: {n_cells_qc} cells retained")

    # --- 3. Deterministic marker cell types ----------------------------------
    obj.meta_data["cell_type"] = assign_cell_types(obj).reindex(obj.cell_names()).values
    ct_counts = obj.meta_data["cell_type"].value_counts().sort_index()
    print("cell_type counts:\n", ct_counts.to_string())

    # --- 4. Deterministic spatial region split (dorsal/ventral by median y) ---
    coords = get_tissue_coordinates(obj)
    ymed = float(coords["y"].median())
    region = np.where(
        coords.set_index("cell").reindex(obj.cell_names())["y"].values >= ymed,
        "ventral", "dorsal")            # image y increases downward
    obj.meta_data["region"] = region

    # --- 5. Standard unsupervised pipeline (structure view) ------------------
    normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(obj, selection_method="vst", nfeatures=n_genes)
    all_genes = obj.assays["Xenium"].features()
    scale_data(obj, features=all_genes)
    run_pca(obj, n_pcs=30, features=all_genes, reduction_name="pca")
    find_neighbors(obj, dims=range(20), k_param=20)
    find_clusters(obj, resolution=0.3, algorithm=1, random_seed=0)
    run_umap(obj, dims=range(20), reduction_name="umap", seed=42)
    n_clusters = obj.meta_data["seurat_clusters"].nunique()
    print(f"clusters (res=0.3): {n_clusters}")

    # --- 6. Spatial statistics on the focal type (deterministic anchors) -----
    nn = nearest_neighbor_distance(obj, "cell_type", FOCAL)      # focal->nearest focal
    nn_median = float(nn["distance"].median())
    nn_mean = float(nn["distance"].mean())
    nb = local_neighborhood(obj, "cell_type", reference=FOCAL, k=10)
    dens_mean = float(nb[f"prop_{FOCAL}"].mean())
    print(f"{FOCAL}: n={len(nn)}  NN median={nn_median:.4f} um  "
          f"local prop mean={dens_mean:.5f}")

    # --- 7. Composition test across the deterministic region split -----------
    comp = composition_test(obj, group_by="cell_type", split_by="region",
                            reference="dorsal")
    print(f"composition_test region (chi2 p={comp.attrs['chisq_p']:.3g}):")
    print(comp[["group", "log2_ratio", "odds_ratio", "padj", "enriched_in"]]
          .to_string(index=False))

    # --- 8. Niche discovery (unsupervised spatial domains) -------------------
    build_niche_assay(obj, "cell_type", k=NICHE_K, niches=NICHES, seed=0)
    n_niche = obj.meta_data["niches"].nunique()
    print(f"niches: {n_niche}")

    # ========================= FIGURES =======================================
    def save(fig, name):
        fig.savefig(FIG / name, dpi=130, bbox_inches="tight")
        print("  wrote", name)

    save(vln_plot(obj, features=["nCount_Xenium", "nFeature_Xenium"],
                  group_by="cell_type", ncol=2, figsize=(12, 4)),
         "01_qc_violin.png")
    save(dim_plot(obj, reduction="umap", group_by="cell_type", label=True,
                  pt_size=2.0, title="UMAP — marker cell types"),
         "02_umap_celltype.png")
    save(image_dim_plot(obj, group_by="cell_type", size=1.5),
         "03_image_celltype.png")
    save(image_dim_plot(obj, group_by="seurat_clusters", size=1.5),
         "04_image_clusters.png")
    save(image_feature_plot(obj, "Slc17a7", size=1.5, cmap="viridis"),
         "05_image_feature_Slc17a7.png")
    save(image_dim_plot(obj, group_by="niches", size=1.5),
         "06_image_niches.png")
    # focal spatial map: highlight the focal type
    obj.meta_data["is_focal"] = np.where(
        obj.meta_data["cell_type"] == FOCAL, FOCAL, "other")
    save(image_dim_plot(obj, group_by="is_focal", size=1.2,
                        cols={FOCAL: "#d62728", "other": "#dddddd"}),
         "07_image_focal.png")

    # ========================= ANCHORS =======================================
    anchors = {
        "n_cells_raw": n_cells_raw,
        "n_genes": n_genes,
        "n_cells_qc": n_cells_qc,
        "celltype_counts": {k: int(v) for k, v in ct_counts.items()},
        "focal_type": FOCAL,
        "n_focal": int(len(nn)),
        "focal_nn_median": nn_median,
        "focal_nn_mean": nn_mean,
        "focal_local_density_mean": dens_mean,
        "region_ymed": ymed,
        "composition_chisq_p": float(comp.attrs["chisq_p"]),
        "composition": {
            r["group"]: {"log2_ratio": float(r["log2_ratio"]),
                         "odds_ratio": float(r["odds_ratio"]),
                         "p": float(r["p"]), "padj": float(r["padj"]),
                         "enriched_in": r["enriched_in"]}
            for _, r in comp.iterrows()
        },
        # structural (stochastic) — reported, not exact-compared
        "n_clusters": int(n_clusters),
        "n_niches": int(n_niche),
    }
    (FIG / "anchors.json").write_text(json.dumps(anchors, indent=2))
    print(f"\nwrote anchors.json  ({time.time() - t0:.1f}s total)")
    return anchors


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None,
                    help="Xenium output folder (default: auto-download cache)")
    main(ap.parse_args().data_dir)
