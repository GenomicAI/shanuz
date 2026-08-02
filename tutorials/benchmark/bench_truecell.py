#!/usr/bin/env python
"""The Truecell arm of the performance benchmark.

Run through :mod:`run_benchmarks`, which measures this process from outside::

    python tutorials/benchmark/run_benchmarks.py run --bench pbmc3k_core

Every bench here has a line-for-line counterpart in ``bench_seurat.R``, with
the same parameters in the same order. Where the two stacks genuinely cannot do
the same thing the difference is a *separate step*, never a silent one — R gets
an extra ``neighbours_annoy`` step because Seurat's default neighbour search is
approximate and truecell's is exact, so the honest comparison needs both
numbers rather than a choice between them.

Each step records an anchor: a scalar the report uses to show the two arms
produced the same result before it compares how long they took.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from steps import StepLog  # noqa: E402

# Pipeline parameters, identical on both arms.
N_HVG = 2000
N_PCS = 50
DIMS = range(10)
K_PARAM = 20
RESOLUTION = 0.5


def _import_truecell(log: StepLog):
    """Time the import, which for a CLI tool is part of what the user waits for."""
    t0 = time.time()
    import numpy as np  # noqa: F401
    from truecell.truecell import create_truecell_object  # noqa: F401
    from truecell.preprocessing import (  # noqa: F401
        normalize_data, find_variable_features, scale_data, percentage_feature_set)
    from truecell.reduction import run_pca  # noqa: F401
    from truecell.neighbors import find_neighbors  # noqa: F401
    from truecell.clustering import find_clusters  # noqa: F401
    from truecell.umap import run_umap  # noqa: F401
    from truecell.markers import find_all_markers  # noqa: F401
    log.mark("import", time.time() - t0)


# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------

def _load(bench: str):
    """Return (counts, genes, cells, meta) for a bench's dataset."""
    import truecell.datasets as ds
    key = bench.split("_")[0]
    if key == "pbmc3k":
        counts, genes, cells = ds.pbmc3k()
        return counts, genes, cells, None
    if key == "pbmc8k":
        counts, genes, cells = ds.pbmc8k()
        return counts, genes, cells, None
    if key == "ifnb":
        counts, genes, cells, meta = ds.ifnb()
        return counts, genes, cells, meta
    if key == "thp1":
        rna, rna_genes, _adt, _adt_names, meta, cells = ds.thp1_eccite()
        return rna, rna_genes, cells, meta
    raise SystemExit(f"unknown dataset for bench {bench!r}")


def _write_idents(bench: str, obj) -> None:
    """Hand this run's cell-to-cluster assignment to the R arm."""
    import pandas as pd
    out = Path(__file__).resolve().parent / "results" / f"{bench}_idents.csv"
    out.parent.mkdir(exist_ok=True)
    pd.DataFrame({"cell": list(obj.cell_names()),
                  "ident": obj.meta_data["seurat_clusters"].astype(str).values}
                 ).to_csv(out, index=False)


# ---------------------------------------------------------------------------
# The standard workflow
# ---------------------------------------------------------------------------

def bench_core(bench: str, log: StepLog) -> None:
    """Counts on disk to a marker table: the pipeline every tutorial starts with."""
    from truecell.truecell import create_truecell_object
    from truecell.preprocessing import (
        normalize_data, find_variable_features, scale_data, percentage_feature_set)
    from truecell.reduction import run_pca
    from truecell.neighbors import find_neighbors
    from truecell.clustering import find_clusters
    from truecell.umap import run_umap
    from truecell.markers import find_all_markers

    with log.step("read_counts") as rec:
        counts, genes, cells, _meta = _load(bench)
        rec["anchor"] = counts.shape[1]

    with log.step("create_object") as rec:
        obj = create_truecell_object(
            counts=counts, assay="RNA", min_cells=3, min_features=200,
            project=bench, feature_names=genes, cell_names=cells)
        rec["anchor"] = len(obj.assays["RNA"]._all_feature_names)
    del counts

    with log.step("qc_metrics") as rec:
        percentage_feature_set(obj, pattern=r"^MT-", col_name="percent.mt")
        rec["anchor"] = round(float(obj.meta_data["percent.mt"].mean()), 6)

    with log.step("normalize") as rec:
        normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
        rec["anchor"] = round(float(obj.assays["RNA"].layers["data"].sum()), 3)

    with log.step("hvg_vst") as rec:
        find_variable_features(obj, selection_method="vst", nfeatures=N_HVG)
        hvg = obj.assays["RNA"].variable_features
        rec["anchor"] = len(hvg)

    with log.step("scale_hvg") as rec:
        scale_data(obj, features=hvg)
        rec["anchor"] = len(hvg)

    # Scaling every gene is what the PBMC 3k vignette does and what dominates
    # its memory. Only measured on the small dataset: at 20k cells the dense
    # result is tens of gigabytes in both languages, which measures the machine
    # rather than either tool.
    if bench.startswith("pbmc3k"):
        with log.step("scale_all_genes") as rec:
            all_genes = list(obj.assays["RNA"]._all_feature_names)
            scale_data(obj, features=all_genes)
            rec["anchor"] = len(all_genes)
        with log.step("rescale_hvg") as rec:
            scale_data(obj, features=hvg)
            rec["anchor"] = len(hvg)

    with log.step("pca") as rec:
        run_pca(obj, n_pcs=N_PCS, features=hvg, reduction_name="pca")
        rec["anchor"] = round(float(obj.reductions["pca"].stdev[0]), 4)

    with log.step("neighbours_exact") as rec:
        find_neighbors(obj, dims=DIMS, k_param=K_PARAM)
        rec["anchor"] = int(obj.graphs["RNA_snn"].nnz)

    with log.step("cluster_louvain") as rec:
        find_clusters(obj, resolution=RESOLUTION, algorithm=1, random_seed=0)
        rec["anchor"] = int(obj.meta_data["seurat_clusters"].nunique())

    # Untimed, and the reason the truecell arm has to run first. The two tools
    # do not always land on the same number of clusters — on PBMC 8k it is 9
    # against Seurat's 12 — and one-vs-rest markers cost one test per cluster.
    # Timing `find_all_markers` against a different number of comparisons would
    # report a clustering difference as a speed difference. R reads this file
    # and adopts the assignment before its own marker step.
    _write_idents(bench, obj)

    # metric="cosine" because that is `RunUMAP`'s default and this benchmark is
    # only worth reading if both arms embed the same distances. umap-learn's
    # own default is euclidean, which is where truecell's default comes from.
    with log.step("umap") as rec:
        run_umap(obj, dims=DIMS, reduction_name="umap", seed=42, metric="cosine")
        rec["anchor"] = int(obj.reductions["umap"].cell_embeddings.shape[0])

    # The same embedding with no seed. umap-learn drops to a single thread the
    # moment a random_state is set — it says so in a warning — so the seeded
    # step above is measuring reproducibility as much as speed, and the gap
    # between these two steps is what that reproducibility costs.
    with log.step("umap_unseeded") as rec:
        run_umap(obj, dims=DIMS, reduction_name="umap_unseeded", seed=None,
                 metric="cosine")
        rec["anchor"] = int(obj.reductions["umap_unseeded"].cell_embeddings.shape[0])

    with log.step("find_all_markers") as rec:
        markers = find_all_markers(obj, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
        rec["anchor"] = int(len(markers))

    print(f"{bench}: {len(obj.cell_names())} cells, "
          f"{obj.meta_data['seurat_clusters'].nunique()} clusters, "
          f"{len(markers)} markers")


# ---------------------------------------------------------------------------
# Named heavy operations
# ---------------------------------------------------------------------------

def bench_sctransform(bench: str, log: StepLog) -> None:
    """Regularized negative-binomial normalization — the expensive alternative
    to LogNormalize, and the step most likely to dominate a real run."""
    from truecell import sctransform
    from truecell.truecell import create_truecell_object
    from truecell.reduction import run_pca

    with log.step("read_counts") as rec:
        counts, genes, cells, _ = _load(bench)
        rec["anchor"] = counts.shape[1]
    with log.step("create_object") as rec:
        obj = create_truecell_object(counts=counts, assay="RNA", min_cells=3,
                                     min_features=200, project=bench,
                                     feature_names=genes, cell_names=cells)
        rec["anchor"] = len(obj.assays["RNA"]._all_feature_names)
    with log.step("sctransform") as rec:
        sctransform(obj, n_features=3000, seed=42)
        rec["anchor"] = len(obj.assays["SCT"].variable_features)
    with log.step("pca_on_sct") as rec:
        run_pca(obj, n_pcs=30, assay="SCT",
                features=obj.assays["SCT"].variable_features)
        rec["anchor"] = round(float(obj.reductions["pca"].stdev[0]), 4)


def bench_integration(bench: str, log: StepLog) -> None:
    """The three batch-correction paths, on a shared PCA so only the
    integration algorithms are being compared."""
    from truecell.truecell import create_truecell_object
    from truecell.preprocessing import normalize_data, find_variable_features, scale_data
    from truecell.reduction import run_pca
    from truecell.integration import run_harmony, integrate_layers

    n_pcs = 30
    with log.step("read_counts") as rec:
        counts, genes, cells, meta = _load(bench)
        rec["anchor"] = counts.shape[1]
    with log.step("prep_to_pca") as rec:
        obj = create_truecell_object(counts=counts, assay="RNA", min_cells=3,
                                     project=bench, feature_names=genes,
                                     cell_names=cells, meta_data=meta)
        normalize_data(obj)
        find_variable_features(obj, selection_method="vst", nfeatures=N_HVG)
        hvg = obj.assays["RNA"].variable_features
        scale_data(obj, features=hvg)
        run_pca(obj, n_pcs=n_pcs, features=hvg, reduction_name="pca")
        rec["anchor"] = round(float(obj.reductions["pca"].stdev[0]), 4)

    with log.step("harmony") as rec:
        run_harmony(obj, group_by="stim", reduction="pca",
                    reduction_name="harmony", seed=0)
        rec["anchor"] = int(obj.reductions["harmony"].cell_embeddings.shape[0])

    # k_weight cannot exceed the smallest batch's anchor count; the tutorial
    # picks it the same way, so both arms weight the same number of anchors.
    smallest = int(meta["stim"].value_counts().min())
    k_weight = int(min(100, smallest // 2))
    for method in ("cca", "rpca"):
        with log.step(f"integrate_{method}") as rec:
            integrate_layers(obj, method=method, group_by="stim",
                             new_reduction=method, seed=42, k_weight=k_weight)
            rec["anchor"] = int(obj.reductions[method].cell_embeddings.shape[0])


def bench_de(bench: str, log: StepLog) -> None:
    """Every differential-expression test, one cluster against the rest, on a
    fixed cell assignment so no clustering difference leaks into the timing."""
    from truecell.truecell import create_truecell_object
    from truecell.preprocessing import normalize_data, find_variable_features, scale_data
    from truecell.reduction import run_pca
    from truecell.neighbors import find_neighbors
    from truecell.clustering import find_clusters
    from truecell.markers import find_markers

    with log.step("prep") as rec:
        counts, genes, cells, _ = _load(bench)
        obj = create_truecell_object(counts=counts, assay="RNA", min_cells=3,
                                     min_features=200, project=bench,
                                     feature_names=genes, cell_names=cells)
        normalize_data(obj)
        find_variable_features(obj, selection_method="vst", nfeatures=N_HVG)
        hvg = obj.assays["RNA"].variable_features
        scale_data(obj, features=hvg)
        run_pca(obj, n_pcs=N_PCS, features=hvg)
        find_neighbors(obj, dims=DIMS, k_param=K_PARAM)
        find_clusters(obj, resolution=RESOLUTION, algorithm=1, random_seed=0)
        rec["anchor"] = int(obj.meta_data["seurat_clusters"].nunique())

    # Written to disk so the R arm tests the same two groups of cells; a DE
    # timing on a different cell split is not a comparison.
    _write_idents(bench, obj)

    # Step names are Seurat's spelling so the two logs line up; the values are
    # truecell's, which lower-cases the two tests named after their R packages.
    tests = (("wilcox", "wilcox"), ("t", "t"), ("bimod", "bimod"), ("LR", "LR"),
             ("negbinom", "negbinom"), ("roc", "roc"), ("MAST", "mast"),
             ("DESeq2", "deseq2"))
    for step_name, test in tests:
        with log.step(f"de_{step_name}") as rec:
            try:
                res = find_markers(obj, ident_1="0", ident_2="1", test_use=test,
                                   logfc_threshold=0.25, min_pct=0.1)
                rec["anchor"] = int(len(res))
            except Exception as exc:  # a missing optional backend, not a timing
                print(f"  de_{step_name} unavailable: {exc}")
                rec["anchor"] = -1


def bench_spatial(bench: str, log: StepLog) -> None:
    """Moran's I on a Xenium slide, at two sizes.

    ``morans_i_2k`` is the like-for-like number, on a 2,000-cell subset both
    arms read from the same file. ``morans_i_full`` is the whole 36,602-cell
    slide, and only truecell has it: Seurat's ``RunMoransI`` computes
    ``as.matrix(dist(pos))``, a dense n x n distance matrix, which is 10.7 GB
    at this n before any statistic is computed. That is not a slow step to
    measure — it is the reason the spatial tutorial subsets in the first place.
    """
    from truecell import find_spatially_variable_features, load_xenium
    from truecell.datasets import xenium_mouse_brain
    from truecell.preprocessing import normalize_data

    n_subset = 2000
    with log.step("read_xenium") as rec:
        path = xenium_mouse_brain()
        obj = load_xenium(path)
        rec["anchor"] = len(obj.cell_names())
    with log.step("normalize") as rec:
        normalize_data(obj, scale_factor=10000)
        rec["anchor"] = len(obj.cell_names())

    # A deterministic subset, written for the R arm. Sorted then evenly spaced
    # rather than random: it has to be reproducible from either language, and
    # spacing it keeps the subset spread over the slide instead of clipping a
    # corner, which would flatter Moran's I on both sides equally but make the
    # statistic meaningless.
    all_cells = sorted(obj.cell_names())
    step_by = max(1, len(all_cells) // n_subset)
    cells = all_cells[::step_by][:n_subset]
    out = Path(__file__).resolve().parent / "results" / "xenium_cells.txt"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(cells) + "\n")
    sub = obj.subset(cells=cells)

    with log.step("morans_i_2k") as rec:
        res = find_spatially_variable_features(sub, method="moransi",
                                               weights="inverse_square")
        rec["anchor"] = int(len(res))
    with log.step("morans_i_full") as rec:
        res = find_spatially_variable_features(obj, method="moransi",
                                               weights="inverse_square")
        rec["anchor"] = int(len(res))


def bench_blas(bench: str, log: StepLog) -> None:
    """Dense linear algebra on identical inputs, touching neither library.

    Not a single-cell benchmark: a control. R here links the reference BLAS it
    ships with, numpy links Accelerate. Every PCA, scaling and SVD number in
    this suite inherits that difference, and without measuring it directly
    there is no way to say how much of a gap belongs to the two projects and
    how much belongs to the two BLAS builds underneath them.
    """
    import numpy as np

    n, k = 2000, 500
    with log.step("blas_setup") as rec:
        rng = np.random.default_rng(0)
        a = rng.standard_normal((n, n))
        rec["anchor"] = n
    with log.step("blas_gemm") as rec:
        g = a @ a.T
        rec["anchor"] = round(float(g[0, 0]), 3)
    with log.step("blas_svd") as rec:
        s = np.linalg.svd(a[:, :k], full_matrices=False, compute_uv=False)
        rec["anchor"] = round(float(s[0]), 3)
    with log.step("blas_crossprod_chol") as rec:
        c = np.linalg.cholesky(a.T @ a + n * np.eye(n))
        rec["anchor"] = round(float(c[0, 0]), 3)


BENCHES = {
    "blas_probe": bench_blas,
    "pbmc3k_core": bench_core,
    "pbmc8k_core": bench_core,
    "ifnb_core": bench_core,
    "thp1_core": bench_core,
    "pbmc3k_sctransform": bench_sctransform,
    "ifnb_integration": bench_integration,
    "pbmc3k_de": bench_de,
    "xenium_spatial": bench_spatial,
}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--bench", required=True, choices=sorted(BENCHES))
    p.add_argument("--steps", required=True)
    args = p.parse_args()

    log = StepLog(args.steps)
    _import_truecell(log)
    try:
        BENCHES[args.bench](args.bench, log)
    finally:
        log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
