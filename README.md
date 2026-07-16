# Shanuz — Python Single-Cell Genomics Toolkit

[![PyPI](https://img.shields.io/pypi/v/shanuz.svg)](https://pypi.org/project/shanuz/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Shanuz** is a Python port of the [Seurat](https://satijalab.org/seurat/) single-cell RNA-seq
analysis framework, implementing Seurat's core data structures, preprocessing pipeline,
dimensionality reduction, clustering, and marker detection — entirely in Python.

> The package is spiritually and algorithmically faithful to Seurat v5 while providing a
> pure-Python, pip-installable alternative that integrates naturally with NumPy, SciPy,
> and AnnData ecosystems.

---

## Features

- **Shanuz object** — mirrors the R `Seurat` S4 class with `__slots__`-based Python classes
- **Assay5** — sparse-matrix-backed multi-layer assay (counts, data, scale.data)
- **Preprocessing** — `normalize_data`, `find_variable_features` (VST), `scale_data`, `percentage_feature_set`
- **SCTransform** — `sctransform` (regularized negative-binomial Pearson residuals)
- **Signature scoring** — `add_module_score`, `cell_cycle_scoring` (S/G2M + Phase)
- **Dimensionality reduction** — `run_pca`, `run_spca` (supervised, off a cell graph), `run_ica`, `run_tsne`, `glm_pca` (Poisson or negative binomial, straight on counts)
- **Batch correction / integration** — `run_harmony` (via harmonypy), CCA/RPCA anchors (`find_integration_anchors` + `integrate_data`), and the `integrate_layers` dispatcher (`method="harmony"|"cca"|"rpca"`)
- **Reference mapping** — `find_transfer_anchors` (project a query into a reference; `pcaproject` or `cca`) + `transfer_data` (annotate the query with reference labels, or impute reference expression onto it); `project_umap` / `map_query` place the query in the reference's own UMAP in one call
- **Scale (sketching)** — `sketch_data` draws a leverage-weighted subset of a huge dataset (rare states kept, not lost), `leverage_score` computes the per-cell scores via a CountSketch (no full SVD), and `project_data` extends the sketch's PCA/UMAP/labels back to every cell
- **Scale (lazy on-disk matrices)** — `LazyMatrix` keeps a matrix out-of-core as memory-mapped compressed-sparse-column arrays (BPCells-style); `write_lazy_matrix` / `open_lazy_matrix` persist and map it, a slice reads only the touched cells off disk, `col_blocks` streams a million cells at bounded RAM, and it drops straight into an `Assay5` layer — no new dependency
- **Cell hashing (demultiplexing)** — `hto_demux` (Seurat's `HTODemux`) demultiplexes pooled samples from hashtag counts: CLR normalize → cluster into `k = n_hashtags + 1` groups (`kfunc="kmeans"`, or `"clara"` for Seurat's k-medoids) → per-hashtag negative-binomial background threshold → singlet / doublet / negative calls, written to `meta_data` (`HTO_maxID`, `HTO_classification`, …) plus a `hash.ID` identity. `multiseq_demux` (Seurat's `MULTIseqDemux`) is the MULTI-seq alternative — a Gaussian-KDE quantile threshold per barcode, with an `autothresh` sweep — writing `MULTI_ID` / `MULTI_classification`
- **Pooled CRISPR screens (Mixscape)** — `calc_perturb_sig` (Seurat's `CalcPerturbSig`) subtracts each cell's nearest non-targeting controls to isolate its perturbation signature, then `run_mixscape` (Seurat's `RunMixscape`) separates true knockouts from non-perturbed escapers per guide — gene-vs-NT DE, then an iterative 2-component Gaussian mixture over the perturbation score — writing `mixscape_class` (`"<gene> KO"` / `NP` / `NT`, also the identity), `mixscape_class.global`, and `mixscape_class_p_ko`. `mixscape_lda` (Seurat's `MixscapeLDA`) adds the supervised map on which each guide population forms its own cloud — per-guide DE-gene PCA subspaces, every cell projected onto each, then one linear discriminant analysis over the concatenation → an `lda` reduction plus `lda_assignments` / `LDAP_<class>`. Two diagnostics complete the workflow: `plot_perturb_score` (Seurat's `PlotPerturbScore`) overlays the NT control density against one guide's own along the perturbation score — the axis mixscape actually splits on, bimodal when the guide has a real effect — and `mixscape_heatmap` (Seurat's `MixscapeHeatmap`) shows the DE genes underneath it with every cell ordered by its knockout probability
- **Nearest-neighbour graph** — `find_neighbors` (KNN + SNN)
- **Multimodal WNN** — `find_multi_modal_neighbors` (full two-stage port: per-cell RNA/protein weights via exponential kernel + softmax, then a joint neighbour search building the `wknn`/`wsnn` graphs)
- **Clustering** — `find_clusters` (Louvain via python-igraph, Leiden via leidenalg)
- **UMAP** — `run_umap` (via umap-learn; embeds a reduction or a precomputed graph)
- **PC significance** — `jack_straw`, `score_jackstraw` (JackStraw permutation test)
- **Differential expression** — `find_markers`, `find_all_markers` (`wilcox` tie-corrected, `t`, `bimod`, `LR`, `negbinom`, `mast` hurdle, `deseq2` pseudobulk, `roc`), `find_conserved_markers` (cross-condition, Fisher-combined)
- **Pseudobulk** — `aggregate_expression` (sum counts per group → matrix or one-cell-per-group object), pseudobulk DESeq2 via `find_markers(test_use="deseq2", sample_col=...)`
- **Plotting** — `dim_plot`, `feature_plot`, `vln_plot`, `dot_plot`, `elbow_plot`, `do_heatmap`, `dim_heatmap`, `feature_scatter`, `variable_feature_plot`, `ridge_plot`, `plot_perturb_score`, `mixscape_heatmap` (matplotlib/seaborn)
- **AnnData interoperability** — `as_anndata`, `from_anndata`
- **Spatial (Xenium / Visium / CosMx / MERSCOPE)** — `load_xenium`/`load_visium`/`load_cosmx`/`load_merscope`, `get_tissue_coordinates`, `nearest_neighbor_distance`, `local_neighborhood`, `build_niche_assay`, `find_spatially_variable_features` (Moran's I + mark variogram), `composition_test`, `image_dim_plot`, `image_feature_plot`
- **Visium tissue images** — `load_visium` reads the H&E PNG + `scalefactors_json.json` into a `VisiumV2` image (Seurat v5's class): `get_image()`, `scale_factors`, `radius()`, `scale_coordinates()`; `spatial_dim_plot` / `spatial_feature_plot` draw spots over that image at their true diameter
- **PBMC 3k tutorial** — end-to-end validated against the official Seurat tutorial
- **PBMC 8k advanced tutorial** — larger dataset + T/NK subclustering workflow
- **CITE-seq multimodal tutorial** — RNA + surface protein (ADT) with CLR normalization and WNN joint clustering
- **Xenium spatial tutorial** — spatial neighbourhood/niche analysis, verified to 8 s.f. against R Seurat

---

## Installation

Shanuz is published on [PyPI](https://pypi.org/project/shanuz/) — `pip install shanuz` just works.
Installing from source (editable install) is only needed if you want to modify shanuz itself.

### From PyPI (recommended)

```bash
pip install shanuz                 # core: object model, preprocessing, PCA, markers
pip install "shanuz[analysis]"     # + clustering, UMAP, plotting (matplotlib/seaborn)
pip install "shanuz[anndata]"      # + AnnData interoperability
pip install "shanuz[integration]"  # + Harmony batch correction (harmonypy)
pip install "shanuz[all]"          # everything (analysis + anndata + integration + dev/test tooling)
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv pip install "shanuz[analysis]"
```

### From source (for development)

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
uv venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -e ".[all]"  # editable install + tests/linting
```

With `pip` instead of `uv`:

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
pip install -e ".[analysis]"
```

### Troubleshooting

**`ModuleNotFoundError: No module named 'numpy._core._multiarray_umath'`** (or a
similar broken NumPy/SciPy import) means the virtual environment has a corrupt or
partially-installed NumPy — usually left over from an interrupted install or from
mixing installers. It is not a shanuz issue. Fix it by reinstalling NumPy, or by
recreating the environment:

```bash
uv pip install --reinstall --no-cache numpy      # quick fix
# or start clean:
rm -rf .venv && uv venv && uv pip install "shanuz[all]"
```

---

## Quick Start

```python
import scipy.sparse as sp
import numpy as np
from shanuz import create_shanuz_object

# Create a Shanuz object from a counts matrix
counts = sp.random(2000, 500, density=0.2, format="csc")
sobj = create_shanuz_object(counts, project="my_project", min_cells=3, min_features=200)
print(sobj)
# Shanuz object — my_project
#   500 cells × 2000 features
#   Active assay: 'RNA'
#   Reductions: []
#   Version: 5.4.0

# Access metadata
print(sobj.meta_data.head())
```

---

## Tutorials

Five end-to-end tutorials — from basic guided clustering through multimodal
CITE-seq to Xenium spatial — each pairing R Seurat code side-by-side with the
Python Shanuz equivalent.
See **[`tutorials/README.md`](https://github.com/GenomicAI/shanuz/blob/main/tutorials/README.md)** for the full index.

| # | Tutorial | Dataset | Complexity |
|---|----------|---------|-----------|
| 1 | [PBMC 3k — Guided Clustering](https://github.com/GenomicAI/shanuz/blob/main/tutorials/pbmc3k_tutorial.md) | 3k PBMCs · 10x Genomics | Beginner |
| 2 | [PBMC 8k — Advanced Subclustering](https://github.com/GenomicAI/shanuz/blob/main/tutorials/advanced_pbmc8k_subclustering.md) | 8k PBMCs · GRCh38 | Intermediate |
| 3 | [CBMC CITE-seq — Multimodal](https://github.com/GenomicAI/shanuz/blob/main/tutorials/multimodal_citeseq.md) | 8,600 CBMCs · RNA + 13 proteins | Advanced |
| 4 | [PBMC 3k — SCTransform](https://github.com/GenomicAI/shanuz/blob/main/tutorials/sctransform_vignette.md) | 3k PBMCs · 10x Genomics | Advanced |
| 5 | [Xenium — Spatial (R vs Python)](https://github.com/GenomicAI/shanuz/blob/main/tutorials/xenium_spatial_tutorial.md) | 36k cells · 10x Xenium mouse brain | Spatial |

```bash
# Tutorial 1 — PBMC 3k
python tutorials/pbmc3k_tutorial.py && python tutorials/generate_plots.py

# Tutorial 2 — PBMC 8k subclustering
python tutorials/pbmc8k_subclustering_tutorial.py && python tutorials/generate_advanced_plots.py

# Tutorial 3 — CITE-seq multimodal
python tutorials/cbmc_citeseq_tutorial.py && python tutorials/generate_multimodal_plots.py

# Tutorial 4 — SCTransform
python tutorials/pbmc3k_sctransform_tutorial.py && python tutorials/generate_sctransform_plots.py

# Tutorial 5 — Xenium spatial (auto-downloads ~20 MB)
python tutorials/generate_spatial_plots.py
```

---

## API Reference

### Object creation

```python
from shanuz import create_shanuz_object

pbmc = create_shanuz_object(
    counts,             # scipy.sparse CSC/CSR or numpy ndarray (genes × cells)
    project="pbmc3k",
    min_cells=3,        # filter genes present in fewer than N cells
    min_features=200,   # filter cells with fewer than N detected genes
)
```

### Preprocessing

```python
from shanuz.preprocessing import (
    normalize_data,
    find_variable_features,
    scale_data,
    percentage_feature_set,
)

percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
scale_data(pbmc)
```

### Dimensionality reduction & clustering

```python
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap

run_pca(pbmc, n_pcs=50)
find_neighbors(pbmc, dims=range(10), k_param=20)
find_clusters(pbmc, resolution=0.5)
run_umap(pbmc, dims=range(10))
```

### Differential expression

```python
from shanuz import (
    find_markers, find_all_markers, find_conserved_markers, aggregate_expression,
)

markers = find_markers(pbmc, ident_1=1)
all_markers = find_all_markers(pbmc, only_pos=True, logfc_threshold=0.25)

# Markers up in cluster 1 across every condition (Fisher-combined p per gene).
conserved = find_conserved_markers(pbmc, ident_1=1, grouping_var="condition")

# Pseudobulk counts summed per (cell type × donor) — input for sample-level DE.
pseudobulk = aggregate_expression(pbmc, group_by=["cell_type", "donor"])

# Pseudobulk DESeq2 between two conditions, one profile per donor (needs
# `pip install shanuz[deseq2]`). pbmc.idents must hold the two conditions.
de = find_markers(pbmc, ident_1="stim", ident_2="ctrl",
                  test_use="deseq2", sample_col="donor")
```

### Plotting

All plotting functions return a `matplotlib.figure.Figure` — save or display as needed.

```python
from shanuz.plotting import (
    dim_plot,            # DimPlot   — cells on UMAP/PCA coloured by ident
    feature_plot,        # FeaturePlot — gene expression on embedding
    vln_plot,            # VlnPlot   — violin plots per cluster
    elbow_plot,          # ElbowPlot — stdev per PC
    feature_scatter,     # FeatureScatter — two features vs each other
    variable_feature_plot, # VariableFeaturePlot — mean-variance HVG plot
    dim_heatmap,         # DimHeatmap — top loading genes per PC
    do_heatmap,          # DoHeatmap  — expression heatmap sorted by cluster
    ridge_plot,          # RidgePlot  — ridgeline plots per cluster
)

# Quick examples
fig = dim_plot(pbmc, reduction="umap", label=True)
fig = feature_plot(pbmc, ["LYZ", "MS4A1", "NKG7"], reduction="umap", ncol=3)
fig = vln_plot(pbmc, ["LYZ", "CD3D", "PPBP"], group_by=None)
fig = elbow_plot(pbmc, ndims=20)
fig = do_heatmap(pbmc, top_marker_genes)
fig.savefig("output.png", dpi=150, bbox_inches="tight")
```

| Shanuz function | R Seurat equivalent |
|-----------------|---------------------|
| `dim_plot` | `DimPlot` |
| `feature_plot` | `FeaturePlot` |
| `vln_plot` | `VlnPlot` |
| `dot_plot` | `DotPlot` |
| `elbow_plot` | `ElbowPlot` |
| `feature_scatter` | `FeatureScatter` |
| `variable_feature_plot` | `VariableFeaturePlot` |
| `dim_heatmap` | `DimHeatmap` |
| `do_heatmap` | `DoHeatmap` |
| `ridge_plot` | `RidgePlot` |

---

## Data Structures

```
Shanuz
├── assays: dict[str, Assay5]
│   └── "RNA"
│       ├── layers["counts"]    # raw integer counts (genes × cells)
│       ├── layers["data"]      # log-normalized (genes × cells)
│       └── layers["scale.data"] # z-scored (genes × cells)
├── meta_data: pd.DataFrame     # per-cell metadata
├── reductions: dict
│   ├── "pca": DimReduc         # PCA embeddings + loadings
│   └── "umap": DimReduc        # UMAP embeddings
├── graphs: dict
│   ├── "RNA_nn": Graph         # KNN graph
│   └── "RNA_snn": Graph        # SNN graph
└── commands: list[ShanuzCommand]  # audit log
```

---

## Roadmap

See **[`ROADMAP.md`](https://github.com/GenomicAI/shanuz/blob/main/ROADMAP.md)** for the full development plan. Milestones:

| Milestone | Focus |
|-----------|-------|
| v0.2.0 | Batch correction — Harmony, CCA/RPCA anchors, `IntegrateLayers` dispatcher ✅ *(complete)* |
| v0.3.0 | Reference mapping — `FindTransferAnchors`, `TransferData`, `MapQuery`/`ProjectUMAP` ✅ *(complete)* |
| v0.4.0 | Multimodal WNN — `FindMultiModalNeighbors`, joint UMAP/clustering ✅ *(delivered — see Tutorial 3)* |
| v0.5.0 | Additional reductions — t-SNE, ICA, `run_spca`, `glm_pca` (Poisson + negative binomial) ✅ *(complete)* |
| v0.6.0 | Pseudobulk & advanced DE — `AggregateExpression`, `FindConservedMarkers`, DESeq2 (`test_use="deseq2"`), MAST (`test_use="mast"`), bimod (`test_use="bimod"`) ✅ *(complete)* |
| v0.7.0 | Spatial — Xenium/Visium/CosMx/MERSCOPE loaders, niche/neighbourhood analysis, `find_spatially_variable_features` (Moran's I + markvariogram), `image_*` plots, `VisiumV2` tissue images, `spatial_*` H&E plots ✅ *(delivered — see Tutorial 5)* |
| v0.8.0 | Scale — `SketchData`/`ProjectData` (leverage-score sketching) ✅; BPCells-style lazy on-disk matrices (`LazyMatrix`) ✅ *(complete)* |
| v0.9.0 | Specialized — `HTODemux` ✅ + `MULTIseqDemux` ✅ (cell hashing); Mixscape ✅ (`CalcPerturbSig` + `RunMixscape` + `MixscapeLDA` + `PlotPerturbScore` + `MixscapeHeatmap`, CRISPR screens) — **complete** |
| v0.10.0 | Infrastructure — PyPI, GitHub Actions CI, type annotations, MkDocs site |

---

## Running Tests

```bash
uv pip install -e ".[dev]"
pytest tests/ -v
```

All 156 unit tests pass.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| numpy, scipy, pandas | Core numerics and data frames |
| statsmodels | LOESS smoothing for VST |
| scikit-learn | PCA |
| umap-learn | UMAP embedding |
| python-igraph | Louvain clustering |
| leidenalg | Leiden clustering |
| packaging | Version handling |

---

## Credits

**Development assistance:** This package was developed with the help of
[Claude](https://claude.ai) (Anthropic's AI assistant, claude-sonnet-4-6),
which assisted in porting the R Seurat codebase to Python, implementing
the VST algorithm, degree-2 LOESS, Louvain clustering, and validating
results against the official PBMC 3k tutorial.

**Original R Seurat package:**  
The algorithms and data structures in Shanuz are direct Python translations of the
R [Seurat](https://satijalab.org/seurat/) package by the Satija Lab.
Please cite the original Seurat papers if you use Shanuz in published work:

> Hao Y, Stuart T, Kowalski MH, et al. (2024).
> **Dictionary learning for integrative, multimodal and scalable single-cell analysis.**
> *Nature Biotechnology*, 42, 293–304.
> https://doi.org/10.1038/s41587-023-01767-y

> Hao Y, Hao S, Andersen-Nissen E, et al. (2021).
> **Integrated analysis of multimodal single-cell data.**
> *Cell*, 184(13), 3573–3587.
> https://doi.org/10.1016/j.cell.2021.04.048

> Stuart T, Butler A, Hoffman P, et al. (2019).
> **Comprehensive Integration of Single-Cell Data.**
> *Cell*, 177(7), 1888–1902.
> https://doi.org/10.1016/j.cell.2019.05.031

> Butler A, Hoffman P, Smibert P, Papalexi E, Satija R. (2018).
> **Integrating single-cell transcriptomic data across different conditions, technologies, and species.**
> *Nature Biotechnology*, 36, 411–420.
> https://doi.org/10.1038/nbt.4096

**PBMC 3k dataset:**  
10x Genomics. (2016). *3k PBMCs from a Healthy Donor*.
https://www.10xgenomics.com/resources/datasets/3-k-pb-mcs-from-a-healthy-donor-1-standard-1-1-0

---

## License

MIT License — see [LICENSE](https://github.com/GenomicAI/shanuz/blob/main/LICENSE) for details.

This software is an independent reimplementation for educational and research purposes.
It is not affiliated with, endorsed by, or maintained by the Satija Lab or 10x Genomics.
