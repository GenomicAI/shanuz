# Shanuz — Python Single-Cell Genomics Toolkit

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Shanuz** is a Python port of the [Seurat](https://satijalab.org/seurat/) single-cell RNA-seq
analysis framework, implementing Seurat's core data structures, preprocessing pipeline,
dimensionality reduction, clustering, and marker detection — entirely in Python.

> **Name origin:** *Shanuz* is derived from the author's name (*Shanika*).
> The package is spiritually and algorithmically faithful to Seurat v5 while providing a
> pure-Python, pip-installable alternative that integrates naturally with NumPy, SciPy,
> and AnnData ecosystems.

---

## Features

- **Shanuz object** — mirrors the R `Seurat` S4 class with `__slots__`-based Python classes
- **Assay5** — sparse-matrix-backed multi-layer assay (counts, data, scale.data)
- **Preprocessing** — `normalize_data`, `find_variable_features` (VST), `scale_data`, `percentage_feature_set`
- **Dimensionality reduction** — `run_pca` (via scikit-learn)
- **Nearest-neighbour graph** — `find_neighbors` (KNN + SNN)
- **Clustering** — `find_clusters` (Louvain via python-igraph, Leiden via leidenalg)
- **UMAP** — `run_umap` (via umap-learn)
- **Differential expression** — `find_markers`, `find_all_markers` (Wilcoxon rank-sum)
- **AnnData interoperability** — `as_anndata`, `from_anndata`
- **PBMC 3k tutorial** — end-to-end validated against the official Seurat tutorial

---

## Installation

Shanuz is not yet published to PyPI. Install directly from the GitHub repository.

### With uv (recommended)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh  # macOS/Linux
# or: powershell -c "irm https://astral.sh/uv/install.ps1 | iex"  # Windows

git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
uv venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -e ".[analysis]"
```

### With pip

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
pip install -e ".[analysis]"
```

### Full development installation (includes tests and linting)

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
uv venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -e ".[all]"
```

---

## Quick Start

```python
import scipy.sparse as sp
import numpy as np
from shanuz import create_shanuz_object

# Create a Shanuz object from a counts matrix
counts = sp.random(2000, 500, density=0.05, format="csc")
sobj = create_shanuz_object(counts, project="my_project", min_cells=3, min_features=200)
print(sobj)
# Shanuz object: my_project
#   500 cells × 2000 features
#   Active assay: 'RNA'

# Access metadata
print(sobj.meta_data.head())
```

---

## PBMC 3k Tutorial

The complete PBMC 3k guided clustering tutorial (mirroring the official
[Seurat PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial))
is included in the `tutorials/` folder.

```bash
python tutorials/pbmc3k_tutorial.py
```

This downloads the PBMC 3k dataset (~24 MB) automatically, runs the full analysis pipeline,
and validates each step against the R Seurat expected outputs:

| Step | Expected | Status |
|------|----------|--------|
| Features after filtering | 13,714 | ✅ |
| Cells after QC | 2,638 | ✅ |
| HVG top-10 overlap (≥50%) | PPBP, LYZ, S100A9 … | ✅ |
| Number of clusters (resolution=0.5) | 9 | ✅ |
| Canonical markers per cell type | all 6 cell types | ✅ |

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
from shanuz.markers import find_markers, find_all_markers

markers = find_markers(pbmc, ident_1=1)
all_markers = find_all_markers(pbmc, only_pos=True, logfc_threshold=0.25)
```

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

## Running Tests

```bash
uv pip install -e ".[dev]"
pytest tests/ -v
```

All 101 unit tests pass.

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

MIT License — see [LICENSE](LICENSE) for details.

This software is an independent reimplementation for educational and research purposes.
It is not affiliated with, endorsed by, or maintained by the Satija Lab or 10x Genomics.
