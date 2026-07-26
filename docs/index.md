---
hide:
  - navigation
---

# shanuz

**Seurat's single-cell pipeline, in Python, checked against R Seurat number by
number.**

`shanuz` ports Seurat v5 — the object model, preprocessing, dimensional
reduction, clustering, differential expression, integration, spatial and the
rest — to pure Python on NumPy, SciPy and pandas. It is not a reimplementation
that borrows the ideas. It follows the same code paths, keeps the same defaults,
and where the two disagree, the disagreement is measured and written down.

```python
import shanuz
from shanuz.datasets import pbmc3k

counts, genes, cells = pbmc3k()
pbmc = shanuz.create_shanuz_object(
    counts=counts, feature_names=genes, cell_names=cells,
    project="pbmc3k", min_cells=3, min_features=200,
)

shanuz.normalize_data(pbmc)
shanuz.find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
shanuz.scale_data(pbmc)
shanuz.run_pca(pbmc, n_pcs=50)
shanuz.find_neighbors(pbmc, dims=range(10), k_param=20)
shanuz.find_clusters(pbmc, resolution=0.5)
shanuz.run_umap(pbmc, dims=range(10), seed=42)

markers = shanuz.find_all_markers(pbmc, only_pos=True, min_pct=0.25)
```

[Install it](installation.md){ .md-button .md-button--primary }
[Run the first tutorial](tutorials/pbmc3k_tutorial.md){ .md-button }

## What "checked against R" means here

Eighteen tutorials, each with a matching R script that runs the same analysis
under real Seurat 5.5.1 on the same data, and a comparison that names specific
numbers rather than declaring success. A few of them:

| Comparison | Result |
|---|---|
| Object model — `Cells`, `Layers`, `FetchData`, `Idents`, `Command` | 91 of 91 anchors exact, no tolerance |
| Differential expression, all eight tests | seven per-cell tests reproduce Seurat's top 50 exactly; `avg_log2FC` to 7.1e-15 |
| Anchors, RPCA | 649 of 649 of Seurat's v4 anchors, 30 of 30 v5 embedding PCs |
| Moran's I on 36,602 Xenium cells | 1.6e-14, on a slide R cannot hold in memory |
| Cell hashing against cross-species ground truth | 99.81 % call-concordant |
| Label transfer, celseq2 → smartseq2 | 98.71 % per-cell concordant |

The tutorials are not a demo. **They are how the defects were found** — dozens of
them in this package, and one in Seurat. Each vignette says which bugs its
comparison caught and what the number was before and after. See
[Fidelity](fidelity.md) for how the checking works and where the two tools
genuinely differ.

## Where to start

<div class="grid cards" markdown>

-   **[Installation](installation.md)**

    Python 3.12+. `pip install shanuz`, plus which extra you need for what.

-   **[Quickstart](quickstart.md)**

    PBMC 3k from counts to labelled clusters, and the same thing in R.

-   **[Tutorials](tutorials/README.md)**

    Eighteen workflows, each side by side with its R original.

-   **[API reference](api/index.md)**

    Every public function, with the Seurat call it mirrors.

-   **[Fidelity](fidelity.md)**

    How the port is verified, and the differences that are real.

-   **[Changelog](CHANGELOG.md)**

    What has shipped, and what is only on `main`.

</div>

!!! warning "The PyPI release lags `main`"
    The newest release is **0.2.0**. Reference mapping, sketching, `LazyMatrix`,
    cell hashing, Mixscape, `run_spca`/`glm_pca`, pseudobulk DE and the
    MERSCOPE/Visium additions are on `main` but **not** in `pip install shanuz`
    yet. [The changelog](CHANGELOG.md) records exactly what is where; these docs
    are built from `main`.
