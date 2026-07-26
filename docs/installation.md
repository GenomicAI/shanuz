# Installation

**Python 3.12 or newer.** CI tests 3.12 and 3.13.

## From PyPI

```bash
pip install shanuz
```

That gets the core: the object model, preprocessing, PCA and marker detection,
on `numpy`, `scipy`, `pandas` and `packaging` alone. Everything heavier is an
extra, and everything heavier is imported lazily — a base install imports and
runs without matplotlib, scikit-learn or umap-learn anywhere on the system.

| Extra | Adds | You need it for |
|---|---|---|
| `analysis` | statsmodels, scikit-learn, numba, umap-learn, igraph, leidenalg, matplotlib, seaborn, scikit-misc | Clustering, UMAP/t-SNE, every plot, and the `LR`/`negbinom`/`mast` DE tests |
| `anndata` | anndata | `as_anndata` / `from_anndata` |
| `integration` | harmonypy | `run_harmony`, and `integrate_layers(method="harmony")` |
| `deseq2` | pydeseq2 | `find_markers(test_use="deseq2", ...)` |
| `all` | all of the above, plus the dev tooling | Running the test suite |

```bash
pip install "shanuz[analysis]"      # what most analyses want
pip install "shanuz[all]"           # everything
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv pip install "shanuz[analysis]"
```

## From source

The PyPI release is **0.2.0** and `main` is well ahead of it — reference
mapping, sketching, `LazyMatrix`, cell hashing, Mixscape, `run_spca`/`glm_pca`,
pseudobulk DE and the MERSCOPE/Visium work are all source-only right now.
[The changelog](CHANGELOG.md) is the authority on which is which.

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
uv venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -e ".[all]"
```

A source checkout is also what the [tutorials](tutorials/README.md) expect: each
one is a script in `tutorials/` next to the R script it is checked against.

## Why the floor is 3.12

`requires-python` tracks [SPEC 0](https://scientific-python.org/specs/spec-0000/)
— three years past each Python release — because numpy, scipy, pandas and
scikit-learn are what actually constrain this package, and that is the calendar
they keep. It is stricter than CPython's own EOL schedule, which would have held
3.11 until October 2027.

Published releases stay installable on the versions they declared: on 3.10 or
3.11, `pip` resolves to **0.2.0**, the last release with `>=3.10`. Nothing breaks
retroactively.

**3.14 is one package away.** Everything in the dependency set has cp314 wheels
except `harmonypy`, which publishes manylinux wheels only through cp313. Without
one, a 3.14 install either builds it from source — needing BLAS and a
CMake-fetched armadillo — or lets the resolver backtrack to harmonypy 0.2.0,
which depends on torch and drags in the whole CUDA stack. Neither is worth
declaring support for, so the CI matrix stops at 3.13 until that wheel exists.

## R, for the comparisons

Nothing in `shanuz` needs R. The R side is only for reproducing the
[fidelity checks](fidelity.md) yourself — each tutorial ships a `*_verify.R`
that runs the same analysis under Seurat and writes the numbers the Python
script compares against.

```r
install.packages("Seurat")            # 5.5.1 is what the references were taken on
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
```

Four tutorials need more: `harmony` for integration, `MAST` and `DESeq2`
(Bioconductor) for two of the DE tests, and `BPCells` for the out-of-core
comparison.
