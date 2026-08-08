"""Emit ``truecell_guided_tour.ipynb`` — the story-shaped tour of the package.

**The notebook is generated, not hand-edited.** A fix goes here and the notebook is
rebuilt; patching the ``.ipynb`` directly means the next regeneration silently
reverts it. Same arrangement as the ``generate_*.py`` scripts next door, which own
the tutorial figures.

    uv run --with nbformat python tutorials/build_guided_tour.py [OUT.ipynb]

The notebook is written **without outputs** so the committed file stays a readable
diff. To check that it still runs end to end, execute it against real data — the
tour downloads PBMC 3k itself, so nothing needs staging:

    uv run --with nbformat,nbclient,ipykernel python -c "
    import nbformat, nbclient
    nb = nbformat.read('tutorials/truecell_guided_tour.ipynb', as_version=4)
    for c in nb.cells:                      # the install cell is for Colab, not CI
        if c.cell_type == 'code' and '%pip install' in c.source:
            c.source = 'pass'
    nbclient.NotebookClient(nb, timeout=1800).execute()
    print('ok')"

Prose in the notebook quotes numbers the run produces (cluster sizes, marker genes,
score distributions). If you change the pipeline, re-run it and reconcile the text —
several claims in the first draft were wrong in exactly that way.
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf

_DEFAULT_OUT = Path(__file__).parent / "truecell_guided_tour.ipynb"

MD: list = []


def md(text: str) -> None:
    MD.append(("md", text.strip("\n")))


def code(text: str) -> None:
    MD.append(("code", text.strip("\n")))


# ===========================================================================
# 0. Title
# ===========================================================================

md(r"""
<img src="https://raw.githubusercontent.com/GenomicAI/truecell/main/docs/assets/logo/truecell-lockup-1200.png" width="360" alt="truecell">

# A vial of blood, and nobody told you what's in it

### A guided tour of **truecell** — single-cell analysis in Python, ported from Seurat

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GenomicAI/truecell/blob/main/tutorials/truecell_guided_tour.ipynb)
[![PyPI](https://img.shields.io/pypi/v/truecell.svg)](https://pypi.org/project/truecell/)
[![Docs](https://img.shields.io/badge/docs-genomicai.github.io%2Ftruecell-17423a.svg)](https://genomicai.github.io/truecell/)

---

Somebody drew blood, separated the white cells, and pushed 2,700 of them through a
droplet machine one at a time. What came back is a spreadsheet: **32,738 rows of genes,
2,700 columns of cells, and not one label anywhere.**

No cell in that file is marked "T cell". Nothing says which of them fight bacteria and
which remember a vaccine you had as a child. Every one of those answers is *in* the
numbers, and this notebook is the walk from the spreadsheet to the answer.

By the end you will have taken an unlabelled matrix and produced a map of human immune
cells with names on it — and, more usefully, you will know **why each step exists**, which
is the part that transfers to your own data.

**Along the way you'll meet the whole package**: the object model, quality control,
normalization, feature selection, PCA, graph clustering, UMAP, differential expression,
annotation, pseudobulk, gene-signature scoring, and a tour of everything else truecell
does that this one dataset can't show — integration, spatial, CITE-seq, CRISPR screens.

|  |  |
| --- | --- |
| **Runtime** | ~3–5 minutes end to end on a default Colab runtime, install and download included |
| **Data** | Downloads ~24 MB automatically. Nothing to upload. |
| **Requires** | Python 3.12+. Colab's default runtime is fine. |
| **Prior knowledge** | None assumed. If you know Seurat, there's a translation table at the end. |

> **Run it top to bottom.** Every cell depends on the one before it. `Runtime → Run all`
> in Colab, or `Kernel → Restart & Run All` in Jupyter.
""")

# ===========================================================================
# 1. Setup
# ===========================================================================

md(r"""
---

## 0 · Setting up

One install, and the only thing it needs from you is a Python new enough to run it.
""")

code(r"""
import sys

# Worth checking before anything installs. On Python 3.11 or older, pip resolves
# `truecell` to 0.2.0 — the last release that supported them — which predates most
# of this notebook, and the failure would show up later as a confusing ImportError
# rather than here as a clear message.
if sys.version_info < (3, 12):
    raise SystemExit(
        f"truecell needs Python 3.12+; this kernel is {sys.version.split()[0]}. "
        "On Colab, use a standard runtime; locally, upgrade Python."
    )
print(f"Python {sys.version.split()[0]} — good")
""")

code(r"""
# `analysis` pulls in the scientific extras: scikit-learn, umap-learn, igraph,
# leidenalg, statsmodels, matplotlib, seaborn. The bare `pip install truecell`
# gives you the object model and preprocessing but not clustering or plots.
%pip install --quiet "truecell[analysis]"
print("installed")
""")

md(r"""
If Colab asks you to restart the session after that install, do it — then carry on from
the next cell. It's the usual NumPy-already-imported dance and it costs nothing.
""")

code(r"""
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import truecell
from truecell.datasets import pbmc3k

%matplotlib inline
plt.rcParams["figure.dpi"] = 110

# Two warnings from dependencies, hidden by name rather than by blanket-silencing
# everything — you want to see any *other* warning this notebook produces.
#   - umap-learn reports "n_jobs value 1 overridden to 1", which is a no-op it
#     announces whenever you pass a random_state, which we do for reproducibility.
#   - numba emits deprecation chatter on some Colab images.
warnings.filterwarnings("ignore", message=".*n_jobs value.*")
warnings.filterwarnings("ignore", category=FutureWarning, module="numba")

print("truecell", truecell.__version__)
print("python  ", sys.version.split()[0])
""")

md(r"""
### One helper, so the plots behave

truecell's plotting functions **return a matplotlib `Figure`** rather than drawing to a
global canvas. That's deliberate — you can save it, embed it in a multi-panel layout, or
restyle it after the fact. In a notebook it means we display it explicitly and then close
it, so nothing renders twice and long notebooks don't accumulate open figures.
""")

code(r"""
from IPython.display import display

def show(fig):
    "Display a truecell figure once, then release it."
    display(fig)
    plt.close(fig)
""")

# ===========================================================================
# 2. The specimen
# ===========================================================================

md(r"""
---

## 1 · The specimen

This is the PBMC 3k dataset from 10x Genomics — **peripheral blood mononuclear cells**
from one healthy donor. It is the fruit fly of single-cell genomics: small, public,
and analysed so many times that we know what the right answer looks like, which is
exactly what you want when you're learning the machinery.

The download is cached in `~/.truecell_data/`, so re-running this notebook won't
re-fetch it.
""")

code(r"""
counts, genes, cells = pbmc3k()

print(f"matrix : {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")
print(f"stored : {counts.nnz:,} non-zero values")
print(f"density: {100 * counts.nnz / (counts.shape[0] * counts.shape[1]):.2f}% "
      "of the matrix is non-zero")
print(f"\nfirst genes : {genes[:4]}")
print(f"first cells : {cells[:2]}")
""")

md(r"""
### That number deserves a second look

**Over 97% of this matrix is zero.** A beginner's instinct is that something went
wrong. It didn't — two true things stack up:

1. **Cells really don't express most genes.** A platelet has no use for the machinery of
   antibody production, so those genes are genuinely off.
2. **Sampling.** We capture only a small share of the mRNA in each cell, so a gene that
   *is* on but quiet often lands on zero purely by chance. This is the *dropout* problem,
   and it's why single-cell statistics look different from bulk RNA-seq statistics.

Because of that, the matrix is stored **sparse** — only the non-zero values are kept in
memory. truecell keeps it sparse for as long as it possibly can, which is what lets a
laptop handle datasets with hundreds of thousands of cells.

### Building the object

Everything from here lives in a `Truecell` object. If you've used Seurat in R, this is
the `Seurat` S4 class, ported: same slots, same vocabulary, Python naming.
""")

code(r"""
from truecell import create_truecell_object

pbmc = create_truecell_object(
    counts=counts,
    assay="RNA",
    project="pbmc3k",
    feature_names=genes,
    cell_names=cells,
    min_cells=3,        # drop genes seen in fewer than 3 cells
    min_features=200,   # drop barcodes with fewer than 200 genes detected
)

pbmc
""")

md(r"""
Those two thresholds already threw away most of the file. **32,738 genes became 13,714**:
the rest were detected in fewer than three cells, so nothing downstream could have used
them anyway — you cannot compare a gene across groups if only two cells in the entire
experiment ever showed it.

### What's inside the object

Four containers, and every step from here writes into one of them:

| Slot | Holds | Seurat equivalent |
| --- | --- | --- |
| `.assays` | the expression matrices, in **layers** | `obj@assays` |
| `.meta_data` | one row per cell — QC stats, cluster IDs, anything you add | `obj@meta.data` |
| `.reductions` | PCA, UMAP, t-SNE embeddings | `obj@reductions` |
| `.graphs` | cell–cell neighbour graphs | `obj@graphs` |

The **layers** idea is worth pausing on. One assay holds several versions of the same
matrix at once — raw `counts`, log-normalized `data`, z-scored `scale.data` — so a
function can ask for the representation it needs and the raw numbers are never destroyed.
Right now only `counts` exists, because that's all we've done.
""")

code(r"""
rna = pbmc.assays["RNA"]

print("assays     :", list(pbmc.assays))
print("layers     :", rna.layers_list())
print("reductions :", list(pbmc.reductions) or "(none yet)")
print("graphs     :", list(pbmc.graphs) or "(none yet)")

print("\nper-cell metadata, populated automatically at creation:")
pbmc.meta_data.head()
""")

md(r"""
`nCount_RNA` is the total molecules captured in that cell; `nFeature_RNA` is how many
distinct genes were detected in it. Both were computed for you, and both are about to
earn their keep.
""")

# ===========================================================================
# 3. QC
# ===========================================================================

md(r"""
---

## 2 · Triage — which of these are even cells?

Here is the uncomfortable truth about droplet sequencing: **a barcode is not a cell.**
It's a droplet, and droplets go wrong in two directions.

- **Empty or nearly empty droplets** catch ambient RNA floating in the buffer instead of a
  cell. They show very few distinct genes.
- **Doublets** are two cells sharing one droplet. They look like a chimera — a T cell and
  a monocyte fused into a creature that exists nowhere in the body — and they show
  unusually *many* genes.

And a third failure isn't about droplets at all. When a cell is **dying**, its membrane
leaks. Cytoplasmic mRNA escapes, but the mitochondria — which carry their own genome —
stay behind. So the share of reads coming from mitochondrial genes climbs. A high
mitochondrial fraction is the classic signature of a cell that was already dying when we
caught it.

Human mitochondrial genes all start with `MT-`, so that fraction is one pattern match away.
""")

code(r"""
from truecell import percentage_feature_set

percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")

pbmc.meta_data[["nFeature_RNA", "nCount_RNA", "percent.mt"]].describe().round(2)
""")

md(r"""
Summary statistics hide the shape of a distribution, and the shape is the whole point
here. Violin plots show it.
""")

code(r"""
from truecell import vln_plot

show(vln_plot(pbmc, ["nFeature_RNA", "nCount_RNA", "percent.mt"],
              ncol=3, figsize=(12, 4), pt_size=1.0))
""")

md(r"""
Read them left to right:

- **`nFeature_RNA`** — the bulk sits around 800 genes per cell, with a thin tail climbing
  past 3,000. That tail is where doublets live.
- **`percent.mt`** — most cells are under 5% mitochondrial, and a few are far above it.
  Those are the dying ones.

The relationship between two QC measures is often clearer than either alone:
""")

code(r"""
from truecell import feature_scatter

show(feature_scatter(pbmc, "nCount_RNA", "nFeature_RNA", figsize=(6, 5)))
""")

md(r"""
That tight upward curve is what a healthy experiment looks like: cells with more molecules
have more distinct genes detected, flattening off as the common genes saturate. Points
falling *below* the curve — many molecules, few genes — are cells whose reads pile onto a
handful of transcripts, which is what ambient RNA contamination does.

### Drawing the line

Now the judgement call. We keep cells with **200–2,500 detected genes** and **under 5%
mitochondrial reads**.

Be clear about what those numbers are: **defaults from the original tutorial, not laws of
nature.** 5% mitochondrial is reasonable for blood; cardiac and muscle tissue are
legitimately mitochondria-rich and a 5% cutoff there would delete your best cells. The
right move on new data is always to plot first, as we just did, and cut where *your*
distribution says to.
""")

code(r"""
md_df = pbmc.meta_data
keep = (
    (md_df["nFeature_RNA"] > 200)
    & (md_df["nFeature_RNA"] < 2500)
    & (md_df["percent.mt"] < 5)
)

before = len(pbmc)
pbmc = pbmc.subset(cells=list(md_df.index[keep]))

print(f"kept {len(pbmc):,} of {before:,} cells "
      f"({before - len(pbmc)} removed, {100 * len(pbmc) / before:.1f}% retained)")
""")

md(r"""
62 cells gone. `.subset()` returns a **new object** — the original is untouched — and it
carries the metadata, layers and any reductions along with it, so nothing falls out of
alignment.
""")

# ===========================================================================
# 4. Normalization
# ===========================================================================

md(r"""
---

## 3 · Making cells comparable

A problem we have to fix before comparing anything: **cells were not sequenced equally
deeply.** One cell yielded 2,000 molecules, its neighbour 8,000. If you compare raw counts
between them, the deeper cell looks like it expresses more of *everything*, which tells you
about the droplet, not the biology.

`normalize_data` fixes this the way Seurat does, in three moves per cell:

1. divide each gene's count by that cell's total (now it's a proportion),
2. multiply by 10,000 (a readable scale — "counts per 10k"),
3. take `log(x + 1)` — which compresses a handful of screaming-loud genes so they don't
   drown out everything else, while `+1` keeps the zeros at zero.

The result goes into a new layer. **The raw counts stay exactly where they were.**
""")

code(r"""
from truecell import normalize_data

normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)

print("layers now:", pbmc.assays["RNA"].layers_list())

# The same gene in the same cells, before and after.
raw  = pbmc.assays["RNA"].layer_data("counts", features=["LYZ"]).toarray().ravel()
norm = pbmc.assays["RNA"].layer_data("data",   features=["LYZ"]).toarray().ravel()

pd.DataFrame({
    "counts": raw[:6],
    "data (log-normalized)": norm[:6].round(3),
    "cell total (nCount_RNA)": pbmc.meta_data["nCount_RNA"].values[:6],
}, index=pbmc.cell_names()[:6])
""")

md(r"""
Look at the relationship between the three columns. The same raw count maps to a *different*
normalized value depending on how deeply its cell was sequenced — which is precisely the
correction we wanted. `counts` and `data` now sit side by side in the same assay, and every
later function names the layer it wants.
""")

# ===========================================================================
# 5. HVG
# ===========================================================================

md(r"""
---

## 4 · Which genes are worth listening to?

We have 13,714 genes. Most of them are **useless for telling cells apart** — not because
they don't matter biologically, but because every cell expresses them at roughly the same
level. Ribosomal proteins, core metabolism, the housekeeping machinery: essential to life,
uninformative about identity. If you're trying to tell a B cell from a monocyte, a gene
both express identically contributes nothing but noise.

What we want are genes that are **high in some cells and low in others**. But there's a
catch: in count data, highly expressed genes are *automatically* more variable — it's a
property of the Poisson-ish sampling, not of the biology. So the honest question isn't
"which genes vary most" but **"which genes vary more than a gene of their brightness
should"**.

That's what `vst` does: fit the mean–variance trend across all genes, then rank each gene by
how far it sits *above its own expected variance*.
""")

code(r"""
from truecell import find_variable_features

find_variable_features(pbmc, selection_method="vst", nfeatures=2000)

hvg = pbmc.assays["RNA"].variable_features
print(f"selected {len(hvg):,} highly variable genes out of {len(pbmc.feature_names()):,}")
print("\ntop 10:")
for i, g in enumerate(hvg[:10], 1):
    print(f"  {i:2d}. {g}")
""")

code(r"""
from truecell import variable_feature_plot

show(variable_feature_plot(pbmc, label=True, n_label=10, figsize=(9, 5)))
""")

md(r"""
Every dot is a gene: mean expression across, variability up. The black cloud is the
expected trend; the coloured points broke away from it.

**The top of that list is already telling you the answer.** `PPBP` and `PF4` are platelet
genes. `LYZ` and `S100A9` are myeloid. `GNLY` is a cytotoxic granule protein. Before any
clustering, the most variable genes in this blood sample are the ones that distinguish
immune lineages from each other — which is a good sign that the variation we selected is
the variation we came for.

From here, the pipeline works with these 2,000 genes instead of all 13,714. That's an
almost 7× reduction in width, with the discriminating signal kept.
""")

# ===========================================================================
# 6. Scale + PCA
# ===========================================================================

md(r"""
---

## 5 · Compression, and what it compresses into

Two steps here, and the first is preparation for the second.

**Scaling** puts every gene on the same footing: subtract its mean, divide by its standard
deviation. Without it, PCA is dominated by whichever genes happen to be expressed loudest,
and "loud" is not the same as "informative". After scaling, a gene's *pattern* is all that
counts.

We scale all 13,714 genes rather than just the variable ones — slightly slower, but it means
any gene can be drawn on a heatmap later without going back a step.
""")

code(r"""
from truecell import scale_data

scale_data(pbmc, features=pbmc.feature_names())
print("layers now:", pbmc.assays["RNA"].layers_list())
""")

md(r"""
### PCA: 2,000 dimensions down to a few dozen

Each cell is currently a point in 2,000-dimensional space. That's both too big to work with
and mostly redundant — genes travel in packs, and all the B-cell genes rise and fall
together, carrying one piece of information between them, not a hundred.

PCA finds those packs. Each **principal component** is a weighted combination of genes that
co-vary, ordered so PC1 captures the most variation, PC2 the next, and so on.

This is where single-cell PCA becomes genuinely interpretable, and it's worth looking at
rather than treating as a black box.
""")

code(r"""
from truecell import run_pca

run_pca(pbmc, n_pcs=50, features=hvg, reduction_name="pca")

pca = pbmc.reductions["pca"]
print("embeddings :", pca.cell_embeddings.shape, "(cells x PCs)")
print("loadings   :", pca.feature_loadings.shape, "(genes x PCs)")
""")

code(r"""
from truecell import viz_dim_loadings

show(viz_dim_loadings(pbmc, reduction="pca", dims=[1, 2], n_features=15, figsize=(10, 6)))
""")

md(r"""
**Read PC1 as a sentence.** At one end: `CST3`, `LYZ`, `TYROBP`, `FCN1` — the myeloid
programme. At the other: `MALAT1`, `LTB`, `IL32`, `CD3D` — lymphocytes. PC1 is not an
abstract axis. It is *myeloid versus lymphoid*, the deepest split in this sample,
discovered without being told that such a division exists.

PC2 does the same job one level down, separating B cells from T and NK cells.

### How many components are real?

We computed 50. Most are noise. The elbow plot ranks them by how much variation each
explains, and the point where the curve flattens is where the signal runs out.
""")

code(r"""
from truecell import elbow_plot

show(elbow_plot(pbmc, ndims=30, figsize=(7, 4)))
""")

md(r"""
The curve drops steeply and levels off somewhere around PC 9–12. We'll take **10** and move
on.

Two honest notes about that choice. First, it's a soft decision — anything from 8 to 15
gives you substantially the same biology here. Second, err **high** rather than low: an
extra noisy PC is diluted by the ones carrying signal, whereas dropping a real one can
delete a rare cell type from your results entirely.

If you'd rather have a statistical test than an eyeball, truecell ports Seurat's
permutation approach too — `jack_straw()` and `score_jackstraw()` give each PC a p-value.
It costs minutes rather than seconds, which is why the elbow remains the working default.
""")

# ===========================================================================
# 7. Clustering + UMAP
# ===========================================================================

md(r"""
---

## 6 · Finding the groups

Now the central act. Nobody told us how many cell types are in this vial, so we're not going
to tell the algorithm either. Instead:

1. **Connect each cell to its 20 nearest neighbours** in the 10-dimensional PCA space —
   a graph, where an edge means "these two cells look alike".
2. **Refine those edges** by shared neighbours (SNN): two cells that merely happen to be
   close get a weak link; two cells that also share most of *their* neighbours get a strong
   one. This is what makes the graph robust to the noise a raw distance would trip over.
3. **Find communities** in that graph — the Louvain algorithm, which searches for groups
   with many internal edges and few external ones.

Community detection doesn't need to be told how many groups to find. It discovers them.
""")

code(r"""
from truecell import find_neighbors, find_clusters

find_neighbors(pbmc, dims=range(10), k_param=20)
print("graphs:", list(pbmc.graphs))

find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)

sizes = pbmc.meta_data["seurat_clusters"].value_counts().sort_index()
print(f"\nfound {len(sizes)} clusters at resolution 0.5:")
for cid, n in sizes.items():
    print(f"  cluster {cid}: {n:>4} cells")
""")

md(r"""
Eight groups, from 692 cells down to 14. **Clusters are numbered by size**, largest first —
which is why cluster 7 here is that 14-cell sliver. Hold on to it; small clusters are
routinely the interesting ones, and this one turns out to be a genuine cell type rather
than debris.

### `resolution` is a dial, not a truth

That parameter controls how hard Louvain looks for splits. It is the single most consequential
knob in the whole pipeline, and there is no correct setting — only settings appropriate to
the question. Watch what it does:
""")

code(r"""
# Several resolutions in one call; each lands in its own metadata column.
resolutions = [0.1, 0.3, 0.5, 1.0, 1.5]
columns = [f"res_{r}" for r in resolutions]

find_clusters(pbmc, resolution=resolutions, cluster_name=columns,
              algorithm=1, random_seed=0)

for r, col in zip(resolutions, columns):
    print(f"  resolution {r:<4} -> {pbmc.meta_data[col].nunique():>2} clusters")
""")

md(r"""
At 0.1 you get the big lineages. At 1.5 you get subtypes — and eventually you get noise
split off into clusters of its own, because Louvain will always find *something* if you push
it hard enough.

The judgement is biological, not statistical: **turn it up until your clusters stop having
distinct marker genes.** That's the point where you've started splitting on noise.

We named those columns ourselves with `cluster_name=`. Left to itself, truecell names them
the way Seurat does — `RNA_snn_res.0.5` — including R's habit of printing `1.0` as `1`, so
`resolution=1.0` lands in `RNA_snn_res.1`. Ported scripts trip over that, which is why the
convention is reproduced rather than tidied up.

A multi-resolution call leaves the object on the **last** resolution given, so we're
currently sitting on 1.5. Back to 0.5, explicitly:
""")

code(r"""
pbmc.idents = pbmc.meta_data["res_0.5"]
print(f"active identity: {len(pd.unique(pbmc.idents))} clusters")
""")

md(r"""
### Drawing it

The clusters were found in **10-dimensional space** — real, but unviewable. UMAP squashes
those 10 dimensions to 2 for display, trying to keep cells that were neighbours as
neighbours.

One warning, because it's the most common way single-cell figures get over-read: **in a UMAP,
distance between clusters is not meaningful.** Two blobs sitting far apart are not more
different than two sitting close. Within-cluster structure is trustworthy; between-cluster
geometry is an artefact of the layout. Use it to look, not to measure.

This is the slowest cell in the notebook, and most of that isn't the embedding — it's numba
compiling UMAP's inner loops the first time they're called. Subsequent runs are much faster.
""")

code(r"""
from truecell import run_umap, dim_plot

run_umap(pbmc, dims=range(10), reduction_name="umap", seed=42)
show(dim_plot(pbmc, reduction="umap", label=True,
              title="Eight clusters, no names yet", figsize=(8, 6.5)))
""")

md(r"""
There's the map. Eight islands, each a group of cells the data insists belong together —
and every one of them still anonymous.
""")

# ===========================================================================
# 8. Markers
# ===========================================================================

md(r"""
---

## 7 · Naming the strangers

A cluster is a number. Biology needs a name, and the way to get one is to ask: **what does
this group express that the others don't?**

That's a differential expression test, run per gene: cells in cluster *X* against all the
rest. Start with one cluster, to see the shape of the output.
""")

code(r"""
from truecell import find_markers

cluster0 = find_markers(pbmc, ident_1="0", only_pos=True, min_pct=0.25)
print(f"cluster 0 vs everything else — {len(cluster0)} genes pass the thresholds\n")
cluster0.head(8).round(4)
""")

md(r"""
Five columns, and each answers a different question:

| Column | What it tells you |
| --- | --- |
| `p_val` | raw Wilcoxon p-value |
| `avg_log2FC` | how many doublings higher in this cluster — **effect size** |
| `pct.1` | fraction of cells *in* the cluster expressing it |
| `pct.2` | fraction expressing it *elsewhere* |
| `p_val_adj` | Bonferroni-corrected across all genes tested |

### Read that table again, because it's a trap

`find_markers` returns results **sorted by p-value**, as Seurat's does, and look what floats
to the top: `RPS6`, `RPS12`, `RPL32` — ribosomal proteins. Their `pct.1` is 1.00 and their
`pct.2` is 0.99. They are expressed in *every cell in the dataset*, cluster or not.

They are not markers. They're an artefact of having thousands of cells: with samples that
large, a 3% difference is detected with overwhelming confidence, and the p-value collapses
to zero. **Statistical significance is nearly free here; it is not evidence of anything
interesting.**

What makes a marker useful is a large `avg_log2FC` *and* a wide gap between `pct.1` and
`pct.2` — present in nearly all of these cells, and in few of the others. So sort by effect
size:
""")

code(r"""
cluster0.sort_values("avg_log2FC", ascending=False).head(8).round(4)
""")

md(r"""
A completely different list — and now it means something. `CCR7`, `LEF1`, `TCF7` and `MAL`
are the naive-T-cell programme, and every one of them has a `pct.1` three to four times its
`pct.2`. Cluster 0 has declared itself.

Now all eight clusters at once — eight tests across thousands of genes, and still a second
or two, because the marker code filters genes before it densifies anything.
""")

code(r"""
from truecell import find_all_markers

all_markers = find_all_markers(pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
print(f"{len(all_markers):,} cluster-marker pairs\n")

top3 = (all_markers.sort_values("avg_log2FC", ascending=False)
        .groupby("cluster", observed=True).head(3)
        .sort_values(["cluster", "avg_log2FC"], ascending=[True, False]))

for cid, grp in top3.groupby("cluster", observed=True):
    print(f"  cluster {cid}: {', '.join(grp['gene'])}")
""")

md(r"""
### Now do the reading

This is the step no algorithm does for you. You take those gene names to the literature —
or to a reference atlas — and match them to known cell types. For blood the canonical panel
is well established:

| Marker | Marks |
| --- | --- |
| `IL7R`, `CCR7` | CD4 T cells (`CCR7` high = naive, not yet activated) |
| `CD8A` | CD8 cytotoxic T cells |
| `MS4A1`, `CD79A` | B cells |
| `GNLY`, `NKG7` | NK cells |
| `LYZ`, `CD14`, `S100A9` | CD14+ monocytes |
| `FCGR3A`, `MS4A7` | FCGR3A+ (non-classical) monocytes |
| `FCER1A`, `CST3` | dendritic cells |
| `PPBP`, `PF4` | platelets |

A dot plot puts the whole panel against every cluster in one panel — **size** is the
fraction of cells expressing the gene, **colour** is the average level.
""")

code(r"""
from truecell import dot_plot

canonical = ["CCR7", "IL7R", "CD8A", "MS4A1", "CD79A", "GNLY", "NKG7",
             "LYZ", "CD14", "S100A9", "FCGR3A", "MS4A7", "FCER1A", "PPBP", "PF4"]

show(dot_plot(pbmc, canonical, figsize=(11, 5)))
""")

md(r"""
Read down each column and the assignment falls out. `MS4A1` and `CD79A` fire in exactly one
cluster — that's B. `PPBP` and `PF4` fire in exactly one — platelets. `LYZ` fires in two,
one of which also has `FCGR3A`: the two monocyte flavours.

The same information, painted on the map:
""")

code(r"""
from truecell import feature_plot

show(feature_plot(pbmc, ["IL7R", "CD8A", "MS4A1", "GNLY", "LYZ", "FCGR3A"],
                  reduction="umap", ncol=3, pt_size=1.5, figsize=(13, 7)))
""")

code(r"""
from truecell import do_heatmap

top5 = [g for _, grp in all_markers.groupby("cluster", observed=True)
        for g in grp.nlargest(5, "avg_log2FC")["gene"]]
top5 = list(dict.fromkeys(top5))

show(do_heatmap(pbmc, top5, figsize=(13, 8)))
""")

md(r"""
Blocks down the diagonal: each cluster's top genes are high in that cluster and low
everywhere else. That block structure is the visual proof that the clustering found real
groups rather than slicing a continuum at arbitrary points.
""")

# ===========================================================================
# 9. Annotation
# ===========================================================================

md(r"""
---

## 8 · The payoff

Numbers become names.
""")

code(r"""
cell_types = {
    "0": "Naive CD4 T",
    "1": "CD14+ Mono",
    "2": "Memory CD4 T",
    "3": "B",
    "4": "CD8 T",
    "5": "FCGR3A+ Mono",
    "6": "NK",
    "7": "Platelet",
}

pbmc.rename_idents(cell_types)
show(dim_plot(pbmc, reduction="umap", label=True,
              title="PBMC 3k — annotated", figsize=(9, 7)))
""")

md(r"""
**That is the whole journey in one picture.** A 32,738 × 2,700 grid of unlabelled numbers,
turned into a map of the human immune system — every label derived from the data, none of
it supplied by us.

### One caveat, stated plainly

The Seurat tutorial this mirrors ends with **nine** cell types. We have eight, and the
missing one is dendritic cells.

They aren't lost — at resolution 0.5 they sit inside the CD14+ monocyte cluster, which is
biologically reasonable, since DCs and monocytes are close relatives sharing much of the
myeloid programme. Check it yourself: `FCER1A`, the DC marker, is faintly on in the CD14+
Mono cluster and absent everywhere else.

You have two ways to pull them out, and both are legitimate:

- **Raise the resolution** — we already computed `res_1.5`.
- **Subcluster** — take the monocytes alone, `.subset()`, and re-run the pipeline on them.
  This is usually the better move, because a rare population competes for attention against
  the whole dataset at high resolution but has the field to itself in a subset.

Don't take that on trust — the columns are sitting in the metadata, so go and look:
""")

code(r"""
from truecell import average_expression

# Where does FCER1A live once the resolution is turned up?
saved = pbmc.idents                       # keep the annotated labels
pbmc.idents = pbmc.meta_data["res_1.5"]

dc_check = average_expression(pbmc, features=["FCER1A", "CD14", "CST3"]).round(2)
winner = dc_check.loc["FCER1A"].idxmax()

print(f"FCER1A is highest in res-1.5 cluster {winner} "
      f"({(pbmc.meta_data['res_1.5'] == winner).sum()} cells)\n")
print(dc_check[[winner] + [c for c in dc_check.columns if c != winner][:4]].to_string())

pbmc.idents = saved                       # put the names back
""")

md(r"""
There they are. A **32-cell cluster** where `FCER1A` reads ~15 against ≤0.3 in every other
cluster in the dataset — dendritic cells, cleanly separated, and the same 32 cells Seurat
finds. They were never missing; they were just below the resolution we chose to look at.

We're leaving the main annotation at eight on purpose. **A tutorial that quietly matched a target number
would be teaching you the wrong lesson.** Cluster counts move with parameters, with package
versions, with the random seed. What's stable is the marker evidence, and that's what you
should trust — the identity of a cluster is the genes it expresses, not its position in a
list.

### One more thing worth knowing about `rename_idents`

It maps **positionally**. Hand it a dictionary whose length doesn't match the number of
clusters and it won't error — it will silently shift every label from the mismatch onward,
and your figure will be confidently, invisibly wrong.

This is not hypothetical. Exactly that bug shipped in this project: a nine-entry map applied
to eight clusters, captioning the platelets as "DC" in a published figure. It's the reason
truecell's own tutorial now has a test asserting the map's keys equal the clusters the
pipeline actually produced. **If you take one operational habit from this notebook: after
renaming, look at the plot and check a marker.**
""")

code(r"""
# Verify rather than trust: does each label lead on the gene it should?
check = {"Naive CD4 T": "CCR7", "CD14+ Mono": "CD14", "Memory CD4 T": "IL7R",
         "B": "MS4A1", "CD8 T": "CD8A", "FCGR3A+ Mono": "FCGR3A",
         "NK": "GNLY", "Platelet": "PPBP"}

avg = average_expression(pbmc, features=list(check.values()))

for label, marker in check.items():
    winner = avg.loc[marker].idxmax()
    flag = "ok " if winner == label else "!! "
    print(f"  {flag}{label:<14} {marker:<8} highest in: {winner}")
""")

md(r"""
Every label leads on its own marker. *Now* the figure can be believed.
""")

# ===========================================================================
# 10. Beyond the tour
# ===========================================================================

md(r"""
---

## 9 · What else is in the box

The pipeline above is the backbone. Here's a quick pass through the rest of what you
actually reach for — all of it running live on the object we just built.

### Provenance — what was done to this object

The preprocessing and reduction steps record themselves, with their parameters, so an object
you were handed six months ago can still tell you how it was made. R's `Command(obj)`, ported.
""")

code(r"""
for cmd in pbmc.commands:
    print(f"  {cmd.name:<24} {cmd.call_string}")
""")

md(r"""
Note what *isn't* in that list: `find_clusters` and `run_umap` don't log themselves yet, though
Seurat's equivalents do. Worth knowing before you rely on the log as a complete record —
the resolution you clustered at is in the metadata column name, not in the command history.
""")

md(r"""
### Pseudobulk

Sum counts within each cell type and you get one profile per group — the format bulk RNA-seq
tools expect, and the statistically correct unit when you have biological replicates. (With
multiple donors, DE *between* pseudobulk samples is far better calibrated than DE between
thousands of non-independent cells.)
""")

code(r"""
from truecell import aggregate_expression

pseudobulk = aggregate_expression(pbmc, group_by="ident", layer="counts")
print(f"pseudobulk: {pseudobulk.shape[0]:,} genes x {pseudobulk.shape[1]} cell types\n")
pseudobulk.loc[["LYZ", "MS4A1", "CD8A", "PPBP"]].round(0)
""")

md(r"""
### Scoring a gene programme

Often the question isn't "which cluster is this" but "how much of *this behaviour* is each
cell doing". `add_module_score` averages a gene set per cell and — importantly — subtracts a
matched control set drawn from genes of similar expression, so a cell doesn't score high just
for being deeply sequenced.

Here's a cytotoxicity signature. It should light up NK and CD8 T cells and nothing else.
""")

code(r"""
from truecell import add_module_score

add_module_score(pbmc, features={"Cytotoxic": ["GNLY", "NKG7", "GZMB", "PRF1", "KLRD1"]})

show(feature_plot(pbmc, "Cytotoxic", reduction="umap", pt_size=1.5, figsize=(6.5, 5.5)))
""")

md(r"""
### Cell-cycle phase — and how to not be fooled by it

The same machinery with a purpose-built gene set: score the S and G2/M programmes, assign a
phase. truecell ships the Tirosh et al. human gene lists as `CC_GENES`.

Resting blood cells shouldn't be proliferating, so this is effectively a negative control.
Watch what it reports.
""")

code(r"""
from truecell import cell_cycle_scoring

cell_cycle_scoring(pbmc)

print(pbmc.meta_data["Phase"].value_counts().to_string())
print()
print(pbmc.meta_data[["S.Score", "G2M.Score"]].describe().loc[
    ["mean", "std", "50%", "max"]].round(3).to_string())
""")

md(r"""
**Taken at face value, that says over half this blood sample is actively replicating.** For
resting PBMCs that is nonsense — and the second table shows why.

Both scores have a median near zero (about −0.01 and −0.04) and a standard deviation around
0.05. Essentially every cell is sitting at zero, which is the correct answer. But `Phase` is
assigned by taking whichever score is *larger*, and when both are noise centred on zero, that
comparison is close to a coin flip. The label is confident; the underlying signal is absent.

**So check the score distribution before you believe the phase.** A genuinely cycling
population looks nothing like this — you'd see a long right tail and a clearly bimodal spread,
which the `max` column hints at here (a handful of cells do reach 0.6–0.8; those few are real).

This matters most on tumours and cell lines, where proliferation is real and strong enough
that cells cluster by *cell cycle* instead of by cell type. That's what
`scale_data(vars_to_regress=["S.Score", "G2M.Score"])` is for.

### The DE test menu

`find_markers` is one function with a choice of statistics behind it. The default `wilcox` is
rank-based and assumption-light. Others suit other questions — and it's instructive to see
how much the answer depends on the test.
""")

code(r"""
results = {t: find_markers(pbmc, ident_1="NK", ident_2="CD8 T", test_use=t, only_pos=True)
           for t in ["wilcox", "t", "bimod", "roc"]}

# Rank each test by its own verdict: p-value for the three that produce one,
# AUC for roc, which doesn't.
pd.DataFrame({
    t: list(r.sort_values("myAUC" if t == "roc" else "p_val",
                          ascending=(t != "roc")).index[:6])
    for t, r in results.items()
})
""")

md(r"""
NK versus CD8 T is a genuinely hard call — closely related killers — and five of the top six
are the same across all four tests, in slightly different orders. That's the usual outcome,
and it's worth knowing: **the choice of statistic matters far less than people expect.**
`wilcox` is the
default because it's rank-based and makes almost no assumptions, not because it's uniquely
correct.

### What *does* change the answer is the column you sort by

Here's a subtlety with real consequences. `avg_log2FC` is computed from the expression values
themselves — **it doesn't depend on the test at all.** Only the p-value does. So ranking by
fold change makes every test produce an identical list, which looks like agreement but is
really just a fact about that column.

And it's a *different* list:
""")

code(r"""
roc = results["roc"]
pd.DataFrame({
    "by avg_log2FC": list(roc.sort_values("avg_log2FC", ascending=False).index[:6]),
    "by myAUC":      list(roc.sort_values("myAUC", ascending=False).index[:6]),
})
""")

md(r"""
**Same test, same cells, two different answers** — and neither is wrong.

Fold change rewards the biggest *relative* jump, which favours lowly-expressed genes where a
small absolute difference is a large ratio. AUC rewards genes that actually **classify**
(1.0 = perfect separator, 0.5 = useless), which favours abundant, reliable ones — `GZMB`,
`GNLY`, `PRF1`, `NKG7`, the genes you'd put on a flow-cytometry panel.

Want a marker to name a cluster? Take the AUC list. Want to know what's biologically different
between two populations? The fold-change list holds things AUC buries. The lesson generalises
well past this function: **when a tool ranks things for you, find out what it ranked by.**

Also available: `LR` (with covariates via `latent_vars`), `negbinom`, `mast`, and `deseq2` for
pseudobulk DE with `sample_col=`.

### Pulling values out for your own analysis

`fetch_data` assembles genes, metadata and embeddings into one tidy frame — the escape hatch
into pandas, seaborn, scikit-learn, or whatever else you use.
""")

code(r"""
frame = pbmc.fetch_data(["LYZ", "MS4A1", "percent.mt", "umap"])
frame.insert(0, "cell_type", list(pbmc.idents))
frame.head().round(3)
""")

md(r"""
### Styling every plot at once

`set_theme` changes the look of everything the plotting module draws, rather than passing the
same overrides to each call. `theme_context` scopes it to a block.
""")

code(r"""
from truecell import theme_context

with theme_context(base_size=13, style="minimal"):
    show(dim_plot(pbmc, reduction="umap", label=True,
                  title="same plot, different theme", figsize=(8, 6.5)))
""")

md(r"""
### Handing off to AnnData

If the rest of your stack is scanpy, the bridge is one call each way.
""")

code(r"""
try:
    from truecell.compat import as_anndata
    adata = as_anndata(pbmc)
    print(adata)
except ImportError:
    print("anndata not installed — `pip install truecell[anndata]` to enable this")
""")

# ===========================================================================
# 11. Catalogue
# ===========================================================================

md(r"""
---

## 10 · Everything this one dataset can't show

PBMC 3k is a single sample, one modality, no tissue coordinates. Most of truecell's surface
needs data this notebook doesn't have. Here's the map of it — the code is real, and each row
links to a worked tutorial.

### Batch correction and integration

Two experiments run on different days will separate by *day* before they separate by cell
type. Integration removes that while keeping the biology.

```python
from truecell import run_harmony, integrate_layers

run_harmony(obj, group_by_vars="batch")                    # fast, embedding-level
integrate_layers(obj, method="cca")                        # Seurat's anchor-based CCA
integrate_layers(obj, method="rpca")                       # faster, more conservative
```

→ [Integration tutorial](https://genomicai.github.io/truecell/tutorials/integration_vignette/) —
Harmony, CCA and RPCA all reach batch mixing 0.991 on the ifnb IFN-β benchmark, above
Seurat's own 0.917.

### Reference mapping — annotate by borrowing labels

Instead of naming clusters by hand, project your cells onto an annotated reference and
transfer its labels.

```python
from truecell import find_transfer_anchors, transfer_data, map_query

anchors = find_transfer_anchors(reference=ref, query=query, dims=range(30))
labels  = transfer_data(anchors, refdata=ref.meta_data["celltype"])
map_query(anchors, query=query, reference=ref)             # + place it in the reference UMAP
```

→ [Reference mapping tutorial](https://genomicai.github.io/truecell/tutorials/refmap_vignette/) —
98.71% per-cell concordant with R Seurat on the panc8 cross-technology benchmark.

### CITE-seq — RNA and surface protein together

Antibody-derived tags measure protein on the same cells. WNN learns, per cell, how much to
trust each modality — RNA for a cell whose transcriptome is decisive, protein for one whose
surface markers are.

```python
from truecell import find_multi_modal_neighbors

find_multi_modal_neighbors(obj, reduction_list=["pca", "apca"],
                           dims_list=[range(30), range(18)])
```

→ [CITE-seq tutorial](https://genomicai.github.io/truecell/tutorials/multimodal_citeseq/)

### Cell hashing — untangling pooled samples

Pool samples, tag each with a hashtag antibody, sequence together, then demultiplex —
including the doublets, which are visible precisely because they carry two tags.

```python
from truecell import hto_demux, multiseq_demux

hto_demux(obj, assay="HTO")        # Seurat's HTODemux
multiseq_demux(obj, assay="HTO")   # the MULTI-seq alternative
```

→ [Hashing tutorial](https://genomicai.github.io/truecell/tutorials/hashing_vignette/) — 99.81%
call-concordant with R Seurat.

### Pooled CRISPR screens — Mixscape

In a CRISPR screen not every cell carrying a guide is actually perturbed. Mixscape separates
true knockouts from escapers.

```python
from truecell import calc_perturb_sig, run_mixscape, mixscape_lda

calc_perturb_sig(obj, assay="RNA", nt_cell_class="NT")
run_mixscape(obj)                  # -> mixscape_class: "<gene> KO" / NP / NT
mixscape_lda(obj)
```

→ [Mixscape tutorial](https://genomicai.github.io/truecell/tutorials/mixscape_vignette/) — 97.45%
per-cell concordant on the THP-1 ECCITE-seq screen.

### Spatial transcriptomics

Xenium, Visium, CosMx and MERSCOPE, with the cells' physical positions kept — so you can ask
which cell types sit next to which, and which genes vary across the tissue.

```python
from truecell import (load_xenium, load_visium, build_niche_assay,
                      find_spatially_variable_features, image_dim_plot)

obj = load_xenium("path/to/xenium/")
build_niche_assay(obj, fov="fov", group_by="celltype")
find_spatially_variable_features(obj, method="moransi")
```

→ [Xenium](https://genomicai.github.io/truecell/tutorials/xenium_spatial_tutorial/) ·
[Visium](https://genomicai.github.io/truecell/tutorials/visium_vignette/) — verified to 8 significant
figures against R Seurat.

### SCTransform — normalization as a model

An alternative to the log-normalize/scale pair: fit a regularized negative binomial per gene
and use its Pearson residuals. Usually better for datasets with wide depth variation.

```python
from truecell import sctransform

sctransform(obj, vst_flavor="v2")
```

→ [SCTransform tutorial](https://genomicai.github.io/truecell/tutorials/sctransform_vignette/)

### Millions of cells

Two strategies, and they compose.

```python
from truecell import sketch_data, project_data, write_lazy_matrix, open_lazy_matrix

sketch_data(obj, ncells=50_000)     # leverage-weighted subset — rare states kept, not lost
project_data(obj)                   # extend the sketch's PCA/UMAP/labels back to every cell

write_lazy_matrix(counts, "store/")             # memory-mapped, BPCells-style
obj.assays["RNA"].layers["counts"] = open_lazy_matrix("store/")
```

→ [Sketching](https://genomicai.github.io/truecell/tutorials/sketch_vignette/) ·
[Out-of-core](https://genomicai.github.io/truecell/tutorials/lazy_vignette/)

### And the rest

`run_tsne`, `run_ica`, `run_spca` (supervised PCA off a cell graph), `glm_pca` (Poisson/NB PCA
straight on counts), `jack_straw` (PC significance), `find_conserved_markers` (markers holding
across conditions), `composition_test` (does cell-type proportion shift between groups),
`ridge_plot`, `dim_heatmap`, `image_feature_plot`.

Full list: **[API reference](https://genomicai.github.io/truecell/api/)**.
""")

# ===========================================================================
# 12. Seurat translation
# ===========================================================================

md(r"""
---

## 11 · Coming from Seurat?

The port keeps the algorithms and the vocabulary, and changes only the naming convention —
R's `CamelCase` becomes Python's `snake_case`, and slots become attributes.

| Seurat (R) | truecell (Python) |
| --- | --- |
| `CreateSeuratObject()` | `create_truecell_object()` |
| `PercentageFeatureSet()` | `percentage_feature_set()` |
| `NormalizeData()` | `normalize_data()` |
| `FindVariableFeatures()` | `find_variable_features()` |
| `ScaleData()` | `scale_data()` |
| `SCTransform()` | `sctransform()` |
| `RunPCA()` / `RunUMAP()` / `RunTSNE()` | `run_pca()` / `run_umap()` / `run_tsne()` |
| `FindNeighbors()` / `FindClusters()` | `find_neighbors()` / `find_clusters()` |
| `FindMarkers()` / `FindAllMarkers()` | `find_markers()` / `find_all_markers()` |
| `AddModuleScore()` / `CellCycleScoring()` | `add_module_score()` / `cell_cycle_scoring()` |
| `AggregateExpression()` / `AverageExpression()` | `aggregate_expression()` / `average_expression()` |
| `FindIntegrationAnchors()` / `IntegrateData()` | `find_integration_anchors()` / `integrate_data()` |
| `FindTransferAnchors()` / `TransferData()` | `find_transfer_anchors()` / `transfer_data()` |
| `HTODemux()` / `MULTIseqDemux()` | `hto_demux()` / `multiseq_demux()` |
| `RunMixscape()` / `MixscapeLDA()` | `run_mixscape()` / `mixscape_lda()` |
| `DimPlot()` / `FeaturePlot()` / `VlnPlot()` / `DotPlot()` | `dim_plot()` / `feature_plot()` / `vln_plot()` / `dot_plot()` |
| `DoHeatmap()` / `ElbowPlot()` / `RidgePlot()` | `do_heatmap()` / `elbow_plot()` / `ridge_plot()` |
| `Idents(obj) <- ...` | `obj.idents = ...` |
| `RenameIdents()` | `obj.rename_idents()` |
| `subset(obj, ...)` | `obj.subset(...)` |
| `obj@meta.data` | `obj.meta_data` |
| `obj[["RNA"]]` | `obj.assays["RNA"]` |
| `Embeddings(obj, "umap")` | `obj.reductions["umap"].cell_embeddings` |
| `FetchData()` | `obj.fetch_data()` |

Two differences that will bite if you don't know them:

1. **Indexing is 0-based.** R's `dims = 1:10` is `dims=range(10)`.
2. **Most pipeline functions modify in place and return `None`.** Write
   `run_pca(obj)`, not `obj = run_pca(obj)`. The exceptions are `subset()` and the
   plotting functions, which return something.

### Is it actually the same answer?

Fair question, and the project treats it as the central one. Each tutorial has a matching R
script that runs the same analysis in Seurat, and the two are compared numerically rather
than by eye. Deterministic steps agree to floating-point precision; steps with an RNG
(Louvain, UMAP, the module-score control sets) agree within measured, declared bands.

Where they disagreed, the cause has more than once turned out to be a defect on the
**truecell** side — found by these comparisons and fixed — and at least once a defect in
Seurat itself.

→ [How the port is verified](https://genomicai.github.io/truecell/fidelity/)
""")

# ===========================================================================
# 13. Close
# ===========================================================================

md(r"""
---

## Where to go next

You started with an unlabelled matrix and finished with a named map of the human immune
system. The same eight steps — QC, normalize, select, scale, reduce, cluster, mark, name —
are the backbone of nearly every scRNA-seq analysis you'll run.

**Try this next, on this notebook, before you leave it:**

- Re-run from §6 with `resolution=1.0`. How many clusters? Do the extra ones have their own
  markers, or did you just split a real one in half?
- Subset the CD14+ monocytes and re-run the pipeline on them alone. The dendritic cells are
  in there.
- Change `dims=range(10)` to `dims=range(5)`. Which cell type disappears first?

**Then:**

- 📖 [Documentation](https://genomicai.github.io/truecell/) — API reference and all
  eighteen tutorials
- 🔬 [Fidelity report](https://genomicai.github.io/truecell/fidelity/) — how the port is
  checked against R Seurat
- 💻 [GitHub](https://github.com/GenomicAI/truecell) — issues and contributions welcome
- 📦 [PyPI](https://pypi.org/project/truecell/)

```bash
pip install "truecell[analysis]"
```

---

<sub>truecell is MIT-licensed. It ports the methods of Seurat
(Hao et al., *Nature Biotechnology* 2024; Stuart et al., *Cell* 2019) — if you use it in
published work, please cite the Seurat papers for the methods alongside truecell for the
implementation. The PBMC 3k dataset is provided by 10x Genomics.</sub>
""")


# ===========================================================================
# Assemble
# ===========================================================================

def build() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(src) if kind == "md" else nbf.v4.new_code_cell(src)
        for kind, src in MD
    ]
    # nbformat 4.5 stamps every cell with a random uuid, so two runs over identical
    # content produce two different files and `git diff` reports all 89 cells as
    # changed. Numbering them by position makes regeneration byte-stable, which is
    # what lets a real edit stand out in review.
    for i, cell in enumerate(nb.cells):
        cell["id"] = f"cell-{i:03d}"
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
        "colab": {"provenance": [], "toc_visible": True},
    }
    return nb


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_OUT
    nb = build()
    nbf.write(nb, str(out))
    n_code = sum(1 for k, _ in MD if k == "code")
    print(f"wrote {out}")
    print(f"  {len(MD)} cells — {n_code} code, {len(MD) - n_code} markdown")
