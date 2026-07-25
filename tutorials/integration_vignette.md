# Batch integration — Harmony, CCA & RPCA (R Seurat vs Shanuz)

A side-by-side port of Seurat's [integration vignette](https://satijalab.org/seurat/articles/integration_introduction)
on the Kang et al. 2018 PBMC dataset (`ifnb`): ~14,000 human blood cells, half
left resting (**CTRL**) and half stimulated for six hours with interferon-β
(**STIM**). Interferon drives a strong, near-global transcriptional response, so
without correction the cells split first by *condition* and only then by *cell
type* — the textbook batch effect.

The task integration solves: **make the same cell type from the two conditions
overlap, without erasing the biology that separates cell types.** Shanuz ships
three integration paths, and this walkthrough runs all three against their Seurat
references on identical counts:

- **`run_harmony`** ↔ `RunHarmony` / `HarmonyIntegration` — iteratively nudges the
  PCA embedding so batches mix while clusters hold.
- **`integrate_layers(method="cca")`** ↔ `FindIntegrationAnchors(reduction="cca")`
  + `IntegrateData` — anchors are mutual nearest neighbours in a shared
  canonical-correlation space.
- **`integrate_layers(method="rpca")`** ↔ `RPCAIntegration` — the reciprocal-PCA
  variant: each dataset is projected into the other's PCA space before the
  mutual-nearest-neighbour search.

> **Why this tutorial exists.** Every integration function landed in v0.2.0 and
> had only ever been checked against synthetic fixtures with *balanced* batches.
> This is the first time they meet a real dataset with a Seurat reference and
> unequal batch sizes (CTRL 6,548 vs STIM 7,451). Integration embeddings are
> *not* expected to be coordinate-identical — harmonypy and R's harmony are
> separate implementations — so the target is *do the two tools recover the same
> structure*: the same collapse in batch separation, the same recovery of cell
> type, cluster partitions that agree. **This tutorial is also the first in the
> initiative to find real defects** — see the concordance section.

---

## The data — ifnb, through the export bridge

`ifnb` is a curated SeuratData object with no clean raw source, so both languages
read the **same counts**, exported once from R by `tutorials/export_seuratdata.R`
into a 10x-style matrix folder. That guarantees byte-identical input and cell
order — the same discipline every other tutorial gets from a shared GEO download.

```
14,053 genes × 13,999 cells   ·   CTRL 6,548 / STIM 7,451   ·   13 cell types
```

The `stim` column is the batch to remove; `seurat_annotations` is the cell-type
label to preserve — both ship with the dataset, so the two tools start from the
identical annotated state.

---

## Step 1 · Load and prep to PCA — one shared variable-feature basis

Standard prep is the uncorrected starting point every method shares. One wrinkle
makes the cross-tool comparison fair, exactly as in the Mixscape tutorial: **both
tools use the same variable features.** The Python run writes the 2,000 HVGs it
selected to `figures_integration/hvg_features.txt`, and the R script reads them
back — so the only divergences left are the genuinely method-level ones (PCA
numerics, the integration algorithms, Louvain ties).

<table>
<tr><th>R (Seurat)</th><th>Python (Shanuz)</th></tr>
<tr><td>

```r
# same exported counts Python reads
counts <- Read10X("~/.shanuz_data/ifnb")
obj <- CreateSeuratObject(counts, min.cells = 3,
                          meta.data = meta)

hvg <- readLines("figures_integration/hvg_features.txt")
obj <- NormalizeData(obj, verbose = FALSE)
VariableFeatures(obj) <- hvg          # Python's HVGs
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 30,
              verbose = FALSE)
```

</td><td>

```python
from shanuz.datasets import ifnb
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data)
from shanuz.reduction import run_pca

counts, genes, cells, meta = ifnb()
obj = create_shanuz_object(counts=counts, assay="RNA",
        min_cells=3, feature_names=genes,
        cell_names=cells, meta_data=meta)

normalize_data(obj, assay="RNA")
find_variable_features(obj, assay="RNA", nfeatures=2000)
scale_data(obj, assay="RNA")
run_pca(obj, assay="RNA", n_pcs=30)   # writes hvg file
```

</td></tr>
</table>

---

## Step 2 · Integrate — three ways

Harmony corrects the existing PCA in place; CCA and RPCA split the object by
condition, anchor the two batches, and rebuild a corrected reduction. Each result
is stored under its own name and clustered identically (Louvain, resolution 0.5,
neighbours on 30 dims) so the comparison is like-for-like.

<table>
<tr><th>R (Seurat)</th><th>Python (Shanuz)</th></tr>
<tr><td>

```r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$stim)

obj <- IntegrateLayers(obj, method = HarmonyIntegration,
        orig.reduction = "pca", new.reduction = "harmony")
obj <- IntegrateLayers(obj, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = "cca")
obj <- IntegrateLayers(obj, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = "rpca")
```

</td><td>

```python
from shanuz.integration import run_harmony, integrate_layers

run_harmony(obj, group_by="stim", reduction="pca",
            reduction_name="harmony")
integrate_layers(obj, method="cca", group_by="stim",
                 new_reduction="cca")
integrate_layers(obj, method="rpca", group_by="stim",
                 new_reduction="rpca")
```

</td></tr>
</table>

---

## Step 3 · Score the integration

Two rotation-invariant summaries tell the story. **Batch separation** (silhouette
by `stim`, *lower is better* — a good integration mixes the conditions) and
**cell-type preservation** (silhouette by cell type, and the adjusted Rand index
of the clusters against the known annotations, *higher is better*). The Shanuz
scoreboard:

| method | sil_batch ↓ | sil_celltype ↑ | n_clusters | ARI→celltype ↑ | batch-mix ↑ |
|--------|---:|---:|---:|---:|---:|
| uncorrected (PCA) | 0.107 | 0.141 | 16 | 0.519 | 0.161 |
| **Harmony** | **0.008** | 0.194 | 12 | **0.911** | **0.991** |
| **CCA** | **0.004** | 0.231 | 14 | **0.918** | **0.991** |
| **RPCA** | **0.005** | 0.220 | 14 | **0.922** | **0.991** |

Uncorrected, the cells separate by condition (batch-mix 0.161 — clusters are
nearly single-condition). All three methods collapse that while *raising*
cell-type recovery, and now land within a point of each other (batch-mix
0.991 across the board). Getting RPCA here took real work — see the
concordance section.

<table>
<tr><th>R — uncorrected, by condition</th><th>Shanuz — uncorrected, by condition</th></tr>
<tr>
<td><img src="figures_integration/r_01_uncorrected_stim.png" width="420"/></td>
<td><img src="figures_integration/py_01_uncorrected_stim.png" width="420"/></td>
</tr>
</table>

The two conditions form two clouds — the interferon shift dominates the embedding.
After Harmony they interleave, while the cell types stay apart:

<table>
<tr><th>R — Harmony, by condition</th><th>Shanuz — Harmony, by condition</th></tr>
<tr>
<td><img src="figures_integration/r_02_harmony_stim.png" width="420"/></td>
<td><img src="figures_integration/py_02_harmony_stim.png" width="420"/></td>
</tr>
</table>

<table>
<tr><th>R — Harmony, by cell type</th><th>Shanuz — Harmony, by cell type</th></tr>
<tr>
<td><img src="figures_integration/r_03_harmony_celltype.png" width="420"/></td>
<td><img src="figures_integration/py_03_harmony_celltype.png" width="420"/></td>
</tr>
</table>

---

## The headline · R-vs-Python concordance, and four RPCA bugs (all fixed)

Because integration embeddings are not coordinate-comparable across tools, the
concordance is **partition-based**: the adjusted Rand index between the two tools'
clusterings (`ARI(py,R)`, 1 = identical), each tool's own cell-type recovery
(`ARI→type` — the biological check, the two columns should track), and each tool's
batch mixing (`mix`, 1 = fully mixed). All are computed from the cluster labels
`report_concordance()` reads out of the verify script's `r_calls.csv`.

| method | ARI(py,R) | py ARI→type | R ARI→type | py mix | R mix |
|--------|---:|---:|---:|---:|---:|
| PCA (baseline) | 0.970 | 0.519 | 0.518 | 0.161 | 0.163 |
| **Harmony** | 0.935 | 0.911 | 0.931 | **0.991** | **0.991** |
| **CCA** | 0.972 | 0.918 | 0.927 | **0.991** | **0.991** |
| **RPCA** | 0.774 | **0.922** | 0.736 | **0.991** | 0.917 |

**Harmony and CCA match Seurat closely on every axis**, partition agreement
included. This is the confirmation the initiative was built to get: shanuz's
two most-used integration paths reproduce Seurat's result on the standard
benchmark, cluster-for-cluster.

**RPCA took four bugs to get here, in two rounds.** The first two were caught
by an earlier pass of this tutorial and are described in
[the anchor-internals vignette](anchors_vignette.md): a crash on unequal batch
sizes (#41), and an under-integration where shanuz scaled batches globally
instead of per-object and searched the raw reciprocal projection instead of
Seurat's SD-standardized, L2-normalized one (#42). Those took RPCA's batch
mixing from 0.222 to **0.867** and cell-type recovery from 0.444 (below the
uncorrected baseline) to **0.677** — real fixes, but still short of Seurat's
0.914 / 0.735. At the time that residual read as "the expected implementation
gap" — exact vs. approximate neighbours, scikit-learn vs. irlba PCA.

**It wasn't a gap. It was two more bugs**, found by
[T-int](anchors_vignette.md) treating "the expected implementation gap" as a
claim to verify rather than a place to stop:

1. **`integrate_layers` was running the wrong algorithm.** Seurat v5's
   `IntegrateLayers(method = CCAIntegration/RPCAIntegration)` does not call
   `IntegrateData` — it calls **`IntegrateEmbeddings`**, which corrects the
   *PCA embedding itself* by transposing it into a fake assay and running the
   same anchor machinery over it. shanuz was running the v4 workflow instead
   (correct expression, re-scale, re-run PCA) — a different object with the
   same shape, which is why the two agreed on only 1 of 30 output dimensions
   above \|r\| = 0.99 on a matched-size probe.
2. **`run_pca` used sklearn's randomized SVD.** It switches solvers once
   `max(shape) > 500`, accurate in the leading components and drifting in the
   trailing ones — only 15 of 30 PCs matched Seurat's irlba above \|r\| = 0.99,
   with one PC down at 0.006. Harmless when only the leading PCs matter
   downstream, but `IntegrateEmbeddings` corrects the embedding directly, so
   the drift landed straight in the output. Swapping in ARPACK (deterministic,
   6× faster on this data, and exact enough to match irlba on all 30 PCs)
   fixed it without slowing the normal path.

Fixing both took embedding agreement from 1/30 to **30/30 PCs above \|r\| = 0.99**
on the full 13,999-cell, unequal-batch dataset (STIM 7,451 is the reference,
Seurat's own `PairwiseIntegrateReference` rule), reference-half cells copied
through at **exactly** zero difference. RPCA's batch mixing rose to **0.991**
— now *higher* than Seurat's own 0.917 — and cell-type recovery to **0.922**,
above Seurat's 0.736.

**What that leaves is not an integration gap — it's a clustering one.**
`ARI(py,R)` for RPCA is still only 0.774, which looks like a leftover
disagreement, but it isn't upstream of `find_clusters`: clustering **Seurat's
own** RPCA embedding with shanuz's `find_neighbors` + `find_clusters` gives
batch-mix 0.990 and ARI→type 0.920 — almost identical to shanuz's own
end-to-end numbers, and nothing like Seurat's 0.917 / 0.736 on that same
embedding. The embeddings agree; the two tools' Louvain implementations, given
an identical input, do not. That's now the open question, not integration.

---

## Running it yourself

```bash
Rscript tutorials/export_seuratdata.R ifnb        # one-time counts export (~394 MB SeuratData pkg)
python  tutorials/ifnb_integration_tutorial.py    # writes HVGs, prints the scoreboard
Rscript tutorials/ifnb_integration_verify.R       # Seurat reference → r_calls.csv + r_*.png
python  tutorials/ifnb_integration_tutorial.py    # re-run: now prints the R-vs-Python concordance
python  tutorials/generate_integration_plots.py   # Shanuz figures → figures_integration/py_*.png
```

The R reference needs the `harmony` package and enough headroom for Seurat v5's
parallel integration (`options(future.globals.maxSize = 3 * 1024^3)`, set in the
script).

**Figures** (`tutorials/figures_integration/`, `r_*` = R Seurat, `py_*` = Shanuz):

| Figure | Description |
|---|---|
| `py_01_uncorrected_stim.png` | UMAP of raw PCA, coloured by condition — the batch effect |
| `py_02_harmony_stim.png` | UMAP after Harmony, by condition — now mixed |
| `py_03_harmony_celltype.png` | Same map by cell type — the biology survived |
| `py_04_scoreboard.png` | Batch mixing vs cell-type recovery, per method |

---

## R Seurat → Shanuz API

| Task | R (Seurat) | Python (Shanuz) |
|------|-----------|-----------------|
| Harmony | `RunHarmony(obj, "stim")` / `IntegrateLayers(method=HarmonyIntegration)` | `run_harmony(obj, group_by="stim")` / `integrate_layers(method="harmony")` |
| CCA anchors | `FindIntegrationAnchors(list, reduction="cca")` + `IntegrateData` / `IntegrateLayers(method=CCAIntegration)` | `find_integration_anchors(objs, reduction="cca")` + `integrate_data` / `integrate_layers(method="cca")` |
| RPCA | `IntegrateLayers(method=RPCAIntegration)` | `integrate_layers(method="rpca")` |
| Neighbours on a reduction | `FindNeighbors(obj, reduction="harmony", dims=1:30)` | `find_neighbors(obj, reduction="harmony", dims=range(30))` |
| Cluster | `FindClusters(obj, resolution=0.5)` | `find_clusters(obj, resolution=0.5)` |
| UMAP on a reduction | `RunUMAP(obj, reduction="harmony", dims=1:30)` | `run_umap(obj, reduction="harmony", dims=range(30))` |

---

## References

Kang HM, Subramaniam M, Targ S, et al. (2018) **Multiplexed droplet single-cell
RNA-sequencing using natural genetic variation.** *Nature Biotechnology* 36,
89-94. <https://doi.org/10.1038/nbt.4042>

Korsunsky I, Millard N, Fan J, et al. (2019) **Fast, sensitive and accurate
integration of single-cell data with Harmony.** *Nature Methods* 16, 1289-1296.
<https://doi.org/10.1038/s41592-019-0619-0>
