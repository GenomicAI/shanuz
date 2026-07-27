---
name: shanuz-multimodal
description: Use for shanuz analyses with more than one measurement per cell — CITE-seq RNA + surface protein (ADT) with CLR normalization and WNN joint clustering, cell-hashing demultiplexing with hto_demux or multiseq_demux, and pooled CRISPR screens with calc_perturb_sig, run_mixscape and mixscape_lda. Includes the CLR margin rule, which is the thing most often got wrong.
---

# Multimodal, hashing and screens

Load the `shanuz` skill first. Three workflows that share one idea: a second
assay on the same cells.

---

# 1 · CITE-seq — RNA + protein

## Attach the protein assay

```python
import shanuz
from shanuz.assay5 import create_assay5_object

obj.assays["ADT"] = create_assay5_object(
    counts=adt_counts,                       # proteins × cells, same barcodes
    feature_names=proteins,
    cell_names=obj.cell_names(),
    key="adt_",
)
shanuz.normalize_data(obj, assay="ADT", normalization_method="CLR", margin=2)
```

The ADT matrix must be aligned to the RNA object's cells first — intersect the
barcodes and reindex before constructing the assay.

## The CLR `margin` rule

`margin` means exactly what it means in Seurat, axis included. Seurat's
`CustomNormalize` runs `apply(data, MARGIN = margin, clr_function)`, and counts
are stored **features × cells**, so:

- **`margin=1`** — normalise each **protein/hashtag across cells**. Seurat's
  default. Correct for **HTO**.
- **`margin=2`** — normalise each **cell across its proteins**. Correct for
  **ADT** panels, which are small (a dozen or so markers).

This was inverted in shanuz once and is now pinned by a test against R ground
truth. `hto_demux` and `multiseq_demux` default to `margin=1` **deliberately** —
that is Seurat's default for hashtags. Do not "correct" them to 2.

## WNN joint clustering

The protein modality gets its own reduction first, then the two are fused.

```python
adt_features = obj.feature_names(assay="ADT")
shanuz.scale_data(obj, assay="ADT", features=adt_features)
shanuz.run_pca(obj, assay="ADT", reduction_name="apca", reduction_key="apca_",
               n_pcs=12, features=adt_features)

shanuz.find_multi_modal_neighbors(
    obj, reduction_list=["pca", "apca"],
    dims_list=[range(15), range(12)], k_nn=20,
)
shanuz.find_clusters(obj, graph_name="wsnn", resolution=0.6, random_seed=0)
shanuz.run_umap(obj, graph="wsnn", reduction_name="wnn_umap", seed=42)
```

- `n_pcs` must stay **below the protein count** — a 13-protein panel supports at
  most 12 components.
- `find_multi_modal_neighbors` writes graphs `wknn` and `wsnn`, plus per-cell
  weights `meta_data["RNA.weight"]` and `["ADT.weight"]`.
- Cluster and embed on `wsnn` / `graph="wsnn"`, not on `"pca"`. That is the whole
  point of WNN, and skipping it is the usual mistake.

The weights are interpretable — group them to see which lineages the protein
panel is actually deciding:

```python
obj.meta_data.groupby("cell_type")["ADT.weight"].mean().sort_values(ascending=False)
```

Plot protein on the RNA embedding with `assay="ADT"`:

```python
fig = shanuz.feature_plot(obj, ["CD3", "CD4", "CD8", "CD19"], assay="ADT",
                          reduction="umap", min_cutoff="q05", max_cutoff="q95", ncol=4)
```

Verified against Seurat: CLR to 4.2e-15, WNN modality weights at Pearson 0.9847,
cell-type labels 99.29 % concordant.

---

# 2 · Cell hashing — demultiplexing pooled samples

Two independent methods. Run both when the call matters; where they disagree is
where the sample is ambiguous.

```python
obj.assays["HTO"] = create_assay5_object(
    counts=hto_counts, feature_names=hashtags, cell_names=obj.cell_names(), key="hto_",
)

shanuz.hto_demux(obj, assay="HTO", positive_quantile=0.99, kfunc="clara")
shanuz.multiseq_demux(obj, assay="HTO", quantile=0.7)
```

**`hto_demux`** (Seurat's `HTODemux`): CLR normalize → cluster into
`k = n_hashtags + 1` groups → fit a negative binomial to each tag's background
cluster → threshold at `positive_quantile`. Cells positive for zero / one / many
tags become `Negative` / `Singlet` / `Doublet`.

Writes: `HTO_maxID`, `HTO_secondID`, `HTO_margin`, `HTO_classification`,
`HTO_classification.global`, and `hash.ID` (also set as the active identity).
Cutoffs are stashed in `obj.misc["hto_demux"]`.

**`multiseq_demux`** (Seurat's `MULTIseqDemux`): a Gaussian-KDE threshold placed
`quantile` of the way between each barcode's background and positive modes.
Writes `MULTI_ID` and `MULTI_classification`; sets the identity to `MULTI_ID`;
thresholds in `obj.misc["multiseq_demux"]`.

```python
shanuz.multiseq_demux(obj, autothresh=True, maxiter=5)   # sweep for the best singlet rate
```

Parameters worth knowing:

- `kfunc="clara"` (default) is Seurat's k-medoids; `"kmeans"` is the alternative.
  They rarely disagree on which cluster is a tag's background.
  **R's `clara` is architecture-dependent** — different answers on arm64 and
  x86_64. shanuz targets IEEE/x86_64 semantics on purpose.
- `normalize=False` uses the assay's existing `data` layer instead of
  CLR-normalizing internally — for when you have already normalized.
- `positive_quantile` up → stricter → more `Negative`, fewer doublets called.

Then filter:

```python
singlets = obj.subset(cells=list(obj.meta_data.index[
    obj.meta_data["HTO_classification.global"] == "Singlet"]))
```

99.81 % call-concordant with R Seurat's `HTODemux` on GSE108313, against
cross-species doublet ground truth.

---

# 3 · Pooled CRISPR screens — Mixscape

Three stages: isolate the perturbation signature, split real knockouts from
escapers, then map the guides.

```python
# 1. Subtract each cell's nearest non-targeting controls
shanuz.calc_perturb_sig(obj, assay="RNA", labels="gene", nt_class="NT",
                        split_by="replicate", num_neighbors=20,
                        reduction="pca", ndims=15, new_assay="PRTB")

# 2. Knockout vs non-perturbed, per guide
shanuz.run_mixscape(obj, assay="PRTB", labels="gene", nt_class="NT",
                    de_assay="RNA", min_de_genes=5, iter_num=10, prtb_type="KO")

# 3. The supervised map on which each guide population separates
shanuz.mixscape_lda(obj, labels="gene", nt_class="NT", assay="PRTB", npcs=10)
```

`labels` is the metadata column holding each cell's target gene / guide;
`nt_class` is the value marking the non-targeting controls.

`run_mixscape` writes `mixscape_class` (`"<gene> KO"` / `"NP"` / `"NT"`, also the
active identity), `mixscape_class.global`, and `mixscape_class_p_ko`; per-gene
bookkeeping goes to `obj.misc["mixscape"]`.

`mixscape_lda` writes an `lda` reduction plus `lda_assignments` and
`LDAP_<class>` columns.

### Diagnostics — run these, they are the point

```python
fig = shanuz.plot_perturb_score(obj, target_gene_ident="IFNGR2", assay="PRTB")
fig = shanuz.mixscape_heatmap(obj, ident_1="IFNGR2 KO", ident_2="NT", max_genes=100)
```

`plot_perturb_score` overlays the NT control density against one guide's own
along the perturbation score — **the axis mixscape actually splits on**. A guide
with a real effect is bimodal there. If it is unimodal, mixscape's KO/NP split
for that guide is not measuring anything, whatever the class column says.

A gene with fewer than `min_de_genes` DE genes against NT is untestable, and all
its cells are called NP. That is a stated outcome, not a failure — check how many
guides fall into it before reading the KO rates.

97.45 % per-cell call-concordant with R Seurat on the THP-1 ECCITE-seq screen.

---

## Reference

- [CITE-seq multimodal](https://genomicai.github.io/shanuz/tutorials/multimodal_citeseq/)
- [Cell hashing](https://genomicai.github.io/shanuz/tutorials/hashing_vignette/)
- [Mixscape](https://genomicai.github.io/shanuz/tutorials/mixscape_vignette/)
