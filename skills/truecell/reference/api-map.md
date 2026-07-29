# truecell API map

Every public function, with its real signature and its Seurat equivalent.
Signatures are from `truecell` 0.9.0. `seurat` / `obj` as the first parameter means
a `Truecell` object.

**Read the return column.** `None` means the function mutates in place — see
contract 1 in the parent skill.

---

## Objects

| Seurat | truecell |
|---|---|
| `CreateSeuratObject` | `create_truecell_object` |
| `CreateAssayObject` | `create_assay_object` / `create_assay5_object` |
| `Seurat`, `Assay`, `Assay5`, `DimReduc`, `Graph`, `Neighbor` | same class names |

```python
create_truecell_object(counts, assay="RNA", min_cells=0, min_features=0,
                     project="SeuratProject", feature_names=None, cell_names=None,
                     meta_data=None, use_v5=True) -> Truecell
create_assay_object(counts=None, data=None, min_cells=0, min_features=0,
                    feature_names=None, cell_names=None, key="rna_") -> Assay
create_assay5_object(...same...) -> Assay5
log_truecell_command(object_, func_name, params=None, assay=None, reduction=None) -> TruecellCommand
as_graph(x, cell_names=None, assay_used=None, weighted=True) -> Graph
```

Classes exported: `Truecell`, `Assay`, `Assay5`, `StdAssay`, `DimReduc`, `Graph`,
`Neighbor`, `JackStrawData`, `LogMap`, `KeyMixin`, `TruecellCommand`.

## Loading data

Not top-level — import from the submodule.

```python
from truecell.io import read_10x
read_10x(data_dir, var_names="gene_symbols", make_unique=True) -> (csc_matrix, genes, cells)

from truecell.datasets import pbmc3k, pbmc8k, cbmc_citeseq, pbmc_hashing, thp1_eccite, \
                            ifnb, panc8, xenium_mouse_brain, visium_mouse_brain
pbmc3k(data_dir=None, force_download=False) -> (counts, genes, cells)
cbmc_citeseq(data_dir=None, force_download=False, species_prefix="HUMAN_")
xenium_mouse_brain(...) -> Path     # loaders that return a directory
visium_mouse_brain(...) -> Path
ifnb(data_dir=None); panc8(data_dir=None)   # need `Rscript tutorials/export_seuratdata.R <name>` once
```

Everything caches to `~/.truecell_data/` (~770 MB for the full set).

```python
from truecell.compat.anndata import as_anndata, from_anndata
as_anndata(seurat, assay=None)
from_anndata(adata, assay="RNA", spatial_key="spatial", fov_key="fov") -> Truecell
```

## Preprocessing

| Seurat | truecell | Returns |
|---|---|---|
| `PercentageFeatureSet` | `percentage_feature_set` | `None` |
| `NormalizeData` | `normalize_data` | `None` |
| `FindVariableFeatures` | `find_variable_features` | `None` |
| `ScaleData` | `scale_data` | `None` |
| `SCTransform` | `sctransform` | the object |

```python
percentage_feature_set(seurat, pattern, col_name=None, assay=None, layer="counts") -> None
normalize_data(seurat, normalization_method="LogNormalize", scale_factor=10000.0,
               assay=None, margin=1) -> None
find_variable_features(seurat, selection_method="vst", nfeatures=2000, assay=None,
                       layer=None, mean_cutoff=(0.1, 8), dispersion_cutoff=(1, inf),
                       num_bin=20, binning_method="equal_width") -> None
scale_data(seurat, features=None, vars_to_regress=None, assay=None, do_scale=True,
           do_center=True, scale_max=10.0, layer="data") -> None
sctransform(seurat, assay=None, new_assay_name="SCT", n_cells=5000, n_genes=2000,
            n_features=3000, min_cells=5, vars_to_regress=None, clip_range=None,
            gene_chunk=500, seed=42, set_default=True, vst_flavor="v2",
            bw_adjust=3.0, verbose=False)
```

- `normalization_method`: `"LogNormalize"`, `"CLR"`, `"RC"`. `margin=1` = per
  feature, `2` = per cell (matters for CLR on protein/HTO assays).
- `selection_method`: `"vst"` (default; honours `nfeatures`), `"mvp"` /
  `"dispersion"` (honours `mean_cutoff` / `dispersion_cutoff` instead).
- `scale_data(features=None)` scales the **variable features**, as Seurat does.
  Pass `features=truecell.generics.features(obj)` for `ScaleData(features = rownames(obj))`.
- `vst_flavor="v2"` is Seurat 5's model; `"v1"` is the 2019 one.

## Dimensional reduction

| Seurat | truecell |
|---|---|
| `RunPCA` | `run_pca` |
| `RunICA` | `run_ica` |
| `RunSPCA` | `run_spca` |
| `RunTSNE` | `run_tsne` |
| `RunUMAP` | `run_umap` |
| `JackStraw` / `ScoreJackStraw` | `jack_straw` / `score_jackstraw` |

```python
run_pca(seurat, n_pcs=50, features=None, assay=None, reduction_name="pca",
        reduction_key="PC_", seed=42, layer="scale.data") -> None
run_ica(seurat, nics=50, ..., reduction_name="ica", reduction_key="ICA_", max_iter=200) -> None
run_spca(seurat, graph, npcs=50, ..., reduction_name="spca") -> None     # graph is required
run_tsne(seurat, dims=None, reduction="pca", n_components=2, perplexity=30.0,
         reduction_name="tsne", seed=42, assay=None) -> None
run_umap(seurat, dims=None, reduction="pca", graph=None, n_components=2,
         n_neighbors=30, min_dist=0.3, metric="euclidean", reduction_name="umap",
         reduction_key="UMAP_", seed=42, assay=None) -> None
glm_pca(seurat, n_components=10, features=None, assay=None, reduction_name="glmpca",
        family="poisson", layer="counts", max_iter=100, tol=1e-4, penalty=1.0,
        learning_rate=0.1, theta=100.0, optimize_theta=True, seed=42) -> None
jack_straw(seurat, reduction="pca", dims=20, num_replicate=100, prop_freq=0.01,
           layer="scale.data", seed=42) -> JackStrawData
score_jackstraw(seurat, reduction="pca", dims=None, score_thresh=1e-5) -> np.ndarray
```

`run_umap` takes **either** `dims=` on a reduction **or** `graph=` (a graph name),
matching `RunUMAP`'s two modes. `glm_pca` runs on `counts`, not scaled data.

## Graphs and clustering

```python
find_neighbors(seurat, dims=None, k_param=20, assay=None, reduction="pca",
               graph_name=None, nn_name=None, prune_snn=1/15, seed=42) -> None
find_clusters(seurat, resolution=0.5, algorithm=1, graph_name=None, random_seed=0,
              n_iterations=-1, group_singletons=True) -> None
find_multi_modal_neighbors(seurat, reduction_list=("pca", "apca"), dims_list=None,
                           k_nn=20, l2_norm=True, knn_graph_name="wknn",
                           snn_graph_name="wsnn", knn_range=200, prune_snn=1/15,
                           sd_scale=1.0, cross_constant=None, smooth=False, seed=42) -> None
```

- `find_neighbors` writes `graphs["{assay}_nn"]` and `graphs["{assay}_snn"]`.
- `find_clusters` reads `{assay}_snn` unless `graph_name=` is given, and writes
  `meta_data["seurat_clusters"]` plus the active identity.
- `algorithm`: **1** = Louvain (default), **2** = Louvain multilevel,
  **4** = Leiden. **3 (SLM) is not implemented.**
- `group_singletons=True` absorbs size-1 clusters into their best-connected
  neighbour, as Seurat's `GroupSingletons` does.

## Differential expression → `truecell-differential-expression`

```python
find_markers(seurat, ident_1, ident_2=None, assay=None, layer=None, test_use="wilcox",
             only_pos=False, min_pct=0.1, logfc_threshold=0.25, features=None,
             latent_vars=None, sample_col=None, max_cells_per_ident=None,
             random_seed=1) -> pd.DataFrame
find_all_markers(seurat, assay=None, layer=None, test_use="wilcox", only_pos=False,
                 min_pct=0.1, logfc_threshold=0.25, sample_col=None,
                 max_cells_per_ident=None, random_seed=1, return_thresh=0.01) -> pd.DataFrame
find_conserved_markers(seurat, ident_1, grouping_var, ident_2=None, assay=None,
                       layer=None, test_use="wilcox", only_pos=False, min_pct=0.1,
                       logfc_threshold=0.25, features=None) -> pd.DataFrame
aggregate_expression(seurat, group_by="ident", assays=None, features=None,
                     layer="counts", return_object=False)
```

`test_use`: `wilcox` · `t` · `bimod` · `LR` · `negbinom` · `mast` · `deseq2` · `roc`.

## Integration and mapping → `truecell-integration`

```python
run_harmony(seurat, group_by, reduction="pca", dims=None, reduction_name="harmony",
            reduction_key="harmony_", theta=None, lambda_=None, sigma=0.1,
            nclust=None, max_iter_harmony=10, assay=None, seed=0) -> None
integrate_layers(seurat, method="harmony", orig_reduction="pca", new_reduction=None,
                 group_by=None, assay=None, **kwargs) -> None
find_integration_anchors(objects, anchor_features=None, reduction="cca", dims=30,
                         k_anchor=5, k_filter=200, k_score=30, reference=0,
                         layer="scale.data", seed=42) -> IntegrationAnchors
integrate_data(anchors, new_assay="integrated", k_weight=100, sd_weight=1.0,
               add_cell_ids=None, seed=42)            # returns a NEW object
integrate_embeddings(anchors, reduction, new_reduction="integrated_dr",
                     dims_to_integrate=None, k_weight=100, sd_weight=1.0) -> DimReduc
find_transfer_anchors(reference, query, anchor_features=None, reduction="pcaproject",
                      dims=30, k_anchor=5, k_filter=200, k_score=30,
                      layer="scale.data", seed=42) -> TransferAnchors
transfer_data(anchors, refdata, k_weight=50, sd_weight=1.0,
              refdata_features=None) -> pd.DataFrame
map_query(anchors, refdata=None, reference_reduction="pca", reduction_model="umap",
          reduction_name="ref.umap", reduction_key="refUMAP_", k_weight=50,
          sd_weight=1.0, refdata_features=None, layer="scale.data")
project_umap(query, reference, reduction="pca", umap_reduction="umap", dims=None,
             reduction_name="ref.umap", reduction_key="refUMAP_", layer="scale.data") -> DimReduc
```

`integrate_layers(method=)`: `"harmony"` · `"cca"` · `"rpca"`. `group_by=` is
required for every method.

## Signature scoring

```python
add_module_score(seurat, features, pool=None, nbin=24, ctrl=100, name="Cluster",
                 assay=None, layer="data", seed=1, search=False)
cell_cycle_scoring(seurat, s_features=None, g2m_features=None, assay=None,
                   layer="data", set_ident=False, nbin=24, ctrl=100, seed=1)
CC_GENES     # the Tirosh 2016 human S / G2M sets, used when the args are None
```

`features` for `add_module_score` is a gene list, a list of gene lists, or a
`{name: genes}` dict. Columns land as `{name}1`, `{name}2`, … or the dict keys.
`cell_cycle_scoring` writes `S.Score`, `G2M.Score`, `Phase`.

## Demultiplexing and screens → `truecell-multimodal`

```python
hto_demux(seurat, assay="HTO", positive_quantile=0.99, init=None, nstarts=10,
          kfunc="clara", nsamples=100, normalize=True, margin=1, seed=42, verbose=False)
multiseq_demux(seurat, assay="HTO", quantile=0.7, autothresh=False, maxiter=5,
               qrange=None, normalize=True, margin=1, verbose=False)
calc_perturb_sig(seurat, assay="RNA", features=None, layer="data", labels="gene",
                 nt_class="NT", split_by=None, num_neighbors=20, reduction="pca",
                 ndims=15, new_assay="PRTB")
run_mixscape(seurat, assay="PRTB", labels="gene", nt_class="NT", de_assay="RNA",
             layer="data", min_de_genes=5, min_cells=5, logfc_threshold=0.25,
             min_pct=0.05, pval_cutoff=0.05, iter_num=10, prtb_type="KO",
             new_class="mixscape_class", de_test="wilcox", seed=0, verbose=False)
mixscape_lda(seurat, labels="gene", nt_class="NT", assay="PRTB", de_assay="RNA",
             layer="data", npcs=10, logfc_threshold=0.25, min_pct=0.1,
             pval_cutoff=0.05, de_test="wilcox", reduction_name="lda",
             reduction_key="LDA_", scale_max=10.0, seed=42, verbose=False)
```

## Spatial → `truecell-spatial`

```python
load_xenium(path, assay="Xenium", fov_column=None, project="Xenium", keep_controls=False)
load_visium(path, assay="Spatial", project="Visium", image=True,
            image_resolution="lowres", filter_by_tissue=True, slice_name="slice1")
load_cosmx(path, expr_file=None, meta_file=None, assay="Nanostring",
           fov_column="fov", project="CosMx")
load_merscope(path, expr_file=None, meta_file=None, assay="Vizgen",
              fov_column="fov", project="MERSCOPE", keep_controls=False)

create_centroids(coords, nsides=0, radius=None, theta=None, assay="", key="centroids_")
create_segmentation(coords, assay="", key="segmentation_")
create_molecules(coords, assay="", key="molecules_")
create_fov(coords, type_="centroids", nsides=0, radius=None, theta=None, assay="", key="fov_")
create_fovs(coords, fov=None, assay="", default_name="fov") -> dict[str, FOV]

get_tissue_coordinates(seurat, image=None) -> pd.DataFrame
spatial_knn(coords, k=10, query=None) -> (distances, indices)
nearest_neighbor_distance(seurat, group_by, reference, target=None, image=None) -> pd.DataFrame
local_neighborhood(seurat, group_by, reference=None, k=10, image=None) -> pd.DataFrame
build_niche_assay(seurat, group_by, image=None, k=20, niches=4, assay_name="niche",
                  cluster=True, seed=0)
find_spatially_variable_features(seurat, features=None, method="moransi", k=10,
                                 weights="inverse_square", assay=None, layer=None,
                                 image=None, r_metric=5.0, bandwidth=1.0) -> pd.DataFrame
composition_test(seurat, group_by, split_by, reference=None) -> pd.DataFrame
```

Spatial classes: `SpatialImage`, `Centroids`, `Segmentation`, `Molecules`, `FOV`,
`VisiumV2`, `ScaleFactors`.

## Scale → `truecell-at-scale`

```python
leverage_score(obj, nsketch=5000, ndims=None, features=None, assay=None, layer="data",
               var_name="leverage.score", eps=0.5, seed=123) -> np.ndarray
sketch_data(obj, ncells=5000, method="LeverageScore", features=None, assay=None,
            layer="data", nsketch=5000, sketched_assay="sketch",
            var_name="leverage.score", seed=123)        # returns a NEW object
project_data(full, sketch, reduction="pca", full_reduction="pca.full",
             umap_reduction="umap", full_umap_reduction="ref.umap", refdata=None,
             project_umap=True, dims=None, k_weight=50, sd_weight=1.0,
             layer="scale.data")
write_lazy_matrix(matrix, path, *, overwrite=False) -> LazyMatrix
open_lazy_matrix(path) -> LazyMatrix
is_lazy(x) -> bool
```

## Plotting → `truecell-plotting`

All 17 return a `matplotlib.figure.Figure`.

`vln_plot` · `feature_plot` · `dim_plot` · `elbow_plot` · `feature_scatter` ·
`variable_feature_plot` · `viz_dim_loadings` · `dim_heatmap` · `do_heatmap` ·
`ridge_plot` · `dot_plot` · `image_dim_plot` · `image_feature_plot` ·
`spatial_dim_plot` · `spatial_feature_plot` · `plot_perturb_score` ·
`mixscape_heatmap`

## Generics (`truecell.generics.*`, not top-level)

`cells` · `features` · `idents` · `set_ident` · `stash_ident` · `rename_idents` ·
`reorder_ident` · `which_cells` · `fetch_data` · `layer_data` · `set_layer_data` ·
`layers` · `split_layers` · `join_layers` · `get_assay_data` · `set_assay_data` ·
`embeddings` · `loadings` · `set_loadings` · `stdev` · `variable_features` ·
`set_variable_features` · `hvf_info` · `default_assay` · `set_default_assay` ·
`set_default_layer` · `key` · `set_key` · `keys` · `assay_names` · `assay_class` ·
`cast_assay` · `add_meta_data` · `rename_cells` · `match_cells` · `calc_n` ·
`command` · `misc` / `set_misc` · `tool` / `set_tool` · `version` · `as_sparse` ·
`as_graph` · `as_neighbor` · `as_seurat` · `check_matrix` · `is_matrix_empty` ·
`simplify` · `stitch_matrix` · `distances` · `indices` ·
spatial: `boundaries` · `crop` · `overlay` · `radius` · `theta` · `get_image` ·
`get_molecules` · `get_tissue_coordinates` · `default_boundary` · `default_fov` ·
`is_global` · `as_centroids` · `as_segmentation` · `create_fov` · `create_centroids` ·
`create_segmentation`

See [`object-model.md`](object-model.md) for what each dispatches on.
