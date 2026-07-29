# API reference

Most of what follows is exported from the top level, so `truecell.find_markers` and
`truecell.markers.find_markers` are the same object. The grouping below is for
reading; it is not a package layout you need to know.

Two pages are the exception, and on those the import path shown is the one to
use. The [generics](generics.md) live on `truecell.generics` —
`truecell.generics.features(obj)`, not `truecell.features(obj)`, which raises
`AttributeError`. Seven of them are re-exported at the top level as well
(`create_truecell_object`, `create_assay_object`, `create_centroids`,
`create_segmentation`, `create_fov`, `get_tissue_coordinates`, `as_graph`); the
other 66 are not. The loaders on [Loading data](io.md) likewise stay on their own
modules: `truecell.io.read_10x`, `truecell.datasets.pbmc3k`,
`truecell.compat.anndata.as_anndata`.

The docstrings are the primary source. Many of them record a specific decision
about matching R — which of Seurat's two code paths a function follows, where a
default was chosen to agree with `Seurat::` rather than with the wider Python
ecosystem, and the handful of places the two genuinely differ. Those notes are
the reason this reference exists rather than a signature dump.

## The map from Seurat

| Seurat | truecell | Page |
|---|---|---|
| `CreateSeuratObject`, `Seurat`, `Assay5` | `create_truecell_object`, `Truecell`, `Assay5` | [Objects](objects.md) |
| `NormalizeData`, `FindVariableFeatures`, `ScaleData`, `SCTransform` | `normalize_data`, `find_variable_features`, `scale_data`, `sctransform` | [Preprocessing](preprocessing.md) |
| `RunPCA`, `RunUMAP`, `RunTSNE`, `JackStraw` | `run_pca`, `run_umap`, `run_tsne`, `jack_straw` | [Dimensional reduction](dimreduc.md) |
| `FindNeighbors`, `FindClusters`, `FindMultiModalNeighbors` | `find_neighbors`, `find_clusters`, `find_multi_modal_neighbors` | [Graphs and clustering](clustering.md) |
| `FindMarkers`, `FindAllMarkers`, `AggregateExpression` | `find_markers`, `find_all_markers`, `aggregate_expression` | [Differential expression](markers.md) |
| `IntegrateLayers`, `FindIntegrationAnchors`, `MapQuery` | `integrate_layers`, `find_integration_anchors`, `map_query` | [Integration and mapping](integration.md) |
| `AddModuleScore`, `CellCycleScoring` | `add_module_score`, `cell_cycle_scoring` | [Signature scoring](scoring.md) |
| `HTODemux`, `MULTIseqDemux`, `RunMixscape` | `hto_demux`, `multiseq_demux`, `run_mixscape` | [Demultiplexing and screens](demux.md) |
| `LoadXenium`, `BuildNicheAssay`, `FindSpatiallyVariableFeatures` | `load_xenium`, `build_niche_assay`, `find_spatially_variable_features` | [Spatial](spatial.md) |
| `SketchData`, `LeverageScore`, BPCells matrices | `sketch_data`, `leverage_score`, `LazyMatrix` | [Working at scale](scale.md) |
| `DimPlot`, `FeaturePlot`, `VlnPlot`, `DoHeatmap` | `dim_plot`, `feature_plot`, `vln_plot`, `do_heatmap` | [Plotting](plotting.md) |
| `Cells`, `Features`, `Idents`, `FetchData`, `LayerData` | `cells`, `features`, `idents`, `fetch_data`, `layer_data` | [Generics](generics.md) |
| `Read10X`, `SeuratData::` | `read_10x`, `truecell.datasets` | [Loading data](io.md) |

## Reading the signatures

Type annotations are resolved statically, straight from the source, so
annotation-only imports guarded by `if TYPE_CHECKING:` still render and still
cross-link. Three of them matter in practice — `matplotlib.figure.Figure` on
every plotting function, `Neighbor` on `as_graph`, and `Truecell` on
`from_anndata` — and all three are deliberate: they keep matplotlib optional and
break two import cycles. See [Fidelity](../fidelity.md#the-annotations-that-only-exist-for-readers).
