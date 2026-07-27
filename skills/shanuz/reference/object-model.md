# The shanuz object model

Ported from Seurat v5's S4 classes onto `__slots__`-based Python classes. Every
container below is the direct analogue of the R one, and the object-model
tutorial pins **91 of 91 anchors exactly against Seurat, no tolerance**.

## Shanuz

```
Shanuz
├── assays: dict[str, Assay5]
│   └── "RNA"
│       ├── layers["counts"]      # raw integer counts   (features × cells)
│       ├── layers["data"]        # log-normalized       (features × cells)
│       └── layers["scale.data"]  # z-scored             (features × cells)
├── meta_data: pd.DataFrame       # per-cell, indexed by cell name
├── reductions: dict[str, DimReduc]     # "pca", "umap", "harmony", …
├── graphs: dict[str, Graph]            # "RNA_nn", "RNA_snn", "wknn", "wsnn"
├── neighbors: dict[str, Neighbor]
├── images: dict[str, FOV | VisiumV2]   # spatial only
├── commands: list[ShanuzCommand]       # audit log of what was run
├── misc: dict                          # stashed fit details (hto_demux, sketch, …)
└── tools: dict
```

Slots: `assays`, `meta_data`, `active_assay`, `_active_ident`, `graphs`,
`neighbors`, `reductions`, `images`, `project_name`, `misc`, `version`,
`commands`, `tools`.

### Methods

```python
obj.cell_names()                     -> list[str]        # Cells()
obj.feature_names(assay=None)        -> list[str]        # Features() / rownames()
obj.assay_names()                    -> list[str]
obj.reduction_names()                -> list[str]
obj.image_names()                    -> list[str]
obj.get_assay(assay=None)            -> Assay5
obj.embeddings(reduction, dims=None) -> np.ndarray       # Embeddings()
obj.fetch_data(vars, cells=None, layer=None) -> pd.DataFrame   # FetchData()
obj.add_meta_data(metadata, col_name=None)   -> Shanuz   # AddMetaData()
obj.subset(cells=None, features=None, idents=None) -> Shanuz   # NEW object
obj.merge(y, add_cell_ids=None, project=None)      -> Shanuz   # NEW object
obj.rename_cells(new_names)          -> Shanuz
obj.which_cells(ident=None, cells=None) -> list[str]     # WhichCells()
obj.set_ident(cells, ident)          -> None
obj.rename_idents(mapping)           -> Shanuz           # RenameIdents()
obj.reorder_ident(ident, order)      -> Shanuz
obj.stash_ident(save_name)           -> Shanuz
obj.get_tissue_coordinates(image=None) -> pd.DataFrame
obj.tool(key) / obj.set_tool(key, value)
len(obj)      # number of cells
obj[key]      # assay by name
```

`obj.idents` is a property — the active identity, a pandas Series over cells.
`obj.default_assay` is a property too.

`fetch_data` accepts metadata columns, feature names, and reduction columns
(`"PC_1"`, `"UMAP_2"`) in one call, and returns plain numbers — not sparse
objects. Mixing them is the point:

```python
df = obj.fetch_data(["UMAP_1", "UMAP_2", "seurat_clusters", "MS4A1"])
```

## Assay5 (the v5 layered assay)

```python
a = obj.get_assay()                  # active assay
a.layers_list(pattern=None)          -> list[str]        # Layers()
a.layer_data(layer=None, cells=None, features=None)      # LayerData()
a.set_layer_data(layer, value, cell_names=None, feature_names=None) -> None
a.cells(layer=None) / a.features(layer=None)
a.variable_features                  # property
a.default_layer                      # property
a.key                                # property, e.g. "rna_"
a.calc_n()                           -> pd.DataFrame     # nCount / nFeature
a.split_layers(f, layer=None)        -> Assay5           # split counts by a factor
a.join_layers(layers=None)           -> Assay5           # JoinLayers()
a.subset(cells=None, features=None)  -> Assay5
a.merge(y, add_cell_ids=None)        -> Assay5
a.cast_assay(to_sparse=True)         -> Assay5
a.rename_cells(new_names)            -> Assay5
```

**Split / join is how v5 holds batches.** `split_layers(f)` turns `counts` into
`counts.<level>` per level of `f`; `join_layers()` puts them back. The
integration functions expect a split object. A split/join round trip preserves
column order — that was a real defect once, and is now a pinned anchor.

`Assay` (v3, `use_v5=False`) is still supported and has `scale_data` as a slot
rather than a layer. Prefer `Assay5`.

## DimReduc

```python
dr = obj.reductions["pca"]
shanuz.generics.embeddings(dr)              # cells × components
shanuz.generics.loadings(dr, projected=False)   # features × components
shanuz.generics.stdev(dr)                   # per-component standard deviation
shanuz.generics.features(dr, projected=False)   # the features it was computed on
shanuz.generics.key(dr)                     # "PC_", "UMAP_", …
```

A reduction records the features it **actually used**, not the ones requested —
those two diverged once and it was a defect.

## Graph

A cell × cell sparse matrix plus `assay_used` and cell names.
`shanuz.generics.cells(graph)`, `as_graph(x, ...)`, `as_neighbor(x)`.

`find_neighbors` stores Seurat's **directed** KNN graph (not a symmetrized one)
and keeps the SNN diagonal — both matched to Seurat deliberately.

## ShanuzCommand

Every analysis function appends one: `name`, `time_stamp`, `assay_used`,
`call_string`, `params`, `key`. `obj.commands` is the provenance record; check it
when reproducing an object whose history you don't know.

## Generics — what dispatches on what

`shanuz.generics` uses `functools.singledispatch`, so the same name works across
containers. The registered pairs that matter:

| Generic | Registered for |
|---|---|
| `cells` | `Assay`, `StdAssay`, `DimReduc`, `Graph`, `Neighbor`, `SpatialImage`, `Shanuz` |
| `features` | `Assay`(`layer=`), `StdAssay`(`layer=`), `DimReduc`(`projected=`), `Shanuz`(`assay=`) |
| `layer_data` | `Assay`(`layer="data"`), `StdAssay`(`layer=None`, `cells=`, `features=`) |
| `layers` | `StdAssay`(`pattern=`), `Shanuz`(`assay=`, `pattern=`) |
| `split_layers` / `join_layers` | `StdAssay` |
| `variable_features` | `Assay`, `StdAssay`, `Shanuz`(`assay=`) |
| `default_assay` | `Assay`, `StdAssay`, `DimReduc`, `Graph`, `Shanuz` |
| `embeddings` | `DimReduc`, `Shanuz`(`reduction`, `dims=`) |
| `loadings` / `stdev` | `DimReduc` |
| `idents` / `which_cells` / `rename_idents` / `fetch_data` / `add_meta_data` | `Shanuz` |
| `radius` | `Centroids`, `SpatialImage` |
| `get_image` | `SpatialImage` |
| `as_sparse` | `np.ndarray` (`fmt="csc"`) |

## AnnData interoperability

```python
from shanuz.compat.anndata import as_anndata, from_anndata

adata = as_anndata(obj, assay=None)          # needs pip install shanuz[anndata]
obj   = from_anndata(adata, assay="RNA", spatial_key="spatial", fov_key="fov")
```

**Orientation flips.** shanuz/Seurat store features × cells; AnnData stores
cells × features. The conversion handles it — don't transpose by hand on top.

## Subsetting

```python
sub = obj.subset(cells=list_of_barcodes)          # by barcode
sub = obj.subset(idents=["0", "3"])               # by active identity
sub = obj.subset(features=gene_list)              # by feature
```

`subset` returns a new object — this is the one place in the API where not
rebinding is the bug. QC filtering is normally expressed as a boolean mask over
`obj.meta_data` turned into a barcode list, since there is no
`subset(subset = expr)` string-expression form.
