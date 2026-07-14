"""Generic functions via @singledispatch, mirroring R's generics.R.

Each public function dispatches on the first argument's type.
The base (object) implementation raises NotImplementedError so
unregistered types get a clear message.
"""
from __future__ import annotations

import functools
from typing import Any, Optional, Union


def _not_implemented(func_name: str):
    def _impl(x, *args, **kwargs):
        raise NotImplementedError(
            f"{func_name}() is not implemented for {type(x).__name__}."
        )
    return _impl


# ------------------------------------------------------------------
# Helper: build a singledispatch generic with a default error impl
# ------------------------------------------------------------------

def _generic(name: str):
    base = _not_implemented(name)
    base.__name__ = name
    return functools.singledispatch(base)


# ======================================================================
# DATA ACCESS
# ======================================================================

cells = _generic("cells")
features = _generic("features")
embeddings = _generic("embeddings")
distances = _generic("distances")
indices = _generic("indices")
loadings = _generic("loadings")
stdev = _generic("stdev")
fetch_data = _generic("fetch_data")
layer_data = _generic("layer_data")
get_assay_data = _generic("get_assay_data")
get_image = _generic("get_image")
get_tissue_coordinates = _generic("get_tissue_coordinates")
which_cells = _generic("which_cells")
tool = _generic("tool")
misc = _generic("misc")
command = _generic("command")
version = _generic("version")
key = _generic("key")
keys = _generic("keys")

# ======================================================================
# DATA ASSIGNMENT
# ======================================================================

set_layer_data = _generic("set_layer_data")
set_assay_data = _generic("set_assay_data")
set_ident = _generic("set_ident")
set_default_assay = _generic("set_default_assay")
set_default_layer = _generic("set_default_layer")
add_meta_data = _generic("add_meta_data")
set_loadings = _generic("set_loadings")
set_variable_features = _generic("set_variable_features")
set_key = _generic("set_key")
set_tool = _generic("set_tool")
set_misc = _generic("set_misc")

# ======================================================================
# IDENTITY MANAGEMENT
# ======================================================================

idents = _generic("idents")
rename_idents = _generic("rename_idents")
reorder_ident = _generic("reorder_ident")
stash_ident = _generic("stash_ident")

# ======================================================================
# ASSAY OPERATIONS
# ======================================================================

assay_names = _generic("assay_names")
default_assay = _generic("default_assay")
cast_assay = _generic("cast_assay")
create_assay_object = _generic("create_assay_object")
calc_n = _generic("calc_n")
assay_class = _generic("assay_class")

# ======================================================================
# REDUCTION / NEIGHBOR
# ======================================================================

default_dim_reduc = _generic("default_dim_reduc")
as_neighbor = _generic("as_neighbor")

# ======================================================================
# LAYER MANAGEMENT
# ======================================================================

layers = _generic("layers")
join_layers = _generic("join_layers")
split_layers = _generic("split_layers")

# ======================================================================
# VARIABLE FEATURES
# ======================================================================

variable_features = _generic("variable_features")
hvf_info = _generic("hvf_info")

# ======================================================================
# OBJECT CREATION / CONVERSION
# ======================================================================

create_shanuz_object = _generic("create_shanuz_object")
as_seurat = _generic("as_seurat")
as_graph = _generic("as_graph")
as_sparse = _generic("as_sparse")
s4_to_list = _generic("s4_to_list")
list_to_s4 = _generic("list_to_s4")

# ======================================================================
# SPATIAL
# ======================================================================

create_fov = _generic("create_fov")
create_segmentation = _generic("create_segmentation")
create_centroids = _generic("create_centroids")
as_centroids = _generic("as_centroids")
as_segmentation = _generic("as_segmentation")
crop = _generic("crop")
overlay = _generic("overlay")
simplify = _generic("simplify")
boundaries = _generic("boundaries")
default_fov = _generic("default_fov")
default_boundary = _generic("default_boundary")
get_molecules = _generic("get_molecules")
radius = _generic("radius")
theta = _generic("theta")

# ======================================================================
# UTILITIES
# ======================================================================

is_global = _generic("is_global")
is_matrix_empty = _generic("is_matrix_empty")
check_matrix = _generic("check_matrix")
match_cells = _generic("match_cells")
stitch_matrix = _generic("stitch_matrix")
rename_cells = _generic("rename_cells")

# ======================================================================
# Register implementations for core classes
# ======================================================================

def _register_all() -> None:
    from .assay import Assay
    from .assay5 import Assay5, StdAssay
    from .dimreduc import DimReduc
    from .graph import Graph
    from .neighbor import Neighbor
    from .shanuz import Shanuz
    from .spatial.base import SpatialImage
    from .spatial.centroids import Centroids
    from .spatial.fov import FOV
    from .spatial.molecules import Molecules
    from .spatial.segmentation import Segmentation
    from ._sparse import is_matrix_empty as _ime
    import numpy as np
    import scipy.sparse as sp

    # --- cells ---
    @cells.register(Assay)
    def _cells_assay(x, *a, **kw): return x.cells()

    @cells.register(StdAssay)
    def _cells_std(x, layer=None, **kw): return x.cells(layer)

    @cells.register(DimReduc)
    def _cells_dr(x, *a, **kw): return x.cells()

    @cells.register(Graph)
    def _cells_g(x, *a, **kw): return x.cells()

    @cells.register(Neighbor)
    def _cells_n(x, *a, **kw): return x.cells()

    @cells.register(SpatialImage)
    def _cells_si(x, *a, **kw): return x.cells()

    @cells.register(Shanuz)
    def _cells_s(x, *a, **kw): return x.cell_names()

    # --- features ---
    @features.register(Assay)
    def _feat_assay(x, layer=None): return x.features(layer)

    @features.register(StdAssay)
    def _feat_std(x, layer=None): return x.features(layer)

    @features.register(DimReduc)
    def _feat_dr(x, projected=False): return x.features(projected)

    @features.register(Shanuz)
    def _feat_s(x, assay=None): return x.feature_names(assay)

    # --- embeddings ---
    @embeddings.register(DimReduc)
    def _emb_dr(x, *a, **kw): return x.embeddings()

    @embeddings.register(Shanuz)
    def _emb_s(x, reduction, dims=None): return x.embeddings(reduction, dims)

    # --- distances / indices ---
    @distances.register(Neighbor)
    def _dist_n(x, *a, **kw): return x.distances()

    @indices.register(Neighbor)
    def _idx_n(x, *a, **kw): return x.indices()

    # --- loadings ---
    @loadings.register(DimReduc)
    def _load_dr(x, projected=False): return x.loadings(projected)

    # --- stdev ---
    @stdev.register(DimReduc)
    def _stdev_dr(x, *a, **kw): return x.stdev

    # --- layer_data ---
    @layer_data.register(Assay)
    def _ld_assay(x, layer="data", cells=None, features=None):
        return x.layer_data(layer, cells, features)

    @layer_data.register(StdAssay)
    def _ld_std(x, layer=None, cells=None, features=None):
        return x.layer_data(layer, cells, features)

    # --- get_assay_data ---
    @get_assay_data.register(Assay)
    def _gad_assay(x, layer="data"): return x.get_assay_data(layer)

    @get_assay_data.register(StdAssay)
    def _gad_std(x, layer=None): return x.layer_data(layer)

    # --- fetch_data ---
    @fetch_data.register(Shanuz)
    def _fd_s(x, vars, cells=None, layer=None): return x.fetch_data(vars, cells, layer)

    # --- which_cells ---
    @which_cells.register(Shanuz)
    def _wc_s(x, ident=None, cells=None): return x.which_cells(ident, cells)

    # --- default_assay ---
    @default_assay.register(Assay)
    def _da_assay(x): return x.assay_orig or ""

    @default_assay.register(StdAssay)
    def _da_std(x): return x.assay_orig or ""

    @default_assay.register(DimReduc)
    def _da_dr(x): return x.default_assay()

    @default_assay.register(Graph)
    def _da_g(x): return x.default_assay()

    @default_assay.register(Shanuz)
    def _da_s(x): return x.default_assay

    # --- variable_features ---
    @variable_features.register(Assay)
    def _vf_assay(x): return x.variable_features

    @variable_features.register(StdAssay)
    def _vf_std(x): return x.variable_features

    @variable_features.register(Shanuz)
    def _vf_s(x, assay=None): return x.get_assay(assay).variable_features

    # --- layers ---
    @layers.register(StdAssay)
    def _layers_std(x, pattern=None): return x.layers_list(pattern)

    @layers.register(Shanuz)
    def _layers_s(x, assay=None, pattern=None):
        a = x.get_assay(assay)
        if isinstance(a, StdAssay):
            return a.layers_list(pattern)
        return []

    # --- join_layers ---
    @join_layers.register(StdAssay)
    def _jl_std(x, layers=None): return x.join_layers(layers)

    # --- is_global ---
    @is_global.register(DimReduc)
    def _ig_dr(x): return x.is_global()

    @is_global.register(SpatialImage)
    def _ig_si(x): return x.is_global()

    # --- is_matrix_empty ---
    @is_matrix_empty.register(np.ndarray)
    def _ime_nd(x): return _ime(x)

    # --- rename_cells ---
    @rename_cells.register(Assay)
    def _rc_assay(x, new_names): return x.rename_cells(new_names)

    @rename_cells.register(DimReduc)
    def _rc_dr(x, new_names): return x.rename_cells(new_names)

    @rename_cells.register(Neighbor)
    def _rc_n(x, new_names): return x.rename_cells(new_names)

    @rename_cells.register(Shanuz)
    def _rc_s(x, new_names): return x.rename_cells(new_names)

    # --- add_meta_data ---
    @add_meta_data.register(Shanuz)
    def _amd_s(x, metadata, col_name=None): return x.add_meta_data(metadata, col_name)

    # --- idents ---
    @idents.register(Shanuz)
    def _idents_s(x): return x.idents

    # --- rename_idents ---
    @rename_idents.register(Shanuz)
    def _ri_s(x, mapping): return x.rename_idents(mapping)

    # --- stash_ident ---
    @stash_ident.register(Shanuz)
    def _si_s(x, save_name): return x.stash_ident(save_name)

    # --- key ---
    @key.register(Assay)
    def _key_assay(x): return x.key

    @key.register(StdAssay)
    def _key_std(x): return x.key

    @key.register(DimReduc)
    def _key_dr(x): return x.key

    @key.register(SpatialImage)
    def _key_si(x): return x.key

    # --- get_tissue_coordinates ---
    @get_tissue_coordinates.register(SpatialImage)
    def _gtc_si(x, cells=None): return x.get_tissue_coordinates(cells)

    @get_tissue_coordinates.register(FOV)
    def _gtc_fov(x, cells=None, boundary=None):
        return x.get_tissue_coordinates(cells, boundary)

    # --- boundaries ---
    @boundaries.register(FOV)
    def _bounds_fov(x, boundary=None): return x.get_boundaries(boundary)

    # --- simplify ---
    @simplify.register(Segmentation)
    def _simp_seg(x, tol=0.5): return x.simplify(tol)

    # --- radius ---
    @radius.register(Centroids)
    def _rad_c(x): return x.radius()

    @radius.register(SpatialImage)
    def _rad_si(x): return x.radius()

    # --- get_image ---
    # None for every image type except VisiumV2, which carries the H&E photo.
    @get_image.register(SpatialImage)
    def _img_si(x): return x.get_image()

    # --- theta ---
    @theta.register(Centroids)
    def _thet_c(x): return x.theta()

    # --- as_graph ---
    @as_graph.register(Neighbor)
    def _ag_n(x, weighted=True): return x.as_graph(weighted)

    # --- as_neighbor ---
    @as_neighbor.register(Graph)
    def _an_g(x): return x.as_neighbor()

    # --- as_sparse ---
    @as_sparse.register(np.ndarray)
    def _asp_nd(x, fmt="csc"):
        import scipy.sparse as sp
        return sp.csc_matrix(x) if fmt == "csc" else sp.csr_matrix(x)


_register_all()
