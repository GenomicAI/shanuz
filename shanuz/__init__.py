"""shanuz — Python port of satijalab/seurat-object (v5.4.0)."""

from .assay import Assay, create_assay_object
from .assay5 import Assay5, StdAssay, create_assay5_object
from .command import ShanuzCommand, log_shanuz_command
from .dimreduc import DimReduc
from .graph import Graph, as_graph
from .jackstraw import JackStrawData, jack_straw, score_jackstraw
from .logmap import LogMap
from .mixins import KeyMixin
from .neighbor import Neighbor
from .shanuz import Shanuz, create_shanuz_object
from .preprocessing import (
    normalize_data,
    find_variable_features,
    scale_data,
    percentage_feature_set,
)
from .reduction import run_pca
from .neighbors import find_neighbors
from .clustering import find_clusters
from .umap import run_umap
from .markers import find_markers, find_all_markers
from .sctransform import sctransform
from .module_score import add_module_score, cell_cycle_scoring, CC_GENES
from .spatial import (
    Centroids,
    FOV,
    Molecules,
    Segmentation,
    SpatialImage,
    create_centroids,
    create_fov,
    create_fovs,
    create_molecules,
    create_segmentation,
    build_niche_assay,
    get_tissue_coordinates,
    local_neighborhood,
    nearest_neighbor_distance,
    spatial_knn,
    load_cosmx,
    load_visium,
    load_xenium,
)
from .composition import composition_test
from . import generics
from . import plotting
from .plotting import (
    vln_plot,
    feature_plot,
    dim_plot,
    elbow_plot,
    feature_scatter,
    variable_feature_plot,
    viz_dim_loadings,
    dim_heatmap,
    do_heatmap,
    ridge_plot,
    dot_plot,
    image_dim_plot,
    image_feature_plot,
)

__version__ = "0.1.1"

__all__ = [
    # Core classes
    "Shanuz",
    "Assay",
    "Assay5",
    "StdAssay",
    "DimReduc",
    "Graph",
    "Neighbor",
    "JackStrawData",
    "LogMap",
    "KeyMixin",
    "ShanuzCommand",
    # Spatial
    "SpatialImage",
    "Centroids",
    "Segmentation",
    "Molecules",
    "FOV",
    # Factories
    "create_shanuz_object",
    "create_assay_object",
    "create_assay5_object",
    "create_centroids",
    "create_segmentation",
    "create_molecules",
    "create_fov",
    "create_fovs",
    # Spatial analysis
    "get_tissue_coordinates",
    "spatial_knn",
    "nearest_neighbor_distance",
    "local_neighborhood",
    "build_niche_assay",
    "composition_test",
    # Spatial loaders
    "load_xenium",
    "load_visium",
    "load_cosmx",
    "as_graph",
    "log_shanuz_command",
    # Analysis pipeline (mirrors Seurat's top-level functions)
    "normalize_data",
    "find_variable_features",
    "scale_data",
    "percentage_feature_set",
    "run_pca",
    "find_neighbors",
    "find_clusters",
    "run_umap",
    "find_markers",
    "find_all_markers",
    "jack_straw",
    "score_jackstraw",
    "sctransform",
    "add_module_score",
    "cell_cycle_scoring",
    "CC_GENES",
    # Generic functions module
    "generics",
    # Plotting module
    "plotting",
    "vln_plot",
    "feature_plot",
    "dim_plot",
    "elbow_plot",
    "feature_scatter",
    "variable_feature_plot",
    "viz_dim_loadings",
    "dim_heatmap",
    "do_heatmap",
    "ridge_plot",
    "dot_plot",
    "image_dim_plot",
    "image_feature_plot",
    "__version__",
]
