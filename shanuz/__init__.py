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
from .reduction import run_pca, run_ica, run_spca, run_tsne
from .glmpca import glm_pca
from .neighbors import find_neighbors
from .multimodal import find_multi_modal_neighbors
from .clustering import find_clusters
from .umap import run_umap
from .integration import run_harmony, integrate_layers
from .anchors import (
    find_integration_anchors,
    integrate_data,
    IntegrationAnchors,
)
from .transfer import (
    find_transfer_anchors,
    transfer_data,
    TransferAnchors,
)
from .mapping import map_query, project_umap
from .sketch import sketch_data, project_data, leverage_score
from .lazy import LazyMatrix, write_lazy_matrix, open_lazy_matrix, is_lazy
from .hto import hto_demux
from .multiseq import multiseq_demux
from .mixscape import calc_perturb_sig, run_mixscape, mixscape_lda
from .markers import find_markers, find_all_markers, find_conserved_markers
from .aggregate import aggregate_expression
from .sctransform import sctransform
from .module_score import add_module_score, cell_cycle_scoring, CC_GENES
from .spatial import (
    Centroids,
    FOV,
    Molecules,
    ScaleFactors,
    Segmentation,
    SpatialImage,
    VisiumV2,
    create_centroids,
    create_fov,
    create_fovs,
    create_molecules,
    create_segmentation,
    build_niche_assay,
    find_spatially_variable_features,
    get_tissue_coordinates,
    local_neighborhood,
    nearest_neighbor_distance,
    spatial_knn,
    load_cosmx,
    load_merscope,
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
    spatial_dim_plot,
    spatial_feature_plot,
    plot_perturb_score,
    mixscape_heatmap,
)

__version__ = "0.2.0"

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
    "VisiumV2",
    "ScaleFactors",
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
    "find_spatially_variable_features",
    "composition_test",
    # Spatial loaders
    "load_xenium",
    "load_visium",
    "load_cosmx",
    "load_merscope",
    "as_graph",
    "log_shanuz_command",
    # Analysis pipeline (mirrors Seurat's top-level functions)
    "normalize_data",
    "find_variable_features",
    "scale_data",
    "percentage_feature_set",
    "run_pca",
    "run_ica",
    "run_spca",
    "run_tsne",
    "glm_pca",
    "find_neighbors",
    "find_multi_modal_neighbors",
    "find_clusters",
    "run_umap",
    "run_harmony",
    "integrate_layers",
    "find_integration_anchors",
    "integrate_data",
    "IntegrationAnchors",
    "find_transfer_anchors",
    "transfer_data",
    "TransferAnchors",
    "map_query",
    "project_umap",
    "sketch_data",
    "project_data",
    "leverage_score",
    "LazyMatrix",
    "write_lazy_matrix",
    "open_lazy_matrix",
    "is_lazy",
    "hto_demux",
    "multiseq_demux",
    "calc_perturb_sig",
    "run_mixscape",
    "mixscape_lda",
    "find_markers",
    "find_all_markers",
    "find_conserved_markers",
    "aggregate_expression",
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
    "spatial_dim_plot",
    "spatial_feature_plot",
    "plot_perturb_score",
    "mixscape_heatmap",
    "__version__",
]
