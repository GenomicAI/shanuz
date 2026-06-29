"""shanuz — Python port of satijalab/seurat-object (v5.4.0)."""

from .assay import Assay, create_assay_object
from .assay5 import Assay5, StdAssay, create_assay5_object
from .command import ShanuzCommand, log_shanuz_command
from .dimreduc import DimReduc
from .graph import Graph, as_graph
from .jackstraw import JackStrawData
from .logmap import LogMap
from .mixins import KeyMixin
from .neighbor import Neighbor
from .shanuz import Shanuz, create_shanuz_object
from .spatial import (
    Centroids,
    FOV,
    Molecules,
    Segmentation,
    SpatialImage,
    create_centroids,
    create_fov,
    create_molecules,
    create_segmentation,
)
from . import generics

__version__ = "5.4.0"

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
    "as_graph",
    "log_shanuz_command",
    # Generic functions module
    "generics",
    "__version__",
]
