"""AnnData ↔ Seurat conversion helpers.

Requires: pip install seurat-object[anndata]
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


def as_anndata(seurat, assay: Optional[str] = None):
    """Convert a Seurat object to anndata.AnnData.

    Mapping
    -------
    active_assay counts/data layer → adata.X  (+ adata.layers for extras)
    meta_data                       → adata.obs
    assay.meta_features / meta_data → adata.var
    reductions["pca"].embeddings    → adata.obsm["X_pca"]
    reductions["pca"].loadings      → adata.varm["PCs"]
    graphs                          → adata.obsp
    misc                            → adata.uns
    """
    try:
        import anndata
    except ImportError:
        raise ImportError(
            "anndata is required. Install with: pip install 'seurat-object[anndata]'"
        )

    from ..assay import Assay
    from ..assay5 import StdAssay
    from ..shanuz import Shanuz

    if not isinstance(seurat, Shanuz):
        raise TypeError(f"Expected Shanuz, got {type(seurat).__name__}.")

    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays.get(assay_name)
    if assay_obj is None:
        raise KeyError(f"Assay '{assay_name}' not found.")

    cells = seurat.cell_names()

    # ---- X ----
    if isinstance(assay_obj, StdAssay):
        default_layer = assay_obj.default_layer
        X = assay_obj.layers.get(default_layer) if default_layer else None
        extra_layers = {
            k: v for k, v in assay_obj.layers.items() if k != default_layer
        }
        var_df = assay_obj.meta_data.copy() if assay_obj.meta_data is not None else pd.DataFrame()
        feature_names = assay_obj._all_feature_names
    else:
        from .._sparse import is_matrix_empty
        X = assay_obj.data if not is_matrix_empty(assay_obj.data) else assay_obj.counts
        extra_layers = {}
        if not is_matrix_empty(assay_obj.counts):
            extra_layers["counts"] = assay_obj.counts
        if not is_matrix_empty(assay_obj.scale_data):
            extra_layers["scale_data"] = assay_obj.scale_data
        var_df = assay_obj.meta_features.copy()
        feature_names = assay_obj._feature_names

    if X is None:
        X = sp.csc_matrix((len(feature_names), len(cells)))

    # anndata wants obs × var (cells × features), so transpose
    X_t = X.T if sp.issparse(X) else X.T

    # ---- obs ----
    obs = seurat.meta_data.copy()
    obs["ident"] = list(seurat.idents)

    # ---- var ----
    var = var_df.copy()
    var.index = feature_names

    # ---- layers ----
    layers_out = {}
    for layer_name, mat in extra_layers.items():
        layers_out[layer_name] = mat.T if sp.issparse(mat) else mat.T

    # ---- obsm ----
    obsm = {}
    for red_name, dr in seurat.reductions.items():
        key = f"X_{red_name.lower()}"
        obsm[key] = dr.cell_embeddings

    # ---- varm ----
    varm = {}
    for red_name, dr in seurat.reductions.items():
        if dr.feature_loadings.shape[0] == len(feature_names):
            varm[red_name.upper() + "s"] = dr.feature_loadings

    # ---- obsp ----
    obsp = {}
    for g_name, g in seurat.graphs.items():
        obsp[g_name] = g._matrix

    # ---- uns ----
    uns = dict(seurat.misc)
    uns["project_name"] = seurat.project_name
    uns["active_assay"] = assay_name

    return anndata.AnnData(
        X=X_t,
        obs=obs,
        var=var,
        layers=layers_out,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        uns=uns,
    )


def from_anndata(
    adata,
    assay: str = "RNA",
    spatial_key: str = "spatial",
    fov_key: str = "fov",
) -> "Shanuz":
    """Convert an anndata.AnnData to a Seurat object.

    Mapping
    -------
    adata.X                   → Assay5 'counts' layer  (transposed → features × cells)
    adata.layers              → additional Assay5 layers
    adata.obs                 → seurat.meta_data
    adata.var                 → assay.meta_data
    adata.obsm["X_pca"]       → seurat.reductions["pca"].cell_embeddings
    adata.obsm[spatial_key]   → seurat.images  (Centroids/FOV, split by obs[fov_key])
    adata.varm["PCs"]         → seurat.reductions["pca"].feature_loadings
    adata.obsp["connectivities"] → seurat.graphs
    adata.uns                 → seurat.misc

    ``spatial_key`` (default ``"spatial"``) is treated as physical coordinates
    and reconstructed into ``seurat.images`` — NOT as a dimensional reduction —
    so ``get_tissue_coordinates`` and the spatial-analysis functions work. If
    ``obs[fov_key]`` exists it splits the cells into one image per FOV.
    """
    try:
        import anndata as _ann
    except ImportError:
        raise ImportError(
            "anndata is required. Install with: pip install 'seurat-object[anndata]'"
        )

    from ..assay5 import Assay5
    from ..dimreduc import DimReduc
    from ..graph import Graph
    from ..shanuz import Shanuz, _VERSION

    cells = list(adata.obs_names)
    features = list(adata.var_names)

    # ---- Build Assay5 layers ----
    X = adata.X
    if sp.issparse(X):
        X_t = X.T.tocsc()
    else:
        X_t = np.asarray(X).T

    layers: dict = {"counts": X_t}
    for layer_name, mat in adata.layers.items():
        if sp.issparse(mat):
            layers[layer_name] = mat.T.tocsc()
        else:
            layers[layer_name] = np.asarray(mat).T

    meta_data_var = adata.var.copy() if adata.var is not None else pd.DataFrame(index=features)
    assay_obj = Assay5(
        layers=layers,
        feature_names=features,
        cell_names=cells,
        meta_data=meta_data_var,
        key=f"{assay.lower()}_",
    )

    # ---- Metadata ----
    meta_data = adata.obs.copy() if adata.obs is not None else pd.DataFrame(index=cells)

    # ---- Spatial images (obsm[spatial_key] → Centroids/FOV) ----
    images: dict = {}
    if spatial_key in adata.obsm:
        from ..spatial.fov import create_fovs
        xy = np.asarray(adata.obsm[spatial_key])[:, :2]
        coords = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "cell": cells})
        fov_labels = (adata.obs[fov_key].astype(str).to_numpy()
                      if fov_key in adata.obs.columns else None)
        images = create_fovs(coords, fov=fov_labels, assay=assay,
                             default_name=assay.lower())

    # ---- Reductions ----
    reductions: dict = {}
    for obsm_key, emb in adata.obsm.items():
        if obsm_key == spatial_key:
            continue                         # handled as images, not a reduction
        if obsm_key.startswith("X_"):
            red_name = obsm_key[2:]
        else:
            red_name = obsm_key

        varm_key = red_name.upper() + "s"
        loadings = adata.varm.get(varm_key) if adata.varm is not None else None

        dr = DimReduc(
            cell_embeddings=np.asarray(emb),
            cell_names=cells,
            feature_loadings=np.asarray(loadings) if loadings is not None else None,
            feature_names=features if loadings is not None else None,
            assay_used=assay,
            key=f"{red_name.upper()}_",
        )
        reductions[red_name] = dr

    # ---- Graphs ----
    graphs: dict = {}
    if adata.obsp is not None:
        for obsp_key, mat in adata.obsp.items():
            g = Graph(matrix=mat if sp.issparse(mat) else sp.csc_matrix(mat), cell_names=cells)
            graphs[obsp_key] = g

    # ---- misc ----
    misc = dict(adata.uns) if adata.uns is not None else {}
    project_name = misc.pop("project_name", "SeuratProject")

    return Shanuz(
        assays={assay: assay_obj},
        meta_data=meta_data,
        active_assay=assay,
        graphs=graphs,
        reductions=reductions,
        images=images,
        project_name=project_name,
        misc=misc,
        version=_VERSION,
    )
