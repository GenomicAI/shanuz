"""Spatial neighbourhood analysis — the FNN / spatial-stats layer.

Mirrors the spatial-analysis idioms a Seurat user reaches for:

  get_tissue_coordinates   GetTissueCoordinates (object-level, all images)
  spatial_knn              FNN::get.knn / get.knnx  (low-level KD-tree)
  nearest_neighbor_distance per-cell distance to nearest cell of a group
  local_neighborhood       per-cell composition of the k nearest neighbours
  build_niche_assay        Seurat v5 BuildNicheAssay (niche clustering)

Coordinates come from each image's ``get_tissue_coordinates`` (Centroids/FOV),
so anything with populated ``seurat.images`` works — including objects built by
``from_anndata`` on a spatial ``obsm``.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree


# ---------------------------------------------------------------------------
# Coordinate access
# ---------------------------------------------------------------------------

def _image_names(seurat, image: Optional[Union[str, Sequence[str]]]) -> list[str]:
    if not getattr(seurat, "images", None):
        raise ValueError(
            "Object has no spatial images. Populate `seurat.images` (e.g. via a "
            "spatial loader or from_anndata with an obsm spatial key)."
        )
    if image is None:
        return list(seurat.images)
    if isinstance(image, str):
        return [image]
    return list(image)


def get_tissue_coordinates(
    seurat,
    image: Optional[Union[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """Return centroid coordinates for the object, concatenated across images.

    Mirrors object-level ``GetTissueCoordinates``. Returns a DataFrame with
    columns ``x, y, cell, image`` (one row per cell per image).
    """
    frames = []
    for nm in _image_names(seurat, image):
        coords = seurat.images[nm].get_tissue_coordinates()
        c = coords.copy()
        c["cell"] = list(coords.index)
        c["image"] = nm
        frames.append(c.reset_index(drop=True))
    if not frames:
        return pd.DataFrame(columns=["x", "y", "cell", "image"])
    return pd.concat(frames, ignore_index=True)[["x", "y", "cell", "image"]]


# ---------------------------------------------------------------------------
# Low-level KD-tree KNN  (FNN::get.knn / get.knnx)
# ---------------------------------------------------------------------------

def spatial_knn(
    coords: np.ndarray,
    k: int = 10,
    query: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """k-nearest-neighbour distances and indices on a set of coordinates.

    Mirrors ``FNN::get.knn`` (``query=None``, self excluded) and
    ``FNN::get.knnx`` (``query`` given, searched against ``coords``).

    Returns ``(distances, indices)`` each shape ``(n_query, k)``; indices point
    into ``coords``.
    """
    coords = np.asarray(coords, dtype=float)
    tree = cKDTree(coords)
    if query is None:
        d, i = tree.query(coords, k=min(k + 1, len(coords)))
        d = np.atleast_2d(d)
        i = np.atleast_2d(i)
        return d[:, 1:], i[:, 1:]           # drop self
    query = np.asarray(query, dtype=float)
    d, i = tree.query(query, k=k)
    if d.ndim == 1:                          # k == 1 → make 2-D (n_query, 1)
        d = d[:, None]
        i = i[:, None]
    return d, i


# ---------------------------------------------------------------------------
# Nearest-neighbour distance between groups
# ---------------------------------------------------------------------------

def nearest_neighbor_distance(
    seurat,
    group_by: str,
    reference,
    target=None,
    image: Optional[Union[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """Distance from each ``reference`` cell to the nearest ``target`` cell.

    Computed per image (so distances never cross FOVs), then concatenated.
    If ``target`` is None it defaults to ``reference`` (nearest same-type cell,
    self excluded) — this is the mast-to-nearest-mast idiom.

    Returns a DataFrame: ``cell, image, reference, target, distance``.
    """
    target = reference if target is None else target
    labels = seurat.meta_data[group_by].astype(str)
    same = str(target) == str(reference)
    out = []
    for nm in _image_names(seurat, image):
        coords = seurat.images[nm].get_tissue_coordinates()
        cells = list(coords.index)
        lab = labels.reindex(cells).astype(str).to_numpy()
        xy = coords[["x", "y"]].to_numpy(dtype=float)
        ref_mask = lab == str(reference)
        tgt_mask = lab == str(target)
        if ref_mask.sum() == 0 or tgt_mask.sum() < (2 if same else 1):
            continue
        tree = cKDTree(xy[tgt_mask])
        if same:
            dist = tree.query(xy[ref_mask], k=2)[0][:, 1]
        else:
            dist = tree.query(xy[ref_mask], k=1)[0]
            dist = np.atleast_1d(dist)
        out.append(pd.DataFrame({
            "cell": [c for c, m in zip(cells, ref_mask) if m],
            "image": nm, "reference": str(reference),
            "target": str(target), "distance": dist,
        }))
    if not out:
        return pd.DataFrame(columns=["cell", "image", "reference", "target", "distance"])
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------------------
# Local neighbourhood composition
# ---------------------------------------------------------------------------

def local_neighborhood(
    seurat,
    group_by: str,
    reference=None,
    k: int = 10,
    image: Optional[Union[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """Composition of each cell's ``k`` nearest neighbours (self excluded).

    For every ``reference`` cell (all cells if ``reference`` is None) returns
    the count and proportion of each ``group_by`` category among its k spatial
    neighbours. Columns: ``cell, image, n_<group>…, prop_<group>…``.

    ``prop_<reference>`` is the "local density" of that type; the ``n_<group>``
    columns are the raw neighbourhood composition.
    """
    labels = seurat.meta_data[group_by].astype(str)
    groups = sorted(labels.dropna().unique())
    out = []
    for nm in _image_names(seurat, image):
        coords = seurat.images[nm].get_tissue_coordinates()
        cells = list(coords.index)
        if len(cells) < k + 1:
            continue
        lab = labels.reindex(cells).astype(str).to_numpy()
        xy = coords[["x", "y"]].to_numpy(dtype=float)
        ref_idx = (np.arange(len(cells)) if reference is None
                   else np.where(lab == str(reference))[0])
        if ref_idx.size == 0:
            continue
        _, idx = cKDTree(xy).query(xy[ref_idx], k=k + 1)
        neigh = lab[idx[:, 1:]]                       # (n_ref, k), self dropped
        df = pd.DataFrame({"cell": [cells[j] for j in ref_idx], "image": nm})
        for g in groups:
            df[f"n_{g}"] = (neigh == g).sum(axis=1)
        tot = df[[f"n_{g}" for g in groups]].sum(axis=1).replace(0, np.nan)
        for g in groups:
            df[f"prop_{g}"] = df[f"n_{g}"] / tot
        out.append(df)
    if not out:
        return pd.DataFrame(columns=["cell", "image"])
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------------------
# Niche assay  (Seurat v5 BuildNicheAssay)
# ---------------------------------------------------------------------------

def build_niche_assay(
    seurat,
    group_by: str,
    image: Optional[Union[str, Sequence[str]]] = None,
    k: int = 20,
    niches: int = 4,
    assay_name: str = "niche",
    cluster: bool = True,
    seed: int = 0,
) -> "object":
    """Build a neighbourhood-composition assay and cluster cells into niches.

    Mirrors Seurat v5's ``BuildNicheAssay``: each cell's feature vector is the
    count of every ``group_by`` category among its ``k`` spatial neighbours.
    The composition matrix is stored as a new assay (``assay_name``); when
    ``cluster`` is True the proportions are k-means clustered into ``niches``
    groups and written to ``meta_data['niches']``.
    """
    from ..assay5 import create_assay5_object

    comp = local_neighborhood(seurat, group_by, reference=None, k=k, image=image)
    if comp.empty:
        raise ValueError("No cells with enough neighbours to build a niche assay.")
    groups = [c[2:] for c in comp.columns if c.startswith("n_")]
    comp = comp.set_index("cell")
    mat = comp[[f"n_{g}" for g in groups]].to_numpy(dtype=float).T   # groups × cells
    cells = list(comp.index)

    assay = create_assay5_object(
        counts=sp.csc_matrix(mat), feature_names=groups, cell_names=cells,
        key=f"{assay_name.lower()}_",
    )
    seurat.assays[assay_name] = assay

    if cluster:
        try:
            from sklearn.cluster import KMeans
        except ImportError as e:  # pragma: no cover
            raise ImportError("scikit-learn is required for niche clustering "
                              "(pip install 'shanuz[analysis]').") from e
        props = comp[[f"prop_{g}" for g in groups]].fillna(0).to_numpy()
        km = KMeans(n_clusters=niches, random_state=seed, n_init=10).fit(props)
        niche_lab = pd.Series([f"niche_{c + 1}" for c in km.labels_], index=cells)
        seurat.meta_data["niches"] = niche_lab.reindex(seurat.meta_data.index).values
    return seurat
