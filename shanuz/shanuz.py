from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from packaging.version import Version

from .assay import Assay, create_assay_object
from .assay5 import Assay5, create_assay5_object
from .command import ShanuzCommand, log_shanuz_command
from .dimreduc import DimReduc
from .graph import Graph
from .neighbor import Neighbor
from .spatial.fov import FOV

_VERSION = Version("5.4.0")

AnyAssay = Union[Assay, Assay5]


class Shanuz:
    """Top-level Shanuz single-cell data object.

    Mirrors R's Seurat class from seurat.R.

    Slots
    -----
    assays        : dict[str, AnyAssay]
    meta_data     : pd.DataFrame           cells × metadata columns
    active_assay  : str
    active_ident  : pd.Categorical
    graphs        : dict[str, Graph]
    neighbors     : dict[str, Neighbor]
    reductions    : dict[str, DimReduc]
    images        : dict[str, FOV]
    project_name  : str
    misc          : dict
    version       : packaging.version.Version
    commands      : list[ShanuzCommand]
    tools         : dict
    """

    __slots__ = (
        "assays",
        "meta_data",
        "active_assay",
        "_active_ident",
        "graphs",
        "neighbors",
        "reductions",
        "images",
        "project_name",
        "misc",
        "version",
        "commands",
        "tools",
    )

    def __init__(
        self,
        assays: dict[str, AnyAssay],
        meta_data: pd.DataFrame,
        active_assay: str,
        active_ident: Optional[pd.Categorical] = None,
        graphs: Optional[dict[str, Graph]] = None,
        neighbors: Optional[dict[str, Neighbor]] = None,
        reductions: Optional[dict[str, DimReduc]] = None,
        images: Optional[dict[str, FOV]] = None,
        project_name: str = "SeuratProject",
        misc: Optional[dict] = None,
        version: Optional[Version] = None,
        commands: Optional[list[ShanuzCommand]] = None,
        tools: Optional[dict] = None,
    ) -> None:
        self.assays = assays
        self.meta_data = meta_data
        self.active_assay = active_assay
        self._active_ident = active_ident if active_ident is not None else pd.Categorical(
            meta_data.index.tolist()
        )
        self.graphs = graphs or {}
        self.neighbors = neighbors or {}
        self.reductions = reductions or {}
        self.images = images or {}
        self.project_name = project_name
        self.misc = misc or {}
        self.version = version or _VERSION
        self.commands = commands or []
        self.tools = tools or {}

    # ------------------------------------------------------------------
    # Cell / feature names
    # ------------------------------------------------------------------

    def cell_names(self) -> list[str]:
        return list(self.meta_data.index)

    def feature_names(self, assay: Optional[str] = None) -> list[str]:
        a = self.assays.get(assay or self.active_assay)
        if a is None:
            return []
        return a.features()

    # ------------------------------------------------------------------
    # Idents
    # ------------------------------------------------------------------

    @property
    def idents(self) -> pd.Categorical:
        return self._active_ident

    @idents.setter
    def idents(self, value) -> None:
        cells = self.cell_names()
        if isinstance(value, pd.Categorical):
            if len(value) != len(cells):
                raise ValueError("Idents length must match number of cells.")
            self._active_ident = value
        elif isinstance(value, (list, np.ndarray, pd.Series)):
            self._active_ident = pd.Categorical(value)
        elif isinstance(value, dict):
            current = list(self._active_ident)
            for cell, new_id in value.items():
                if cell in cells:
                    idx = cells.index(cell)
                    current[idx] = new_id
            self._active_ident = pd.Categorical(current)
        else:
            raise TypeError(f"Cannot assign idents from {type(value).__name__}.")

    def set_ident(self, cells: list[str], ident: str) -> None:
        current = list(self._active_ident)
        all_cells = self.cell_names()
        for c in cells:
            idx = all_cells.index(c)
            current[idx] = ident
        self._active_ident = pd.Categorical(current)

    def stash_ident(self, save_name: str) -> "Shanuz":
        self.meta_data[save_name] = list(self._active_ident)
        return self

    def rename_idents(self, mapping: dict[str, str]) -> "Shanuz":
        new_idents = [mapping.get(str(x), str(x)) for x in self._active_ident]
        self._active_ident = pd.Categorical(new_idents)
        return self

    def reorder_ident(self, ident: str, order: list[str]) -> "Shanuz":
        self._active_ident = pd.Categorical(list(self._active_ident), categories=order)
        return self

    # ------------------------------------------------------------------
    # Active assay
    # ------------------------------------------------------------------

    @property
    def default_assay(self) -> str:
        return self.active_assay

    @default_assay.setter
    def default_assay(self, value: str) -> None:
        if value not in self.assays:
            raise KeyError(f"Assay '{value}' not found.")
        self.active_assay = value

    def assay_names(self) -> list[str]:
        return list(self.assays)

    def get_assay(self, assay: Optional[str] = None) -> AnyAssay:
        return self.assays[assay or self.active_assay]

    # ------------------------------------------------------------------
    # Reductions
    # ------------------------------------------------------------------

    def reduction_names(self) -> list[str]:
        return list(self.reductions)

    def embeddings(
        self,
        reduction: str,
        dims: Optional[list[int]] = None,
    ) -> np.ndarray:
        dr = self.reductions.get(reduction)
        if dr is None:
            raise KeyError(f"Reduction '{reduction}' not found.")
        emb = dr.cell_embeddings
        if dims is not None:
            emb = emb[:, dims]
        return emb

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------

    def add_meta_data(
        self,
        metadata: Union[pd.DataFrame, pd.Series, dict],
        col_name: Optional[str] = None,
    ) -> "Shanuz":
        if isinstance(metadata, pd.Series):
            if col_name is None:
                col_name = metadata.name or "metadata"
            self.meta_data[col_name] = metadata.reindex(self.meta_data.index)
        elif isinstance(metadata, dict):
            if col_name is None:
                col_name = "metadata"
            self.meta_data[col_name] = pd.Series(metadata).reindex(self.meta_data.index)
        elif isinstance(metadata, pd.DataFrame):
            for col in metadata.columns:
                self.meta_data[col] = metadata[col].reindex(self.meta_data.index)
        else:
            raise TypeError(f"Cannot add metadata from {type(metadata).__name__}.")
        return self

    # ------------------------------------------------------------------
    # FetchData — mirrors R FetchData()
    # ------------------------------------------------------------------

    def fetch_data(
        self,
        vars: list[str],
        cells: Optional[list[str]] = None,
        layer: Optional[str] = None,
    ) -> pd.DataFrame:
        cells = cells or self.cell_names()
        result = {}

        assay = self.get_assay()
        all_features = set(assay.features())

        for v in vars:
            if v in self.meta_data.columns:
                result[v] = self.meta_data.loc[cells, v]
            elif v in all_features:
                if isinstance(assay, Assay5):
                    mat = assay.layer_data(layer=layer, features=[v], cells=cells)
                else:
                    mat = assay.layer_data(layer=layer or "data", features=[v], cells=cells)
                result[v] = np.asarray(mat).flatten()
            elif v in self.reductions:
                emb = self.reductions[v].cell_embeddings
                all_c = self.cell_names()
                cell_idx = [all_c.index(c) for c in cells]
                n_dims = emb.shape[1]
                for d in range(n_dims):
                    col = f"{v}_{d + 1}"
                    result[col] = emb[cell_idx, d]
            else:
                raise KeyError(f"Variable '{v}' not found in metadata, features, or reductions.")

        return pd.DataFrame(result, index=cells)

    # ------------------------------------------------------------------
    # WhichCells
    # ------------------------------------------------------------------

    def which_cells(
        self,
        ident: Optional[Union[str, list[str]]] = None,
        cells: Optional[list[str]] = None,
    ) -> list[str]:
        all_cells = self.cell_names()
        result = cells or all_cells

        if ident is not None:
            ident_set = {ident} if isinstance(ident, str) else set(ident)
            ident_series = pd.Series(list(self._active_ident), index=all_cells)
            result = [c for c in result if str(ident_series.get(c, "")) in ident_set]

        return result

    # ------------------------------------------------------------------
    # Rename cells
    # ------------------------------------------------------------------

    def rename_cells(self, new_names: list[str]) -> "Shanuz":
        old_names = self.cell_names()
        if len(new_names) != len(old_names):
            raise ValueError("new_names must match number of cells.")
        mapping = dict(zip(old_names, new_names))

        new_meta = self.meta_data.copy()
        new_meta.index = new_names

        new_assays = {
            name: a.rename_cells(new_names) for name, a in self.assays.items()
        }
        new_graphs = {
            name: Graph(g._matrix, new_names, g.assay_used)
            for name, g in self.graphs.items()
        }
        new_neighbors = {
            name: n.rename_cells(new_names=new_names)
            for name, n in self.neighbors.items()
        }
        new_reductions = {
            name: r.rename_cells(new_names)
            for name, r in self.reductions.items()
        }
        new_ident = pd.Categorical(list(self._active_ident))

        return Shanuz(
            assays=new_assays,
            meta_data=new_meta,
            active_assay=self.active_assay,
            active_ident=new_ident,
            graphs=new_graphs,
            neighbors=new_neighbors,
            reductions=new_reductions,
            images=self.images,
            project_name=self.project_name,
            misc=dict(self.misc),
            version=self.version,
            commands=list(self.commands),
            tools=dict(self.tools),
        )

    # ------------------------------------------------------------------
    # Subset
    # ------------------------------------------------------------------

    def subset(
        self,
        cells: Optional[list[str]] = None,
        features: Optional[list[str]] = None,
        idents: Optional[Union[str, list[str]]] = None,
    ) -> "Shanuz":
        if idents is not None:
            cells = self.which_cells(ident=idents, cells=cells)
        if cells is None:
            cells = self.cell_names()

        cell_set = set(cells)
        new_meta = self.meta_data.loc[cells].copy()
        new_assays = {name: a.subset(cells=cells, features=features) for name, a in self.assays.items()}
        new_ident_vals = [
            str(i) for c, i in zip(self.cell_names(), self._active_ident) if c in cell_set
        ]
        new_ident = pd.Categorical(new_ident_vals)
        new_reductions = {name: r.subset(cells=cells) for name, r in self.reductions.items()}
        new_images = {name: img.subset(cells) for name, img in self.images.items()}

        # Subset each cell×cell graph to the retained cells. Mirrors Seurat,
        # which subsets graphs rather than carrying the full-size matrix.
        new_graphs = {name: g.subset(cells) for name, g in self.graphs.items()}
        # Neighbor objects store integer KNN indices into the original cell
        # ordering; those indices are invalidated by subsetting, so (as Seurat
        # does) drop them — re-run find_neighbors() on the subset.
        new_neighbors: dict[str, Neighbor] = {}

        return Shanuz(
            assays=new_assays,
            meta_data=new_meta,
            active_assay=self.active_assay,
            active_ident=new_ident,
            graphs=new_graphs,
            neighbors=new_neighbors,
            reductions=new_reductions,
            images=new_images,
            project_name=self.project_name,
            misc=dict(self.misc),
            version=self.version,
            commands=list(self.commands),
            tools=dict(self.tools),
        )

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def merge(
        self,
        y: Union["Shanuz", list["Shanuz"]],
        add_cell_ids: Optional[list[str]] = None,
        project: Optional[str] = None,
    ) -> "Shanuz":
        others = [y] if isinstance(y, Shanuz) else y
        all_objects = [self] + others

        if add_cell_ids is not None and len(add_cell_ids) != len(all_objects):
            raise ValueError("add_cell_ids must have one entry per Shanuz object.")

        # Merge cell names
        new_cell_names = []
        for idx, obj in enumerate(all_objects):
            prefix = add_cell_ids[idx] if add_cell_ids else None
            for c in obj.cell_names():
                new_cell_names.append(f"{prefix}_{c}" if prefix else c)

        # Merge metadata
        meta_frames = []
        for idx, obj in enumerate(all_objects):
            meta = obj.meta_data.copy()
            if add_cell_ids:
                meta.index = [f"{add_cell_ids[idx]}_{c}" for c in meta.index]
            meta_frames.append(meta)
        new_meta = pd.concat(meta_frames, axis=0, join="outer")

        # Merge assays (only shared assay names)
        shared_assay_names = set(all_objects[0].assays)
        for obj in all_objects[1:]:
            shared_assay_names &= set(obj.assays)

        new_assays: dict[str, AnyAssay] = {}
        for aname in shared_assay_names:
            base = all_objects[0].assays[aname]
            rest = [obj.assays[aname] for obj in all_objects[1:]]
            new_assays[aname] = base.merge(
                rest, add_cell_ids=add_cell_ids
            )

        # Merged ident
        ident_vals = []
        for obj in all_objects:
            ident_vals.extend([str(i) for i in obj._active_ident])
        new_ident = pd.Categorical(ident_vals)

        return Shanuz(
            assays=new_assays,
            meta_data=new_meta,
            active_assay=self.active_assay,
            active_ident=new_ident,
            project_name=project or self.project_name,
            version=self.version,
        )

    # ------------------------------------------------------------------
    # Tool storage (mirrors R Tool() / Tool<-())
    # ------------------------------------------------------------------

    def tool(self, key: str) -> object:
        return self.tools.get(key)

    def set_tool(self, key: str, value: object) -> None:
        self.tools[key] = value

    # ------------------------------------------------------------------
    # Operators
    # ------------------------------------------------------------------

    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            cells, features = key
        else:
            cells, features = key, None
        cells = list(cells) if not isinstance(cells, (list, type(None))) else cells
        features = list(features) if not isinstance(features, (list, type(None))) else features
        return self.subset(cells=cells, features=features)

    def __getattr__(self, name: str):
        # Mirrors R $ accessor — check metadata columns
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            meta = object.__getattribute__(self, "meta_data")
            if name in meta.columns:
                return meta[name]
        except AttributeError:
            pass
        raise AttributeError(f"'Shanuz' object has no attribute '{name}'.")

    def __repr__(self) -> str:
        n_cells = len(self.meta_data)
        assay = self.assays.get(self.active_assay)
        n_feat = len(assay.features()) if assay is not None else 0
        reds = list(self.reductions)
        return (
            f"Shanuz object — {self.project_name}\n"
            f"  {n_cells} cells × {n_feat} features\n"
            f"  Active assay: {self.active_assay!r}\n"
            f"  Reductions: {reds}\n"
            f"  Version: {self.version}"
        )

    def __len__(self) -> int:
        return len(self.meta_data)


# ------------------------------------------------------------------
# Factory
# ------------------------------------------------------------------

def create_shanuz_object(
    counts,
    assay: str = "RNA",
    min_cells: int = 0,
    min_features: int = 0,
    project: str = "SeuratProject",
    feature_names: Optional[list[str]] = None,
    cell_names: Optional[list[str]] = None,
    meta_data: Optional[pd.DataFrame] = None,
    use_v5: bool = True,
) -> Shanuz:
    """Create a Shanuz object from a counts matrix.

    Mirrors R's CreateSeuratObject().

    Parameters
    ----------
    counts       : sparse or dense matrix (features × cells)
    assay        : assay name (default "RNA")
    min_cells    : min cells a feature must be detected in to be kept
    min_features : min features a cell must have to be kept
    project      : project name
    feature_names: optional list of feature (gene) names
    cell_names   : optional list of cell barcodes
    meta_data    : optional per-cell metadata DataFrame
    use_v5       : if True, create Assay5 (v5); else Assay (v3)
    """
    key = f"{assay.lower()}_"

    if use_v5:
        assay_obj = create_assay5_object(
            counts=counts,
            min_cells=min_cells,
            min_features=min_features,
            feature_names=feature_names,
            cell_names=cell_names,
            key=key,
        )
        cells = assay_obj._all_cell_names
    else:
        assay_obj = create_assay_object(
            counts=counts,
            min_cells=min_cells,
            min_features=min_features,
            feature_names=feature_names,
            cell_names=cell_names,
            key=key,
        )
        cells = assay_obj._cell_names

    # Build metadata
    raw_meta = assay_obj.calc_n() if hasattr(assay_obj, "calc_n") else _calc_n_for_assay5(assay_obj)
    base_meta = raw_meta.rename(columns={"nCount": f"nCount_{assay}", "nFeature": f"nFeature_{assay}"})

    if meta_data is not None:
        # Align user-supplied metadata to filtered cells
        supplied = meta_data.reindex(cells)
        for col in supplied.columns:
            base_meta[col] = supplied[col].values

    active_ident = pd.Categorical([project] * len(cells), categories=[project])

    obj = Shanuz(
        assays={assay: assay_obj},
        meta_data=base_meta,
        active_assay=assay,
        active_ident=active_ident,
        project_name=project,
    )
    return obj


def _calc_n_for_assay5(assay: Assay5) -> pd.DataFrame:
    from ._utils import calc_n as _calc_n

    default_layer = assay.default_layer
    if default_layer is None:
        n = len(assay._all_cell_names)
        return pd.DataFrame(
            {"nCount_RNA": np.zeros(n), "nFeature_RNA": np.zeros(n)},
            index=assay._all_cell_names,
        )
    mat = assay.layers[default_layer]
    ncount, nfeature = _calc_n(mat)
    key_base = assay._key.rstrip("_").upper()
    return pd.DataFrame(
        {f"nCount_{key_base}": ncount, f"nFeature_{key_base}": nfeature},
        index=assay._all_cell_names,
    )
