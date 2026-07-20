from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from .base import SpatialImage
from .centroids import Centroids, create_centroids
from .molecules import Molecules, create_molecules
from .segmentation import Segmentation, create_segmentation

BoundaryType = Union[Centroids, Segmentation]


class FOV(SpatialImage):
    """Field of view — modern container for spatially-resolved single-cell coordinates.

    Mirrors R's FOV class from fov.R.
    Can hold multiple segmentation boundaries (Centroids, Segmentation) and
    molecule-level FISH data (Molecules).

    Slots
    -----
    molecules             : dict[str, Molecules]
    boundaries            : dict[str, BoundaryType]
    coords_x_orientation  : str   which axis x maps to in visualisation
    """

    __slots__ = (
        "molecules",
        "boundaries",
        "coords_x_orientation",
        "_default_boundary",
        "assay",
        "misc",
        "_key",
    )

    def __init__(
        self,
        boundaries: Optional[dict[str, BoundaryType]] = None,
        molecules: Optional[dict[str, Molecules]] = None,
        coords_x_orientation: str = "horizontal",
        assay: str = "",
        key: str = "fov_",
        misc: Optional[dict] = None,
    ) -> None:
        super().__init__(assay=assay, key=key, misc=misc)
        self.boundaries: dict[str, BoundaryType] = boundaries or {}
        self.molecules: dict[str, Molecules] = molecules or {}
        self.coords_x_orientation = coords_x_orientation
        # Default boundary is first entry by insertion order
        self._default_boundary: Optional[str] = (
            next(iter(self.boundaries)) if self.boundaries else None
        )

    # ------------------------------------------------------------------
    # SpatialImage interface
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        if self._default_boundary and self._default_boundary in self.boundaries:
            return self.boundaries[self._default_boundary].cells()
        # Fall back to union of all boundaries
        seen: set[str] = set()
        result = []
        for b in self.boundaries.values():
            for c in b.cells():
                if c not in seen:
                    seen.add(c)
                    result.append(c)
        return result

    def dim(self) -> tuple[int, int]:
        return (len(self.cells()), 2)

    def get_tissue_coordinates(
        self,
        cells: Optional[list[str]] = None,
        boundary: Optional[str] = None,
    ) -> pd.DataFrame:
        bname = boundary or self._default_boundary
        if bname and bname in self.boundaries:
            return self.boundaries[bname].get_tissue_coordinates(cells=cells)
        return pd.DataFrame(columns=["x", "y", "cell"])

    def rename_cells(self, new_names: list[str]) -> "FOV":
        old_names = self.cells()
        if len(new_names) != len(old_names):
            raise ValueError("new_names length must match number of cells.")
        mapping = dict(zip(old_names, new_names))
        new_boundaries = {
            k: b.rename_cells([mapping.get(c, c) for c in b.cells()])
            for k, b in self.boundaries.items()
        }
        new_molecules = {
            k: m.rename_cells([mapping.get(c, c) for c in m.cells()] if m.cells() else [])
            if m.cells()
            else m
            for k, m in self.molecules.items()
        }
        return FOV(
            boundaries=new_boundaries,
            molecules=new_molecules,
            coords_x_orientation=self.coords_x_orientation,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    def subset(self, cells: list[str]) -> "FOV":
        cell_set = set(cells)
        new_boundaries = {k: b.subset(cells) for k, b in self.boundaries.items()}
        new_molecules = {k: m.subset(cells) for k, m in self.molecules.items()}
        return FOV(
            boundaries=new_boundaries,
            molecules=new_molecules,
            coords_x_orientation=self.coords_x_orientation,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    # ------------------------------------------------------------------
    # Boundary / molecule accessors
    # ------------------------------------------------------------------

    def get_boundaries(self, boundary: Optional[str] = None) -> Union[dict, BoundaryType]:
        if boundary is not None:
            return self.boundaries[boundary]
        return dict(self.boundaries)

    def default_boundary(self) -> Optional[str]:
        return self._default_boundary

    def set_default_boundary(self, value: str) -> None:
        if value not in self.boundaries:
            raise KeyError(f"Boundary '{value}' not found.")
        self._default_boundary = value

    def get_molecules(self, molecule: Optional[str] = None) -> Union[dict, Molecules]:
        if molecule is not None:
            return self.molecules[molecule]
        return dict(self.molecules)

    # ------------------------------------------------------------------
    # Crop / overlay
    # ------------------------------------------------------------------

    def crop(
        self,
        x_range: tuple[float, float],
        y_range: tuple[float, float],
    ) -> "FOV":
        all_cells: set[str] = set()
        for b in self.boundaries.values():
            coords = b.get_tissue_coordinates()
            mask = (
                (coords["x"] >= x_range[0])
                & (coords["x"] <= x_range[1])
                & (coords["y"] >= y_range[0])
                & (coords["y"] <= y_range[1])
            )
            all_cells.update(coords.index[mask])
        return self.subset(list(all_cells))

    def overlay(self, query: "FOV") -> list[str]:
        """Return cells from self whose centroids fall within any boundary of query."""
        query_coords = query.get_tissue_coordinates()
        if query_coords.empty:
            return []
        x_min, x_max = query_coords["x"].min(), query_coords["x"].max()
        y_min, y_max = query_coords["y"].min(), query_coords["y"].max()
        self_coords = self.get_tissue_coordinates()
        mask = (
            (self_coords["x"] >= x_min)
            & (self_coords["x"] <= x_max)
            & (self_coords["y"] >= y_min)
            & (self_coords["y"] <= y_max)
        )
        return list(self_coords.index[mask])

    # ------------------------------------------------------------------
    # Operators
    # ------------------------------------------------------------------

    def __getitem__(self, key: str):
        if key in self.boundaries:
            return self.boundaries[key]
        if key in self.molecules:
            return self.molecules[key]
        raise KeyError(f"'{key}' not found in boundaries or molecules.")

    def __repr__(self) -> str:
        n = len(self.cells())
        bounds = list(self.boundaries)
        mols = list(self.molecules)
        return (
            f"FOV: {n} cells  Key: {self._key!r}  Assay: {self.assay!r}\n"
            f"  Boundaries: {bounds}\n"
            f"  Molecules:  {mols}"
        )


# ------------------------------------------------------------------
# Factory
# ------------------------------------------------------------------

def create_fov(
    coords: pd.DataFrame,
    type_: str = "centroids",
    nsides: int = 0,
    radius: Optional[float] = None,
    theta: Optional[float] = None,
    assay: str = "",
    key: str = "fov_",
) -> FOV:
    """Create an FOV from a coordinate DataFrame.

    type_ : 'centroids' | 'segmentation' | 'molecules'
    """
    if type_ == "centroids":
        boundary = create_centroids(coords, nsides=nsides, radius=radius, theta=theta, assay=assay)
        return FOV(boundaries={"centroids": boundary}, assay=assay, key=key)
    elif type_ == "segmentation":
        boundary = create_segmentation(coords, assay=assay)
        return FOV(boundaries={"segmentation": boundary}, assay=assay, key=key)
    elif type_ == "molecules":
        mol = create_molecules(coords, assay=assay)
        return FOV(molecules={"molecules": mol}, assay=assay, key=key)
    else:
        raise ValueError(f"Unknown type_ '{type_}'. Choose centroids, segmentation, or molecules.")


def create_fovs(
    coords: pd.DataFrame,
    fov: Optional[Union[str, "pd.Series", Sequence]] = None,
    assay: str = "",
    default_name: str = "fov",
) -> dict[str, FOV]:
    """Build a ``{name: FOV}`` dict of centroid FOVs from a coordinate frame.

    ``coords`` must have columns ``x, y, cell``. When ``fov`` is given (a column
    name in ``coords`` or a per-row array of labels) the cells are split into one
    FOV per distinct label — matching a multi-FOV Xenium/CosMx run. Otherwise a
    single FOV named ``default_name`` is returned.

    Shared by the spatial loaders and ``from_anndata`` so both build identical,
    accessor-ready ``seurat.images`` structures.
    """
    coords = coords.copy()
    for col in ("x", "y", "cell"):
        if col not in coords.columns:
            raise ValueError(f"coords must have a '{col}' column.")

    if fov is None:
        return {default_name: create_fov(coords, type_="centroids", assay=assay,
                                         key=f"{default_name}_")}

    if isinstance(fov, str):
        if fov not in coords.columns:
            raise ValueError(f"fov column '{fov}' not in coords.")
        labels = coords[fov].astype(str).to_numpy()
    else:
        labels = pd.Series(fov).astype(str).to_numpy()
        if len(labels) != len(coords):
            raise ValueError("fov label length must match number of rows in coords.")

    images: dict[str, FOV] = {}
    for name in pd.unique(labels):
        sub = coords.loc[labels == name, ["x", "y", "cell"]]
        safe = str(name).replace(" ", "_")
        images[safe] = create_fov(sub, type_="centroids", assay=assay, key=f"{safe}_")
    return images
