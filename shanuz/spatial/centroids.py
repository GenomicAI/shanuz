from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .base import SpatialImage


def _auto_radius(coords: pd.DataFrame) -> Optional[float]:
    """A default spot radius, matching SeuratObject's ``.AutoRadius``.

    ``0.01 * mean(width, height)`` of the bounding box — one percent of the
    slide's mean dimension. It is a drawing hint, not a measurement: nothing
    about the data says how big a cell is, so R picks a size that looks right at
    slide scale and lets the caller override it.

    Without this the radius stays ``None`` and every true-to-scale spot renderer
    silently falls back to a fixed-size scatter, which is what shanuz did for
    every FOV that did not come from a Visium ``scalefactors_json.json``.
    """
    if len(coords) == 0:
        return None
    spans = [
        float(np.ptp(coords[axis].to_numpy(dtype=float)))
        for axis in ("x", "y")
    ]
    radius = 0.01 * float(np.mean(spans))
    return radius if np.isfinite(radius) and radius > 0 else None


class Centroids(SpatialImage):
    """Cell centroid coordinates.

    Mirrors R's Centroids class from centroids.R.
    Stores one (x, y) coordinate per cell representing the center of that cell.

    Slots
    -----
    _coords  : pd.DataFrame   columns: x, y, cell
    nsides   : int            number of polygon sides (0 = circle)
    radius_  : Optional[float] spot radius (for spot-based technologies)
    theta_   : Optional[float] angle offset
    """

    __slots__ = ("_coords", "nsides", "radius_", "theta_", "assay", "misc", "_key")

    def __init__(
        self,
        coords: pd.DataFrame,
        nsides: int = 0,
        radius: Optional[float] = None,  # None → SeuratObject's .AutoRadius
        theta: Optional[float] = None,
        assay: str = "",
        key: str = "centroids_",
        misc: Optional[dict] = None,
    ) -> None:
        super().__init__(assay=assay, key=key, misc=misc)
        coords = coords.copy()
        for col in ("x", "y", "cell"):
            if col not in coords.columns:
                raise ValueError(f"coords must have a '{col}' column.")
        self._coords = coords[["x", "y", "cell"]].copy()
        self.nsides = nsides
        self.radius_ = radius if radius is not None else _auto_radius(self._coords)
        self.theta_ = theta

    # ------------------------------------------------------------------
    # SpatialImage interface
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        return list(self._coords["cell"])

    def dim(self) -> tuple[int, int]:
        return (len(self._coords), 2)

    def get_tissue_coordinates(
        self,
        cells: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        df = self._coords.set_index("cell")[["x", "y"]]
        if cells is not None:
            df = df.loc[cells]
        return df.copy()

    def rename_cells(self, new_names: list[str]) -> "Centroids":
        if len(new_names) != len(self._coords):
            raise ValueError("new_names length must match number of cells.")
        new_coords = self._coords.copy()
        new_coords["cell"] = new_names
        return Centroids(
            coords=new_coords,
            nsides=self.nsides,
            radius=self.radius_,
            theta=self.theta_,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    def subset(self, cells: list[str]) -> "Centroids":
        mask = self._coords["cell"].isin(cells)
        return Centroids(
            coords=self._coords[mask].copy(),
            nsides=self.nsides,
            radius=self.radius_,
            theta=self.theta_,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    # ------------------------------------------------------------------
    # Provided overrides
    # ------------------------------------------------------------------

    def radius(self) -> Optional[float]:
        return self.radius_

    def theta(self) -> Optional[float]:
        return self.theta_


def create_centroids(
    coords: pd.DataFrame,
    nsides: int = 0,
    radius: Optional[float] = None,
    theta: Optional[float] = None,
    assay: str = "",
    key: str = "centroids_",
) -> Centroids:
    return Centroids(
        coords=coords,
        nsides=nsides,
        radius=radius,
        theta=theta,
        assay=assay,
        key=key,
    )
