from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .base import SpatialImage


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
        radius: Optional[float] = None,
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
        self.radius_ = radius
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
