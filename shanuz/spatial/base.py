from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional

import pandas as pd

from ..mixins import KeyMixin


class SpatialImage(KeyMixin, ABC):
    """Abstract base class for all spatial image objects.

    Mirrors R's SpatialImage virtual class from spatial.R.
    Subclasses must implement: cells, dim, get_tissue_coordinates,
    rename_cells, subset.

    Slots
    -----
    assay : str            associated assay name
    misc  : dict           miscellaneous storage
    _key  : str            (from KeyMixin)
    """

    __slots__ = ("assay", "misc", "_key")

    def __init__(
        self,
        assay: str = "",
        key: str = "image_",
        misc: Optional[dict] = None,
    ) -> None:
        self.assay = assay
        self._key = key
        self.misc = misc or {}

    # ------------------------------------------------------------------
    # Abstract methods — subclasses MUST implement these
    # ------------------------------------------------------------------

    @abstractmethod
    def cells(self) -> list[str]:
        ...

    @abstractmethod
    def dim(self) -> tuple[int, int]:
        ...

    @abstractmethod
    def get_tissue_coordinates(
        self,
        cells: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        ...

    @abstractmethod
    def rename_cells(self, new_names: list[str]) -> "SpatialImage":
        ...

    @abstractmethod
    def subset(self, cells: list[str]) -> "SpatialImage":
        ...

    # ------------------------------------------------------------------
    # Provided methods — subclasses MAY override
    # ------------------------------------------------------------------

    def default_assay(self) -> str:
        return self.assay

    def set_default_assay(self, value: str) -> None:
        self.assay = value

    def get_image(self):
        return None

    def is_global(self) -> bool:
        return True

    def radius(self) -> Optional[float]:
        return None

    def theta(self) -> Optional[float]:
        return None

    def __getitem__(self, cells):
        if isinstance(cells, list):
            return self.subset(cells)
        return self.subset([cells])

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}\n"
            f"  Cells: {len(self.cells())}  Key: {self._key!r}  Assay: {self.assay!r}"
        )
