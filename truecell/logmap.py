from __future__ import annotations

from typing import Iterator

import numpy as np


class LogMap:
    """Logical mapping: named boolean arrays indicating membership.

    Mirrors R's LogMap class from logmap.R.
    Keys are cell or feature names; values are boolean numpy arrays.
    """

    __slots__ = ("_map",)

    def __init__(self, data: dict[str, np.ndarray] | None = None) -> None:
        self._map: dict[str, np.ndarray] = {}
        if data:
            for k, v in data.items():
                self[k] = v

    # ------------------------------------------------------------------
    # Mapping protocol
    # ------------------------------------------------------------------

    def __getitem__(self, key: str) -> np.ndarray:
        return self._map[key]

    def __setitem__(self, key: str, value) -> None:
        self._map[key] = np.asarray(value, dtype=bool)

    def __delitem__(self, key: str) -> None:
        del self._map[key]

    def __contains__(self, key: object) -> bool:
        return key in self._map

    def __len__(self) -> int:
        return len(self._map)

    def __iter__(self) -> Iterator[str]:
        return iter(self._map)

    def keys(self):
        return self._map.keys()

    def values(self):
        return self._map.values()

    def items(self):
        return self._map.items()

    def get(self, key: str, default=None):
        return self._map.get(key, default)

    # ------------------------------------------------------------------
    # Bitwise ops (element-wise across all stored arrays)
    # ------------------------------------------------------------------

    def __and__(self, other: "LogMap") -> "LogMap":
        result = LogMap()
        for k in self._map:
            if k in other:
                result[k] = self._map[k] & other[k]
        return result

    def __or__(self, other: "LogMap") -> "LogMap":
        result = LogMap()
        all_keys = set(self._map) | set(other._map)
        for k in all_keys:
            a = self._map.get(k, np.zeros(0, dtype=bool))
            b = other._map.get(k, np.zeros(0, dtype=bool))
            if a.shape == b.shape:
                result[k] = a | b
        return result

    def __invert__(self) -> "LogMap":
        result = LogMap()
        for k, v in self._map.items():
            result[k] = ~v
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def copy(self) -> "LogMap":
        return LogMap({k: v.copy() for k, v in self._map.items()})

    def __repr__(self) -> str:
        keys = list(self._map)
        n = len(keys)
        preview = keys[:4]
        suffix = ", ..." if n > 4 else ""
        return f"LogMap({n} entries: {preview}{suffix})"
