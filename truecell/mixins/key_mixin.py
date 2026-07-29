from __future__ import annotations

import re


_KEY_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9]*_$")


def _validate_key(key: str) -> None:
    """Keys must be non-empty, start with a letter, contain only alphanumerics, and end with '_'."""
    if not _KEY_RE.match(key):
        raise ValueError(
            f"Invalid key '{key}'. Keys must start with a letter, contain only "
            "alphanumeric characters, and end with '_'."
        )


class KeyMixin:
    """Mixin providing a validated 'key' slot, mirroring R's KeyMixin from keymixin.R."""

    __slots__ = ("_key",)

    @property
    def key(self) -> str:
        return self._key

    @key.setter
    def key(self, value: str) -> None:
        _validate_key(value)
        self._key = value
