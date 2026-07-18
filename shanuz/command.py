from __future__ import annotations

import inspect
from datetime import datetime
from typing import Any, Optional

import pandas as pd


class ShanuzCommand:
    """Logs commands executed on a Shanuz object.

    Mirrors R's ShanuzCommand class from command.R.

    Slots
    -----
    name        : str        function/method name
    time_stamp  : datetime   when the command ran
    assay_used  : Optional[str]
    call_string : str        human-readable call representation
    params      : dict       non-function parameters passed to the command
    """

    __slots__ = ("name", "time_stamp", "assay_used", "call_string", "params", "key")

    def __init__(
        self,
        name: str,
        time_stamp: Optional[datetime] = None,
        assay_used: Optional[str] = None,
        call_string: str = "",
        params: Optional[dict] = None,
        key: Optional[str] = None,
    ) -> None:
        self.name = name
        self.time_stamp = time_stamp or datetime.now()
        self.assay_used = assay_used
        self.call_string = call_string
        self.params = params or {}
        #: How R indexes this entry in ``obj@commands`` — the command name, the
        #: assay, and for reduction-consuming steps the reduction, joined with
        #: dots: ``FindNeighbors.RNA.pca``. Scripts look entries up by this, so
        #: it is data rather than decoration and matches Seurat's spelling.
        self.key = key or name

    def default_assay(self) -> Optional[str]:
        return self.assay_used

    def as_list(self, include_meta: bool = True) -> dict:
        d = dict(self.params)
        if include_meta:
            d["call_string"] = self.call_string
            d["time_stamp"] = self.time_stamp
        return d

    # ------------------------------------------------------------------
    # Mirrors R $-accessor for params
    # ------------------------------------------------------------------

    def __getattr__(self, name: str) -> Any:
        if name in ("name", "time_stamp", "assay_used", "call_string", "params"):
            raise AttributeError(name)
        try:
            return self.params[name]
        except KeyError:
            raise AttributeError(f"ShanuzCommand has no parameter '{name}'.")

    def __getitem__(self, key: str) -> Any:
        slot_map = {
            "name": self.name,
            "time_stamp": self.time_stamp,
            "assay_used": self.assay_used,
            "call_string": self.call_string,
            "params": self.params,
        }
        if key in slot_map:
            return slot_map[key]
        return self.params[key]

    def __repr__(self) -> str:
        ts = self.time_stamp.strftime("%Y-%m-%d %H:%M:%S")
        non_fn = {k: v for k, v in self.params.items() if not callable(v)}
        return f"ShanuzCommand: {self.call_string}\n  Time: {ts}\n  Params: {non_fn}"


# ------------------------------------------------------------------
# Factory
# ------------------------------------------------------------------

def log_shanuz_command(
    object_,
    func_name: str,
    params: Optional[dict] = None,
    assay: Optional[str] = None,
    reduction: Optional[str] = None,
) -> "ShanuzCommand":
    """Capture a command log entry, typically called at the end of a function.

    Mirrors R's LogSeuratCommand, including how it names the entry: Seurat's
    function name, the assay, and the reduction where one was consumed. The
    names are R's (``RunPCA``, not ``run_pca``) because the log is a lookup
    table users query — the same reasoning that keeps layer names ``scale.data``
    and reduction keys ``PC_``.
    """
    params = params or {}
    # Remove any callable values (functions/lambdas) from params log
    safe_params = {k: v for k, v in params.items() if not callable(v)}

    call_parts = [f"{func_name}("]
    parts = [f"{k}={v!r}" for k, v in safe_params.items()]
    call_parts.append(", ".join(parts))
    call_parts.append(")")
    call_string = "".join(call_parts)

    key = ".".join(p for p in (func_name, assay, reduction) if p)
    cmd = ShanuzCommand(
        name=func_name,
        assay_used=assay,
        call_string=call_string,
        params=safe_params,
        key=key,
    )
    # Append to object's command log if it has one
    if hasattr(object_, "commands") and isinstance(object_.commands, list):
        object_.commands.append(cmd)
    return cmd
