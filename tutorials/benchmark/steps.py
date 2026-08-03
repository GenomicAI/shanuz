"""Step recorder for the Python side of the benchmark.

The child process writes one JSON line per timed step to a file the parent
(:mod:`run_benchmarks`) reads afterwards. Timestamps are epoch seconds from
:func:`time.time`, which is the same clock ``as.numeric(Sys.time())`` reads in
R — that is what lets the parent line one language's step boundaries up against
RSS samples it took of the other.

Each step also carries an ``anchor``: some cheap scalar summarising what the
step produced (cells kept, clusters found, markers returned). Anchors are not
timing data. They exist so the report can show the two arms did the *same*
work before it claims one of them did it faster.
"""
from __future__ import annotations

import json
import time
from contextlib import contextmanager
from pathlib import Path


class StepLog:
    """Append-only JSONL log of timed steps."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self.path.open("w")
        self._t_origin = time.time()

    @contextmanager
    def step(self, name: str):
        """Time a block, recording whatever the block puts in ``anchor``."""
        rec: dict = {}
        t0 = time.time()
        try:
            yield rec
        finally:
            t1 = time.time()
            self._emit({"step": name, "t0": t0, "t1": t1,
                        "seconds": t1 - t0, "anchor": rec.get("anchor")})

    def mark(self, name: str, seconds: float = 0.0, anchor=None) -> None:
        """Record a step that was timed by hand (an import, say)."""
        t1 = time.time()
        self._emit({"step": name, "t0": t1 - seconds, "t1": t1,
                    "seconds": seconds, "anchor": anchor})

    def _emit(self, rec: dict) -> None:
        self._fh.write(json.dumps(rec, default=str) + "\n")
        self._fh.flush()

    def close(self) -> None:
        self._fh.close()
