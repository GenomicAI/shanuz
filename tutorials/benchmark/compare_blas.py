#!/usr/bin/env python
"""Compare two benchmark sweeps taken with different BLAS builds.

    python tutorials/benchmark/compare_blas.py results_refblas results

Two things get compared, and the second matters more than the first:

**Anchors.** Swapping the BLAS underneath R must not change what Seurat
computes. Every step in the suite records a scalar summarising its result, so
the two sweeps can be diffed on those directly. Anything that moves is either a
last-place floating-point difference — worth seeing, not worth worrying about —
or evidence that the faster BLAS is not computing the same thing, which would
make every timing in the report irrelevant.

**Timings.** How much the swap actually bought, per step.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(root: Path, bench: str, arm: str) -> dict | None:
    p = root / f"{bench}.{arm}.json"
    if not p.exists():
        return None
    doc = json.loads(p.read_text())
    steps: dict[str, dict] = {}
    for run in doc["runs"]:
        for s in run["steps"]:
            d = steps.setdefault(s["step"], {"sec": [], "anchor": s.get("anchor")})
            d["sec"].append(s["seconds"])
    return {k: {"sec": statistics.median(v["sec"]), "anchor": v["anchor"]}
            for k, v in steps.items()}


def relative(a, b) -> float | None:
    """Relative difference between two anchors, if both are numeric."""
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return None
    if fa == fb:
        return 0.0
    denom = max(abs(fa), abs(fb))
    return abs(fa - fb) / denom if denom else 0.0


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    before, after = (HERE / sys.argv[1]), (HERE / sys.argv[2])
    benches = sorted({p.stem.split(".")[0] for p in before.glob("*.seurat.json")})

    print(f"# Anchors: {before.name} vs {after.name}  (R arm)\n")
    identical = moved = 0
    worst: list[tuple] = []
    for bench in benches:
        a, b = load(before, bench, "seurat"), load(after, bench, "seurat")
        if not (a and b):
            continue
        for step in a:
            if step not in b:
                continue
            va, vb = a[step]["anchor"], b[step]["anchor"]
            if va is None and vb is None:
                continue
            rel = relative(va, vb)
            if str(va) == str(vb):
                identical += 1
            else:
                moved += 1
                worst.append((rel if rel is not None else 1.0,
                              bench, step, va, vb))
    print(f"identical: {identical}    moved: {moved}\n")
    for rel, bench, step, va, vb in sorted(worst, reverse=True):
        print(f"  {bench:20s} {step:20s} {va}  ->  {vb}"
              f"   (relative {rel:.2e})")

    print(f"\n\n# Timings: {before.name} vs {after.name}  (R arm)\n")
    for bench in benches:
        a, b = load(before, bench, "seurat"), load(after, bench, "seurat")
        if not (a and b):
            continue
        print(f"\n## {bench}")
        for step in a:
            if step not in b:
                continue
            ta, tb = a[step]["sec"], b[step]["sec"]
            speedup = ta / tb if tb > 0 else float("inf")
            flag = ""
            if speedup >= 1.10:
                flag = f"  {speedup:.1f}x faster on {after.name}"
            elif speedup <= 1 / 1.10:
                flag = f"  {1 / speedup:.1f}x SLOWER on {after.name}"
            print(f"  {step:22s} {ta:8.2f}s -> {tb:8.2f}s{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
