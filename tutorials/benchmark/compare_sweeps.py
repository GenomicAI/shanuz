#!/usr/bin/env python
"""Diff two benchmark sweeps, step by step.

    python tutorials/benchmark/compare_sweeps.py results_refblas results
    python tutorials/benchmark/compare_sweeps.py old new --arm truecell

Written for the reference-BLAS vs Accelerate comparison in `PERFORMANCE.md`
section 2.1, and used since for any change that claims to make one arm faster —
the question is the same either way, so the arm is a flag.

Two things get compared, and the second matters more than the first:

**Anchors.** A change that is supposed to be a speed-up must not change what the
code computes. Every step in the suite records a scalar summarising its result,
so the two sweeps can be diffed on those directly. Anything that moves is either
a last-place floating-point difference — worth seeing, not worth worrying about
— or evidence that the faster version is not computing the same thing, which
would make every timing in the report irrelevant.

**Timings.** How much the change actually bought, per step.
"""
from __future__ import annotations

import argparse
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
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--arm", default="seurat", choices=("seurat", "truecell"),
                    help="which arm's runs to diff (default: seurat)")
    args = ap.parse_args()

    before, after = (HERE / args.before), (HERE / args.after)
    arm = args.arm
    if not (before.is_dir() and after.is_dir()):
        print(f"no such sweep: {before if not before.is_dir() else after}")
        return 2
    benches = sorted({p.stem.split(".")[0] for p in before.glob(f"*.{arm}.json")})

    print(f"# Anchors: {before.name} vs {after.name}  ({arm} arm)\n")
    identical = moved = 0
    worst: list[tuple] = []
    for bench in benches:
        a, b = load(before, bench, arm), load(after, bench, arm)
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

    print(f"\n\n# Timings: {before.name} vs {after.name}  ({arm} arm)\n")
    for bench in benches:
        a, b = load(before, bench, arm), load(after, bench, arm)
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
