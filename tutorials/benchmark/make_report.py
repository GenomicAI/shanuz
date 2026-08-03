#!/usr/bin/env python
"""Turn the collected benchmark JSON into the tables PERFORMANCE.md is built on.

Separate from :mod:`run_benchmarks` on purpose: that module is what the sweep
executes, and this one is edited while a sweep is in flight.

    python tutorials/benchmark/make_report.py > tutorials/benchmark/tables.md
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

# Steps that are measurement scaffolding rather than analysis. They are real
# costs a user pays, and they are reported — but a pipeline total that included
# them would be comparing Python's import time against R's library() time and
# calling the difference an algorithm.
NON_PIPELINE = {"import", "library_load"}

# Steps deliberately outside the like-for-like pipeline: an extra variant one
# arm can do and the other cannot, measured so the report can say so.
ASIDES = {"neighbours_annoy", "umap_unseeded", "scale_all_genes",
          "rescale_hvg", "morans_i_full"}

# Excluded from the headline tally on top of ASIDES, because the two arms are
# not doing the same arithmetic: `blas_gemm` writes R's `a %*% t(a)`, which
# materialises the transpose, against numpy's `a @ a.T`, which passes a flag.
# It stays in the tables — see section 2.1 — but it does not get a vote.
NOT_LIKE_FOR_LIKE = {"blas_gemm"}


def load(bench: str, arm: str) -> dict | None:
    p = RESULTS / f"{bench}.{arm}.json"
    if not p.exists():
        return None
    doc = json.loads(p.read_text())
    steps: dict[str, dict] = {}
    for run in doc["runs"]:
        for s in run["steps"]:
            d = steps.setdefault(s["step"], {"sec": [], "rss": [], "anchor": []})
            d["sec"].append(s["seconds"])
            d["rss"].append(s["peak_rss_mb"])
            d["anchor"].append(s.get("anchor"))
    return {
        "steps": {k: {"sec": statistics.median(v["sec"]),
                      "spread": max(v["sec"]) - min(v["sec"]),
                      "rss": statistics.median(v["rss"]),
                      "anchor": v["anchor"][0]}
                  for k, v in steps.items()},
        "wall": statistics.median(r["wall_seconds"] for r in doc["runs"]),
        "peak": statistics.median(r["process_peak_rss_mb"] for r in doc["runs"]),
        "n": len(doc["runs"]),
        "machine": doc["machine"],
    }


def pipeline_seconds(arm: dict) -> float:
    """Sum of the steps both arms genuinely share."""
    return sum(v["sec"] for k, v in arm["steps"].items()
               if k not in NON_PIPELINE and k not in ASIDES)


def ratio(seurat_s: float, truecell_s: float) -> str:
    if seurat_s <= 0 or truecell_s <= 0:
        return "-"
    if seurat_s >= truecell_s:
        r = seurat_s / truecell_s
        return f"**{r:.1f}x**" if r >= 1.05 else "~equal"
    r = truecell_s / seurat_s
    return f"{r:.1f}x slower" if r >= 1.05 else "~equal"


def step_table(bench: str) -> list[str]:
    tc, sr, t1 = load(bench, "truecell"), load(bench, "seurat"), \
        load(bench, "truecell.t1")
    if not (tc and sr):
        return []
    out = [f"\n#### `{bench}`\n"]
    head = ["Step", "Truecell", "Seurat", "Faster by", "Truecell peak RSS",
            "Seurat peak RSS", "Truecell result", "Seurat result"]
    if t1:
        head.insert(2, "Truecell (1 thread)")
    out.append("| " + " | ".join(head) + " |")
    out.append("|---" + "|---:" * (len(head) - 1) + "|")
    names = list(tc["steps"]) + [n for n in sr["steps"] if n not in tc["steps"]]
    for name in names:
        a, b = tc["steps"].get(name), sr["steps"].get(name)
        if (a and a["anchor"] == -1) or (b and b["anchor"] == -1):
            who = " and ".join(w for w, v in (("Truecell", a), ("Seurat", b))
                               if v and v["anchor"] == -1)
            out.append(f"| `{name}` | _not available in {who}_ "
                       + "| " * (len(head) - 2) + "|")
            continue
        row = [f"`{name}`", f"{a['sec']:.2f}s" if a else "—"]
        if t1:
            c = t1["steps"].get(name)
            row.append(f"{c['sec']:.2f}s" if c else "—")
        row += [f"{b['sec']:.2f}s" if b else "—",
                ratio(b["sec"], a["sec"]) if a and b else "—",
                f"{a['rss']:.0f} MB" if a else "—",
                f"{b['rss']:.0f} MB" if b else "—",
                "" if not a or a["anchor"] is None else str(a["anchor"]),
                "" if not b or b["anchor"] is None else str(b["anchor"])]
        out.append("| " + " | ".join(row) + " |")

    pa, pb = pipeline_seconds(tc), pipeline_seconds(sr)
    tot = ["**shared pipeline**", f"**{pa:.1f}s**"]
    if t1:
        tot.append(f"**{pipeline_seconds(t1):.1f}s**")
    tot += [f"**{pb:.1f}s**", ratio(pb, pa), f"**{tc['peak']:.0f} MB**",
            f"**{sr['peak']:.0f} MB**", "", ""]
    out.append("| " + " | ".join(tot) + " |")
    asides = sorted(ASIDES & (set(tc["steps"]) | set(sr["steps"])))
    note = (f"\nMedian of {tc['n']} timed repeats per arm, warm-up discarded. "
            "*Shared pipeline* excludes interpreter start-up")
    note += (" and the arm-specific asides ("
             + ", ".join(f"`{a}`" for a in asides) + ").") if asides else "."
    out.append(note)
    return out


def scaling_table(benches: list[str]) -> list[str]:
    out = ["\n| Dataset | Cells | Truecell | Seurat | Faster by "
           "| Truecell peak RSS | Seurat peak RSS |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for bench in benches:
        tc, sr = load(bench, "truecell"), load(bench, "seurat")
        if not (tc and sr):
            continue
        cells = tc["steps"].get("read_counts", {}).get("anchor", "?")
        pa, pb = pipeline_seconds(tc), pipeline_seconds(sr)
        out.append(f"| `{bench.replace('_core', '')}` | {cells} | {pa:.1f}s "
                   f"| {pb:.1f}s | {ratio(pb, pa)} | {tc['peak']:.0f} MB "
                   f"| {sr['peak']:.0f} MB |")
    return out


def script_table() -> list[str]:
    p = RESULTS / "tutorial_scripts.json"
    if not p.exists():
        return []
    doc = json.loads(p.read_text())
    out = ["\n| Tutorial | Python script | R script | Faster by "
           "| Python peak RSS | R peak RSS |", "|---|---:|---:|---:|---:|---:|"]
    for name, pair in doc["tutorials"].items():
        tc, sr = pair["truecell"], pair["seurat"]
        def cell(r):
            return (f"{r['wall_seconds']:.1f}s" if r["exit_code"] == 0
                    else f"failed ({r['wall_seconds']:.0f}s)")
        both_ok = tc["exit_code"] == 0 and sr["exit_code"] == 0
        out.append(
            f"| `{name}` | {cell(tc)} | {cell(sr)} "
            f"| {ratio(sr['wall_seconds'], tc['wall_seconds']) if both_ok else '—'} "
            f"| {tc['peak_rss_mb']:.0f} MB | {sr['peak_rss_mb']:.0f} MB |")
    return out


def tally(benches: list[str]) -> list[str]:
    """The report's headline count, computed rather than eyeballed.

    An earlier draft of `PERFORMANCE.md` asserted "faster at 24 of the 31
    operations measured" from memory. The real figure was 50 of 68. Counting
    is the one part of this that should never be done by hand.
    """
    faster = slower = tied = 0
    losses: list[tuple[float, str, str]] = []
    for bench in benches:
        tc, sr = load(bench, "truecell"), load(bench, "seurat")
        if not (tc and sr):
            continue
        for step, rec in tc["steps"].items():
            if step in NON_PIPELINE or step in ASIDES or step in NOT_LIKE_FOR_LIKE:
                continue
            other = sr["steps"].get(step)
            if other is None or rec["sec"] <= 0 or other["sec"] <= 0:
                continue
            r = other["sec"] / rec["sec"]
            if r >= 1.05:
                faster += 1
            elif r <= 1 / 1.05:
                slower += 1
                losses.append((1 / r, bench, step))
            else:
                tied += 1
    out = [f"\nAcross **{faster + slower + tied}** like-for-like step "
           f"comparisons Truecell is faster in **{faster}**, slower in "
           f"**{slower}**, and within 5% in **{tied}**.\n",
           "Where it loses, worst first:\n",
           "| Step | Bench | Slower by |", "|---|---|---:|"]
    for gap, bench, step in sorted(losses, reverse=True):
        out.append(f"| `{step}` | `{bench}` | {gap:.1f}x |")
    return out


def main() -> int:
    core = ["pbmc3k_core", "pbmc8k_core", "ifnb_core", "thp1_core"]
    heavy = ["blas_probe", "pbmc3k_sctransform", "pbmc3k_de",
             "ifnb_integration", "xenium_spatial"]
    lines = ["## Tally"]
    lines += tally(core + heavy)
    lines.append("\n## Scaling: the standard workflow, 2.7k to 20.7k cells")
    lines += scaling_table(core)
    lines.append("\n## Step by step")
    for b in core:
        lines += step_table(b)
    lines.append("\n## Named operations")
    for b in heavy:
        lines += step_table(b)
    lines.append("\n## The tutorial scripts, end to end")
    lines += script_table()
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
