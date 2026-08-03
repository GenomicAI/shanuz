#!/usr/bin/env python
"""Run the Truecell-vs-Seurat performance benchmarks and collect the results.

The problem this solves: timing two languages fairly. Timing *inside* each
process is easy and each side already does it, but memory is not — R's
``gc()`` accounting and Python's ``tracemalloc`` measure different things and
neither sees what the other's allocator is holding. So memory is measured from
*outside*: this parent spawns the child, samples the resident set size of the
whole process tree on a fixed cadence, and afterwards intersects those samples
with the step boundaries the child wrote to its own log. Both arms are measured
by the same instrument, in the same units, on the same clock.

Usage
-----
    python tutorials/benchmark/run_benchmarks.py run --bench pbmc3k_core
    python tutorials/benchmark/run_benchmarks.py run --bench pbmc3k_core \\
        --arm truecell --repeats 5
    python tutorials/benchmark/run_benchmarks.py report

Results land in ``tutorials/benchmark/results/<bench>.<arm>.<n>.json``. The
``report`` subcommand reads whatever is there and writes ``PERFORMANCE.md``.

Notes on fairness
-----------------
* The **first** repeat of every (bench, arm) pair is a warm-up and is dropped:
  R's lazy package loading and umap-learn's numba JIT both cost seconds the
  first time and nothing afterwards, and neither is what the benchmark is
  asking about. ``--keep-warmup`` keeps it if you want to see that cost.
* Thread counts are pinned with ``--threads``. R's shipped BLAS is
  single-threaded reference BLAS; numpy here is on Accelerate, which is not.
  Running the truecell arm a second time at ``--threads 1`` is what separates
  "faster implementation" from "more cores".
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS = HERE / "results"
LOGS = HERE / "logs"

SAMPLE_INTERVAL = 0.05  # seconds between RSS samples

# The environment variables every thread-pool in either stack reads. Setting
# all of them is cheaper than working out which one each library honours.
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
               "RAYON_NUM_THREADS")


# ---------------------------------------------------------------------------
# RSS sampling
# ---------------------------------------------------------------------------

def _tree_rss_kb(root_pid: int) -> int:
    """Resident set of ``root_pid`` and every descendant, in kilobytes.

    Uses one ``ps`` call for the whole process table rather than one per pid:
    the children fork (BPCells, harmony, R's parallel), and walking the tree
    from a single snapshot is the only way to avoid missing a child that comes
    and goes between two calls.
    """
    try:
        out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="],
                             capture_output=True, text=True, timeout=5).stdout
    except (subprocess.SubprocessError, OSError):
        return 0
    kids: dict[int, list[int]] = {}
    rss: dict[int, int] = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            pid, ppid, kb = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            continue
        rss[pid] = kb
        kids.setdefault(ppid, []).append(pid)
    total, stack = 0, [root_pid]
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += rss.get(pid, 0)
        stack.extend(kids.get(pid, ()))
    return total


class Sampler(threading.Thread):
    """Polls the child's tree RSS until told to stop."""

    def __init__(self, pid: int) -> None:
        super().__init__(daemon=True)
        self.pid = pid
        self.samples: list[tuple[float, int]] = []
        self._halt = threading.Event()  # not _stop: Thread already has one

    def run(self) -> None:
        while not self._halt.is_set():
            self.samples.append((time.time(), _tree_rss_kb(self.pid)))
            self._halt.wait(SAMPLE_INTERVAL)

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=5)

    def nearest(self, t: float) -> int:
        """The sample taken closest to ``t``."""
        if not self.samples:
            return 0
        return min(self.samples, key=lambda s: abs(s[0] - t))[1]

    def peak_between(self, t0: float, t1: float) -> int:
        """Peak RSS over [t0, t1], or the nearest sample if the step was
        shorter than the sampling interval."""
        inside = [kb for t, kb in self.samples if t0 <= t <= t1]
        if inside:
            return max(inside)
        if not self.samples:
            return 0
        mid = (t0 + t1) / 2
        return min(self.samples, key=lambda s: abs(s[0] - mid))[1]


# ---------------------------------------------------------------------------
# Running one child
# ---------------------------------------------------------------------------

ARMS = {
    "truecell": lambda: [str(ROOT / ".venv" / "bin" / "python"),
                         str(HERE / "bench_truecell.py")],
    "seurat": lambda: ["/usr/local/bin/Rscript", str(HERE / "bench_seurat.R")],
}


def run_one(bench: str, arm: str, rep: int, threads: int | None) -> dict:
    steps_path = HERE / f".steps.{arm}.{bench}.{rep}.jsonl"
    if steps_path.exists():
        steps_path.unlink()
    LOGS.mkdir(exist_ok=True)
    log_path = LOGS / f"{bench}.{arm}.{rep}.log"

    env = dict(os.environ)
    if threads is not None:
        for var in THREAD_VARS:
            env[var] = str(threads)
    env["TRUECELL_BENCH_STEPS"] = str(steps_path)

    cmd = ARMS[arm]() + ["--bench", bench, "--steps", str(steps_path)]
    t_start = time.time()
    with log_path.open("w") as log:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env,
                                stdout=log, stderr=subprocess.STDOUT)
        sampler = Sampler(proc.pid)
        sampler.start()
        code = proc.wait()
        sampler.stop()
    wall = time.time() - t_start

    steps = []
    if steps_path.exists():
        for line in steps_path.read_text().splitlines():
            if line.strip():
                steps.append(json.loads(line))
        steps_path.unlink()
    for s in steps:
        s["peak_rss_mb"] = sampler.peak_between(s["t0"], s["t1"]) / 1024
        s["rss_start_mb"] = sampler.nearest(s["t0"]) / 1024
        # What this step *added*. Resident set is a high-water mark — neither
        # allocator returns freed pages promptly — so the raw peak of a late
        # step mostly reports what earlier steps left behind. The delta is the
        # step's own footprint; the process peak below is the number that
        # decides whether a machine can run the pipeline at all.
        s["delta_rss_mb"] = max(0.0, s["peak_rss_mb"] - s["rss_start_mb"])

    return {
        "bench": bench, "arm": arm, "repeat": rep, "threads": threads,
        "exit_code": code, "wall_seconds": wall,
        "process_peak_rss_mb": (max((kb for _, kb in sampler.samples),
                                    default=0) / 1024),
        "baseline_rss_mb": (min((kb for _, kb in sampler.samples if kb > 0),
                                default=0) / 1024),
        "n_samples": len(sampler.samples),
        "steps": steps,
        "log": str(log_path.relative_to(ROOT)),
    }


def machine() -> dict:
    def sysctl(key: str) -> str:
        try:
            return subprocess.run(["sysctl", "-n", key], capture_output=True,
                                  text=True).stdout.strip()
        except OSError:
            return "?"
    return {
        "cpu": sysctl("machdep.cpu.brand_string"),
        "cores": sysctl("hw.ncpu"),
        "memory_gb": round(int(sysctl("hw.memsize") or 0) / 1e9, 1),
        "platform": platform.platform(),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    RESULTS.mkdir(exist_ok=True)
    arms = args.arm.split(",")
    failures = 0
    for bench in args.bench.split(","):
        for arm in arms:
            kept = []
            total = args.repeats + (0 if args.keep_warmup else 1)
            for rep in range(total):
                warm = (rep == 0 and not args.keep_warmup)
                tag = " (warm-up, discarded)" if warm else ""
                print(f"  {bench:22s} {arm:9s} rep {rep}{tag} ...",
                      end="", flush=True)
                res = run_one(bench, arm, rep, args.threads)
                ok = res["exit_code"] == 0
                print(f" {res['wall_seconds']:7.1f}s  "
                      f"peak {res['process_peak_rss_mb']:6.0f} MB"
                      f"{'' if ok else '   FAILED (see ' + res['log'] + ')'}")
                if not ok:
                    failures += 1
                    break
                if not warm:
                    kept.append(res)
            if kept:
                suffix = "" if args.threads is None else f".t{args.threads}"
                out = RESULTS / f"{bench}.{arm}{suffix}.json"
                out.write_text(json.dumps(
                    {"machine": machine(), "runs": kept}, indent=2))
                print(f"  -> {out.relative_to(ROOT)}")
    return 1 if failures else 0


def summarise(path: Path) -> dict:
    """Median seconds and peak RSS per step across a result file's repeats."""
    doc = json.loads(path.read_text())
    per_step: dict[str, dict[str, list]] = {}
    order: list[str] = []
    for run in doc["runs"]:
        for s in run["steps"]:
            if s["step"] not in per_step:
                per_step[s["step"]] = {"seconds": [], "rss": [],
                                       "delta": [], "anchor": []}
                order.append(s["step"])
            per_step[s["step"]]["seconds"].append(s["seconds"])
            per_step[s["step"]]["rss"].append(s["peak_rss_mb"])
            per_step[s["step"]]["delta"].append(s.get("delta_rss_mb", 0.0))
            per_step[s["step"]]["anchor"].append(s.get("anchor"))
    out = {}
    for name in order:
        v = per_step[name]
        out[name] = {
            "seconds": statistics.median(v["seconds"]),
            "seconds_min": min(v["seconds"]),
            "seconds_max": max(v["seconds"]),
            "peak_rss_mb": statistics.median(v["rss"]),
            "delta_rss_mb": statistics.median(v["delta"]),
            "anchor": v["anchor"][0],
            "anchors_agree": len(set(map(str, v["anchor"]))) == 1,
        }
    return {
        "steps": out,
        "n_repeats": len(doc["runs"]),
        "wall_seconds": statistics.median(r["wall_seconds"] for r in doc["runs"]),
        "process_peak_rss_mb": statistics.median(
            r["process_peak_rss_mb"] for r in doc["runs"]),
        "baseline_rss_mb": statistics.median(
            r["baseline_rss_mb"] for r in doc["runs"]),
        "machine": doc["machine"],
    }


def cmd_report(args: argparse.Namespace) -> int:
    files = sorted(RESULTS.glob("*.json"))
    if not files:
        print("No results in", RESULTS, "- run some benchmarks first.")
        return 1
    for f in files:
        s = summarise(f)
        print(f"\n=== {f.stem}  ({s['n_repeats']} repeats, "
              f"wall {s['wall_seconds']:.1f}s, "
              f"peak {s['process_peak_rss_mb']:.0f} MB) ===")
        for name, v in s["steps"].items():
            flag = "" if v["anchors_agree"] else "  [anchor varies across repeats]"
            print(f"  {name:24s} {v['seconds']:8.2f}s  "
                  f"{v['peak_rss_mb']:7.0f} MB peak  "
                  f"{v['delta_rss_mb']:+7.0f} MB  anchor={v['anchor']}{flag}")
    return 0


# ---------------------------------------------------------------------------
# End-to-end tutorial scripts
# ---------------------------------------------------------------------------

# The tutorials as they actually ship, each a Python script and the R script
# that reproduces it. Measuring these answers a different question from the
# benches above: not "how fast is this operation" but "how long does running
# the tutorial take". The two are not interchangeable — every one of these
# scripts also prints validation, writes CSVs or draws figures, and the two
# sides do not do equal amounts of that. Read these as script wall-clock, and
# the benches for anything attributable to an algorithm.
#
# Order matters: several R scripts read a file the Python run writes (shared
# variable features, gene lists, cell assignments), so Python goes first.
TUTORIALS = {
    "pbmc3k":       ("pbmc3k_tutorial.py", "pbmc3k_verify.R"),
    "pbmc8k":       ("pbmc8k_subclustering_tutorial.py",
                     "pbmc8k_subclustering_verify.R"),
    "sctransform":  ("pbmc3k_sctransform_tutorial.py",
                     "pbmc3k_sctransform_verify.R"),
    "de":           ("pbmc3k_de_tutorial.py", "pbmc3k_de_verify.R"),
    "dimreduc":     ("pbmc3k_dimreduc_tutorial.py", "pbmc3k_dimreduc_verify.R"),
    "objects":      ("pbmc3k_objects_tutorial.py", "pbmc3k_objects_verify.R"),
    "integration":  ("ifnb_integration_tutorial.py", "ifnb_integration_verify.R"),
    "sketch":       ("ifnb_sketch_tutorial.py", "ifnb_sketch_verify.R"),
    "anchors":      ("anchors_tutorial.py", "anchors_verify.R"),
    "citeseq":      ("cbmc_citeseq_tutorial.py", "cbmc_citeseq_verify.R"),
    "hashing":      ("pbmc_hashing_tutorial.py", "pbmc_hashing_verify.R"),
    "mixscape":     ("thp1_mixscape_tutorial.py", "thp1_mixscape_verify.R"),
    "cellcycle":    ("thp1_cellcycle_tutorial.py", "thp1_cellcycle_verify.R"),
    "refmap":       ("panc8_reference_mapping_tutorial.py",
                     "panc8_reference_mapping_verify.R"),
    "svf":          ("xenium_svf_tutorial.py", "xenium_svf_verify.R"),
    "visium":       ("visium_tutorial.py", "visium_verify.R"),
    "lazy":         ("lazy_bpcells_tutorial.py", "lazy_bpcells_verify.R"),
}

TUTORIAL_DIR = ROOT / "tutorials"


def run_script(cmd: list[str], log_path: Path) -> dict:
    """Wall clock and peak RSS for one script, measured from outside."""
    t_start = time.time()
    with log_path.open("w") as log:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=log,
                                stderr=subprocess.STDOUT)
        sampler = Sampler(proc.pid)
        sampler.start()
        code = proc.wait()
        sampler.stop()
    return {
        "cmd": " ".join(Path(c).name if "/" in c else c for c in cmd),
        "exit_code": code,
        "wall_seconds": time.time() - t_start,
        "peak_rss_mb": max((kb for _, kb in sampler.samples), default=0) / 1024,
        "log": str(log_path.relative_to(ROOT)),
    }


def cmd_scripts(args: argparse.Namespace) -> int:
    RESULTS.mkdir(exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    names = (list(TUTORIALS) if args.tutorial == "all"
             else args.tutorial.split(","))
    # Merge into whatever is already there. `run` writes one file per bench, so
    # re-running one bench leaves the rest alone; this subcommand writes them
    # all to a single file, and rewriting it from scratch would silently drop
    # every tutorial not named on this invocation. Re-running one row of the
    # table has to stay a safe thing to do.
    dest = RESULTS / "tutorial_scripts.json"
    out = {"machine": machine(), "tutorials": {}}
    if dest.exists():
        out["tutorials"] = json.loads(dest.read_text()).get("tutorials", {})
    failures = 0
    for name in names:
        if name not in TUTORIALS:
            print(f"  unknown tutorial {name!r}; known: {', '.join(TUTORIALS)}")
            failures += 1
            continue
        py, r = TUTORIALS[name]
        pair = {}
        for arm, cmd in (
            ("truecell", [str(ROOT / ".venv" / "bin" / "python"),
                          str(TUTORIAL_DIR / py)]),
            ("seurat", ["/usr/local/bin/Rscript", str(TUTORIAL_DIR / r)]),
        ):
            print(f"  {name:12s} {arm:9s} ...", end="", flush=True)
            res = run_script(cmd, LOGS / f"tutorial.{name}.{arm}.log")
            ok = res["exit_code"] == 0
            print(f" {res['wall_seconds']:7.1f}s  peak {res['peak_rss_mb']:6.0f} MB"
                  f"{'' if ok else '   FAILED (see ' + res['log'] + ')'}")
            failures += 0 if ok else 1
            pair[arm] = res
        out["tutorials"][name] = pair
        # Canonical order, not measurement order — otherwise re-running one row
        # moves it to the end and the report's table silently reshuffles.
        out["tutorials"] = {k: out["tutorials"][k]
                            for k in TUTORIALS if k in out["tutorials"]}
        dest.write_text(json.dumps(out, indent=2))
    print(f"  -> {dest.relative_to(ROOT)}")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# Markdown tables
# ---------------------------------------------------------------------------

def _fmt_ratio(a: float, b: float) -> str:
    """``b`` relative to ``a``, written as the reader would say it out loud."""
    if a <= 0 or b <= 0:
        return "-"
    if a >= b:
        return f"**{a / b:.1f}x**" if a / b >= 1.05 else "~1x"
    return f"{b / a:.1f}x slower"


def cmd_tables(args: argparse.Namespace) -> int:
    """Emit the per-bench comparison tables the report is written around."""
    benches: dict[str, dict[str, dict]] = {}
    for f in sorted(RESULTS.glob("*.json")):
        parts = f.stem.split(".")
        if len(parts) < 2:
            continue  # tutorial_scripts.json — a different shape, see make_report
        bench, arm = parts[0], ".".join(parts[1:])
        benches.setdefault(bench, {})[arm] = summarise(f)

    lines: list[str] = []
    for bench, arms in benches.items():
        tc, sr = arms.get("truecell"), arms.get("seurat")
        if not (tc and sr):
            continue
        # The single-threaded truecell arm, when it was run. It is the control
        # for the BLAS gap: R's shipped BLAS is single-threaded, so a truecell
        # number measured on twelve cores answers a different question.
        t1 = arms.get("truecell.t1")
        lines.append(f"\n### `{bench}`\n")
        head = ["Step", "Truecell", "Truecell 1-thread", "Seurat", "Faster by",
                "Truecell peak", "Seurat peak", "Truecell anchor", "Seurat anchor"]
        if t1 is None:
            head.pop(2)
        lines.append("| " + " | ".join(head) + " |")
        lines.append("|---" + "|---:" * (len(head) - 1) + "|")
        names = list(tc["steps"]) + [n for n in sr["steps"] if n not in tc["steps"]]
        for name in names:
            a, b = tc["steps"].get(name), sr["steps"].get(name)
            # anchor == -1 is the benches' marker for "this call raised". The
            # step still has a duration, and it is the duration of a failure —
            # printing it in a speed column would read as a 200x win.
            failed = [arm for arm, v in (("truecell", a), ("Seurat", b))
                      if v is not None and v["anchor"] == -1]
            if failed:
                lines.append(f"| `{name}` | _unavailable in "
                             f"{' and '.join(failed)} — not timed_ "
                             + "| " * (len(head) - 2) + "|")
                continue
            row = [f"`{name}`", f"{a['seconds']:.2f}s" if a else "-"]
            if t1 is not None:
                c = t1["steps"].get(name)
                row.append(f"{c['seconds']:.2f}s" if c else "-")
            row += [
                f"{b['seconds']:.2f}s" if b else "-",
                _fmt_ratio(b["seconds"], a["seconds"]) if a and b else "-",
                f"{a['peak_rss_mb']:.0f} MB" if a else "-",
                f"{b['peak_rss_mb']:.0f} MB" if b else "-",
                "" if a is None or a["anchor"] is None else str(a["anchor"]),
                "" if b is None or b["anchor"] is None else str(b["anchor"]),
            ]
            lines.append("| " + " | ".join(row) + " |")
        total = ["**total (process)**", f"**{tc['wall_seconds']:.1f}s**"]
        if t1 is not None:
            total.append(f"**{t1['wall_seconds']:.1f}s**")
        total += [f"**{sr['wall_seconds']:.1f}s**",
                  _fmt_ratio(sr["wall_seconds"], tc["wall_seconds"]),
                  f"**{tc['process_peak_rss_mb']:.0f} MB**",
                  f"**{sr['process_peak_rss_mb']:.0f} MB**", "", ""]
        lines.append("| " + " | ".join(total) + " |")
    out = "\n".join(lines)
    (HERE / "tables.md").write_text(out)
    print(out)
    print(f"\n-> {(HERE / 'tables.md').relative_to(ROOT)}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run one or more benchmarks")
    r.add_argument("--bench", required=True,
                   help="comma-separated bench names (see bench_truecell.py)")
    r.add_argument("--arm", default="truecell,seurat",
                   help="comma-separated: truecell, seurat")
    r.add_argument("--repeats", type=int, default=3,
                   help="timed repeats, on top of the discarded warm-up")
    r.add_argument("--threads", type=int, default=None,
                   help="pin every thread pool to this many threads")
    r.add_argument("--keep-warmup", action="store_true",
                   help="keep the first repeat instead of discarding it")
    r.set_defaults(func=cmd_run)

    rep = sub.add_parser("report", help="summarise collected results")
    rep.set_defaults(func=cmd_report)

    sc = sub.add_parser("scripts",
                        help="time the tutorial scripts themselves, end to end")
    sc.add_argument("--tutorial", default="all",
                    help="comma-separated tutorial keys, or 'all'")
    sc.set_defaults(func=cmd_scripts)

    tab = sub.add_parser("tables", help="write the side-by-side markdown tables")
    tab.set_defaults(func=cmd_tables)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
