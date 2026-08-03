# Benchmarks

The performance counterpart to the accuracy work in [`../README.md`](../README.md).
That asks whether Truecell and R Seurat produce the same answers. This asks what
each answer costs, in seconds and in bytes.

The results and the discussion are in [`PERFORMANCE.md`](PERFORMANCE.md).

## Running it

```bash
bash tutorials/benchmark/sweep.sh
```

About 40 minutes on an M4 Pro. Run it on an otherwise idle machine — memory is
sampled from outside the process, so anything else competing for cores or RAM
lands in the numbers.

> **The last phase rewrites tutorial figures.** `run_benchmarks.py scripts`
> runs the tutorial scripts as they ship, and those scripts draw figures into
> `tutorials/figures*/` — which are committed. A sweep will therefore leave a
> dirty working tree with a dozen or so modified PNGs that have nothing to do
> with what you were measuring. Check `git status` afterwards and
> `git checkout -- tutorials/figures*` unless you actually meant to regenerate
> them. Needs the datasets cached (any tutorial run downloads
them) and an R with Seurat, plus `presto`, `harmony`, `RANN`, `Rfast2` and
`data.table`.

Individual pieces:

```bash
python tutorials/benchmark/run_benchmarks.py run --bench pbmc3k_core
python tutorials/benchmark/run_benchmarks.py run --bench ifnb_core --arm truecell --threads 1
python tutorials/benchmark/run_benchmarks.py scripts --tutorial pbmc3k,de
python tutorials/benchmark/make_report.py
python tutorials/benchmark/compare_blas.py results_refblas results
```

`results_refblas/` holds a full sweep taken before R's BLAS was switched from
the reference build it ships with to Accelerate. Keeping it is what let the
report show that the swap changed no result anywhere while taking up to 12.9x
off individual steps.

## How it works

| File | Role |
|---|---|
| `run_benchmarks.py` | Spawns each child, samples its process-tree RSS, joins the samples to the child's step log, writes `results/*.json` |
| `bench_truecell.py` | The Truecell arm — eight benches |
| `bench_seurat.R` | The Seurat arm — the same eight, step for step |
| `steps.py` / `steps.R` | The step recorder each arm writes its boundaries with |
| `make_report.py` | Turns `results/*.json` into the tables in `PERFORMANCE.md` |
| `sweep.sh` | The full run, in the order the report is written from |
| `compare_blas.py` | Diffs two sweeps — anchors first, then timings. Used for the reference-BLAS vs Accelerate comparison in `PERFORMANCE.md` section 2.1 |

**Time is measured inside each process; memory from outside.** R's `gc()` and
Python's `tracemalloc` measure different things and neither sees the other's
allocator, so a memory figure from either would not be comparable across the
two. The parent samples resident set size every 50 ms instead, and both arms
timestamp their steps with the same epoch clock (`time.time()` /
`as.numeric(Sys.time())`) so the samples can be attributed afterwards.

**Every step reports an anchor** — cells kept, clusters found, markers
returned. They are printed beside the timings so a reader can check that the
two arms did the same work before believing that one did it faster.

**The Truecell arm must run first.** It writes `results/<bench>_idents.csv` and
`results/xenium_cells.txt`, which the R arm reads back. Without that the two
sides would cluster differently and then time a different number of one-vs-rest
marker tests against each other.

## Adding a bench

Add a function to `BENCHES` in *both* arms under the same key, with the same
steps in the same order. Where the two stacks genuinely cannot do the same
thing, add a **separate step** rather than substituting quietly — `bench_seurat.R`
has `neighbours_annoy` because Seurat's default neighbour search is approximate
and Truecell's is exact, and reporting only one of the two would be a choice
about which tool to flatter.
