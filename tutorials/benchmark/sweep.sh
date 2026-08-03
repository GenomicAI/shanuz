#!/usr/bin/env bash
# The full benchmark sweep, in the order the report is written from.
#
#   bash tutorials/benchmark/sweep.sh
#
# Takes about an hour on an M4 Pro. Run it on an otherwise idle machine — the
# parent samples resident set size from outside the child, so anything else
# competing for cores or memory lands in the numbers.
#
# Order is not arbitrary. Within every bench the truecell arm runs first,
# because it writes the cell-to-cluster assignment and the Xenium cell subset
# that the R arm reads back; that is what keeps the two arms doing the same
# work rather than merely the same named step.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
RUN="$PY tutorials/benchmark/run_benchmarks.py"

rm -rf tutorials/benchmark/results tutorials/benchmark/logs
find tutorials/benchmark -maxdepth 1 -name '.steps.*.jsonl' -delete

echo "== 1/5  BLAS control"
$RUN run --bench blas_probe --repeats 3
$RUN run --bench blas_probe --arm truecell --threads 1 --repeats 3

echo "== 2/5  standard workflow, 2.7k to 20.7k cells"
$RUN run --bench pbmc3k_core,pbmc8k_core,ifnb_core,thp1_core --repeats 3

echo "== 3/5  standard workflow, truecell pinned to one thread"
$RUN run --bench pbmc3k_core,pbmc8k_core,ifnb_core --arm truecell \
    --threads 1 --repeats 3

echo "== 4/5  named heavy operations"
$RUN run --bench pbmc3k_sctransform,pbmc3k_de,ifnb_integration,xenium_spatial \
    --repeats 2

echo "== 5/5  the tutorial scripts themselves, end to end"
$RUN scripts --tutorial pbmc3k,sctransform,de,dimreduc,objects,integration,cellcycle,svf,visium,lazy

$RUN tables
echo "done"
