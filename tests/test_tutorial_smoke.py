"""End-to-end tutorial runs.

**These do not run in CI, and are not meant to be mistaken for coverage.** They
need the cached datasets (~200 MB under ``~/.shanuz_data``) and take minutes, so
they are gated behind an explicit opt-in:

    SHANUZ_TUTORIAL_SMOKE=1 pytest tests/test_tutorial_smoke.py -v

The gate is an env var rather than a bare "skip if the data is missing" check on
purpose. A test that silently skips wherever the data happens to be absent reads
as green in CI while proving nothing — which is the failure mode that let
``pbmc3k_tutorial.py`` ship broken in the first place. Requiring the opt-in means
a skip here always means "nobody asked for this", never "this passed".

What they cover that ``test_tutorial_marker_tables.py`` cannot: the actual script
path, top to bottom, against real data. The pandas 3 regression lived in the
middle of ``run_tutorial()`` — importable helpers and unit fixtures would not
have reached it. Worth running before cutting a release.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = Path(os.path.expanduser("~/.shanuz_data"))

OPT_IN = os.environ.get("SHANUZ_TUTORIAL_SMOKE") == "1"

# script -> the dataset directory it needs cached
TUTORIALS = [
    ("pbmc3k_tutorial.py", "pbmc3k"),
    ("pbmc8k_subclustering_tutorial.py", "pbmc8k"),
    ("cbmc_citeseq_tutorial.py", "cbmc"),
    ("pbmc3k_sctransform_tutorial.py", "pbmc3k"),
]

pytestmark = pytest.mark.skipif(
    not OPT_IN,
    reason="set SHANUZ_TUTORIAL_SMOKE=1 to run the tutorials end-to-end "
           "(needs ~200MB of cached datasets, takes minutes)",
)


@pytest.mark.parametrize("script,dataset", TUTORIALS)
def test_tutorial_runs_to_completion(script, dataset):
    """The script exits 0 and reaches its final section.

    Deliberately asserts on the exit code and the tail of stdout rather than on
    any number the tutorial prints: cluster counts and marker rankings are
    library-version dependent (see the R-comparison notes), and pinning them here
    would turn a dependency bump into a test failure that says nothing about the
    tutorial. What this pins is that the script *runs*.
    """
    if not (DATA_ROOT / dataset).is_dir():
        pytest.skip(f"dataset {dataset!r} not cached under {DATA_ROOT}")

    proc = subprocess.run(
        [sys.executable, f"tutorials/{script}"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, (
        f"tutorials/{script} exited {proc.returncode}\n"
        f"--- last 40 lines of stdout ---\n{os.linesep.join(proc.stdout.splitlines()[-40:])}\n"
        f"--- stderr ---\n{proc.stderr[-3000:]}"
    )
    assert proc.stdout.strip(), f"tutorials/{script} produced no output"


@pytest.mark.parametrize("script,dataset", [TUTORIALS[0]])
def test_pbmc3k_prints_the_marker_table(script, dataset):
    """The exact block that the pandas 3 regression killed.

    ``run_tutorial`` crashed with KeyError: 'cluster' immediately after printing
    the "Top 3 markers per cluster:" header, so a run that reaches the header but
    dies right after still looks half-right in a log. Assert a cluster line
    actually follows it.
    """
    if not (DATA_ROOT / dataset).is_dir():
        pytest.skip(f"dataset {dataset!r} not cached under {DATA_ROOT}")

    proc = subprocess.run(
        [sys.executable, f"tutorials/{script}"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    lines = proc.stdout.splitlines()
    header = next(
        (i for i, ln in enumerate(lines) if "Top 3 markers per cluster" in ln), None
    )
    assert header is not None, "tutorial never reached the marker table"
    cluster_lines = [
        ln for ln in lines[header + 1: header + 12] if ln.strip().startswith("Cluster ")
    ]
    assert len(cluster_lines) >= 8, (
        f"expected a line per cluster after the header, got {cluster_lines!r}"
    )
    # "Cluster 0: GENE, GENE, GENE" — genes present, not an empty table
    assert all(":" in ln and ln.split(":", 1)[1].strip() for ln in cluster_lines)
