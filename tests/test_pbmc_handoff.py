"""Guards for the PBMC 3k / PBMC 8k numeric handoffs against R Seurat.

Tutorials 1 and 2 were the last two in the suite compared entirely by eye:
every figure in `pbmc3k_tutorial.md` links the canonical satijalab.org image,
so nothing failed if the numbers behind them drifted. These cover the pieces of
that handoff that can be exercised without R or a dataset download — the
cluster-matching arithmetic, the gene-symbol normalisation, and the contract
between each `--report` and the R script that feeds it.
"""
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tutorials.pbmc3k_tutorial import _r_symbols, match_partitions  # noqa: E402

TUTORIALS = Path(__file__).resolve().parent.parent / "tutorials"


# ---------------------------------------------------------------------------
# match_partitions — the arithmetic both reports lean on
# ---------------------------------------------------------------------------

def test_match_partitions_is_blind_to_relabelling():
    """Cluster ids are arbitrary, so a pure renaming must score perfectly.

    This is the whole reason the matching exists. Comparing the label columns
    directly, shanuz's cluster 3 against Seurat's cluster 3, answers a question
    neither tool makes a promise about; on PBMC 3k the two runs agree about
    2,554 of 2,638 cells while numbering one of the clusters differently.
    """
    a = ["0"] * 10 + ["1"] * 10 + ["2"] * 5
    b = ["2"] * 10 + ["0"] * 10 + ["1"] * 5      # same partition, renamed
    got = match_partitions(a, b)
    assert got["concordance"] == pytest.approx(1.0)
    assert got["ari"] == pytest.approx(1.0)
    assert got["mapping"] == {"0": "2", "1": "0", "2": "1"}


def test_match_partitions_scores_a_genuine_split_below_one():
    """One cluster split in two must not read as agreement.

    Concordance alone is forgiving here — the split half still lands inside its
    matched pair — which is exactly why the report prints ARI beside it.
    """
    a = ["0"] * 20 + ["1"] * 20
    b = ["0"] * 20 + ["1"] * 10 + ["2"] * 10     # cluster 1 split in half
    got = match_partitions(a, b)
    assert got["n_a"] == 2 and got["n_b"] == 3
    assert got["concordance"] == pytest.approx(0.75)
    assert got["ari"] < 0.75, (
        "ARI must penalise the split harder than best-match concordance does; "
        "if it does not, printing both tells the reader nothing"
    )


def test_match_partitions_counts_only_cells_inside_matched_pairs():
    """Concordance is the matched-pair cell count over the total."""
    a = ["0"] * 6 + ["1"] * 4
    b = ["x"] * 5 + ["y"] * 5                    # one cell crosses the boundary
    got = match_partitions(a, b)
    assert got["concordance"] == pytest.approx(0.9)
    assert got["table"].loc["0", "x"] == 5
    assert got["table"].loc["1", "y"] == 4


def test_match_partitions_rejects_mismatched_inputs():
    with pytest.raises(ValueError):
        match_partitions(["0", "1"], ["0"])
    with pytest.raises(ValueError):
        match_partitions([], [])


# ---------------------------------------------------------------------------
# Gene symbols
# ---------------------------------------------------------------------------

def test_gene_symbols_are_mapped_to_reads10x_spelling():
    """R's Read10X() rewrites underscores to dashes; shanuz's loader does not.

    Without this the two per-gene tables join on ~30 fewer genes than they
    have, and every gene whose symbol contains an underscore silently drops out
    of the VST and marker comparisons. The SCTransform handoff hit the same
    thing as an off-by-one in the shared-gene count.
    """
    assert _r_symbols(["RP11-34P13_3", "MIR1302_2", "LYZ"]) == \
        ["RP11-34P13-3", "MIR1302-2", "LYZ"]
    # idempotent: running it on already-mapped symbols is a no-op
    assert _r_symbols(_r_symbols(["A_B_C"])) == ["A-B-C"]


# ---------------------------------------------------------------------------
# The R references
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("script", [
    "pbmc3k_objects_verify.R",
    "pbmc3k_verify.R",
    "pbmc8k_subclustering_verify.R",
])
def test_r_references_pin_exact_neighbours(script):
    """No R reference may fall back to Seurat's approximate `annoy`.

    `FindNeighbors` defaults to `nn.method = "annoy"`, which is approximate,
    while this port's neighbour search is exact. With the default, the two
    sides build their graphs from different neighbour tables and report a
    difference that belongs to annoy — 182 SNN edges on pbmc3k — which reads as
    a shanuz defect and cost a real investigation once already. On PBMC 8k it
    is worse than cosmetic: the graph decides the global clusters, which decide
    which cells enter the subclustering stage.

    This asserts the script text rather than running R, which CI has no Seurat
    for. It therefore proves only that the pin is present, not that Seurat
    honours it — but "the pin was silently dropped" is the regression that
    actually happened, and it is the one this catches.
    """
    text = (TUTORIALS / script).read_text()
    calls = [ln for ln in text.splitlines() if "FindNeighbors(" in ln]
    assert calls, f"FindNeighbors call vanished from {script}"
    for call in calls:
        assert 'nn.method = "rann"' in call, (
            f"{script} uses Seurat's approximate default: {call.strip()!r}"
        )


@pytest.mark.parametrize("module,script", [
    ("tutorials.pbmc3k_tutorial", "pbmc3k_verify.R"),
    ("tutorials.pbmc8k_subclustering_tutorial", "pbmc8k_subclustering_verify.R"),
])
def test_every_file_the_report_reads_is_written_by_one_of_the_two_sides(module, script):
    """The report's inputs and the two producers must not drift apart.

    A `--report` that names a file nobody writes degrades to the "missing ...
    run the tutorial first" branch and prints nothing — a silent no-op that
    looks exactly like the handoff not having been run yet. Rather than trust
    the file list by eye, read the names out of `report()` itself and check
    each one against the side that is supposed to produce it.
    """
    import importlib

    mod = importlib.import_module(module)
    src = Path(mod.__file__).read_text()
    body = src[src.index("def report("):]
    names = set(re.findall(r'"((?:py|r)_[A-Za-z0-9_.]+\.(?:csv|json))"', body))
    assert names, f"no handoff filenames found in {module}.report()"
    assert any(n.startswith("r_") for n in names), "report reads no R-side files"
    assert any(n.startswith("py_") for n in names), "report reads no Python-side files"

    r_text = (TUTORIALS / script).read_text()
    for name in sorted(n for n in names if n.startswith("r_")):
        assert name in r_text, f"{script} never writes {name}, which report() reads"
    for name in sorted(n for n in names if n.startswith("py_")):
        assert name in src, f"{module} never writes {name}, which report() reads"


@pytest.mark.parametrize("script", ["pbmc3k_verify.R", "pbmc8k_subclustering_verify.R"])
def test_r_references_serialise_json_at_full_precision(script):
    """`jsonlite::toJSON` defaults to 4 significant digits.

    The anchors carry sums over millions of values — `data_sum` is ~4.5e6 on
    PBMC 3k — and at the default the R side would round them to four digits
    while Python wrote all seventeen, so a relative difference of ~1e-13 would
    print as ~1e-5 and every scalar would look like a near-miss. `digits = 22`
    is what round-trips; `NA` does not.
    """
    text = (TUTORIALS / script).read_text()
    calls = [ln for ln in text.splitlines() if "toJSON(" in ln]
    assert calls, f"{script} no longer writes an anchors JSON"
    for call in calls:
        assert "digits = 22" in call, (
            f"{script} serialises JSON at reduced precision: {call.strip()!r}"
        )


def test_reports_are_reachable_without_a_dataset():
    """`--report` must not need the 24 MB download to tell you what is missing."""
    import importlib

    for module in ("tutorials.pbmc3k_tutorial",
                   "tutorials.pbmc8k_subclustering_tutorial"):
        mod = importlib.import_module(module)
        assert callable(mod.report)
        assert callable(mod.write_anchors)
        assert isinstance(mod.FIGURES, Path)


def test_numeric_anchor_names_agree_across_the_two_sides():
    """The scalars the Python side writes and the R side writes must line up.

    Both scripts hand-build their anchors dict, so a key renamed on one side
    would just vanish from the comparison loop — which skips anything missing
    from either — rather than failing. Checking the literal key lists keeps the
    two in step.
    """
    for module, script in (
        ("tutorials.pbmc3k_tutorial", "pbmc3k_verify.R"),
        ("tutorials.pbmc8k_subclustering_tutorial", "pbmc8k_subclustering_verify.R"),
    ):
        import importlib

        src = Path(importlib.import_module(module).__file__).read_text()
        block = src[src.index("    anchors = {"):]
        block = block[:block.index("\n    }")]
        py_keys = set(re.findall(r'^\s+"([a-z0-9_]+)":', block, re.M))

        r_text = (TUTORIALS / script).read_text()
        r_block = r_text[r_text.index("anchors <- list("):]
        r_block = r_block[:r_block.index("\n)")]
        r_keys = set(re.findall(r'^\s*([a-z0-9_.]+)\s*=', r_block, re.M))

        assert py_keys, f"no anchors parsed from {module}"
        assert py_keys == r_keys, (
            f"{module} and {script} disagree about the anchor set: "
            f"only python {sorted(py_keys - r_keys)}, only R {sorted(r_keys - py_keys)}"
        )
