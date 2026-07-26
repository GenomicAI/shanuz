"""Guards for the declared bands and the handoff check they rest on.

The bands exist because prose does not fail. These exist because a band that
cannot fail is prose with extra steps — every assertion here is written so that
removing the thing it guards makes it go red, and the ones that matter were
checked that way rather than assumed.
"""
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tutorials.bands import (  # noqa: E402
    Band, StaleReferenceError, check_bands, check_same_cells,
    check_same_features, check_shared_groups, render_verdicts,
    summarise_spread, PCT_ROUNDING_TOL,
)


# ---------------------------------------------------------------------------
# Band
# ---------------------------------------------------------------------------

def test_band_is_inclusive_at_both_ends():
    """The measured values sit *on* two of the boundaries.

    `roc max |dAUC|` is bounded by R's own rounding at 5e-4 and measures
    0.00049986; an exclusive comparison would fail the run for being correct.
    """
    band = Band(0, 5e-4, "R rounds myAUC to three decimals")
    assert band.holds(0)
    assert band.holds(5e-4)
    assert not band.holds(5e-4 + 1e-9)


def test_band_rejects_nan():
    """A missing measurement must fail, not pass.

    This is the whole failure mode the module was written for: the dim-reduction
    report reindexed R's embedding onto a cell set it did not cover, got NaN
    correlations, and printed them without comment.
    """
    assert not Band(0, 1, "any").holds(float("nan"))


def test_band_refuses_to_be_empty_or_unexplained():
    with pytest.raises(ValueError, match="empty band"):
        Band(1, 0, "backwards")
    with pytest.raises(ValueError, match="why it is where it is"):
        Band(0, 1, "   ")


def test_exact_band_admits_one_value():
    band = Band(50, 50, "the same 50 genes")
    assert band.exact and band.describe() == "= 50"
    assert band.holds(50) and not band.holds(49)


def test_describe_reads_as_a_bound_when_one_side_is_open():
    assert Band(0.8, math.inf, "floor").describe().startswith(">=")
    assert Band(-math.inf, 2, "ceiling").describe().startswith("<=")


# ---------------------------------------------------------------------------
# check_bands / render_verdicts
# ---------------------------------------------------------------------------

def test_a_band_with_no_measurement_fails():
    """Deleting the line that computes a number must not retire its band.

    Otherwise the cheapest way to make `--report` green is to stop measuring.
    """
    verdicts = check_bands({"never computed": Band(0, 1, "x")}, {})
    assert len(verdicts) == 1
    assert not verdicts[0].ok
    assert math.isnan(verdicts[0].value)


def test_check_bands_judges_each_declared_name():
    bands = {"a": Band(0, 1, "x"), "b": Band(0, 1, "y")}
    verdicts = check_bands(bands, {"a": 0.5, "b": 7.0, "c": 0.0})
    assert [v.name for v in verdicts] == ["a", "b"]      # 'c' has no band
    assert [v.ok for v in verdicts] == [True, False]


def test_render_returns_false_and_names_the_reason(capsys):
    """A failure has to print *why* the band is where it is.

    A bare "OUT" tells the next reader nothing about whether to widen the band
    or fix the port.
    """
    bands = {"deseq2 top50": Band(15, 32, "pseudobulk against per-cell")}
    ok = render_verdicts(check_bands(bands, {"deseq2 top50": 48}), "T")
    out = capsys.readouterr().out
    assert ok is False
    assert "OUT" in out and "pseudobulk against per-cell" in out


def test_render_is_true_when_everything_holds(capsys):
    ok = render_verdicts(check_bands({"a": Band(0, 1, "x")}, {"a": 0.5}), "T")
    assert ok is True
    assert "OUT" not in capsys.readouterr().out


# ---------------------------------------------------------------------------
# check_shared_groups — the DE handoff
# ---------------------------------------------------------------------------

def _marker_table(pct1, pct2, genes=("A", "B", "C")):
    return pd.DataFrame({"pct.1": pct1, "pct.2": pct2}, index=list(genes))


def test_shared_groups_accepts_r_rounding():
    """Seurat rounds pct to three decimals; that is not staleness."""
    py = _marker_table([0.1234, 0.5006, 0.9], [0.2, 0.3, 0.4])
    r = _marker_table([0.123, 0.501, 0.9], [0.2, 0.3, 0.4])
    worst = check_shared_groups(py, r, source="r_wilcox.csv")
    assert worst <= PCT_ROUNDING_TOL


def test_shared_groups_rejects_a_different_cluster_assignment():
    """The real failure: 25 July Python against a 19 July R run.

    On the working copy this was reached with pct.1 differing by up to 0.0174
    across 12,491 of 13,712 genes, and the report printed a full parity table
    with wilcox at 48/50 instead of 50/50 — indistinguishable from a regression.
    """
    py = _marker_table([0.1234, 0.5, 0.9], [0.2, 0.3, 0.4])
    r = _marker_table([0.141, 0.5, 0.9], [0.2, 0.3, 0.4])
    with pytest.raises(StaleReferenceError, match="pct.1 differs by up to"):
        check_shared_groups(py, r, source="figures_de/r_wilcox.csv")


def test_shared_groups_names_the_command_that_fixes_it():
    py = _marker_table([0.1], [0.9], genes=("A",))
    r = _marker_table([0.9], [0.1], genes=("A",))
    with pytest.raises(StaleReferenceError) as exc:
        check_shared_groups(py, r, source="figures_de/r_mast.csv")
    assert "pbmc3k_de_verify.R" in str(exc.value)


def test_shared_groups_rejects_a_table_with_no_genes_in_common():
    py = _marker_table([0.1], [0.9], genes=("A",))
    r = _marker_table([0.1], [0.9], genes=("Z",))
    with pytest.raises(StaleReferenceError, match="shares no genes"):
        check_shared_groups(py, r, source="r_t.csv")


def test_shared_groups_checks_both_pct_columns():
    """pct.2 alone moving is still a different group assignment.

    Guarding only pct.1 would miss a change confined to the second group, which
    is exactly half the ways the handoff can go stale.
    """
    py = _marker_table([0.5, 0.5], [0.10, 0.5], genes=("A", "B"))
    r = _marker_table([0.5, 0.5], [0.31, 0.5], genes=("A", "B"))
    with pytest.raises(StaleReferenceError, match="pct.2"):
        check_shared_groups(py, r, source="r_bimod.csv")


# ---------------------------------------------------------------------------
# check_same_cells / check_same_features — the dim-reduction handoff
# ---------------------------------------------------------------------------

def test_same_cells_passes_when_the_sets_match_in_any_order():
    frame = pd.DataFrame({"PC_1": [1.0, 2.0]}, index=["c2", "c1"])
    check_same_cells(["c1", "c2"], frame, source="r_pca.csv")


def test_same_cells_rejects_a_reindex_that_would_have_gone_to_nan():
    """`reindex` fills a missing barcode with NaN instead of raising.

    Without this the correlations downstream come out NaN and print as `nan`,
    which is why `Band.holds` also treats NaN as a failure — two guards, because
    either one alone leaves the number looking merely unavailable.
    """
    frame = pd.DataFrame({"PC_1": [1.0, 2.0]}, index=["c1", "c9"])
    with pytest.raises(StaleReferenceError, match="different cell set"):
        check_same_cells(["c1", "c2"], frame, source="figures_dimreduc/r_pca.csv")


def test_same_features_rejects_extras_as_well_as_absences():
    """Set equality, not one-sided containment.

    The check this replaced reindexed the Python features onto R's table and
    complained only about the ones that came back empty, so an R run holding
    *additional* features — a different HVG selection of the same size — passed.
    """
    with pytest.raises(StaleReferenceError, match="unexpected"):
        check_same_features(["A", "B"], ["A", "B", "C"], source="r_jackstraw_p.csv")


def test_same_features_applies_the_r_renaming_rule():
    """Read10X turns `Y_RNA` into `Y-RNA`; that is not a different feature set."""
    check_same_features(["Y_RNA", "CD14"], ["Y-RNA", "CD14"],
                        source="r_jackstraw_p.csv",
                        key=lambda s: s.replace("_", "-"))


def test_same_features_without_a_key_treats_the_renaming_as_drift():
    with pytest.raises(StaleReferenceError):
        check_same_features(["Y_RNA"], ["Y-RNA"], source="r_jackstraw_p.csv")


# ---------------------------------------------------------------------------
# summarise_spread
# ---------------------------------------------------------------------------

def test_summarise_spread_reports_the_gap_to_the_reference():
    """The JackStraw band is |shanuz - R|, so the gap is the number that matters."""
    got = summarise_spread([12, 13, 13, 14, 15], reference=13)
    assert got["n"] == 5 and got["min"] == 12 and got["max"] == 15
    assert got["median"] == 13
    assert got["max_abs_gap"] == 2


def test_summarise_spread_refuses_an_empty_sweep():
    with pytest.raises(ValueError, match="no values"):
        summarise_spread([])
