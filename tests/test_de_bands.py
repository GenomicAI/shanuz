"""The DE tutorial's declared bands, and the guard that has to pass first.

`de_vignette.md` carries a parity table — top-50 overlap, p-value Spearman,
max |Δlog2FC| — that nothing checked. `deseq2`'s row had drifted from 25/50 to
22/50 unnoticed, and the same silence covered every other row. These tests drive
`measure_bands` and `BANDS` on hand-built tables, so each assertion can be shown
to fail rather than assumed to.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.pbmc3k_de_tutorial as de  # noqa: E402
from tutorials.bands import check_bands  # noqa: E402


def _clean_table() -> pd.DataFrame:
    """The concordance table as it reads on a good run, measured 2026-07-26."""
    rows = {
        "wilcox":   (50, 1.000000, 6.44e-15, np.nan),
        "t":        (50, 1.000000, 6.44e-15, np.nan),
        "bimod":    (50, 1.000000, 6.44e-15, np.nan),
        "LR":       (50, 1.000000, 6.44e-15, np.nan),
        "negbinom": (50, 0.916463, 6.44e-15, np.nan),
        "roc":      (np.nan, np.nan, 6.44e-15, 4.9986e-4),
        "mast":     (50, 0.997925, 6.44e-15, np.nan),
        "deseq2":   (22, 0.195879, 3.47, np.nan),
    }
    return pd.DataFrame(
        [{"test": k, f"top{de.TOP_N}_overlap": v[0], "p_spearman_expressed": v[1],
          "log2fc_max_abs_diff": v[2], "auc_max_abs_diff": v[3]}
         for k, v in rows.items()]).set_index("test")


def _holds(table) -> list[str]:
    """Names of the bands that fail on ``table``."""
    return [v.name for v in check_bands(de.BANDS, de.measure_bands(table))
            if not v.ok]


# ---------------------------------------------------------------------------
# The bands, against the numbers actually measured
# ---------------------------------------------------------------------------

def test_every_band_holds_on_the_measured_run():
    """The declared bands must admit the run they were derived from.

    A band that fails on a good run is worse than no band: it trains the reader
    to ignore the section.
    """
    assert _holds(_clean_table()) == []


def test_every_band_is_fed_by_measure_bands():
    """No band may be declared without something computing it.

    `check_bands` scores an absent name NaN and NaN fails, so this would show up
    as a mystery failure at `--report` time; here it names the band.
    """
    measured = de.measure_bands(_clean_table())
    assert set(de.BANDS) <= set(measured), (
        f"declared but never measured: {sorted(set(de.BANDS) - set(measured))}")


def test_a_single_dropped_gene_fails_the_parity_band():
    """49/50 is a regression, not scatter — the band has to say so.

    The exact bands are the ones most likely to be dismissed as brittle, so this
    pins the intent: for a test that shares Seurat's statistic and cells, one
    gene is a failure.
    """
    table = _clean_table()
    table.loc["wilcox", f"top{de.TOP_N}_overlap"] = 49
    assert _holds(table) == ["wilcox top50"]


def test_deseq2_reaching_parity_fails_too():
    """The upper bound is the load-bearing half of the deseq2 band.

    A jump to 50/50 would not read as a problem anywhere else in the report, but
    it would mean `sample_col` had stopped being honoured and the pseudobulk
    aggregation was no longer happening — a silent revert to the per-cell test
    Squair et al. warn against.
    """
    table = _clean_table()
    table.loc["deseq2", f"top{de.TOP_N}_overlap"] = 50
    assert "deseq2 top50" in _holds(table)


def test_the_deseq2_band_admits_both_recorded_measurements():
    """22 today, 25 on the previous cluster assignment; both are good runs."""
    band = de.BANDS["deseq2 top50"]
    assert band.holds(22) and band.holds(25)
    assert band.holds(20) and band.holds(26)     # ends of the resampling sweep
    assert not band.holds(5)                     # a collapse


def test_the_fold_change_band_excludes_deseq2_by_name():
    """deseq2's 3.47 must not set the tolerance for the other seven.

    Taking a max over every row would let a real cell-level fold-change
    regression hide under the pseudobulk row, which is correct at ~3.5.
    """
    measured = de.measure_bands(_clean_table())
    assert measured["max |dlog2FC| (parity tests)"] == pytest.approx(6.44e-15)

    table = _clean_table()
    table.loc["mast", "log2fc_max_abs_diff"] = 1e-6
    assert "max |dlog2FC| (parity tests)" in _holds(table)


def test_the_auc_band_is_r_s_rounding_and_not_a_free_tolerance():
    """5e-4 is half a unit in Seurat's third decimal, not a chosen slack."""
    assert de.BANDS["roc max |dAUC|"].high == de.AUC_TOLERANCE
    table = _clean_table()
    table.loc["roc", "auc_max_abs_diff"] = 1e-3
    assert "roc max |dAUC|" in _holds(table)


def test_a_missing_column_fails_rather_than_disappearing():
    """Dropping the measurement must not be a way to pass.

    `compare` only emits `top50_overlap` when enough genes survive the NaN and
    underflow filters, so this is reachable without anyone editing a band.
    """
    table = _clean_table().drop(columns=[f"top{de.TOP_N}_overlap"])
    failed = _holds(table)
    assert "wilcox top50" in failed and "deseq2 top50" in failed


# ---------------------------------------------------------------------------
# The handoff guard, in the report path
# ---------------------------------------------------------------------------

def _marker_csv(path, pct1, pct2, genes=("A", "B")):
    # myAUC rides along for the roc tables; the other tests ignore the column.
    pd.DataFrame({"p_val": [0.1] * len(genes), "avg_log2FC": [0.0] * len(genes),
                  "myAUC": [0.5] * len(genes),
                  "pct.1": pct1, "pct.2": pct2,
                  "p_val_adj": [0.1] * len(genes)},
                 index=list(genes)).rename_axis("gene").to_csv(path)


def test_report_concordance_refuses_a_stale_r_run(tmp_path, monkeypatch):
    """End to end: a mismatched pct column stops the report, not just a helper.

    The helper being right is not enough — the defect was that nothing called
    it. Damaging only the R side of one test is exactly the 19-July-vs-25-July
    situation, in miniature.
    """
    from tutorials.bands import StaleReferenceError

    monkeypatch.setattr(de, "FIGURES", tmp_path)
    for test in de.TEST_MAP:
        _marker_csv(tmp_path / f"py_{test}.csv", [0.500, 0.400], [0.300, 0.200])
        _marker_csv(tmp_path / f"r_{test.lower()}.csv", [0.500, 0.400], [0.300, 0.200])
    _marker_csv(tmp_path / "r_mast.csv", [0.500, 0.470], [0.300, 0.200])

    with pytest.raises(StaleReferenceError, match="not computed on the current"):
        de.report_concordance()


def test_report_concordance_accepts_a_reference_within_r_s_rounding(tmp_path,
                                                                    monkeypatch):
    """The guard must not fire on R's three-decimal pct columns.

    Seurat rounds; truecell does not. If this were exact, every real run would be
    rejected as stale.
    """
    monkeypatch.setattr(de, "FIGURES", tmp_path)
    for test in de.TEST_MAP:
        _marker_csv(tmp_path / f"py_{test}.csv", [0.50049, 0.4], [0.3, 0.2])
        _marker_csv(tmp_path / f"r_{test.lower()}.csv", [0.500, 0.4], [0.3, 0.2])

    table = de.report_concordance()
    assert list(table.index) == list(de.TEST_MAP)
