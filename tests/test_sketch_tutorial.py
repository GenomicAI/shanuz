"""Tests for the leverage-score sketching tutorial's helpers and pipeline.

Covers ``tutorials/ifnb_sketch_tutorial.py``:

  * the pure summary helpers (``group_enrichment``, ``rarity_correlation``,
    ``sketch_fold_change``, ``agreement``, ``label_accuracy``) on synthetic
    frames, where the right answer can be written down by hand;
  * the pipeline itself on a small synthetic object, so the wiring is exercised
    without the 14,000-cell download.

The real-data numbers live in ``tests/test_tutorial_smoke.py`` and only run when
``SHANUZ_TUTORIAL_SMOKE=1`` is set. Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.preprocessing import normalize_data, find_variable_features  # noqa: E402
from tutorials.ifnb_sketch_tutorial import (  # noqa: E402
    agreement,
    composition,
    group_enrichment,
    group_means,
    label_accuracy,
    rarity_correlation,
    sketch_fold_change,
    _r_feature_key,
    _read_r_scores,
)


# ----------------------------------------------------------------------
# group summaries
# ----------------------------------------------------------------------


def test_group_means_averages_within_each_group():
    values = [1.0, 3.0, 10.0, 20.0]
    groups = ["a", "a", "b", "b"]

    got = group_means(values, groups)

    assert got["a"] == pytest.approx(2.0)
    assert got["b"] == pytest.approx(15.0)


def test_group_enrichment_is_relative_to_the_overall_mean():
    # Overall mean 2.0; group "hi" averages 4.0 and "lo" averages 1.0.
    values = [4.0, 4.0, 1.0, 1.0, 1.0, 1.0]
    groups = ["hi", "hi", "lo", "lo", "lo", "lo"]

    table = group_enrichment(values, groups)

    assert table.loc["hi", "vs_overall"] == pytest.approx(2.0)
    assert table.loc["lo", "vs_overall"] == pytest.approx(0.5)
    assert table.loc["hi", "n"] == 2
    # Sorted commonest-first, which is how the report reads it.
    assert list(table.index) == ["lo", "hi"]


def test_rarity_correlation_is_negative_when_rare_groups_score_high():
    # Three groups: sizes 100 / 10 / 2, mean leverage rising as size falls.
    groups = ["big"] * 100 + ["mid"] * 10 + ["rare"] * 2
    values = [1.0] * 100 + [2.0] * 10 + [5.0] * 2

    assert rarity_correlation(values, groups) == pytest.approx(-1.0)


def test_rarity_correlation_is_positive_when_the_relationship_inverts():
    groups = ["big"] * 100 + ["mid"] * 10 + ["rare"] * 2
    values = [5.0] * 100 + [2.0] * 10 + [1.0] * 2

    assert rarity_correlation(values, groups) == pytest.approx(1.0)


def test_rarity_correlation_needs_three_groups_to_rank():
    groups = ["a"] * 5 + ["b"] * 5
    assert np.isnan(rarity_correlation([1.0] * 5 + [2.0] * 5, groups))


# ----------------------------------------------------------------------
# sketch composition
# ----------------------------------------------------------------------


def test_composition_returns_fractions():
    frac = composition(["a", "a", "a", "b"])

    assert frac["a"] == pytest.approx(0.75)
    assert frac["b"] == pytest.approx(0.25)


def test_sketch_fold_change_flags_over_and_under_representation():
    full = ["common"] * 90 + ["rare"] * 10          # 90 % / 10 %
    sketch = ["common"] * 60 + ["rare"] * 40        # 60 % / 40 %

    table = sketch_fold_change(full, sketch)

    assert table.loc["rare", "fold"] == pytest.approx(4.0)      # 0.40 / 0.10
    assert table.loc["common", "fold"] == pytest.approx(0.6 / 0.9)
    assert table.loc["rare", "n_full"] == 10


def test_sketch_fold_change_reports_zero_for_a_missed_group():
    full = ["common"] * 95 + ["rare"] * 5
    sketch = ["common"] * 50                        # rare cells all missed

    table = sketch_fold_change(full, sketch)

    assert table.loc["rare", "frac_sketch"] == 0.0
    assert table.loc["rare", "fold"] == 0.0


# ----------------------------------------------------------------------
# agreement between two score vectors
# ----------------------------------------------------------------------


def test_agreement_is_perfect_for_identical_vectors():
    a = np.array([3.0, 1.0, 2.0, 5.0, 4.0])

    got = agreement(a, a, top_k=(2,))

    assert got["spearman"] == pytest.approx(1.0)
    assert got["pearson"] == pytest.approx(1.0)
    assert got["max_abs_diff"] == 0.0
    assert got["top2_overlap"] == pytest.approx(1.0)


def test_agreement_top_k_catches_a_reordered_head():
    # Same values, but the two highest are swapped with two low ones: rank
    # correlation stays high while the *set of cells you would draw* changes,
    # which is the thing that actually matters for a sampling weight.
    a = np.array([1.0, 2.0, 3.0, 9.0, 10.0])
    b = np.array([9.0, 10.0, 3.0, 1.0, 2.0])

    got = agreement(a, b, top_k=(2,))

    assert got["top2_overlap"] == 0.0


def test_agreement_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="shape"):
        agreement(np.zeros(4), np.zeros(5))


# ----------------------------------------------------------------------
# label accuracy and the R-side readers
# ----------------------------------------------------------------------


def test_label_accuracy_counts_exact_matches():
    assert label_accuracy(["a", "b", "c"], ["a", "b", "x"]) == pytest.approx(2 / 3)


def test_label_accuracy_skips_unannotated_cells():
    # Cells with no ground-truth annotation must not count as errors.
    assert label_accuracy(["a", "b"], ["a", "nan"]) == pytest.approx(1.0)


def test_label_accuracy_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        label_accuracy(["a"], ["a", "b"])


def test_r_feature_key_maps_underscores_to_dashes():
    assert _r_feature_key("Y_RNA") == "Y-RNA"
    assert _r_feature_key("CD14") == "CD14"


def test_read_r_scores_aligns_by_name_not_position(tmp_path):
    # The R file deliberately lists the cells in a different order.
    path = tmp_path / "r_leverage.csv"
    pd.DataFrame({"cell": ["c2", "c0", "c1"],
                  "r_leverage_exact": [20.0, 0.0, 10.0]}).to_csv(path, index=False)

    got = _read_r_scores(path, "r_leverage_exact", ["c0", "c1", "c2"])

    assert list(got) == [0.0, 10.0, 20.0]


def test_read_r_scores_raises_when_a_cell_is_missing(tmp_path):
    path = tmp_path / "r_leverage.csv"
    pd.DataFrame({"cell": ["c0"], "r_leverage_exact": [1.0]}).to_csv(path, index=False)

    with pytest.raises(KeyError, match="absent"):
        _read_r_scores(path, "r_leverage_exact", ["c0", "c1"])


# ----------------------------------------------------------------------
# the pipeline, on synthetic data
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_object():
    """A tall, clustered object: many cells, few features, one rare state."""
    rng = np.random.default_rng(0)
    G, sizes = 150, (500, 400, 60)
    baseline = rng.gamma(0.6, size=G) + 0.05
    block = G // (len(sizes) + 1)
    cols, ct, cells = [], [], []
    for k, nk in enumerate(sizes):
        prof = baseline.copy()
        prof[k * block:(k + 1) * block] *= 5.0
        prof = prof / prof.sum() * 3000.0
        for _ in range(nk):
            cols.append(rng.poisson(prof))
            ct.append(f"T{k}")
            cells.append(f"c{len(cells)}")
    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.array(cols).T), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
        meta_data=pd.DataFrame({"seurat_annotations": ct}, index=cells),
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=100)
    return obj, np.array(ct)


def test_run_leverage_covers_both_regimes(small_object):
    from tutorials.ifnb_sketch_tutorial import run_leverage

    obj, _ = small_object
    hvg = list(obj.assays[obj.active_assay].variable_features)

    scores = run_leverage(obj, hvg)

    assert set(scores) >= {"exact", "sketched", "exact_seconds", "sketched_seconds"}
    for regime in ("exact", "sketched"):
        s = scores[regime]
        assert s.shape == (len(obj),)
        assert np.all(np.isfinite(s)) and s.min() >= 0
    # Only the exact regime is normalised to the truncation.
    assert scores["exact"].sum() == pytest.approx(50.0, abs=1e-6)


def test_run_sketch_returns_both_methods_at_the_requested_size(small_object):
    from tutorials.ifnb_sketch_tutorial import run_sketch

    obj, _ = small_object
    hvg = list(obj.assays[obj.active_assay].variable_features)

    sketches = run_sketch(obj, hvg, ncells=200, seed=0)

    assert set(sketches) == {"LeverageScore", "Uniform"}
    for sk in sketches.values():
        assert len(sk) == 200
        assert set(sk.cell_names()) <= set(obj.cell_names())


def test_summarize_reports_every_number_the_vignette_quotes(small_object):
    from tutorials.ifnb_sketch_tutorial import run_leverage, run_sketch, summarize

    obj, ct = small_object
    hvg = list(obj.assays[obj.active_assay].variable_features)
    scores = run_leverage(obj, hvg)
    sketches = run_sketch(obj, hvg, ncells=200, seed=0)

    summary = summarize(obj, scores, sketches, hvg)

    assert summary["n_cells"] == len(obj)
    assert summary["n_hvg"] == len(hvg)
    assert set(summary["leverage"]) == {"exact", "sketched"}
    assert set(summary["sketch_composition"]) == {"LeverageScore", "Uniform"}
    assert set(summary["rarity_spearman"]) == {"exact", "sketched"}
    # Every annotated type appears in the enrichment table.
    assert set(summary["enrichment"].index) == set(np.unique(ct))


def test_report_concordance_is_silent_without_the_r_reference(tmp_path, small_object):
    from tutorials.ifnb_sketch_tutorial import run_leverage, run_sketch, summarize
    from tutorials.ifnb_sketch_tutorial import report_concordance

    obj, _ = small_object
    hvg = list(obj.assays[obj.active_assay].variable_features)
    summary = summarize(obj, run_leverage(obj, hvg),
                        run_sketch(obj, hvg, ncells=200, seed=0), hvg)

    # No r_leverage.csv in an empty directory → returns None rather than raising,
    # so the Python half of the tutorial runs standalone.
    assert report_concordance(summary, figures=tmp_path) is None
