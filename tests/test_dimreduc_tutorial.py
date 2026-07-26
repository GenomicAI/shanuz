"""Network-free tests for the dim-reduction tutorial (tutorials/pbmc3k_dimreduc_tutorial.py).

Covers the pure metric helpers directly — the PC-cutoff readings, the
sign/order-free ICA matcher, the kNN structure comparison and the per-PC
JackStraw table — and drives the whole pipeline on a small synthetic dataset
with planted structure, never touching the network or the real PBMC download.

The helpers here carry more weight than usual: JackStraw's per-PC scores are not
comparable across the two tools by construction (R aggregates with ``prop.test``,
shanuz with a KS test), so the *derived* readings — which PCs clear alpha, how
many to keep — are what the comparison actually rests on. They need to be right.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.pbmc3k_dimreduc_tutorial as tut  # noqa: E402
from tutorials.pbmc3k_dimreduc_tutorial import (  # noqa: E402
    significant_dims,
    n_leading_significant,
    matched_component_correlation,
    knn_overlap,
    basis_agreement,
    pc_table,
    build_scoreboard,
    _r_feature_key,
)


# ---------------------------------------------------------------------------
# PC cutoff readings
# ---------------------------------------------------------------------------

def test_significant_dims_is_one_based():
    p = np.array([1e-8, 1e-4, 0.5, 0.01])
    assert significant_dims(p, alpha=0.05).tolist() == [1, 2, 4]


def test_significant_dims_can_be_empty_or_all():
    assert significant_dims(np.array([0.9, 0.5])).tolist() == []
    assert significant_dims(np.array([0.0, 0.0])).tolist() == [1, 2]


def test_significant_dims_rejects_bad_shapes():
    with pytest.raises(ValueError):
        significant_dims(np.array([]))
    with pytest.raises(ValueError):
        significant_dims(np.array([[0.1, 0.2]]))


def test_n_leading_significant_stops_at_the_first_failure():
    """A significant PC *after* a gap must not extend the run.

    This is the whole convention: Seurat's guidance is to cut at the drop-off,
    so PC 5 clearing alpha after PC 4 failed does not make the answer 5. Getting
    this wrong would silently inflate the recommended dimensionality.
    """
    p = np.array([1e-9, 1e-9, 1e-9, 0.4, 1e-9, 0.9])
    assert n_leading_significant(p, alpha=0.05) == 3
    assert significant_dims(p, alpha=0.05).tolist() == [1, 2, 3, 5]   # differs, by design


def test_n_leading_significant_edges():
    assert n_leading_significant(np.array([0.9, 0.9])) == 0           # none survive
    assert n_leading_significant(np.array([0.0, 0.0, 0.0])) == 3      # all survive
    with pytest.raises(ValueError):
        n_leading_significant(np.array([]))


# ---------------------------------------------------------------------------
# ICA component matching
# ---------------------------------------------------------------------------

def test_matched_components_ignore_order_and_sign():
    """The point of the matcher: ICA names components arbitrarily.

    A permuted, negated, rescaled copy of the same components is the *same*
    answer, and must score 1.0 — a naive column-wise correlation would score it
    near zero and report a perfectly faithful port as broken.
    """
    rng = np.random.default_rng(0)
    a = rng.normal(size=(300, 4))
    b = a[:, [2, 0, 3, 1]] * np.array([-1.0, 1.0, -2.0, 0.5])
    out = matched_component_correlation(a, b)
    assert out["mean_abs_r"] == pytest.approx(1.0, abs=1e-9)
    assert out["min_abs_r"] == pytest.approx(1.0, abs=1e-9)
    assert sorted(c for _, c in out["pairs"]) == [1, 2, 3, 4]         # a bijection


def test_matched_components_score_low_for_unrelated():
    rng = np.random.default_rng(1)
    out = matched_component_correlation(rng.normal(size=(400, 5)),
                                        rng.normal(size=(400, 5)))
    assert out["mean_abs_r"] < 0.3


def test_matched_components_handle_constant_column():
    """A zero-variance component must not produce NaN or a divide-by-zero."""
    rng = np.random.default_rng(2)
    a = rng.normal(size=(100, 3))
    b = a.copy()
    b[:, 1] = 5.0
    out = matched_component_correlation(a, b)
    assert np.isfinite(out["mean_abs_r"])


def test_matched_components_truncate_to_the_shorter_set():
    rng = np.random.default_rng(3)
    a = rng.normal(size=(120, 6))
    out = matched_component_correlation(a, a[:, :3])
    assert len(out["pairs"]) == 3
    assert out["mean_abs_r"] == pytest.approx(1.0, abs=1e-9)


def test_matched_components_reject_bad_shapes():
    with pytest.raises(ValueError):
        matched_component_correlation(np.zeros((10, 2)), np.zeros((9, 2)))
    with pytest.raises(ValueError):
        matched_component_correlation(np.zeros(10), np.zeros(10))


# ---------------------------------------------------------------------------
# kNN structure comparison
# ---------------------------------------------------------------------------

def test_knn_overlap_is_one_for_identical_neighbourhoods():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(120, 5))
    assert knn_overlap(x, x, k=10) == pytest.approx(1.0)
    # Rotation + rescaling preserves every neighbourhood, so it must still be 1.
    assert knn_overlap(x, x * 3.0, k=10) == pytest.approx(1.0)


def test_knn_overlap_is_low_for_unrelated_embeddings():
    rng = np.random.default_rng(5)
    a = rng.normal(size=(300, 4))
    b = rng.normal(size=(300, 4))
    assert knn_overlap(a, b, k=10) < 0.2


def test_knn_overlap_clamps_k_to_the_sample_size():
    rng = np.random.default_rng(6)
    x = rng.normal(size=(6, 3))
    assert knn_overlap(x, x, k=50) == pytest.approx(1.0)   # k clamped to n-1


def test_knn_overlap_rejects_bad_shapes():
    with pytest.raises(ValueError):
        knn_overlap(np.zeros((10, 2)), np.zeros((9, 2)))
    with pytest.raises(ValueError):
        knn_overlap(np.zeros(10), np.zeros(10))
    with pytest.raises(ValueError):
        knn_overlap(np.zeros((1, 2)), np.zeros((1, 2)))


# ---------------------------------------------------------------------------
# The shared-basis check
# ---------------------------------------------------------------------------

def test_basis_agreement_is_perfect_for_the_same_embedding():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(200, 6))
    out = basis_agreement(x, x * -1.0)                 # sign flips are free
    assert out["pca_basis_aligned_through"] == 6
    assert out["pca_basis_min_aligned_r"] == pytest.approx(1.0, abs=1e-9)


def test_basis_agreement_localises_a_reordered_tail():
    """A permuted noise tail must not read as a disagreeing basis.

    The real pbmc3k run swaps PC 16-19 between the tools while PC 1-15 match
    one-to-one. Reporting that as "min |r| = 0.15 over 20 PCs" would wrongly
    suggest the two PCAs disagree, and would undercut a JackStraw finding that
    is actually decisive over the aligned range.
    """
    rng = np.random.default_rng(8)
    x = rng.normal(size=(300, 6))
    swapped = x[:, [0, 1, 2, 3, 5, 4]]                 # tail two components swap
    out = basis_agreement(x, swapped)
    assert out["pca_basis_aligned_through"] == 4       # PC 1-4 still in order
    assert out["pca_basis_min_aligned_r"] == pytest.approx(1.0, abs=1e-9)
    assert out["pca_basis_best_match"][4] == 6         # py PC5 matches R PC6
    assert out["pca_basis_best_match"][5] == 5


def test_basis_agreement_rejects_mismatched_cells():
    with pytest.raises(ValueError):
        basis_agreement(np.zeros((10, 3)), np.zeros((9, 3)))


# ---------------------------------------------------------------------------
# The per-PC JackStraw table
# ---------------------------------------------------------------------------

def test_pc_table_counts_significant_features():
    scores = np.array([0.0, 0.5])
    emp = np.array([                 # 4 features x 2 PCs
        [0.0, 0.4],
        [0.0, 0.6],
        [1e-6, 0.2],
        [0.5, 0.9],
    ])
    table = pc_table(scores, emp, score_thresh=1e-5)
    assert table["PC"].tolist() == [1, 2]
    assert table["py_n_sig_features"].tolist() == [3, 0]
    assert table["py_median_p"].iloc[0] == pytest.approx(5e-7)
    assert "r_score" not in table.columns          # absent until R has run


def test_pc_table_adds_the_r_columns_when_given():
    scores = np.array([0.0, 0.5])
    emp = np.array([[0.0, 0.4], [1e-6, 0.6]])
    table = pc_table(scores, emp, r_pvals=np.array([1e-3, 0.8]),
                     r_emp=np.array([[0.0, 0.5], [0.3, 0.7]]), score_thresh=1e-5)
    assert table["r_score"].tolist() == [1e-3, 0.8]
    assert table["r_n_sig_features"].tolist() == [1, 0]


def test_pc_table_rejects_a_short_pvalue_matrix():
    with pytest.raises(ValueError):
        pc_table(np.zeros(5), np.zeros((10, 2)))


def test_build_scoreboard_orders_and_drops_absent_columns():
    df = build_scoreboard([{"metric": "m", "method": "ICA", "agreement": 0.9}])
    assert df.columns.tolist() == ["method", "metric", "agreement"]


def test_r_feature_key_matches_read10x_renaming():
    assert _r_feature_key("Y_RNA") == "Y-RNA"
    assert _r_feature_key("CD8A") == "CD8A"


# ---------------------------------------------------------------------------
# Pipeline on synthetic data (no network)
# ---------------------------------------------------------------------------

def _synthetic_object(n_cells=300, n_genes=400, n_programs=3, seed=0):
    """Counts with a few planted co-expression programs.

    Real structure in the leading PCs and noise after it, so JackStraw has
    something to find *and* something to reject — the shape the tutorial's real
    run has, at a size that runs in a second.
    """
    from shanuz.shanuz import create_shanuz_object

    rng = np.random.default_rng(seed)
    counts = rng.poisson(0.3, size=(n_genes, n_cells)).astype(float)
    block = n_genes // (n_programs * 4)
    for p in range(n_programs):
        cells = slice(p * (n_cells // n_programs), (p + 1) * (n_cells // n_programs))
        genes = slice(p * block, (p + 1) * block)
        counts[genes, cells] += rng.poisson(8.0, size=(block, n_cells // n_programs))
    genes = [f"GENE{i}" for i in range(n_genes)]
    cells = [f"CELL{i}" for i in range(n_cells)]
    return create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="RNA", min_cells=0, min_features=0,
        project="synthetic", feature_names=genes, cell_names=cells,
    )


@pytest.fixture(scope="module")
def scored(tmp_path_factory):
    """Run the tutorial's own pipeline once on synthetic data.

    ``out_dir`` is not optional here. Without it ``prep`` writes its handoff to
    the real ``tutorials/figures_dimreduc/``, so running the test suite left
    ``hvg_features.txt`` holding 100 synthetic ``GENE*`` names and ``cells.txt``
    holding ``CELL0..CELL119`` — and ``pbmc3k_dimreduc_verify.R`` reads both.
    """
    obj = _synthetic_object()
    tut.prep(obj, n_hvg=200, n_pcs=20, out_dir=tmp_path_factory.mktemp("handoff"))
    js, scores = tut.run_jackstraw(obj, dims=10, num_replicate=20, seed=0)
    tut.run_reductions(obj, n_ics=5, tsne_dims=5, seed=0)
    return obj, js, scores


def test_pipeline_populates_every_reduction(scored):
    obj, js, scores = scored
    assert {"pca", "ica", "tsne"} <= set(obj.reductions)
    assert obj.reductions["ica"].cell_embeddings.shape == (300, 5)
    assert obj.reductions["tsne"].cell_embeddings.shape == (300, 2)
    assert js.empirical_p_values.shape == (200, 10)
    assert scores.shape == (10,)


def test_jackstraw_p_values_are_probabilities(scored):
    _obj, js, _scores = scored
    emp = js.empirical_p_values
    assert np.isfinite(emp).all()
    assert emp.min() >= 0.0 and emp.max() <= 1.0


def test_leading_pcs_beat_the_trailing_ones(scored):
    """The planted programs must load more significantly than the noise PCs.

    Deliberately *not* an assertion about the absolute scores: shanuz's
    aggregation saturates at 0.0 on real data, which is the finding the tutorial
    reports rather than something to bake in here. What must hold either way is
    the ordering — the leading PCs carry more sub-threshold features than the
    trailing ones.
    """
    _obj, js, _scores = scored
    n_sig = (js.empirical_p_values <= 1e-5).sum(axis=0)
    assert n_sig[:3].mean() > n_sig[-3:].mean()


def test_summarize_reports_the_cutoff_and_structure(scored):
    obj, js, scores = scored
    summary = tut.summarize(obj, scores, js, verbose=False)
    assert summary["n_cells"] == 300
    assert summary["n_features"] == 200
    assert 0 <= summary["n_leading_significant"] <= 10
    assert summary["n_ics"] == 5
    # t-SNE of a 5-PC space keeps a good share of its neighbourhoods.
    assert 0.0 <= summary["tsne_knn_vs_pca"] <= 1.0
    assert isinstance(summary["pc_table"], pd.DataFrame)


def test_summarize_survives_without_the_optional_reductions(tmp_path):
    """JackStraw half of the tutorial must stand alone (run_full(do_reductions=False))."""
    obj = _synthetic_object(n_cells=120, n_genes=200, seed=1)
    tut.prep(obj, n_hvg=100, n_pcs=10, out_dir=tmp_path)
    js, scores = tut.run_jackstraw(obj, dims=5, num_replicate=10, seed=1)
    summary = tut.summarize(obj, scores, js, verbose=False)
    assert "n_ics" not in summary and "tsne_knn_vs_pca" not in summary
    assert summary["pc_scores"].shape == (5,)


def test_report_concordance_returns_none_without_the_r_run(scored, tmp_path, monkeypatch):
    obj, js, scores = scored
    monkeypatch.setattr(tut, "FIGURES", tmp_path)
    summary = tut.summarize(obj, scores, js, verbose=False)
    assert tut.report_concordance(obj, summary, verbose=False) is None


# ---------------------------------------------------------------------------
# The declared bands
# ---------------------------------------------------------------------------

def _write_fake_r_reference(obj, scores, js, dest):
    """Write the Python run back out in the R script's file format.

    Comparing a run against itself is the one case where every band's answer is
    known in advance, which makes it the right fixture for testing the *wiring*:
    each band has to be populated and hold. No network, no R.
    """
    dest.mkdir(parents=True, exist_ok=True)
    emb = obj.reductions["pca"].cell_embeddings
    pd.DataFrame(emb, index=obj.cell_names(),
                 columns=[f"PC_{i + 1}" for i in range(emb.shape[1])]).rename_axis(
        "cell").to_csv(dest / "r_pca.csv")
    pd.DataFrame({"PC": np.arange(1, len(scores) + 1), "R_Score": scores}).to_csv(
        dest / "r_jackstraw_pcs.csv", index=False)
    emp = js.empirical_p_values
    pd.DataFrame(emp, index=list(obj.reductions["pca"].features()),
                 columns=[f"PC{i + 1}" for i in range(emp.shape[1])]).rename_axis(
        "feature").to_csv(dest / "r_jackstraw_p.csv")


def test_every_band_is_populated_and_holds_against_an_identical_run(
        scored, tmp_path, monkeypatch):
    """A declared band with nothing feeding it must not pass silently.

    `check_bands` scores an absent name as NaN, which fails — but only if
    something hands it the dict `report_concordance` actually builds. This ties
    the two together: adding a name to `BANDS` without computing it, or renaming
    a key in `out`, goes red here rather than on someone else's `--report`.
    """
    obj, js, scores = scored
    monkeypatch.setattr(tut, "FIGURES", tmp_path)
    _write_fake_r_reference(obj, scores, js, tmp_path)
    summary = tut.summarize(obj, scores, js, verbose=False)

    out = tut.report_concordance(obj, summary, verbose=False)

    missing = [name for name in tut.BANDS if name not in out]
    assert not missing, f"declared but never measured: {missing}"
    assert out["bands_hold"], [
        (v.name, v.value) for v in out["band_verdicts"] if not v.ok]
    # Against itself the two "differences" are zero by construction, which is
    # what makes the previous assertion meaningful rather than lucky.
    assert out["jackstraw_keep_gap"] == 0
    assert out["significant_agreement"] == pytest.approx(1.0)


def test_the_report_refuses_an_r_reference_from_a_different_run(
        scored, tmp_path, monkeypatch):
    """The failure that started this: a reference older than its handoff.

    `_read` reindexes R's embedding onto the Python cell order, which fills an
    absent barcode with NaN rather than raising, so the whole comparison came
    out `nan` and printed without comment.
    """
    from tutorials.bands import StaleReferenceError

    obj, js, scores = scored
    monkeypatch.setattr(tut, "FIGURES", tmp_path)
    _write_fake_r_reference(obj, scores, js, tmp_path)
    stale = pd.read_csv(tmp_path / "r_pca.csv")
    stale.loc[0, "cell"] = "CELL_FROM_A_PREVIOUS_RUN"
    stale.to_csv(tmp_path / "r_pca.csv", index=False)
    summary = tut.summarize(obj, scores, js, verbose=False)

    with pytest.raises(StaleReferenceError, match="different cell set"):
        tut.report_concordance(obj, summary, verbose=False)


def test_prep_leaves_the_real_handoff_alone_when_given_an_out_dir(tmp_path):
    """Running the test suite must not rewrite what the R script reads.

    It did: `prep` wrote unconditionally to `tutorials/figures_dimreduc/`, so a
    `pytest` run replaced the 2,000 real HVGs with 100 synthetic `GENE*` names
    and the 2,700 barcodes with `CELL0..CELL119`.
    """
    obj = _synthetic_object(n_cells=60, n_genes=120, seed=3)
    handoff, sentinel = tmp_path / "handoff", tmp_path / "figures"
    sentinel.mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(tut, "FIGURES", sentinel)
        tut.prep(obj, n_hvg=50, n_pcs=5, out_dir=handoff)
    assert len((handoff / "hvg_features.txt").read_text().split()) == 50
    assert (handoff / "cells.txt").read_text().split()[0] == "CELL0"
    assert list(sentinel.iterdir()) == [], "prep wrote to FIGURES anyway"


def test_the_jackstraw_band_admits_the_measured_spread():
    """R is deterministic at 13; shanuz spans 12-15 over a 60-seed sweep.

    Both ends of that measurement have to sit inside the band, or the tutorial
    fails on a good run. A gap of 3 has to sit outside, or it fails on nothing.
    """
    band = tut.BANDS["jackstraw_keep_gap"]
    assert band.holds(abs(12 - 13)) and band.holds(abs(15 - 13))
    assert not band.holds(abs(16 - 13))


def test_the_basis_band_reaches_at_least_r_s_cutoff():
    """`aligned_through` is a precondition, and its floor is not arbitrary.

    The per-PC JackStraw comparison is only like-for-like while the two bases
    are matched in order, so it has to reach the PC where the cutoff decision is
    actually made. R's drop-off is at 14, i.e. it keeps 13.
    """
    band = tut.BANDS["pca_basis_aligned_through"]
    assert band.low >= 13
    assert not band.holds(12)


def test_bands_hold_is_false_when_a_band_is_violated(monkeypatch):
    """`run_full` has to *act* on a failure, not just print one.

    Without this the `--report` exit code is decorative: the tutorial would
    print OUT-of-band lines and still return 0.
    """
    from tutorials.bands import Band
    monkeypatch.setitem(tut.BANDS, "significant_agreement",
                        Band(1.0, 1.0, "deliberately impossible"))
    verdicts = tut.check_bands(tut.BANDS, {"significant_agreement": 0.9})
    assert not all(v.ok for v in verdicts)
