"""Network-free tests for the reference-mapping tutorial (tutorials/panc8_reference_mapping_tutorial.py).

Covers the pure metric helpers (transfer accuracy, label concordance, per-class
recall, macro recall, scoreboard) directly, and drives the whole pipeline on a
small synthetic two-technology dataset — an annotated reference and a query that
carries a technology/batch gene block the reference never sees — never touching
the network or the real panc8 export. The synthetic techs are named the same as
the real reference/query (``celseq2`` / ``smartseq2``) so ``run_full`` runs on
its defaults.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.panc8_reference_mapping_tutorial as tut  # noqa: E402
from tutorials.panc8_reference_mapping_tutorial import (  # noqa: E402
    transfer_accuracy,
    label_concordance,
    per_class_recall,
    macro_recall,
    build_scoreboard,
)


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_transfer_accuracy_perfect_and_partial():
    truth = np.array(["a", "b", "c", "a"])
    assert transfer_accuracy(truth, truth) == pytest.approx(1.0)
    assert transfer_accuracy(np.array(["a", "b", "x", "x"]), truth) == pytest.approx(0.5)


def test_transfer_accuracy_rejects_bad_shapes():
    with pytest.raises(ValueError):
        transfer_accuracy(np.array(["a"]), np.array(["a", "b"]))
    with pytest.raises(ValueError):
        transfer_accuracy(np.array([]), np.array([]))


def test_label_concordance_agreement():
    a = np.array(["a", "b", "c", "d"])
    assert label_concordance(a, a) == pytest.approx(1.0)
    # Two predictors that agree on 3 of 4 cells.
    assert label_concordance(a, np.array(["a", "b", "c", "z"])) == pytest.approx(0.75)
    with pytest.raises(ValueError):
        label_concordance(np.array(["a", "b"]), np.array(["a"]))


def test_per_class_recall_orders_by_support_and_scores():
    # type A: 3 cells, 2 correct -> 2/3; type B: 1 cell, correct -> 1.0.
    truth = np.array(["A", "A", "A", "B"])
    pred = np.array(["A", "A", "X", "B"])
    df = per_class_recall(pred, truth)
    # Support-descending: A (3) before B (1).
    assert list(df["celltype"]) == ["A", "B"]
    assert list(df["support"]) == [3, 1]
    assert df.loc[df["celltype"] == "A", "recall"].iloc[0] == pytest.approx(2 / 3)
    assert df.loc[df["celltype"] == "B", "recall"].iloc[0] == pytest.approx(1.0)


def test_macro_recall_is_unweighted_mean_of_classes():
    # A recall 2/3, B recall 1.0 -> macro = (2/3 + 1)/2, unaffected by A's larger support.
    truth = np.array(["A", "A", "A", "B"])
    pred = np.array(["A", "A", "X", "B"])
    assert macro_recall(pred, truth) == pytest.approx((2 / 3 + 1.0) / 2)


def test_build_scoreboard_column_order():
    rows = [{"tool": "shanuz (pcaproject)", "n_query": 100, "n_anchors": 300,
             "accuracy": 0.98, "macro_recall": 0.8, "mean_score": 0.97}]
    board = build_scoreboard(rows)
    assert list(board.columns) == [
        "tool", "n_query", "n_anchors", "accuracy", "macro_recall", "mean_score"]
    assert list(board["tool"]) == ["shanuz (pcaproject)"]


# ---------------------------------------------------------------------------
# Synthetic end-to-end pipeline (no network)
# ---------------------------------------------------------------------------

def _synthetic_panc8(seed=0):
    """A tiny panc8 stand-in: 3 cell types x 2 technologies + a query batch block.

    Returns the loader's 4-tuple (counts, genes, cells, meta) so it can be
    substituted for shanuz.datasets.panc8. The reference tech is ``celseq2`` and
    the query tech ``smartseq2`` (the tutorial's defaults); the query carries a
    technology-specific gene block the reference never sees, so a naive nearest
    neighbour would fail and the projection has to earn the transfer.
    """
    rng = np.random.default_rng(seed)
    G = 300
    types = ["A", "B", "C"]
    plan = {"celseq2": [40, 40, 40], "smartseq2": [45, 45, 45]}
    cols, tech, celltype, cells = [], [], [], []
    c = 0
    for t in ("celseq2", "smartseq2"):
        for ti, ct in enumerate(types):
            for _ in range(plan[t][ti]):
                base = rng.gamma(0.3, size=G) + 0.05
                base[ti * 50:(ti + 1) * 50] += 5.0        # cell-type marker block
                if t == "smartseq2":
                    base[200:250] += 4.0                   # technology/batch block
                cols.append(rng.poisson(base * 3000.0 / base.sum()))
                tech.append(t)
                celltype.append(ct)
                cells.append(f"{t}_c{c}")
                c += 1
    counts = sp.csc_matrix(np.asarray(cols).T)             # G x n_cells
    genes = [f"g{i}" for i in range(G)]
    meta = pd.DataFrame({"tech": tech, "celltype": celltype}, index=cells)
    return counts, genes, cells, meta


@pytest.fixture(scope="module")
def mapped():
    """Run the full pipeline once on the synthetic dataset, figures in a tmp dir."""
    data = _synthetic_panc8()
    mp = pytest.MonkeyPatch()
    mp.setattr(tut, "panc8", lambda data_dir=None: data)
    tmp = Path(tempfile.mkdtemp())
    mp.setattr(tut, "FIGURES", tmp)
    reference, query, anchors, predictions, summary = tut.run_full(
        verbose=False, do_umap=False)
    yield tut, reference, query, anchors, predictions, summary, tmp
    mp.undo()


def test_load_split_by_technology():
    data = _synthetic_panc8()
    mp = pytest.MonkeyPatch()
    mp.setattr(tut, "panc8", lambda data_dir=None: data)
    reference, query = tut.load_panc8_split()
    mp.undo()
    assert set(np.asarray(reference.meta_data[tut.TECH])) == {"celseq2"}
    assert set(np.asarray(query.meta_data[tut.TECH])) == {"smartseq2"}
    assert tut.CELLTYPE in query.meta_data.columns
    # Reference and query share an identical gene universe (split from one object).
    assert reference.assays["RNA"].features() == query.assays["RNA"].features()


def test_prep_writes_shared_hvgs(mapped):
    _tut, _ref, _query, _anchors, _pred, _summary, tmp = mapped
    assert (tmp / "hvg_features.txt").exists()
    hvg = (tmp / "hvg_features.txt").read_text().split()
    assert len(hvg) > 0


def test_transfer_predicts_query_celltypes(mapped):
    _tut, _ref, query, anchors, predictions, _summary, _tmp = mapped
    assert list(predictions.index) == query.cell_names()
    assert "predicted.id" in predictions.columns
    assert "prediction.score.max" in predictions.columns
    assert len(anchors.anchors) > 0
    truth = np.asarray(query.meta_data[tut.CELLTYPE]).astype(str)
    pred = predictions["predicted.id"].to_numpy().astype(str)
    # The projection should annotate the query despite its technology block.
    assert transfer_accuracy(pred, truth) > 0.8


def test_summary_scoreboard_and_perclass(mapped):
    _tut, _ref, query, _anchors, _pred, summary, _tmp = mapped
    board = summary["scoreboard"]
    assert list(board["tool"]) == ["shanuz (pcaproject)"]
    assert board["n_anchors"].iloc[0] > 0
    assert summary["n_celltypes"] == 3
    # Every true cell type appears in the per-class recall table.
    assert set(summary["per_class"]["celltype"]) == {"A", "B", "C"}


def test_report_concordance_reads_fabricated_r_calls(mapped):
    _tut, _ref, query, _anchors, predictions, _summary, tmp = mapped
    cells = query.cell_names()
    py_pred = predictions["predicted.id"].reindex(cells).to_numpy()
    # Fabricate an R calls file whose predicted.id equals Python's -> concordance 1.
    df = pd.DataFrame({
        "cell": cells,
        "tech": "smartseq2",
        "celltype": np.asarray(query.meta_data[tut.CELLTYPE]),
        "R_predicted": py_pred,
        "R_score_max": 1.0,
    })
    df.to_csv(tmp / "r_calls.csv", index=False)
    out = tut.report_concordance(query, predictions, verbose=False)
    assert out is not None
    assert out["concordance"] == pytest.approx(1.0)
    # Identical predictions -> identical accuracy against the truth.
    assert out["r_accuracy"] == pytest.approx(out["py_accuracy"])


def test_report_concordance_absent_returns_none(mapped):
    _tut, _ref, query, _anchors, predictions, _summary, tmp = mapped
    missing = tmp / "does_not_exist.csv"
    assert tut.report_concordance(
        query, predictions, r_calls_path=missing, verbose=False) is None
