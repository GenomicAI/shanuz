"""Tests for reference mapping (v0.3.0).

  * find_transfer_anchors  (transfer.py)
  * transfer_data          (transfer.py)

Builds an annotated reference and an unlabelled query (same two cell types, the
query carrying a batch-effect gene block the reference never sees), then checks
that projecting the query into the reference and transferring the ``celltype``
label recovers the query's true cell types. Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data,
    find_variable_features,
    scale_data,
)
from shanuz.transfer import (  # noqa: E402
    find_transfer_anchors,
    transfer_data,
    TransferAnchors,
)


def _labelled_object(batch, seed=0, n_per=60, G=200):
    """Two cell types (A/B); batch '2' carries a batch-effect gene block."""
    rng = np.random.default_rng(seed)
    mat = np.zeros((G, 2 * n_per))
    celltype, cells = [], []
    c = 0
    for t in ("A", "B"):
        for _ in range(n_per):
            base = rng.gamma(0.3, size=G) + 0.05
            if t == "A":
                base[0:50] += 5.0
            else:
                base[50:100] += 5.0
            if batch == "2":
                base[100:150] += 4.0          # batch-specific gene block
            mat[:, c] = rng.poisson(base * 3000.0 / base.sum())
            celltype.append(t)
            cells.append(f"b{batch}_c{c}")
            c += 1

    meta = pd.DataFrame({"celltype": celltype, "batch": batch}, index=cells)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(mat), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
        meta_data=meta,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=150)
    scale_data(obj)
    return obj, np.array(celltype)


def _ref_query(seed=0):
    reference, ct_ref = _labelled_object("1", seed=seed)
    query, ct_query = _labelled_object("2", seed=seed + 1)
    return reference, query, ct_ref, ct_query


# ----------------------------------------------------------------------
# find_transfer_anchors
# ----------------------------------------------------------------------


def test_find_transfer_anchors_returns_scored_anchors():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query, k_filter=0)

    assert isinstance(anchors, TransferAnchors)
    df = anchors.anchors
    assert list(df.columns) == ["cell1", "cell2", "score"]
    assert len(df) > 0
    assert df["cell1"].between(0, len(reference) - 1).all()
    assert df["cell2"].between(0, len(query) - 1).all()
    assert df["score"].between(0.0, 1.0).all()
    assert anchors.query_embedding.shape[0] == len(query)


def test_find_transfer_anchors_rejects_unknown_reduction():
    reference, query, _, _ = _ref_query()
    with pytest.raises(ValueError):
        find_transfer_anchors(reference, query, reduction="mnn")


def test_find_transfer_anchors_is_deterministic():
    reference, query, _, _ = _ref_query()
    a1 = find_transfer_anchors(reference, query)
    a2 = find_transfer_anchors(reference, query)
    pd.testing.assert_frame_equal(a1.anchors, a2.anchors)


def test_find_transfer_anchors_cca_reduction():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query, reduction="cca")
    assert anchors.reduction == "cca"
    assert len(anchors.anchors) > 0


# ----------------------------------------------------------------------
# transfer_data — classification
# ----------------------------------------------------------------------


def _accuracy(predicted, truth):
    return float(np.mean(np.asarray(predicted) == np.asarray(truth)))


def test_transfer_data_predicts_query_celltypes():
    reference, query, _, ct_query = _ref_query()
    anchors = find_transfer_anchors(reference, query, reduction="pcaproject")

    pred = transfer_data(anchors, refdata="celltype")

    assert list(pred.index) == query.cell_names()
    assert "predicted.id" in pred.columns
    assert "prediction.score.max" in pred.columns
    # Projection should annotate the query despite its batch block.
    assert _accuracy(pred["predicted.id"], ct_query) > 0.85


def test_transfer_data_scores_form_a_distribution():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)
    pred = transfer_data(anchors, refdata="celltype")

    score_cols = [c for c in pred.columns if c.startswith("prediction.score.")
                  and c != "prediction.score.max"]
    assert set(score_cols) == {"prediction.score.A", "prediction.score.B"}
    # Per-class scores are a probability distribution over the classes.
    row_sums = pred[score_cols].sum(axis=1).to_numpy()
    assert np.allclose(row_sums, 1.0)
    assert (pred["prediction.score.max"] >= 0.5).all()
    assert pred[score_cols].to_numpy().min() >= 0.0


def test_transfer_data_accepts_label_array():
    reference, query, ct_ref, ct_query = _ref_query()
    anchors = find_transfer_anchors(reference, query)

    # Same labels, passed as a raw per-reference-cell array instead of a column.
    pred = transfer_data(anchors, refdata=ct_ref)
    assert _accuracy(pred["predicted.id"], ct_query) > 0.85


def test_transfer_data_unknown_column_raises():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)
    with pytest.raises(KeyError):
        transfer_data(anchors, refdata="no_such_column")


def test_transfer_data_wrong_length_array_raises():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)
    with pytest.raises(ValueError):
        transfer_data(anchors, refdata=np.array(["A", "B", "A"]))


# ----------------------------------------------------------------------
# transfer_data — continuous imputation
# ----------------------------------------------------------------------


def test_transfer_data_imputes_expression():
    reference, query, _, ct_query = _ref_query()
    anchors = find_transfer_anchors(reference, query)

    # Impute the two cell-type marker blocks from the reference onto the query.
    feats = ["g10", "g60"]     # g10 ↑ in type A, g60 ↑ in type B
    ref_expr = reference.get_assay().layer_data("data", features=feats)
    ref_expr = ref_expr.toarray() if sp.issparse(ref_expr) else np.asarray(ref_expr)

    imputed = transfer_data(anchors, refdata=ref_expr, refdata_features=feats)

    assert list(imputed.index) == feats
    assert list(imputed.columns) == query.cell_names()
    assert imputed.shape == (2, len(query))

    # g10 should read higher in the query cells the reference calls type A.
    a_mask = ct_query == "A"
    assert imputed.loc["g10"].to_numpy()[a_mask].mean() > \
        imputed.loc["g10"].to_numpy()[~a_mask].mean()
    assert imputed.loc["g60"].to_numpy()[~a_mask].mean() > \
        imputed.loc["g60"].to_numpy()[a_mask].mean()


def test_transfer_data_no_anchors_raises():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)
    anchors.anchors = anchors.anchors.iloc[0:0]  # drop every anchor
    with pytest.raises(ValueError):
        transfer_data(anchors, refdata="celltype")
