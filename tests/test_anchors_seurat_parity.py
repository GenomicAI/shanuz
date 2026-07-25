"""Guards for the v4 anchor path against Seurat 5.5.1's own algorithm.

Every number quoted here came out of a Seurat 5.5.1 run or out of Seurat's
compiled helpers called directly — see ``tutorials/anchors_verify.R``. The unit
guards below are deliberately *properties* rather than re-computations of the
implementation: a test that recomputes what the code does can only ever prove
the code equals itself (see the ``avg_log2FC`` case in the DE tutorial, which
was green for months while the formula was wrong).

Reference run: 2,400-cell ifnb subsample (CTRL 1,200 / STIM 1,200), 2,000
anchor features, ``reduction = "cca"``, ``dims = 1:30``, ``nn.method = "rann"``.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.anchors import (  # noqa: E402
    _anchor_feature_matrix,
    _anchor_weights,
    _cca,
    _check_features,
    _filter_anchors,
    _mutual_nn_from,
    _neighbor_sets,
    _pca_loadings,
    _score_anchors,
    _standardize_cols,
    _top_dim_features,
    find_integration_anchors,
    integrate_data,
)
from shanuz.integration import integrate_layers  # noqa: E402
from shanuz.preprocessing import normalize_data, scale_data  # noqa: E402
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

# ---------------------------------------------------------------------------
# Transcribed from Seurat 5.5.1 (see the module docstring for the run)
# ---------------------------------------------------------------------------
SEURAT_N_ANCHORS = 2768
SEURAT_N_ANCHOR_FEATURES = 2000
SEURAT_N_AFTER_CHECKFEATURES = 1917      # 83 constant in one object or the other
SEURAT_TOP_FEATURE_CAP = 200             # TopDimFeatures max.features
SEURAT_CORRECTION_MEAN_ABS = 0.103181
SEURAT_CORRECTION_FRAC_NONZERO = 0.6222

# Seurat:::FindWeightsC called directly with four anchors at one query cell,
# scaled distances (0.8, 0.6, 0.4, 0.2), all scores 1, sd = 1.
SEURAT_FINDWEIGHTSC_4 = np.array([0.39025, 0.29988, 0.20487, 0.105])


def _pair(n=260, n_genes=150, seed=0):
    """Two batches with shared cell types and a batch-specific gene block.

    ``n`` must stay above the default ``k_filter`` of 200, or the retain-all
    guard fires, the anchor filter never runs, and every test that depends on
    filtering silently exercises nothing.
    """
    rng = np.random.default_rng(seed)
    objs = []
    for b, shift in (("1", 0.0), ("2", 4.0)):
        mat = np.zeros((n_genes, n))
        cells = []
        for c in range(n):
            base = rng.gamma(0.3, size=n_genes) + 0.05
            base[0:40] += 5.0 if c < n // 2 else 0.0
            base[40:80] += 0.0 if c < n // 2 else 5.0
            base[80:110] += shift
            mat[:, c] = rng.poisson(base * 2000.0 / base.sum())
            cells.append(f"b{b}_c{c}")
        obj = create_shanuz_object(
            counts=sp.csc_matrix(mat), assay="RNA",
            feature_names=[f"g{i}" for i in range(n_genes)], cell_names=cells,
        )
        normalize_data(obj)
        scale_data(obj)
        objs.append(obj)
    return objs


# ---------------------------------------------------------------------------
# CCA input transform — Standardize, not per-cell L2
# ---------------------------------------------------------------------------


def test_standardize_is_a_per_cell_zscore_not_a_per_cell_l2():
    rng = np.random.default_rng(0)
    mat = rng.normal(size=(30, 8)) + 5.0
    out = _standardize_cols(mat)
    assert np.allclose(out.mean(axis=0), 0.0, atol=1e-12)
    assert np.allclose(out.std(axis=0, ddof=1), 1.0, atol=1e-12)
    # An L2-normalised column keeps its mean; a standardised one does not. This
    # is the whole difference, and it is what moved anchor recall 70% -> 100%.
    l2 = mat / np.linalg.norm(mat, axis=0, keepdims=True)
    assert not np.allclose(l2.mean(axis=0), 0.0, atol=1e-6)


def test_cca_ignores_a_constant_shift_in_one_cells_profile():
    """A property of centring that per-cell L2 normalisation does not have.

    Standardising subtracts each cell's own mean, so adding a constant to one
    cell's whole expression profile cannot move the CCA. Dividing by the L2 norm
    does not centre, so under the old transform the same shift moves that cell.
    """
    rng = np.random.default_rng(1)
    A = rng.normal(size=(40, 25))
    B = rng.normal(size=(40, 25))
    ea, eb, _ = _cca(A, B, 5)

    shifted = A.copy()
    shifted[:, 3] += 7.5
    ea2, eb2, _ = _cca(shifted, B, 5)
    assert np.abs(ea - ea2).max() < 1e-9
    assert np.abs(eb - eb2).max() < 1e-9


def test_check_features_drops_features_constant_in_either_object():
    """``RunCCA`` runs ``CheckFeatures`` on *both* objects' scale.data.

    A constant feature carries no signal, but ``Standardize`` works down each
    *cell*, so leaving it in still shifts that cell's mean and sd and therefore
    every standardised value. The drop is not hygiene, it changes the answer:
    on the ifnb pair it is 83 of 2,000 features and it is most of the gap
    between 70% and 100% anchor recall.
    """
    rng = np.random.default_rng(2)
    A = rng.normal(size=(40, 25))
    B = rng.normal(size=(40, 25))
    A[7, :] = 3.0                       # constant in A only
    B[19, :] = -1.0                     # constant in B only

    keep = _check_features([A, B])
    assert 7 not in keep and 19 not in keep
    assert len(keep) == 38

    kept_a, _, _ = _cca(A[keep], B[keep], 5)
    naive_a, _, _ = _cca(A, B, 5)
    assert np.abs(kept_a - naive_a).max() > 1e-3, (
        "if dropping them changed nothing, the drop would not be worth doing"
    )


# ---------------------------------------------------------------------------
# The filter runs on TopDimFeatures of the log-normalised data
# ---------------------------------------------------------------------------


def test_top_dim_features_respects_seurats_cap():
    rng = np.random.default_rng(3)
    loadings = rng.normal(size=(2000, 30))
    top = _top_dim_features(loadings)
    assert len(top) <= SEURAT_TOP_FEATURE_CAP
    # and it is a real subspace, not "everything"
    assert len(top) < loadings.shape[0] / 2
    assert len(set(top)) == len(top)


def test_top_dim_features_never_exceeds_the_cap_for_any_dim_count():
    rng = np.random.default_rng(4)
    for n_dims in (2, 10, 30, 50):
        loadings = rng.normal(size=(2000, n_dims))
        assert len(_top_dim_features(loadings)) <= max(n_dims * 2, 200)


def test_filter_keeps_every_anchor_when_either_dataset_is_smaller_than_k_filter():
    """Seurat's guard is ``min(len(cells1), len(cells2)) < k.filter``.

    The asymmetric case is the one that separates it from just clamping k to the
    query size: a 30-cell reference against a 500-cell query still trips
    Seurat's guard and keeps everything, while ``min(k_filter, n_query)`` leaves
    k at 200 and filters away most of the anchors.
    """
    rng = np.random.default_rng(5)
    small_ref = rng.normal(size=(20, 30))
    big_query = rng.normal(size=(20, 500))
    pairs = [(i, i * 3) for i in range(25)]
    assert _filter_anchors(pairs, small_ref, big_query, k_filter=200) == pairs

    # With both sides above k_filter the filter does engage.
    big_a, big_b = rng.normal(size=(20, 300)), rng.normal(size=(20, 300))
    kept = _filter_anchors([(i, i) for i in range(50)], big_a, big_b, k_filter=5)
    assert len(kept) < 50


# ---------------------------------------------------------------------------
# ScoreAnchors — four neighbour tables, and it runs AFTER the filter
# ---------------------------------------------------------------------------


def test_score_uses_within_and_across_dataset_neighbours():
    """Each member contributes k_score own-dataset plus k_score other-dataset cells.

    With a batch effect present, a single kNN over the pooled stack returns
    almost only same-batch neighbours, so the two members of an anchor share
    nearly nothing and every score collapses towards the low end.
    """
    objs = _pair()
    anchors = find_integration_anchors(objs, reduction="cca", dims=10)
    scores = anchors.anchors["score"].to_numpy()
    assert len(scores) > 0
    assert scores.max() == pytest.approx(1.0)
    # A pooled-kNN score would pile up at zero; the four-way one spreads out.
    assert scores.std() > 0.1
    assert 0.2 < scores.mean() < 0.8


def test_returned_scores_are_consistent_with_the_returned_anchor_set():
    """FindAnchors filters, then scores — and the score is rescaled against the
    percentiles of whatever set it is handed.

    So re-scoring the anchors that came back must reproduce the scores that came
    back. If scoring ran first, the percentiles were taken from anchors that
    were then discarded and this identity fails.
    """
    objs = _pair()
    dims, k_anchor, k_score, k_filter = 10, 5, 30, 200
    anchors = find_integration_anchors(
        objs, reduction="cca", dims=dims, k_anchor=k_anchor,
        k_score=k_score, k_filter=k_filter,
    )
    from shanuz.anchors import _anchor_feature_matrix

    feats = anchors.anchor_features
    ref = _anchor_feature_matrix(objs[0], feats, "scale.data")
    qry = _anchor_feature_matrix(objs[1], feats, "scale.data")
    keep = _check_features([ref, qry])
    ea, eb, _ = _cca(ref[keep], qry[keep], min(dims, ref.shape[1] - 1))
    nbrs = _neighbor_sets(ea, eb, max(k_anchor, k_score))

    # The filter has to have actually removed something, or both orders agree
    # trivially and this proves nothing.
    unfiltered = _mutual_nn_from(nbrs, k_anchor)
    assert len(anchors.anchors) < len(unfiltered), (
        "fixture too small for the anchor filter to engage"
    )

    pairs = list(zip(anchors.anchors["cell1"], anchors.anchors["cell2"]))
    rescored = _score_anchors([(int(i), int(j)) for i, j in pairs],
                              nbrs, ref.shape[1], k_score)
    assert np.allclose(rescored, anchors.anchors["score"].to_numpy(), atol=1e-12)


def test_the_filter_sees_the_data_layer_not_scale_data():
    """``FilterAnchors`` uses ``slot = "data"``.

    Log-normalised data is non-negative; scale.data is centred and therefore has
    negatives. Recording what the filter was handed distinguishes the two
    without asserting anything about the filtering itself.
    """
    import shanuz.anchors as A

    seen = []
    original = A._filter_anchors

    def spy(pairs, a_feat, b_feat, k_filter):
        seen.append((a_feat.min(), b_feat.min()))
        return original(pairs, a_feat, b_feat, k_filter)

    A._filter_anchors = spy
    try:
        find_integration_anchors(_pair(), reduction="cca", dims=10)
    finally:
        A._filter_anchors = original

    assert seen, "the CCA path must run the anchor filter"
    for a_min, b_min in seen:
        assert a_min >= 0.0 and b_min >= 0.0


def test_pca_loadings_are_exact_not_randomized():
    """The reciprocal-PCA anchors need every PC right, not just the leading ones.

    ``_pca_loadings`` must use an exact SVD. sklearn's default ``PCA`` picks a
    *randomized* solver for a matrix this shape, whose trailing components are
    visibly wrong — on ifnb only ~12 of 30 PCs matched Seurat. Because
    reciprocal PCA standardizes each projected dimension before the neighbour
    search, a wrong trailing axis becomes a different anchor. This pins the
    helper against the reference NumPy SVD directly, and separately against
    sklearn's *randomized* result to prove the distinction is real.
    """
    from sklearn.decomposition import PCA

    rng = np.random.default_rng(0)
    # A shape where sklearn defaults to randomized SVD (min dim > 500,
    # n_components well below it) and the trailing PCs actually differ.
    mat = rng.normal(size=(900, 600))          # features × cells
    got = _pca_loadings(mat, 30)

    _u, _s, vt = np.linalg.svd(mat.T, full_matrices=False)
    ref = vt[:30, :].T
    per_pc = np.array([abs(np.corrcoef(got[:, d], ref[:, d])[0, 1]) for d in range(30)])
    assert per_pc.min() > 0.999, "loadings must match the exact SVD on every PC"

    randomized = PCA(n_components=30, svd_solver="randomized",
                     random_state=0).fit(mat.T).components_.T
    rand_pc = np.array(
        [abs(np.corrcoef(randomized[:, d], ref[:, d])[0, 1]) for d in range(30)]
    )
    assert rand_pc.min() < 0.9, (
        "if randomized SVD agreed on every PC the exact solver would not matter"
    )


def test_rpca_and_cca_find_different_anchors():
    """The reciprocal-PCA path is wired to its own spaces, not aliased to CCA.

    A guard against the whole rpca branch silently degrading to something
    CCA-like: on a batched pair the two reductions must genuinely disagree on
    which pairs are anchors, while both still recovering real structure.
    """
    objs = _pair()
    cca = find_integration_anchors(objs, reduction="cca", dims=10)
    rpca = find_integration_anchors(objs, reduction="rpca", dims=10)
    cca_pairs = set(zip(cca.anchors["cell1"], cca.anchors["cell2"]))
    rpca_pairs = set(zip(rpca.anchors["cell1"], rpca.anchors["cell2"]))
    assert len(rpca_pairs) > 0
    overlap = len(cca_pairs & rpca_pairs) / max(len(rpca_pairs), 1)
    assert overlap < 0.95, "rpca must not be a rename of the cca path"


def test_mutual_nn_only_looks_at_the_first_k_anchor_neighbours():
    """FindNN searches at max(k_anchor, k_score) so one set of neighbour tables
    serves both steps; only the leading k_anchor columns make an anchor."""
    nbrs = {
        # A-cell 0's nearest B-cell is 0, and B-cell 0's nearest A-cell is 0.
        # A-cell 1's nearest B-cell is 2, and B-cell 2's nearest A-cell is 1.
        "ab": np.array([[0, 1, 2], [2, 0, 1]]),
        "ba": np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]]),
    }
    narrow = _mutual_nn_from(nbrs, k_anchor=1)
    assert narrow == [(0, 0), (1, 2)]
    # Widening the window admits pairs the narrow one rejected — proving the
    # slice is load-bearing rather than decorative.
    assert set(narrow) < set(_mutual_nn_from(nbrs, k_anchor=3))


# ---------------------------------------------------------------------------
# FindWeights — the kernel, the cap, and the guard
# ---------------------------------------------------------------------------


def _weight_setup(n_query=60, n_anchor_cells=40, dups=2, seed=0):
    rng = np.random.default_rng(seed)
    emb = rng.normal(size=(n_query, 8))
    cells = np.repeat(np.arange(n_anchor_cells), dups)
    scores = rng.uniform(0.2, 1.0, size=len(cells))
    return emb, cells, scores


def test_the_furthest_neighbour_contributes_exactly_zero():
    """``d̃ = 1 - d/dₖ`` is zero at the k-th neighbour, so ``1 - exp(0) = 0``.

    A Gaussian in the raw distance gives it ``exp(-sd²)`` — small but non-zero —
    so this single assertion separates the two kernels.
    """
    from sklearn.neighbors import NearestNeighbors

    emb, cells, scores = _weight_setup()
    k = 20
    W = _anchor_weights(emb, cells, scores, k_weight=k, sd_weight=1.0)

    uniq = np.unique(cells)
    _, idx = NearestNeighbors(n_neighbors=k).fit(emb[uniq]).kneighbors(emb)
    for c in range(emb.shape[0]):
        furthest_cell = uniq[idx[c, k - 1]]
        rows = np.where(cells == furthest_cell)[0]
        assert np.allclose(W[rows, c], 0.0), "the k-th neighbour must weigh nothing"


def test_k_weight_caps_anchors_not_anchor_cells():
    """FindWeightsC expands each neighbour cell into all of its anchor rows and
    stops at k entries — so a duplicated anchor set uses far fewer cells."""
    emb, cells, scores = _weight_setup(dups=3)
    k = 30
    W = _anchor_weights(emb, cells, scores, k_weight=k, sd_weight=1.0)
    for c in range(emb.shape[0]):
        nz = np.nonzero(W[:, c])[0]
        assert len(nz) <= k
        # 3 anchors per cell means well under k distinct cells contribute
        assert len(set(cells[nz])) < k


def test_weights_are_column_stochastic():
    emb, cells, scores = _weight_setup()
    W = _anchor_weights(emb, cells, scores, k_weight=20, sd_weight=1.0)
    assert np.allclose(W.sum(axis=0), 1.0)
    assert (W >= 0).all()


def test_a_zero_score_anchor_gets_no_weight_anywhere():
    emb, cells, scores = _weight_setup()
    scores = scores.copy()
    scores[5] = 0.0
    W = _anchor_weights(emb, cells, scores, k_weight=20, sd_weight=1.0)
    assert np.allclose(W[5, :], 0.0)


def test_sd_weight_enters_inside_the_exponent_and_flattens_the_kernel():
    """``sd_weight`` divides through ``(2/sd)²``, so raising it *widens* the kernel.

    Larger sd → smaller divisor → a larger exponent → every weight saturates
    towards 1 → after column normalisation they even out. A Gaussian multiplied
    by the score has the opposite response to its bandwidth, so the direction of
    this trend is itself a discriminator between the two kernels.
    """
    emb, cells, scores = _weight_setup()
    tops = [
        _anchor_weights(emb, cells, scores, k_weight=20, sd_weight=sd)[:, 0].max()
        for sd in (0.5, 1.0, 2.0, 4.0)
    ]
    assert tops == sorted(tops, reverse=True), (
        "a larger sd_weight must spread the weight, not concentrate it"
    )
    assert tops[0] > tops[-1] * 1.5


def test_pca_loadings_are_exact_not_randomized():
    """The anchor path needs the *exact* SVD, not sklearn's randomized default.

    Checked against ``np.linalg.svd`` — the mathematically exact answer, not
    against the implementation. Randomized SVD is accurate in the leading
    components and drifts in the trailing ones; reciprocal PCA standardizes each
    projected dimension by its own SD, which is not rotation-invariant, so a
    drifted trailing axis silently becomes a different anchor set. Against
    Seurat this was the difference between 45% and 100% RPCA anchor recall.
    """
    rng = np.random.default_rng(11)
    # The matrix has to be big enough that sklearn would pick its *randomized*
    # solver (it switches once max(shape) > 500). On a smaller fixture sklearn
    # runs exact LAPACK anyway, so swapping the implementation back would change
    # nothing and this guard would pass while proving nothing.
    n_feat, n_cells, k = 900, 700, 30
    mat = (rng.normal(size=(n_feat, 60)) @ rng.normal(size=(60, n_cells))
           + 0.6 * rng.normal(size=(n_feat, n_cells)))
    mat = mat - mat.mean(axis=1, keepdims=True)      # gene-centred, like scale.data

    got = _pca_loadings(mat, k)
    _u, _s, vt = np.linalg.svd(mat.T, full_matrices=False)
    want = vt[:k, :].T

    assert got.shape == (n_feat, k)
    # Compare subspaces axis by axis; a whole-column sign flip is not an error.
    per_axis = np.abs((got * want).sum(axis=0))
    assert per_axis.min() > 1 - 1e-8, (
        f"trailing axes drifted from the exact SVD: min |cos| = {per_axis.min()}"
    )
    # Orthonormal, and deterministic regardless of the (unused) seed.
    assert np.allclose(got.T @ got, np.eye(k), atol=1e-10)
    assert np.array_equal(got, _pca_loadings(mat, k, seed=999))


def test_integrate_layers_makes_the_largest_batch_the_reference():
    """Seurat corrects the smaller dataset onto the larger one.

    ``PairwiseIntegrateReference`` reverses the merge pair when the second
    object is bigger, so the reference is whichever batch has more cells. Taking
    the first batch is invisible on an even split and pulls the wrong way on a
    real one.
    """
    import shanuz.anchors as A
    from shanuz.integration import _integrate_anchor_reduction

    captured = {}
    # _integrate_anchor_reduction imports from .anchors at call time, so the
    # patch has to land on the source module, not on shanuz.integration.
    original = A.find_integration_anchors

    def spy(objects, **kw):
        captured["reference"] = kw.get("reference")
        captured["sizes"] = [len(o.cell_names()) for o in objects]
        return original(objects, **kw)

    objs = _pair(n=90)
    small = objs[0].subset(cells=objs[0].cell_names()[:40])
    merged = small.merge([objs[1]])
    merged.meta_data["batch"] = (["a"] * 40) + (["b"] * 90)
    scale_data(merged)
    run_pca(merged, n_pcs=10)

    A.find_integration_anchors = spy
    try:
        _integrate_anchor_reduction(
            merged, group_by="batch", reduction="cca",
            new_reduction="int", k_weight=15, dims=10,
        )
    finally:
        A.find_integration_anchors = original

    assert captured["sizes"] == [40, 90], "the spy must actually have run"
    assert captured["reference"] == 1, (
        "the larger batch (index 1) must be the reference"
    )


def test_integrate_data_leaves_the_reference_untouched():
    """Only the query half is corrected — so any comparison that looks at
    reference cells is comparing the log-normalised data to itself."""
    objs = _pair()
    anchors = find_integration_anchors(objs, reduction="cca", dims=10)
    merged = integrate_data(anchors, k_weight=20)

    ref_cells = objs[0].cell_names()
    ia = merged.get_assay("integrated")
    feats = anchors.anchor_features
    got = np.asarray(ia.layer_data("data", features=feats))
    names = merged.cell_names()
    cols = [names.index(c) for c in ref_cells]

    raw = objs[0].get_assay().layer_data("data", features=feats)
    raw = raw.toarray() if sp.issparse(raw) else np.asarray(raw)
    assert np.allclose(got[:, cols], raw)

    query_cols = [names.index(c) for c in objs[1].cell_names()]
    qraw = objs[1].get_assay().layer_data("data", features=feats)
    qraw = qraw.toarray() if sp.issparse(qraw) else np.asarray(qraw)
    assert not np.allclose(got[:, query_cols], qraw), "the query must move"


# ----------------------------------------------------------------------
# IntegrateEmbeddings — the Seurat v5 IntegrateLayers path
#
# These pin the v5 path as a *different algorithm* from v4's IntegrateData,
# not a wrapper over it. Running v4 behind this API produced an embedding that
# agreed with Seurat's on 1 of 30 dimensions, and the suite could not tell:
# both routes return a (cells x dims) array of the right shape.
# ----------------------------------------------------------------------


def _batched(n_a=90, n_b=140, n_pcs=10):
    """One object with two batches, scaled and reduced — the v5 entry state."""
    objs = _pair(n=max(n_a, n_b))
    a = objs[0].subset(cells=objs[0].cell_names()[:n_a])
    b = objs[1].subset(cells=objs[1].cell_names()[:n_b])
    merged = a.merge([b])
    merged.meta_data["batch"] = (["a"] * n_a) + (["b"] * n_b)
    scale_data(merged)
    run_pca(merged, n_pcs=n_pcs)
    return merged, n_a, n_b


@pytest.mark.parametrize("method", ["cca", "rpca"])
def test_integrate_layers_corrects_the_input_reduction_not_a_fresh_pca(method):
    """The v5 output must live in ``orig_reduction``'s basis.

    The sharpest observable: ``IntegrateEmbeddings`` copies the *reference*
    batch through untouched, so those rows are bit-identical to the input
    reduction. The v4 route (correct expression, re-scale, re-run PCA) moves
    every cell, reference included, and lands in a new basis entirely.
    """
    merged, n_a, n_b = _batched()
    before = np.asarray(merged.reductions["pca"].cell_embeddings).copy()
    cells = merged.cell_names()

    integrate_layers(merged, method=method, group_by="batch",
                     new_reduction="out", k_weight=15, dims=10)
    after = np.asarray(merged.reductions["out"].cell_embeddings)

    assert after.shape == before.shape
    # batch "b" is larger, so it is the reference and must be untouched
    ref = np.array([c.startswith("b2_") for c in cells])
    assert ref.sum() == n_b
    assert np.abs(after[ref] - before[ref]).max() == 0.0, (
        "reference cells moved — this is not IntegrateEmbeddings"
    )
    assert np.abs(after[~ref] - before[~ref]).max() > 0.0, "the query must move"


def test_integrate_embeddings_keeps_the_input_loadings():
    """``CreateDimReducObject(..., loadings = Loadings(reductions)[, dims])`` —
    the corrected reduction still maps back to genes. A re-run PCA would have
    its own, different loadings."""
    merged, _a, _b = _batched()
    integrate_layers(merged, method="rpca", group_by="batch",
                     new_reduction="out", k_weight=15, dims=10)
    got = merged.reductions["out"].feature_loadings
    want = merged.reductions["pca"].feature_loadings
    assert got is not None and np.array_equal(np.asarray(got), np.asarray(want))


@pytest.mark.parametrize("method", ["cca", "rpca"])
def test_v5_path_disables_the_anchor_filter(method):
    """``CCAIntegration``/``RPCAIntegration`` both call FindIntegrationAnchors
    with ``k.filter = NA``. v4's default of 200 belongs to the object-list API
    only, and applying it here drops ~15% of Seurat's CCA anchors."""
    import shanuz.anchors as A

    captured = {}
    original = A.find_integration_anchors

    def spy(objects, **kw):
        captured.update(kw)
        return original(objects, **kw)

    merged, _a, _b = _batched()
    A.find_integration_anchors = spy
    try:
        integrate_layers(merged, method=method, group_by="batch",
                         new_reduction="out", k_weight=15, dims=10)
    finally:
        A.find_integration_anchors = original

    assert captured, "the spy must actually have run"
    assert captured["k_filter"] is None, (
        f"v5 must pass k_filter=None, got {captured['k_filter']!r}"
    )


def test_rpca_rescales_each_batch_but_cca_does_not():
    """The two v5 methods differ here and it is deliberate: ``RPCAIntegration``
    runs ``ScaleData`` per object, ``CCAIntegration`` slices the object's
    existing ``scale.data``. Reciprocal PCA needs each batch centred on its own
    mean; CCA is given the pooled centring on purpose."""
    import shanuz.anchors as A

    seen = {}
    original = A.find_integration_anchors

    def spy(objects, **kw):
        seen[kw["reduction"]] = [
            _anchor_feature_matrix(o, o.feature_names(), "scale.data")
            for o in objects
        ]
        return original(objects, **kw)

    A.find_integration_anchors = spy
    try:
        for method in ("cca", "rpca"):
            merged, _a, _b = _batched()
            integrate_layers(merged, method=method, group_by="batch",
                             new_reduction="out", k_weight=15, dims=10)
    finally:
        A.find_integration_anchors = original

    # per-batch scaling drives each batch's own gene means to ~0; pooled
    # scaling leaves each batch offset by its share of the batch effect.
    rpca_off = max(abs(m.mean(axis=1)).max() for m in seen["rpca"])
    cca_off = max(abs(m.mean(axis=1)).max() for m in seen["cca"])
    assert rpca_off < 1e-8, f"rpca batches are not re-scaled (offset {rpca_off:.3g})"
    assert cca_off > 1e-3, (
        f"cca batches were re-scaled (offset {cca_off:.3g}); it must inherit "
        "the pooled scale.data"
    )


def test_integrate_layers_requires_the_reduction_to_exist():
    """v5 corrects an existing reduction, so a missing one is a user error to
    report, not something to silently compute."""
    objs = _pair(n=90)
    merged = objs[0].merge([objs[1]])
    merged.meta_data["batch"] = (["a"] * 90) + (["b"] * 90)
    with pytest.raises(KeyError, match="orig_reduction"):
        integrate_layers(merged, method="cca", group_by="batch", k_weight=15)


def test_integrate_embeddings_rejects_a_reduction_missing_cells():
    """A reduction that does not span every dataset would silently correct a
    subset; Seurat errors on the same condition in ValidateParams."""
    from shanuz.anchors import integrate_embeddings

    merged, n_a, _b = _batched()
    objs = _pair(n=90)
    anchors = find_integration_anchors(objs, reduction="cca", dims=10)
    # a reduction covering only the first batch's cells
    partial = merged.reductions["pca"].subset(cells=merged.cell_names()[:n_a])
    with pytest.raises(ValueError, match="absent from reduction"):
        integrate_embeddings(anchors, partial, k_weight=15)
