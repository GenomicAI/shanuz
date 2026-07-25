"""Regression guards for the DE defects T-de found against Seurat 5.5.1.

Both defects here were invisible to the suite: `find_markers` had tests, they
passed throughout, and neither the fold-change formula nor the negative-binomial
statistic was pinned to anything outside shanuz itself.

Constants marked "R:" were read off a live Seurat 5.5.1 / MASS session on pbmc3k
clusters 0 vs 1 (695 vs 477 cells), not derived from shanuz's own output.
"""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from shanuz import create_shanuz_object
from shanuz.markers import PSEUDOCOUNT, find_markers
from shanuz.preprocessing import normalize_data


# ---------------------------------------------------------------------------
# avg_log2FC — where the pseudocount goes
# ---------------------------------------------------------------------------

def _seurat_log2fc(mat1, mat2, pseudocount=1.0):
    """Seurat 5's `log1pdata.mean.fxn`, transcribed.

        log2((rowSums(expm1(x)) + pseudocount) / NCOL(x))

    Written out here rather than imported so the test compares against the R
    formula, not against whatever shanuz currently does.
    """
    m1 = (np.expm1(mat1).sum(axis=1) + pseudocount) / mat1.shape[1]
    m2 = (np.expm1(mat2).sum(axis=1) + pseudocount) / mat2.shape[1]
    return np.log2(m1) - np.log2(m2)


@pytest.fixture
def two_group_object():
    """40 genes x 80 cells, with a block of genes silent in group A.

    The silent block is the point: that is where adding the pseudocount to the
    mean rather than the sum does its damage.
    """
    rng = np.random.default_rng(11)
    n1 = n2 = 40
    g = 40
    counts = rng.poisson(3.0, size=(g, n1 + n2)).astype(float)
    counts[:10, :n1] = 0.0          # silent in group A only
    counts[10:15, :] = 0.0          # silent everywhere
    obj = create_shanuz_object(
        sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(g)],
        cell_names=[f"c{i}" for i in range(n1 + n2)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n1 + ["B"] * n2
    return obj


def test_log2fc_matches_seurats_formula(two_group_object):
    obj = two_group_object
    res = find_markers(obj, "A", "B", test_use="wilcox",
                       logfc_threshold=0, min_pct=0)
    assay = obj.assays["RNA"]
    data = assay.layer_data("data")
    X = np.asarray(data.todense() if hasattr(data, "todense") else data, dtype=float)
    names = list(assay.features())
    cells = list(assay.cells())
    idents = np.asarray([str(i) for i in obj.idents])
    i1 = [k for k, c in enumerate(cells) if idents[k] == "A"]
    i2 = [k for k, c in enumerate(cells) if idents[k] == "B"]
    want = pd.Series(_seurat_log2fc(X[:, i1], X[:, i2]), index=names)
    got = res["avg_log2FC"]
    shared = got.index.intersection(want.index)
    assert np.abs(got[shared] - want[shared]).max() < 1e-12


def test_pseudocount_goes_on_the_sum_not_the_mean(two_group_object):
    """The specific error, pinned.

    Adding the pseudocount to the mean is Seurat 4's formula; it floors the fold
    change of a gene that is silent in one group at roughly -1 instead of the
    large negative value the group size implies.
    """
    obj = two_group_object
    res = find_markers(obj, "A", "B", test_use="wilcox",
                       logfc_threshold=0, min_pct=0)
    assay = obj.assays["RNA"]
    data = assay.layer_data("data")
    X = np.asarray(data.todense() if hasattr(data, "todense") else data, dtype=float)
    names = list(assay.features())
    idents = np.asarray([str(i) for i in obj.idents])
    i1 = np.where(idents == "A")[0]
    i2 = np.where(idents == "B")[0]

    # Both candidate formulas, computed from the same matrix. An earlier version
    # of this test asserted only `fc < -4`, which BOTH formulas satisfy on this
    # fixture — it carried the name of the defect while proving nothing. Assert
    # the value, and assert it is not the other one.
    sum_form = _seurat_log2fc(X[:, i1], X[:, i2])
    e1, e2 = np.expm1(X[:, i1]).mean(axis=1), np.expm1(X[:, i2]).mean(axis=1)
    mean_form = np.log2(e1 + 1.0) - np.log2(e2 + 1.0)

    silent_in_a = [f"g{i}" for i in range(10)]
    rows = [names.index(g) for g in silent_in_a]
    got = res.loc[silent_in_a, "avg_log2FC"].to_numpy()
    assert np.abs(got - sum_form[rows]).max() < 1e-12
    assert np.abs(got - mean_form[rows]).min() > 1.0, (
        "silent-gene fold changes match the mean-pseudocount formula"
    )


def test_a_gene_silent_everywhere_has_zero_fold_change(two_group_object):
    res = find_markers(two_group_object, "A", "B", test_use="wilcox",
                       logfc_threshold=0, min_pct=0)
    # Both groups reduce to log2(pseudocount / n) with equal n, so the difference
    # is exactly 0 — a real constraint, since unequal group sizes would not
    # cancel and would show up as spurious signal.
    silent = [f"g{i}" for i in range(10, 15)]
    assert np.abs(res.loc[silent, "avg_log2FC"]).max() < 1e-12


def test_unequal_group_sizes_still_cancel_for_a_silent_gene():
    """With n1 != n2 a silent gene is log2(1/n1) - log2(1/n2), not 0."""
    rng = np.random.default_rng(3)
    counts = rng.poisson(3.0, size=(5, 90)).astype(float)
    counts[0, :] = 0.0
    obj = create_shanuz_object(
        sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(5)],
        cell_names=[f"c{i}" for i in range(90)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * 60 + ["B"] * 30
    res = find_markers(obj, "A", "B", test_use="wilcox",
                       logfc_threshold=0, min_pct=0)
    want = np.log2(PSEUDOCOUNT / 60) - np.log2(PSEUDOCOUNT / 30)
    assert res.loc["g0", "avg_log2FC"] == pytest.approx(want, abs=1e-12)


def test_logfc_threshold_filters_on_the_corrected_value(two_group_object):
    """The fold change is not just reported — it gates what comes back.

    This is why the formula mattered enough to be a defect rather than a display
    quirk: on pbmc3k the wrong formula returned 2,298 genes at a 0.25 threshold
    where Seurat returned 11,931.
    """
    obj = two_group_object
    loose = find_markers(obj, "A", "B", test_use="wilcox",
                         logfc_threshold=0, min_pct=0)
    strict = find_markers(obj, "A", "B", test_use="wilcox",
                          logfc_threshold=4.0, min_pct=0)
    assert set(strict.index) == set(loose.index[np.abs(loose["avg_log2FC"]) >= 4.0])
    assert len(strict) < len(loose)


# ---------------------------------------------------------------------------
# negbinom — ML dispersion + Wald, not moment dispersion + LRT
# ---------------------------------------------------------------------------

def test_negbinom_matches_glm_nb_wald():
    """Against statsmodels' ML fit directly, which is what MASS::glm.nb does.

    Seurat reads `summary(glm.nb(...))$coef[2, 4]` — a Wald p-value on an
    ML-estimated dispersion. The previous implementation fixed the dispersion by
    method of moments and ran a likelihood-ratio test instead, which on pbmc3k
    put HLA-DRA at 5.5e-128 against R's 1.1e-321.
    """
    import statsmodels.api as sm

    rng = np.random.default_rng(5)
    n1 = n2 = 60
    counts = np.vstack([
        np.r_[rng.poisson(8.0, n1), rng.poisson(2.0, n2)],   # clearly different
        np.r_[rng.poisson(3.0, n1), rng.poisson(3.0, n2)],   # null
    ]).astype(float)
    obj = create_shanuz_object(
        sp.csc_matrix(counts), assay="RNA", feature_names=["hit", "null"],
        cell_names=[f"c{i}" for i in range(n1 + n2)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n1 + ["B"] * n2
    res = find_markers(obj, "A", "B", test_use="negbinom",
                       logfc_threshold=0, min_pct=0)

    grp = np.r_[np.zeros(n1), np.ones(n2)]
    X = sm.add_constant(grp)
    for gene, row in (("hit", counts[0]), ("null", counts[1])):
        want = sm.NegativeBinomial(row, X).fit(disp=0, maxiter=200).pvalues[1]
        got = res.loc[gene, "p_val"]
        # `abs=0` matters. `pytest.approx` carries a default *absolute* tolerance
        # of 1e-12, which is larger than these p-values, so the plain form would
        # call any two tiny numbers equal and prove nothing. Setting abs=0 leaves
        # a pure relative comparison, which works at both ends of the range.
        #
        # 1e-3 rather than something tighter because both sides are iterative ML
        # fits and agree to ~5 significant figures, not to machine precision. It
        # is still 26 orders of magnitude away from the moment-dispersion LRT
        # this replaced, which is the thing being guarded against.
        assert got == pytest.approx(want, rel=1e-3, abs=0), gene

    assert res.loc["hit", "p_val"] < 1e-6
    assert res.loc["null", "p_val"] > 0.01


def test_negbinom_is_not_a_likelihood_ratio_test():
    """Guards the statistic, not just the plumbing.

    A moment-dispersion LRT and an ML-dispersion Wald test agree in direction and
    roughly in ranking, so a test asserting only "the hit is significant" passes
    under both. This asserts the actual number.
    """
    import statsmodels.api as sm
    from scipy.stats import chi2

    rng = np.random.default_rng(9)
    n = 60
    row = np.r_[rng.poisson(9.0, n), rng.poisson(2.0, n)].astype(float)
    obj = create_shanuz_object(
        sp.csc_matrix(row.reshape(1, -1)), assay="RNA", feature_names=["hit"],
        cell_names=[f"c{i}" for i in range(2 * n)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n + ["B"] * n
    got = find_markers(obj, "A", "B", test_use="negbinom",
                       logfc_threshold=0, min_pct=0).loc["hit", "p_val"]

    grp = np.r_[np.zeros(n), np.ones(n)]
    m = row.mean()
    alpha = max((row.var() - m) / (m * m), 1e-6) if row.var() > m else 1e-6
    fam = sm.families.NegativeBinomial(alpha=alpha)
    full = sm.GLM(row, np.column_stack([np.ones(2 * n), grp]), family=fam).fit()
    red = sm.GLM(row, np.ones((2 * n, 1)), family=fam).fit()
    old = float(chi2.sf(max(red.deviance - full.deviance, 0.0), df=1))
    # Orders of magnitude, for the same reason as above: `!= approx(old)` with
    # approx's default 1e-12 absolute tolerance is trivially false for any pair
    # of small p-values, so that form of the assertion could never fail.
    assert abs(np.log10(got) - np.log10(old)) > 1.0, (
        f"negbinom still matches the old moment-dispersion LRT "
        f"(got {got:.3e}, old LRT {old:.3e})"
    )


# ---------------------------------------------------------------------------
# find_all_markers — Seurat's return.thresh, and the order rows come back in
# ---------------------------------------------------------------------------
#
# T1's handoff found this. On PBMC 3k, two clusters ended up with *identical*
# cell membership on both sides; on those two, shanuz returned 190 and 383
# genes against Seurat's 151 and 242. Every extra row was a gene shanuz agreed
# with Seurat about numerically (max |avg_log2FC diff| 4.9e-15) and that Seurat
# simply does not return: `FindAllMarkers(return.thresh = 1e-2)` drops
# everything with p_val >= 0.01. Applying the same filter reproduced Seurat's
# gene sets exactly, 151/151 and 242/242.

def _three_cluster_object():
    """3 clusters x 30 cells, each with one strong private marker.

    The rest of the genes are pure noise, which is what makes the threshold
    test meaningful: without a filter they still pass the pct and logfc
    pre-filters often enough to reach the output table.
    """
    rng = np.random.default_rng(4)
    g, n = 30, 30
    counts = rng.negative_binomial(4, 0.4, size=(g, 3 * n)).astype(float)
    counts[0, :n] += 30          # g0 marks cluster A
    counts[1, n:2 * n] += 30     # g1 marks cluster B
    counts[2, 2 * n:] += 30      # g2 marks cluster C
    obj = create_shanuz_object(
        sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(g)],
        cell_names=[f"c{i}" for i in range(3 * n)],
    )
    normalize_data(obj)
    obj.idents = ["0"] * n + ["1"] * n + ["2"] * n
    return obj


def test_find_all_markers_applies_seurats_return_thresh():
    from shanuz.markers import find_all_markers

    obj = _three_cluster_object()
    unfiltered = find_all_markers(obj, return_thresh=None)
    default = find_all_markers(obj)

    # Anti-vacuity: the fixture must actually produce rows the filter removes,
    # or "all p_val < 0.01" below would hold no matter what the code did.
    assert (unfiltered["p_val"] >= 1e-2).any(), (
        "fixture yields no non-significant rows — the threshold assertions "
        "below would be vacuous"
    )
    assert len(default) < len(unfiltered)
    assert (default["p_val"] < 1e-2).all()
    # and the filter is the *only* difference: same rows otherwise
    kept = unfiltered[unfiltered["p_val"] < 1e-2]
    assert set(zip(default["cluster"], default["gene"])) == \
        set(zip(kept["cluster"], kept["gene"]))


def test_find_all_markers_return_thresh_is_configurable():
    from shanuz.markers import find_all_markers

    obj = _three_cluster_object()
    loose = find_all_markers(obj, return_thresh=0.5)
    tight = find_all_markers(obj, return_thresh=1e-6)
    assert len(tight) < len(loose)
    assert (loose["p_val"] < 0.5).all()
    assert (tight["p_val"] < 1e-6).all()


def _tied_marker_object():
    """Six markers that separate the two groups *perfectly*, at six magnitudes.

    A Wilcoxon test on a perfect separation gives the maximal U statistic no
    matter how large the gap is, so all six genes come back with the same
    p-value while their fold changes differ by an order of magnitude. That is
    the real situation this ordering exists for: on PBMC 3k, Seurat reported
    between 40 and 302 genes per cluster tied at exactly p = 0.

    The `arange` matters. Every marker must have the *same* pattern of ties
    within itself — 25 distinct positive values against 25 zeros — because the
    normal approximation carries a tie correction that depends on that pattern.
    Draw the group-A values randomly and the six p-values land within a few
    parts in a thousand of each other instead of on top of each other, which is
    close enough to look tied and not close enough to be.
    """
    rng = np.random.default_rng(7)
    g, n = 20, 25
    counts = rng.integers(0, 3, size=(g, 2 * n)).astype(float)
    for i, lift in enumerate([10, 20, 40, 80, 160, 320]):
        counts[i, :n] = lift + np.arange(n)                 # group A, well clear
        counts[i, n:] = 0.0
    obj = create_shanuz_object(
        sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(g)],
        cell_names=[f"c{i}" for i in range(2 * n)],
    )
    normalize_data(obj)
    obj.idents = ["0"] * n + ["1"] * n
    return obj


def test_find_all_markers_breaks_p_value_ties_on_fold_change():
    """Seurat orders each cluster by `order(gde$p_val, -gde[, 2])`.

    This is not cosmetic. Without the tie-break, "the top N markers" — which is
    what every tutorial in the suite prints — is decided by incoming row order
    among the genes that tie, and the strongest markers are exactly the ones
    that tie.
    """
    from shanuz.markers import find_all_markers

    got = find_all_markers(_tied_marker_object())
    checked = 0
    for cluster, block in got.groupby("cluster", sort=False):
        p = block["p_val"].to_numpy()
        fc = block["avg_log2FC"].to_numpy()
        assert np.all(np.diff(p) >= 0), f"cluster {cluster} not sorted by p_val"
        ties = p[:-1] == p[1:]
        checked += int(ties.sum())
        assert np.all(np.diff(fc)[ties] <= 0), (
            f"cluster {cluster} does not break p_val ties on descending log2FC"
        )
    # Anti-vacuity: with no tied pairs anywhere, the loop above asserts nothing.
    assert checked >= 4, (
        f"fixture produced only {checked} tied adjacent pairs — the tie-break "
        f"assertion would be vacuous"
    )


def test_cluster_labels_are_ordered_numerically_past_ten():
    """`sorted()` on strings puts "10" before "2".

    Seurat's identities are a factor ordered 0, 1, ... 10, 11, so a dataset
    with eleven or more clusters came back in a different cluster order —
    silently, and only past ten clusters, which is why every tutorial in the
    suite (eight and nine clusters on PBMC 3k) missed it.
    """
    from shanuz.markers import _ident_sort_key

    labels = [str(i) for i in range(12)]
    assert sorted(labels, key=_ident_sort_key) == labels
    assert sorted(labels) != labels          # the bug this replaces
    # non-numeric labels still sort, and sort after the numeric ones
    mixed = sorted(["B", "2", "10", "CD4 T"], key=_ident_sort_key)
    assert mixed == ["2", "10", "B", "CD4 T"]
