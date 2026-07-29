"""What `selection_method="mvp"` computes, and which features it hands back.

PR #64 gave the dispersion path Seurat's column names — `mvp.mean`,
`mvp.dispersion`, `mvp.dispersion.scaled`, replacing three that had been
borrowed from the vst path. The *values* underneath were a different set of
quantities entirely, so the rename left three columns making a promise the
numbers did not keep. Four separate divergences:

**The statistics were computed on the wrong scale.** `CalcDispersion` calls
`FastExpMean` and `FastLogVMR`, and both undo the log first::

    mvp.mean       = log1p(mean(expm1(x)))
    mvp.dispersion = log(var(expm1(x)) / mean(expm1(x)))

truecell took the mean and variance of the log-normalized values directly and
added an `eps` inside each logarithm. On PBMC 3k the resulting `mvp.mean`
column ran 0-2 where Seurat's ran 1-7: not a rounding difference, a different
quantity.

**The bins were equal-frequency, not equal-width.** R's
`cut(x, breaks = 20)` lays 20 bins of equal width across the range of the mean;
truecell used percentiles of `log(mean)`. Since the scaled dispersion is a
z-score *within a bin*, changing which genes share a bin changes every value.

**Bins holding one gene were scored 0 rather than NaN.** R's `sd` of a single
value is NA and Seurat propagates it. A gene alone in a bin has no
within-bin spread to be measured against, and 0 is not that fact — it is a
dead-centre score the gene never earned, and dead centre is a real rank.

**`mean.cutoff` and `dispersion.cutoff` were accepted and discarded.** Seurat
has two selectors here, not one: `MVP` (`"mvp"` / `"mean.var.plot"`) keeps
every gene inside both cutoffs and ignores `nfeatures`; `DISP` (`"dispersion"`
/ `"disp"`) takes the top `nfeatures` and ignores the cutoffs. truecell ran the
second for all four spellings — so `mean.var.plot` returned exactly
`nfeatures` genes under a name that promises a cutoff. On PBMC 3k Seurat
returns 1,006 features for `mvp`; truecell returned 2,000. Both also rank by the
*raw* dispersion; truecell ranked by the scaled one, which reorders the list.

The golden values below come from Seurat 5.5.1 / SeuratObject 5.4.0, via
`Seurat:::CalcDispersion` and `FindVariableFeatures` on the fixture in
`_counts()`. To regenerate, dump `_counts()` to JSON and run it through
`CreateSeuratObject` |> `NormalizeData` |> `CalcDispersion`.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from truecell import create_truecell_object
from truecell.preprocessing import (
    _dispersion_hvg,
    find_variable_features,
    normalize_data,
)

N_GENES, N_CELLS = 60, 40


def _counts() -> np.ndarray:
    """Per-gene rates spanning three orders of magnitude.

    The spread is the point. Equal-width bins over a range like this are wildly
    uneven — the top bins hold one or two genes each — which is what makes the
    binning method, the singleton-bin NaN and the raw-vs-scaled ranking all
    visible at once. A flat Poisson fixture puts every gene in the middle bins,
    where all three choices give the same answer and none of them is tested.
    """
    rng = np.random.default_rng(0)
    rate = np.geomspace(0.05, 30.0, N_GENES)
    return rng.poisson(rate[:, None], size=(N_GENES, N_CELLS)).astype(float)


@pytest.fixture
def obj():
    o = create_truecell_object(
        counts=sp.csc_matrix(_counts()),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(N_GENES)],
        cell_names=[f"c{j}" for j in range(N_CELLS)],
    )
    normalize_data(o)
    return o


# ---------------------------------------------------------------------------
# Seurat 5.5.1, verbatim
# ---------------------------------------------------------------------------

R_MEAN = np.array([
    1.015434782189202, 1.0083342792361392, 1.5368899715599122,
    1.2306123980604995, 1.7098037573975582, 1.280121943378925,
    1.669756240086109, 1.259396387735282, 1.6645755935503108,
    1.9548715343306553, 1.295016555784717, 1.6790156596154713,
    1.6787610357053442, 1.48377692434249, 1.9422278270882576,
    2.212452669310992, 2.8657256808435316, 2.2620643079419085,
    2.830577179007586, 2.6535022404905515, 2.8720402800109306,
    2.175163253691229, 3.0030358518484777, 3.3646063929194847,
    3.0065043784104786, 3.3157239238412486, 2.947430064665999,
    3.7965272377897095, 3.614916531387552, 3.695211657615681,
    3.7000958511275766, 3.9475650018865394, 3.800685659475411,
    3.964829071166174, 4.053595182119726, 4.349248909458298,
    4.390888548676176, 4.614795635877549, 4.708787346715267,
    4.807144890760756, 4.994054490675863, 5.043562186568431,
    5.094110575719659, 5.171222889906178, 5.324168000325411,
    5.432462108898828, 5.61167847855422, 5.63422020727049,
    5.835839466733609, 5.756986637064261, 6.037307233113902,
    6.077312227044667, 6.190275747956429, 6.260725704469873,
    6.382110623825269, 6.526747966177989, 6.6014570724060215,
    6.694252609531169, 6.828310350730889, 6.930409847567599,
])

R_DISPERSION = np.array([
    4.254513314374922, 3.5253831717075803, 3.518454987653382,
    3.423371827616801, 3.485287423878462, 3.4928258932570024,
    3.825498210229419, 3.4635586810828207, 3.767638628427326,
    3.380950411536381, 3.515473654388456, 3.8050244511721036,
    3.4459963362252988, 3.449911690450033, 3.36486871744154,
    3.36310601706432, 3.4066005635067365, 3.5098393445077085,
    3.6079240220736764, 3.5181987953168967, 3.113677535705483,
    3.6119646717222897, 3.3895073024904487, 3.4017356887040098,
    3.255930567033634, 3.4775455957847243, 3.5551303296911403,
    3.5166464027846622, 3.345701488016725, 3.471472381116591,
    3.6407402603855323, 3.4589102875185915, 3.438755499091012,
    3.2091563129466802, 3.209210841675657, 3.256606950735305,
    3.2755779523139665, 3.7435578494715904, 3.643700701541898,
    3.5109363527417927, 3.1716502609917567, 2.9182715544860467,
    3.517170760397887, 3.6041040234354376, 3.581028901925284,
    3.624173687627673, 3.42475242431699, 3.3746953699174314,
    3.58436348931286, 3.510703396864211, 3.6165840888510856,
    3.206917129623096, 3.4779722821847727, 3.394937122335303,
    3.6331174552293573, 3.272190796713545, 3.4355953530945604,
    3.7542775091055036, 3.6096717685775976, 3.29760089341375,
])

R_SCALED = np.array([
    2.027201284331052, -0.27515302226572097, 0.7071067811865476,
    -0.597271453043979, -0.9788866563693474, -0.37795824984490406,
    0.8651048959120003, -0.4703745251408382, 0.5514976137738961,
    -0.5182123887812119, -0.30644403403560544, 0.754134167587746,
    -1.1918500209042926, -0.7071067811865476, -0.6345335205310447,
    -0.7071067811865496, 0.10037375904413952, 0.7071067811865454,
    1.1943212989823115, np.nan, -1.491305759753095, 1.1527459093122598,
    0.00749272584814182, -0.7071067811865434, -0.7183339720742351,
    0.7071067811865517, 0.9074519479527348, 0.43075552911295767, np.nan,
    0.11017718783018302, 1.3113900523213577, 0.02103004130698518,
    -0.12199881475754865, -1.7513539958139286, np.nan, -0.707106781186564,
    0.7071067811865309, 0.9496836721124282, 0.09399530425047713,
    -1.0436789763629015, -0.10216568810373171, -0.9449952799578517,
    1.0471609680615834, 0.04640235671587958, -1.022393410097353,
    0.975991053381494, 0.7071067811865476, -0.7071067811865476,
    0.2545029447977612, -1.1026596886432096, 0.8481567438454566,
    -1.1018921688220074, 0.849902998056682, 0.2519891707653254,
    1.0299764661167345, -0.9670516793070792, -0.06292478680965537,
    0.8587304644975832, 0.23916785743923705, -1.0978983219368221,
])

# `CalcDispersion(binning.method = "equal_frequency")` on the same fixture. Not
# the default, but pinned because it is the only path where a gene can fall
# *outside* every bin: the breaks are quantiles of the positive means, so the
# smallest one sits exactly on the first break and stays in range only because
# R passes `include.lowest = TRUE`. Without that the gene scores NaN and drops
# out of the selection, and no equal-width fixture can show it.
R_SCALED_EQUAL_FREQUENCY = np.array([
    1.4914754622563091, -0.3585884534518545, -1.1364031943024369,
    -0.6174277136212774, -0.47470853329323426, 0.2028736721798044,
    0.745517233755587, -0.5154592951831773, 0.3908859605468526,
    -0.5182123887812119, 0.8830080353235383, 1.1489401684210687,
    -0.6742316351278344, -1.0858817075033294, -0.6345335205310447,
    -1.1533742018060194, 0.12285321479498881, 0.5287708932458641,
    0.9328974422517646, 0.6246033085601603, -1.0557506570467534,
    1.1527459093122598, -0.07126831690858536, -0.09962254730692333,
    -0.962459334287628, 1.0460825775918283, 1.0337276511962104,
    -0.002905972658949763, -0.9464600302849118, -0.5132387890198477,
    1.3989875383476449, 1.1547005177339429, -0.8828427766688374,
    -0.577539358637862, -0.5771611590960808, -0.6113979730565118,
    -0.5426195791690849, 1.154017552225595, 0.8281563928044748,
    0.2827845890819698, -1.1109409818864464, -1.1468506864270551,
    0.4570203292263913, 0.6898303572006625, 0.35941279084690847,
    0.7706181987714569, -1.1300309896183696, -1.0832726774056252,
    0.8878880243158007, 0.1953846530898286, 0.8770869357377304,
    -1.088957859136962, 0.21187092339923147, -0.2096732572835613,
    1.0882123675559956, -0.8785391102724319, -0.44449537950937756,
    1.152653161305256, 0.4279284061807997, -1.1360861879766757,
])

# `VariableFeatures()` after FindVariableFeatures(nfeatures = 10) — in Seurat's
# order, which is by descending mvp.dispersion for both methods.
R_MVP_VARIABLE = ["g0", "g30", "g54", "g21", "g18", "g42"]
R_DISP_VARIABLE = ["g0", "g6", "g11", "g8", "g57", "g37", "g38", "g30", "g54",
                   "g45"]

# 8.9e-16 on PBMC 3k's 13,714 genes; the fixture is smaller and lands tighter.
# Loose enough for the C++ loop's accumulation order, far too tight for any of
# the four divergences above — the smallest of them moves `mvp.mean` by 1.8.
TOL = 1e-12


# ---------------------------------------------------------------------------
# 1. The three columns
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("column,expected", [
    ("mvp.mean", R_MEAN),
    ("mvp.dispersion", R_DISPERSION),
    ("mvp.dispersion.scaled", R_SCALED),
])
def test_columns_match_seurat(obj, column, expected):
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    got = obj.assays["RNA"].meta_data[column].to_numpy(dtype=float)
    np.testing.assert_allclose(got, expected, rtol=TOL, atol=TOL)


def test_mean_is_measured_after_undoing_the_log(obj):
    """Stated as an identity rather than a golden vector, so it says *why*.

    `log1p(mean(expm1(x)))` is not `mean(x)` and is not close to it: the mean
    of the un-logged values is dominated by the cells with the highest
    expression, which is exactly the property the mvp method is built on. On
    this fixture the two disagree by more than 1.8 in absolute terms and put
    the genes in a different order.
    """
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    got = obj.assays["RNA"].meta_data["mvp.mean"].to_numpy(dtype=float)

    dense = np.asarray(obj.assays["RNA"].layers["data"].todense())
    np.testing.assert_allclose(got, np.log1p(np.expm1(dense).mean(axis=1)),
                               rtol=TOL, atol=TOL)

    naive = dense.mean(axis=1)
    assert np.abs(got - naive).max() > 1.0, (
        "fixture cannot tell log1p(mean(expm1(x))) from mean(x)"
    )


def test_dispersion_is_the_variance_to_mean_ratio_of_the_unlogged_values(obj):
    """`log(v/rm)` with a *sample* variance, and with no epsilon anywhere.

    The epsilons that used to sit inside both logarithms were not a guard
    against log(0) — `FastLogVMR` has none and lets -Inf through — they shifted
    every value on a scale where the smallest variance here is order 10^4.
    """
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    got = obj.assays["RNA"].meta_data["mvp.dispersion"].to_numpy(dtype=float)

    e = np.expm1(np.asarray(obj.assays["RNA"].layers["data"].todense()))
    np.testing.assert_allclose(got, np.log(e.var(axis=1, ddof=1) / e.mean(axis=1)),
                               rtol=TOL, atol=TOL)

    # ddof is not cosmetic at 40 cells: N vs N-1 shifts every value by
    # log(40/39), which is larger than the gap between adjacent genes here.
    ddof0 = np.log(e.var(axis=1, ddof=0) / e.mean(axis=1))
    assert np.abs(got - ddof0).max() > 1e-3


# ---------------------------------------------------------------------------
# 2. Binning
# ---------------------------------------------------------------------------

def test_bins_are_equal_width_over_the_mean_not_percentiles_of_it(obj):
    """The two binnings give different z-scores for the same dispersions.

    Guarded by running both and requiring they disagree, because that is the
    thing at issue: if the fixture's genes happen to fall into the same groups
    either way, the golden-vector test above passes under the wrong binning and
    proves nothing about it.
    """
    data = obj.assays["RNA"].layers["data"]
    _, _, _, wide = _dispersion_hvg(data, 10, binning_method="equal_width")
    _, _, _, freq = _dispersion_hvg(data, 10, binning_method="equal_frequency")

    assert not np.allclose(wide, freq, equal_nan=True), (
        "fixture cannot distinguish equal-width from equal-frequency bins"
    )
    np.testing.assert_allclose(wide, R_SCALED, rtol=TOL, atol=TOL)
    np.testing.assert_allclose(freq, R_SCALED_EQUAL_FREQUENCY, rtol=TOL, atol=TOL)
    # Seurat scores every gene here, including the one sitting exactly on the
    # first break; a bare `searchsorted` would put it out of range.
    assert not np.isnan(freq).any()


def test_a_gene_detected_in_no_cell_scores_nan_under_either_binning():
    """The only way a gene falls outside *every* bin, and it happens for free.

    An undetected gene has `mvp.mean` 0, and `equal_frequency` builds its
    breaks from `feature.mean[feature.mean > 0]` — so the gene sits below the
    first break, `cut` returns NA, and Seurat scores it NA. Clamping it into
    the first bin instead (the obvious way to keep the index in range) gives it
    a real z-score against genes it has nothing in common with.

    `equal_width` reaches the same NaN by a different route: mean 0 becomes the
    bottom of the range, the undetected genes are alone down there with
    identical dispersions, and their within-bin standard deviation is 0.
    Verified both ways against Seurat 5.5.1: 2 NaN for equal_frequency, 3 for
    equal_width (the extra one is the singleton bin at the top).
    """
    counts = _counts()
    counts[[3, 17]] = 0.0
    obj = create_truecell_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(N_GENES)],
        cell_names=[f"c{j}" for j in range(N_CELLS)],
    )
    normalize_data(obj)
    data = obj.assays["RNA"].layers["data"]

    for binning, n_nan in [("equal_frequency", 2), ("equal_width", 3)]:
        top, means, dispersion, scaled = _dispersion_hvg(
            data, 10, binning_method=binning
        )
        assert means[3] == 0.0 and means[17] == 0.0
        # `feature.dispersion[is.na(...)] <- 0`: log(0/0) is NaN, and Seurat
        # replaces it. Leaving it NaN would be invisible in the scaled column
        # (NaN in, NaN out) but changes the sort, where NaN goes last and 0
        # goes wherever 0 belongs — below any gene with a negative dispersion.
        assert (dispersion[[3, 17]] == 0.0).all(), (
            f"{binning}: an undetected gene's dispersion is not 0"
        )
        assert 3 not in top and 17 not in top
        assert np.isnan(scaled[[3, 17]]).all(), (
            f"{binning}: an undetected gene was given a scaled dispersion"
        )
        assert int(np.isnan(scaled).sum()) == n_nan


def test_ties_in_dispersion_keep_assay_order():
    """R's `order()` is a radix sort, so tied genes come back in input order.

    Ties are not hypothetical: two genes with the same expression pattern have
    the same dispersion to the last bit, and every gene whose dispersion came
    back NA is set to exactly 0 before the sort. numpy's default quicksort is
    free to return either, and which one it picks decides which gene makes the
    `nfeatures` cut.
    """
    counts = _counts()
    counts[41] = counts[5]  # an exact duplicate, so the dispersions tie
    obj = create_truecell_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(N_GENES)],
        cell_names=[f"c{j}" for j in range(N_CELLS)],
    )
    normalize_data(obj)

    top, _, dispersion, _ = _dispersion_hvg(
        obj.assays["RNA"].layers["data"], N_GENES, select="disp"
    )
    assert dispersion[5] == dispersion[41], "fixture no longer produces a tie"
    order = list(top)
    assert order.index(5) < order.index(41)


def test_a_gene_alone_in_its_bin_scores_nan_not_zero(obj):
    """R's `sd` of one value is NA, and Seurat carries the NA through.

    Three genes here land in a bin of their own. Scoring them 0 puts them at
    the exact centre of a distribution that does not exist — and 0 is the
    median of this column, so the fabricated value is indistinguishable from a
    genuine middling gene by inspection.

    The difference shows up at the cutoff, which is why that is what is
    asserted: NaN fails every comparison and drops out, while 0 passes any
    threshold below it and is selected. The default `dispersion_cutoff` of
    (1, Inf) happens to reject both, so the distinction is only observable
    with a cutoff the fabricated zero would clear.
    """
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    got = obj.assays["RNA"].meta_data["mvp.dispersion.scaled"].to_numpy(dtype=float)

    assert np.isnan(R_SCALED).sum() == 3, "golden vector no longer has singletons"
    np.testing.assert_array_equal(np.isnan(got), np.isnan(R_SCALED))

    singletons = {f"g{i}" for i in np.where(np.isnan(R_SCALED))[0]}
    find_variable_features(obj, selection_method="mvp", nfeatures=10,
                           dispersion_cutoff=(-1.0, float("inf")))
    selected = set(obj.assays["RNA"].variable_features)
    assert not (selected & singletons), (
        f"{selected & singletons} passed a cutoff of -1, so their scaled "
        f"dispersion is a filled-in number rather than NaN"
    )
    assert len(selected) > 6, "cutoff did not widen; the check above is vacuous"


# ---------------------------------------------------------------------------
# 3. Two selectors, not one
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method,expected", [
    ("mvp", R_MVP_VARIABLE),
    ("mean.var.plot", R_MVP_VARIABLE),
    ("dispersion", R_DISP_VARIABLE),
    ("disp", R_DISP_VARIABLE),
])
def test_each_spelling_selects_what_seurat_selects(obj, method, expected):
    find_variable_features(obj, selection_method=method, nfeatures=10)
    assert list(obj.assays["RNA"].variable_features) == expected


def test_mvp_ignores_nfeatures_and_disp_honours_it(obj):
    """The two return different counts *and* different genes.

    `MVP` is cutoff-based, so `nfeatures` does not reach it: asking for 10 gets
    6 here and 1,006 on PBMC 3k. Routing all four spellings to the top-N
    selector — which is what truecell did — made `nfeatures` look respected and
    the cutoffs look applied, and neither was true of both.
    """
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    mvp = list(obj.assays["RNA"].variable_features)
    find_variable_features(obj, selection_method="disp", nfeatures=10)
    disp = list(obj.assays["RNA"].variable_features)

    assert len(mvp) == 6 and len(disp) == 10
    assert set(mvp) != set(disp), (
        "fixture cannot distinguish the two selectors"
    )


def test_mvp_applies_both_cutoffs(obj):
    """Widen one cutoff, get more genes; narrow it, get fewer.

    Pins that the parameters reach the computation at all. They were in the
    signature, documented, and threaded down into `_dispersion_hvg` — which
    then never read them.
    """
    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    base = len(obj.assays["RNA"].variable_features)

    find_variable_features(obj, selection_method="mvp", nfeatures=10,
                           dispersion_cutoff=(0.5, float("inf")))
    assert len(obj.assays["RNA"].variable_features) > base

    find_variable_features(obj, selection_method="mvp", nfeatures=10,
                           mean_cutoff=(0.1, 3.0))
    assert len(obj.assays["RNA"].variable_features) < base


def test_ranking_is_by_raw_dispersion_not_the_scaled_one(obj):
    """Both selectors `order(hvf.info$mvp.dispersion, decreasing = TRUE)`.

    The scaled value is what the *cutoff* tests, and ranking by it instead is
    an easy thing to assume — it is the more principled quantity of the two.
    Seurat does not, and on this fixture the two orders differ from the second
    gene onward, so the returned list is a different list, not a permutation of
    the same one.
    """
    data = obj.assays["RNA"].layers["data"]
    top, _, dispersion, scaled = _dispersion_hvg(data, 10, select="disp")

    by_raw = np.argsort(-dispersion, kind="stable")[:10]
    by_scaled = np.argsort(-scaled, kind="stable")[:10]
    assert list(by_raw) != list(by_scaled), (
        "fixture cannot distinguish the two rankings"
    )
    assert list(top) == list(by_raw)


# ---------------------------------------------------------------------------
# 4. The streaming path agrees with the dense one
# ---------------------------------------------------------------------------

def test_variable_feature_plot_draws_the_stored_mvp_dispersion(obj):
    """The plot has to read the fixed columns, or the fix stops at the table.

    `variable_feature_plot` looked for the vst columns and, not finding them
    after an mvp run, fell through to recomputing `E[x^2] - E[x]^2` off the
    data matrix — a third quantity, on a third scale, drawn without complaint
    under a y axis reading "Dispersion". Seurat plots `mvp.mean` against
    `mvp.dispersion`, so those are the numbers that have to be on the figure;
    the axis label alone cannot tell the two branches apart, which is why the
    scatter's own data is what gets checked.
    """
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from truecell.plotting import variable_feature_plot

    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    fig = variable_feature_plot(obj, label=False)
    try:
        ax = fig.axes[0]
        assert ax.get_ylabel() == "Dispersion"
        plotted = np.vstack([c.get_offsets() for c in ax.collections])
        order = np.lexsort((plotted[:, 1], plotted[:, 0]))
        want = np.column_stack([R_MEAN, R_DISPERSION])
        want = want[np.lexsort((want[:, 1], want[:, 0]))]
        np.testing.assert_allclose(plotted[order], want, rtol=1e-9, atol=1e-9)
    finally:
        plt.close(fig)


def test_rerunning_with_another_method_retires_the_previous_columns(obj):
    """One `variable_features` list, one set of statistics describing it.

    SeuratObject can hold both methods at once because it namespaces the
    columns by method and layer and makes you name one to read them back.
    truecell's `meta_data` has a single flat name per statistic, so a vst run
    followed by an mvp run left `variance.standardized` sitting beside an mvp
    selection — and `variable_feature_plot`, which picks its axis off whichever
    columns it finds, drew the vst figure over the mvp genes.
    """
    find_variable_features(obj, selection_method="vst", nfeatures=10)
    assert "variance.standardized" in obj.assays["RNA"].meta_data.columns

    find_variable_features(obj, selection_method="mvp", nfeatures=10)
    columns = set(obj.assays["RNA"].meta_data.columns)
    assert {"mvp.mean", "mvp.dispersion", "mvp.dispersion.scaled"} <= columns
    assert not columns & {"mean", "variance", "variance.expected",
                          "variance.standardized"}

    find_variable_features(obj, selection_method="vst", nfeatures=10)
    columns = set(obj.assays["RNA"].meta_data.columns)
    assert "variance.standardized" in columns
    assert not columns & {"mvp.mean", "mvp.dispersion", "mvp.dispersion.scaled"}


def test_sparse_and_dense_layers_give_the_same_statistics():
    """One layer format must not select different genes from another.

    The sparse path reconstructs the variance from the non-zeros plus a
    `(n_cells - nnz) * rm^2` term for the zeros it never visits, which is
    Seurat's own trick and not obviously equal to the dense two-pass sum. Where
    the two disagree even slightly, a tie or a near-tie in the dispersion
    ranking sends the two formats home with different variable features.
    """
    counts = _counts()
    args = dict(assay="RNA", feature_names=[f"g{i}" for i in range(N_GENES)],
                cell_names=[f"c{j}" for j in range(N_CELLS)])

    dense_obj = create_truecell_object(counts=counts.copy(), **args)
    sparse_obj = create_truecell_object(counts=sp.csc_matrix(counts), **args)
    for o in (dense_obj, sparse_obj):
        normalize_data(o)

    dense_layer = dense_obj.assays["RNA"].layers["data"]
    if sp.issparse(dense_layer):
        dense_layer = np.asarray(dense_layer.todense())

    a = _dispersion_hvg(dense_layer, 10)
    b = _dispersion_hvg(sparse_obj.assays["RNA"].layers["data"], 10)
    assert list(a[0]) == list(b[0])
    for dense_col, sparse_col in zip(a[1:], b[1:]):
        np.testing.assert_allclose(dense_col, sparse_col, rtol=1e-12, atol=1e-12)
