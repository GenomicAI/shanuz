"""SCTransform fidelity against R's sctransform / MASS.

Every ``R = ...`` constant below was produced by running R 4.6.1 (sctransform
0.4.3, MASS 7.3) on the *exact* arrays ``_fixtures()`` builds. The inputs come
from ``numpy.random.default_rng``, whose PCG64 stream is stable across NumPy
versions and platforms, so the fixtures regenerate identically without shipping
data files and without needing R at test time. To refresh them, dump the arrays
to text and feed them to the matching R function.

Why this file exists: shanuz's SCT model was silently wrong. Its Poisson GLM was
faithful (intercept and slope matched R at Spearman 1.0), but a moment estimator
stood in for ``theta.ml``, the regularization was smoothed against the arithmetic
rather than the geometric gene mean and targeted log(theta) rather than the
overdispersion factor, and residual variance was computed from residuals clipped
at sqrt(N/30) instead of sqrt(N). Together those flattened every residual: on
PBMC 3k the regularized theta came out *anti*-correlated with R's (Spearman
-0.89), residual variance ranked essentially at random against R (-0.07), and
only 414 of 3,000 variable features agreed. The tutorial's SCT run resolved 9
clusters where R found 12 and — the tell — fewer than the 11 from plain
log-normalization, inverting the vignette's whole point.

The end-to-end numbers are not asserted here: they need the cached PBMC 3k
download and take minutes (see tests/test_tutorial_smoke.py). What these pin are
the numerical primitives that made the model wrong, each against R directly.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from shanuz.preprocessing import percentage_feature_set
from shanuz.sctransform import (
    _bw_nrd,
    _bw_sj,
    _dispersion_par,
    _fit_nb_offset,
    _ksmooth_normal,
    _row_gmean,
    _row_var,
    _theta_from_dispersion_par,
    _theta_ml,
    sctransform,
)
from shanuz.shanuz import create_shanuz_object


def _fixtures():
    """Regenerate the exact arrays R was run on.

    The draw *order* is load-bearing: these all share one PCG64 stream, so
    reproducing R's inputs means reproducing this sequence. Keep it in sync with
    whatever script last generated the constants below.
    """
    rng = np.random.default_rng(0)
    f = {}
    f["x_norm"] = rng.normal(0, 1, 200)
    f["x_bimodal"] = np.concatenate([rng.normal(0, 1, 150), rng.normal(4, 0.5, 50)])
    f["x_gamma"] = rng.gamma(2.0, 1.0, 300)
    f["kx"] = np.sort(rng.uniform(0, 10, 100))
    f["ky"] = np.sin(f["kx"]) + rng.normal(0, 0.2, 100)
    f["mu1"] = np.repeat([0.5, 2.0, 10.0], 100)
    f["y1"] = rng.negative_binomial(1.5, 1.5 / (1.5 + f["mu1"])).astype(float)
    f["mu2"] = np.repeat([0.2, 1.0, 5.0], 100)
    f["y2"] = rng.negative_binomial(0.4, 0.4 / (0.4 + f["mu2"])).astype(float)
    f["mat"] = rng.poisson(1.5, size=(50, 20)).astype(float)
    return f


FIX = _fixtures()


# ----------------------------------------------------------------------
# bw.SJ / bw.nrd — bandwidth selection
# ----------------------------------------------------------------------


def test_bw_sj_matches_r():
    """R's bw.SJ, ported because SciPy has no equivalent.

    The bandwidth sets how hard each gene's fitted parameters are pulled toward
    the trend, so a "close enough" substitute would shift every residual;
    Silverman's rule on the normal sample gives 0.354 against bw.SJ's 0.391.

    Tolerance: R stops uniroot at tol=0.1*lower, so R's own value is only defined
    to a few tenths of a percent. We solve the same equation tightly, which lands
    inside that band (and in fact agrees to ~1e-12 on the two samples where R's
    uniroot happens to converge fully).
    """
    assert _bw_sj(FIX["x_norm"]) == pytest.approx(0.391029768240, rel=5e-3)
    # bimodal: the solve-the-equation step must find the *right* root, so a
    # bracketing bug shows here and not on the unimodal samples.
    assert _bw_sj(FIX["x_bimodal"]) == pytest.approx(0.366574506835, rel=5e-3)
    assert _bw_sj(FIX["x_gamma"]) == pytest.approx(0.260695523555, rel=5e-3)


def test_bw_nrd_matches_r():
    """R's bw.nrd — weights the step-1 gene sample toward the sparse tails.

    The gamma sample is here on purpose: bw.nrd divides the IQR by 1.34 while
    bw.SJ uses 1.349, and the two only diverge when the IQR term wins the min().
    On the normal sample sd() wins either way, so it cannot see the difference —
    this pair caught exactly that bug during the port.
    """
    assert _bw_nrd(FIX["x_norm"]) == pytest.approx(0.353991260460, rel=1e-9)
    assert _bw_nrd(FIX["x_gamma"]) == pytest.approx(0.384197300833, rel=1e-9)


# ----------------------------------------------------------------------
# ksmooth — the regularizing smoother
# ----------------------------------------------------------------------


def test_ksmooth_normal_matches_r():
    """R's ksmooth(kernel="normal"), which replaced a lowess call.

    R scales the kernel so its quartiles land at +/-0.25*bandwidth (sd =
    0.3706506*bandwidth) and truncates at 4 sd. Passing the bandwidth straight in
    as a standard deviation would silently over-smooth by ~2.7x.
    """
    got = _ksmooth_normal(FIX["kx"], FIX["ky"], np.array([0.5, 2.5, 5.0, 7.5, 9.5]),
                          bandwidth=1.5)
    assert got == pytest.approx([0.765303001477, 0.512248191167, -0.898102318480,
                                 0.892481226545, 0.403320441329], rel=1e-9)


# ----------------------------------------------------------------------
# theta.ml — the NB dispersion MLE
# ----------------------------------------------------------------------


def test_theta_ml_matches_r():
    """MASS::theta.ml, which replaced a method-of-moments estimator.

    This is the defect that mattered most. The moment estimator ranked genes
    almost independently of the MLE (Spearman 0.18 against R on PBMC 3k) and its
    denominator went negative on a third of genes, pinning them at the theta
    ceiling.
    """
    assert _theta_ml(FIX["y1"][None, :], FIX["mu1"][None, :])[0] == pytest.approx(
        1.266822522075, rel=1e-6)
    assert _theta_ml(FIX["y2"][None, :], FIX["mu2"][None, :])[0] == pytest.approx(
        0.348714279929, rel=1e-6)


def test_theta_ml_is_vectorised_over_genes():
    """Rows must stay independent — the batched form is not a shared fit."""
    rng = np.random.default_rng(3)
    mu = np.repeat([1.0, 4.0], 150)
    Y = np.stack([rng.negative_binomial(s, s / (s + mu)) for s in (0.5, 2.0, 8.0)])
    Y = Y.astype(float)
    MU = np.stack([mu, mu, mu])

    together = _theta_ml(Y, MU)
    apart = np.array([_theta_ml(Y[i:i + 1], MU[i:i + 1])[0] for i in range(3)])
    assert together == pytest.approx(apart, rel=1e-12)
    # and it recovers the ordering of the true dispersions
    assert together[0] < together[1] < together[2]


# ----------------------------------------------------------------------
# od_factor — what gets smoothed
# ----------------------------------------------------------------------


def test_dispersion_par_round_trips_through_theta():
    """The smoothed quantity must invert back to theta exactly.

    R smooths log10(1 + gmean/theta) and then recovers theta as
    gmean/(10^disp - 1). If the pair ever drift apart, every regularized theta is
    silently wrong while everything still runs.
    """
    log10_gmean = np.log10(np.array([0.001, 0.01, 0.1, 1.0, 10.0, 100.0]))
    theta = np.array([0.05, 0.3, 1.0, 4.0, 25.0, 500.0])
    disp = _dispersion_par(log10_gmean, theta)
    assert _theta_from_dispersion_par(log10_gmean, disp) == pytest.approx(theta, rel=1e-9)


def test_dispersion_par_is_the_od_factor_not_log_theta():
    """Pins R's theta_regularization="od_factor" against the "log_theta" option.

    shanuz used to smooth log10(theta). The two agree only in the limit where
    gmean/theta dominates, so a fixture with a *low* gmean separates them — which
    is most genes in a real count matrix.
    """
    log10_gmean = np.log10(np.array([0.01, 0.01]))
    theta = np.array([0.5, 50.0])
    got = _dispersion_par(log10_gmean, theta)
    assert got == pytest.approx(np.log10(1.0 + 0.01 / theta), rel=1e-12)
    assert not np.allclose(got, np.log10(theta)), "this is log_theta, not od_factor"


def test_dispersion_par_is_monotone_decreasing_in_theta():
    """More overdispersion (smaller theta) => larger od factor.

    The property the smoother relies on; inverting it is what turned R's
    Spearman +0.96 into shanuz's -0.89.
    """
    theta = np.array([0.01, 0.1, 1.0, 10.0, 100.0])
    disp = _dispersion_par(np.full(5, np.log10(0.5)), theta)
    assert np.all(np.diff(disp) < 0)


# ----------------------------------------------------------------------
# gene summaries
# ----------------------------------------------------------------------


def test_row_gmean_matches_r():
    """sctransform's row_gmean. The regularization x-axis is log10 of *this*.

    shanuz previously smoothed against the arithmetic mean, which on sparse
    counts is a different quantity entirely (median log10 -1.59 vs -1.77 on
    PBMC 3k) and so put every gene in the wrong neighbourhood.
    """
    got = _row_gmean(sp.csr_matrix(FIX["mat"]), eps=1.0)
    assert got[:5] == pytest.approx(
        [1.069360948862, 1.391111754785, 0.988368502418,
         0.898148223767, 1.062694036697], rel=1e-9)
    assert got.sum() == pytest.approx(59.655085846863, rel=1e-9)


def test_row_gmean_differs_from_arithmetic_mean_on_sparse_counts():
    """Why it matters: on sparse counts the two are not interchangeable."""
    counts = np.zeros((2, 100))
    counts[0, :5] = 40.0     # rare and bright
    counts[1, :] = 2.0       # ubiquitous and dim
    gm = _row_gmean(sp.csr_matrix(counts), eps=1.0)
    am = counts.mean(axis=1)
    assert gm[0] < am[0] / 5          # sparsity crushes the geometric mean
    assert gm[1] == pytest.approx(am[1])   # a constant row: the two agree


def test_row_var_matches_dense_ddof1():
    rng = np.random.default_rng(4)
    m = rng.poisson(0.7, size=(30, 80)).astype(float)
    assert _row_var(sp.csr_matrix(m)) == pytest.approx(m.var(axis=1, ddof=1), rel=1e-10)


# ----------------------------------------------------------------------
# v2 offset fit
# ----------------------------------------------------------------------


def test_fit_nb_offset_matches_r_glm_nb():
    """v2's per-gene fit is R's glm.nb(y ~ 1 + offset(log_umi)).

    Reference generated with the same arrays this rebuilds; tolerance reflects
    that R alternates IRLS and theta.ml to its own convergence threshold rather
    than to machine precision.
    """
    rng = np.random.default_rng(7)
    n = 400
    umi = np.round(np.exp(rng.normal(np.log(3000), 0.4, n)))
    mu = np.exp(np.log(2e-4) + np.log(umi))
    y = rng.negative_binomial(1.2, 1.2 / (1.2 + mu)).astype(float)

    b0, theta = _fit_nb_offset(y[None, :], np.log10(umi), gene_chunk=500)
    assert b0[0] == pytest.approx(-8.707863138517, rel=1e-4)
    assert theta[0] == pytest.approx(2.215043455562, rel=1e-3)


# ----------------------------------------------------------------------
# the two clip ranges — the defect that cost 3 clusters
# ----------------------------------------------------------------------


def _object(n_cells=300, n_genes=60, seed=5):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(1.0, size=(n_genes, n_cells)).astype(float)
    counts[0, :] = 0.0
    counts[0, :10] = 3000.0          # a rare, bright marker
    return create_shanuz_object(
        counts=sp.csr_matrix(counts), assay="RNA", min_cells=0, min_features=0,
        project="spike", feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"c{i}" for i in range(n_cells)],
    )


def _scale_data(obj):
    sd = obj.assays["SCT"].layers["scale.data"]
    return sd.toarray() if sp.issparse(sd) else np.asarray(sd)


@pytest.mark.parametrize("flavor", ["v1", "v2"])
def test_residual_variance_ignores_the_scale_data_clip(flavor):
    """R uses two different clips; using the tighter one for both is the bug.

    Residual *variance* — which ranks the variable features — is computed from
    residuals clipped at sqrt(N). Only the stored scale.data is clipped to the
    much tighter sqrt(N/30) (or a caller's clip_range). So residual variance must
    be *invariant* to clip_range.

    shanuz previously fed clip_range into both, capping residual variance at 8.05
    on PBMC 3k where R reached 71.75 and leaving the ranking uncorrelated with
    R's (Spearman -0.07). Under that code this test fails: res_var tracks the
    clip. Asserting invariance pins the two clips apart without depending on any
    particular dataset.
    """
    tight, wide = _object(), _object()
    sctransform(tight, n_features=20, min_cells=0, seed=0, vst_flavor=flavor,
                clip_range=(-1.0, 1.0))
    sctransform(wide, n_features=20, min_cells=0, seed=0, vst_flavor=flavor,
                clip_range=(-40.0, 40.0))

    rv_tight = tight.assays["SCT"].meta_data["residual_variance"]
    rv_wide = wide.assays["SCT"].meta_data["residual_variance"]
    assert rv_tight.values == pytest.approx(rv_wide.values, rel=1e-12), (
        "residual_variance changed with clip_range — it must come from the "
        "sqrt(N) clip, not the scale.data clip"
    )
    # and the clip really was in force for scale.data, so the above is not
    # vacuously true because clip_range was ignored altogether
    assert np.abs(_scale_data(tight)).max() <= 1.0 + 1e-9
    assert np.abs(_scale_data(wide)).max() > 1.0


@pytest.mark.parametrize("flavor", ["v1", "v2"])
def test_scale_data_defaults_to_the_tight_seurat_clip(flavor):
    """The other half: scale.data defaults to sqrt(N/30), as Seurat does."""
    obj = _object()
    n = len(obj.cell_names())
    sctransform(obj, n_features=20, min_cells=0, seed=0, vst_flavor=flavor)
    assert np.abs(_scale_data(obj)).max() <= np.sqrt(n / 30.0) + 1e-9


# ----------------------------------------------------------------------
# flavors
# ----------------------------------------------------------------------


def test_v2_marks_non_overdispersed_genes_poisson():
    """v2 excludes genes with variance <= mean and models them as pure Poisson.

    v1 regularizes every gene, so the same gene keeps a finite theta there. This
    is what R's "excluding poisson genes" v2 message refers to.
    """
    rng = np.random.default_rng(9)
    n_cells = 300
    counts = rng.poisson(3.0, size=(40, n_cells)).astype(float)  # var ~= mean
    counts[1] = rng.negative_binomial(0.3, 0.3 / (0.3 + 5.0), n_cells)

    def build():
        return create_shanuz_object(
            counts=sp.csr_matrix(counts), assay="RNA", min_cells=0, min_features=0,
            project="p", feature_names=[f"g{i}" for i in range(40)],
            cell_names=[f"c{i}" for i in range(n_cells)])

    v2 = build()
    sctransform(v2, n_features=10, min_cells=0, seed=0, vst_flavor="v2")
    assert np.isinf(v2.assays["SCT"].meta_data["theta"]["g0"]), \
        "a Poisson gene should get theta=inf under v2"

    v1 = build()
    sctransform(v1, n_features=10, min_cells=0, seed=0, vst_flavor="v1")
    assert np.isfinite(v1.assays["SCT"].meta_data["theta"]["g0"]), \
        "v1 regularizes every gene"


def test_default_flavor_is_v2_matching_seurat5():
    """Seurat 5's SCTransform defaults to vst.flavor="v2"; so do we.

    Pinned because the default silently decides which algorithm users get, and
    the two give materially different models.
    """
    default, v2 = _object(), _object()
    sctransform(default, n_features=20, min_cells=0, seed=0)
    sctransform(v2, n_features=20, min_cells=0, seed=0, vst_flavor="v2")
    assert default.assays["SCT"].meta_data["theta"].equals(
        v2.assays["SCT"].meta_data["theta"])


def test_unknown_flavor_raises():
    with pytest.raises(ValueError, match="vst_flavor"):
        sctransform(_object(), min_cells=0, vst_flavor="v3")


def test_vars_to_regress_still_removes_the_covariate():
    """The regression step survives the model rewrite."""
    rng = np.random.default_rng(12)
    n_cells, n_genes = 200, 40
    depth = rng.uniform(0.5, 2.0, n_cells)
    counts = rng.poisson(2.0 * depth[None, :], size=(n_genes, n_cells)).astype(float)
    obj = create_shanuz_object(
        counts=sp.csr_matrix(counts), assay="RNA", min_cells=0, min_features=0,
        project="p", feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"c{i}" for i in range(n_cells)])
    percentage_feature_set(obj, pattern=r"^g0$", col_name="pct")
    sctransform(obj, vars_to_regress=["pct"], n_features=10, min_cells=0, seed=0)

    sd = _scale_data(obj)
    pct = obj.meta_data["pct"].values
    assert pct.std() > 1e-9, "fixture must vary the covariate for this to mean anything"
    corr = [abs(np.corrcoef(row, pct)[0, 1]) for row in sd if row.std() > 1e-9]
    assert np.nanmax(corr) < 0.05
