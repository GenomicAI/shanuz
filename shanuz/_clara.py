"""CLARA — Clustering LARge Applications (Kaufman & Rousseeuw), ported from R.

A faithful port of the ``clara`` C routine in R's **cluster** package (2.1.8.2,
``src/clara.c``), which is what Seurat's ``HTODemux`` reaches for by default.
CLARA is k-medoids for datasets too big to run PAM on directly: rather than
building the full n x n dissimilarity matrix, it repeatedly draws a small
sub-sample, runs PAM on that, assigns *every* object to the resulting medoids,
and keeps whichever sub-sample gave the lowest total dissimilarity.

Why a port and not a library
----------------------------
The output is fed to ``HTODemux``'s background fit, so it has to agree with what
R produces, and the details that decide agreement are all non-obvious:

* **The RNG is clara's own**, not R's. ``clara(rngR = FALSE)`` — the default, and
  what Seurat uses — draws from a 16-bit LCG seeded to 0, so the result is
  deterministic given the data alone and ``set.seed`` has no effect on it. R's
  ``HTODemux`` calls ``set.seed(seed)`` regardless; for the clara path that call
  does nothing.
* **The swap rule is not PAM's.** At ``pamLike = FALSE`` (again the default)
  ``bswap2`` uses the pre-2011 clara update, which the C source itself flags as
  "seems a bit illogical". A textbook PAM, or a k-medoids library, silently
  disagrees here.
* **Ties break inconsistently, and on purpose.** BUILD takes the *last* candidate
  attaining the maximum (the C carries a comment that ``<`` instead of ``<=``
  does *not* work); SWAP and the final assignment take the *first*. Getting these
  backwards perturbs cluster membership on data with duplicate points — exactly
  what hashtag counts are full of.
* **Cluster numbering is by first appearance**, not by medoid index: ``selec``
  permutes the medoids into the order their clusters are first encountered while
  scanning objects in order.

Rounding is part of the algorithm
---------------------------------
clara takes a swap on *any* improvement below zero, and R really does accept
swaps worth ``-2.2e-16``. A one-ulp difference in a single distance can therefore
flip a swap, and with it the winning sub-sample and the entire clustering. Sums
here are accumulated in the same sequence as the C loops (see ``_pairwise_dys``,
``_bswap2`` and ``_selec``) rather than handed to ``np.sum``, whose pairwise
reordering is enough to change the result.

That sensitivity has a consequence worth stating plainly: **R's clara is not
reproducible across CPU architectures.** ``clara.c`` built for arm64 contracts
``clk += d*d`` into a fused multiply-add — one rounding — while the same source
built for baseline x86_64, whose ISA has no FMA, rounds twice. The two binaries
return materially different clusterings for a few percent of inputs, so there is
no single "R answer" to match. This port follows plain IEEE double arithmetic,
which is what numpy gives everywhere and what ``clara.c`` gives on x86_64;
against that reference it is exact. Where the two R builds agree — the overwhelming
majority of inputs, and every realistic hashtag panel tested — shanuz agrees too.
Only the Euclidean metric is ported, the only one ``HTODemux`` asks for.

``tests/test_hto.py`` pins this module's output to clustering vectors captured
from real ``cluster::clara`` runs rather than to a re-derivation of the algorithm
in Python. Those fixtures use values on a 1/4 grid so their arithmetic is exact
and they cannot drift with any of the above.
"""
from __future__ import annotations

import bisect

import numpy as np

__all__ = ["clara", "clara_sampsize"]


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def clara_sampsize(n: int, k: int) -> int:
    """R's default ``sampsize`` for :func:`clara`: ``min(n, 40 + 2 * k)``."""
    return int(min(n, 40 + 2 * k))


def clara(
    x: np.ndarray,
    k: int,
    samples: int = 5,
    sampsize: int | None = None,
) -> np.ndarray:
    """Cluster ``x`` into ``k`` medoid-based groups, as R's ``cluster::clara``.

    Mirrors ``clara(x, k, samples = samples, sampsize = sampsize)`` at the R
    defaults Seurat relies on — ``metric = "euclidean"``, ``rngR = FALSE``,
    ``pamLike = FALSE``, ``stand = FALSE`` — and reproduces its cluster
    assignments exactly.

    There is no ``seed`` argument because clara at ``rngR = FALSE`` does not take
    one: its sampling is driven by a built-in generator that always starts from
    the same state, so the result is a deterministic function of ``x``, ``k``,
    ``samples`` and ``sampsize``.

    Parameters
    ----------
    x        : ``(n_observations, n_features)`` float array. Note this is
               observations-by-features, the orientation R's ``clara`` expects —
               callers holding a features-by-cells matrix must transpose first.
    k        : number of clusters; ``1 <= k <= n - 1``.
    samples  : number of sub-samples to draw (R's default is 5; Seurat passes
               100). Ignored when ``sampsize >= n``, where one sample is the whole
               dataset and further draws could not differ.
    sampsize : observations per sub-sample; defaults to ``min(n, 40 + 2 * k)``.

    Returns
    -------
    numpy.ndarray
        ``(n_observations,)`` array of **0-based** cluster labels. R's ``clara``
        returns 1-based labels; the shift is the only intentional difference from
        it, for consistency with the rest of shanuz (and scikit-learn).
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"clara() needs a 2-D array; got shape {x.shape}.")
    n, jpp = x.shape

    k = int(k)
    if k < 1 or k > n - 1:
        raise ValueError(
            f"The number of clusters should be at least 1 and at most n-1; "
            f"got k={k} with n={n}."
        )

    nsam = clara_sampsize(n, k) if sampsize is None else int(sampsize)
    if nsam < max(2, k + 1):
        raise ValueError(
            f"'sampsize' should be at least {max(2, k + 1)} "
            f"= max(2, 1 + number of clusters); got {nsam}."
        )
    if nsam > n:
        raise ValueError(
            f"'sampsize' = {nsam} should not be larger than the number of "
            f"objects, {n}."
        )
    samples = int(samples)
    if samples < 1:
        raise ValueError(f"'samples' should be at least 1; got {samples}.")

    return _cl_clara(x, k, samples, nsam)


# ----------------------------------------------------------------------
# The built-in RNG (clara.c :: randm)
# ----------------------------------------------------------------------


class _Randm:
    """clara's own generator — a 16-bit LCG, always started from 0.

    Verbatim from ``clara.c``::

        *nrun = (*nrun * 5761 + 999) & 0177777;
        return ((double) (*nrun) / 65536.);

    ``0177777`` is octal for 65535, so the mask is ``% 65536``. The period is
    65536, which the C source acknowledges is short and deems "good enough". This
    is *not* R's RNG: ``set.seed`` cannot reach it.
    """

    __slots__ = ("nrun",)

    def __init__(self) -> None:
        self.nrun = 0

    def __call__(self) -> float:
        self.nrun = (self.nrun * 5761 + 999) & 0o177777
        return self.nrun / 65536.0


# ----------------------------------------------------------------------
# Distances
# ----------------------------------------------------------------------


def _pairwise_dys(x: np.ndarray, nsel: np.ndarray) -> np.ndarray:
    """Full ``(nsam, nsam)`` Euclidean distance matrix for the selected rows.

    Stands in for ``dysta2`` plus C's ``ind_2`` condensed indexing. ``ind_2(i, i)``
    returns 0 and ``dys[0]`` is held at 0. permanently, so a dense matrix with a
    zero diagonal is exactly equivalent and far easier to read.

    The accumulation loops over features rather than calling a vectorized norm:
    C sums ``clk += (x[lj] - x[kj])^2`` one feature at a time, and numpy's
    pairwise summation would otherwise add them in a different order. The values
    agree to well under a rounding error either way, but ``bswap2`` compares
    distances with ``==``, so matching the exact float is worth the loop — it runs
    once per feature, not per pair.
    """
    sub = x[nsel]
    nsam = sub.shape[0]
    acc = np.zeros((nsam, nsam), dtype=float)
    for j in range(sub.shape[1]):
        diff = sub[:, j][:, None] - sub[:, j][None, :]
        acc += diff * diff
    # dysta2 scales by jpp/npres, which is exactly 1.0 without missing data.
    return np.sqrt(acc)


def _dist_to_medoids(x: np.ndarray, medoids: np.ndarray) -> np.ndarray:
    """``(n, k)`` *squared* Euclidean distances from every row to each medoid.

    Squared, because both ``selec`` and ``resul`` compare sums-of-squares and only
    take the square root once a winner is chosen. Feature-wise accumulation for
    the same reason as :func:`_pairwise_dys`.
    """
    med = x[medoids]
    acc = np.zeros((x.shape[0], med.shape[0]), dtype=float)
    for j in range(x.shape[1]):
        diff = x[:, j][:, None] - med[:, j][None, :]
        acc += diff * diff
    return acc


# ----------------------------------------------------------------------
# PAM on a sub-sample (clara.c :: bswap2)
# ----------------------------------------------------------------------


def _bswap2(kk: int, dys: np.ndarray, s: float) -> np.ndarray:
    """PAM's BUILD then SWAP over one sub-sample; returns a medoid mask.

    ``dys`` is the ``(nsam, nsam)`` distance matrix from :func:`_pairwise_dys` and
    ``s`` its maximum. Ported from ``bswap2`` at ``pam_like = FALSE``, i.e. the
    swap clara has used since before 2011, which is *not* the one ``pam()`` uses.
    """
    n = dys.shape[0]
    s = s * 1.1 + 1.0  # strictly larger than every dissimilarity

    # ---- BUILD: greedily seed kk medoids ----
    nrepr = np.zeros(n, dtype=bool)
    dysma = np.full(n, s)

    for _ in range(kk):
        # beter[i] = sum_j max(0, dysma[j] - dys[i, j]), the gain from adding i.
        beter = np.zeros(n)
        for j in range(n):
            cmd = dysma[j] - dys[:, j]
            beter += np.where(cmd > 0.0, cmd, 0.0)

        # C scans i ascending under `if (ammax <= beter[i])`, so among candidates
        # tied at the maximum the LAST one wins. The source explicitly notes that
        # tightening `<=` to `<` breaks the algorithm.
        cand = np.flatnonzero(~nrepr)
        nmax = cand[np.flatnonzero(beter[cand] == beter[cand].max())[-1]]

        nrepr[nmax] = True
        dysma = np.minimum(dysma, dys[nmax])

    if kk == 1:
        return nrepr

    # ---- SWAP: exchange a medoid for a non-medoid while that helps ----
    while True:
        med = np.flatnonzero(nrepr)
        # dysma[j] = d(j, closest medoid); dysmb[j] = d(j, 2nd closest).
        # Sorting reproduces C's running two-smallest scan: with kk >= 2 both
        # entries are always overwritten (s exceeds every distance), and ties
        # yield identical values whatever the order.
        near = np.sort(dys[med], axis=0)
        dysma, dysmb = near[0], near[1]

        non_med = np.flatnonzero(~nrepr)

        # dz[h, i] = the change in total dissimilarity from swapping medoid i out
        # for non-medoid h. Laid out h-major because that is the order C scans in,
        # which decides who wins a tie; see the argmin below.
        #
        # Accumulated one j at a time rather than summed along an axis, because C
        # adds the terms in that order and the residue is load-bearing: the terms
        # routinely cancel to a mathematical zero, and whether the leftover lands
        # a hair below zero or exactly on it decides whether the swap happens at
        # all. R really does accept swaps worth -2.2e-16.
        dz = np.zeros((len(non_med), len(med)))
        for j in range(n):
            dj_i = dys[med, j][None, :]       # (1, n_med) -- d(i, j)
            dj_h = dys[non_med, j][:, None]   # (n_non_med, 1) -- d(h, j)

            # The pam_like = FALSE branch. Note `dysmb[j] > dys[i, j]` tests the
            # *removed* medoid's distance where pam() tests the candidate's --
            # and inside this branch dys[i, j] == dysma[j], so it is really asking
            # whether the 2nd-closest is further than the closest, which it all
            # but always is. The C comment calls this "a bit illogical"; it is
            # also what R computes, so it is what we compute.
            on_i = dj_i == dysma[j]
            small = np.where(dysmb[j] > dj_i, dj_h, dysmb[j])
            dz += np.where(
                on_i,
                -dysma[j] + small,
                np.where(dj_h < dysma[j], -dysma[j] + dj_h, 0.0),
            )

        # C nests `for h { for i { if (dzsky > dz) ... } }` from dzsky = 1, so the
        # winner is the first pair attaining the minimum scanning h-major -- which
        # is exactly argmin on the row-major flat view. Note this is the opposite
        # tie-break to BUILD's, and h-major rather than i-major: both matter only
        # on ties, and hashtag counts tie often.
        flat = int(np.argmin(dz))
        if dz.flat[flat] >= 0.0:  # no improving swap left
            return nrepr
        h_pos, i_pos = divmod(flat, len(med))
        nrepr[non_med[h_pos]] = True
        nrepr[med[i_pos]] = False


# ----------------------------------------------------------------------
# Whole-dataset assignment (clara.c :: selec / resul)
# ----------------------------------------------------------------------


def _assign(x: np.ndarray, medoids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Assign every row of ``x`` to its nearest medoid.

    Returns ``(labels, dist)`` with ``labels`` indexing into ``medoids`` and
    ``dist`` the Euclidean distance to the chosen one. Shared by ``selec`` (which
    sums ``dist`` into the sample's score) and ``resul`` (which keeps the labels),
    since the two do the same work in the no-missing-data case.
    """
    dsq = _dist_to_medoids(x, medoids)
    # Both C loops replace the incumbent only on a strictly smaller distance, so
    # the first medoid attaining the minimum wins -- which is np.argmin's rule.
    labels = np.argmin(dsq, axis=1)
    # A medoid is assigned to its own cluster even where some other medoid sits
    # at distance 0 from it (duplicate rows): C skips the comparison entirely
    # when the candidate *is* the object, so this assignment cannot be beaten.
    labels[medoids] = np.arange(len(medoids))
    dist = np.sqrt(dsq[np.arange(x.shape[0]), labels])
    return labels, dist


def _selec(x: np.ndarray, medoids: np.ndarray) -> tuple[float, np.ndarray]:
    """Score one sub-sample's medoids over the whole dataset (``selec``).

    Returns ``(zb, medoids_reordered)`` where ``zb`` is the summed distance —
    lower is better — and the medoids come back permuted into the order their
    clusters are *first encountered* scanning objects in order, which is the
    ordering that ends up deciding R's cluster numbering.
    """
    labels, dist = _assign(x, medoids)
    # C accumulates `*zb += dnull` one object at a time, and zb is compared across
    # sub-samples with `>`, so the rounding has to match or near-ties between two
    # samples resolve the wrong way. np.sum() reorders (pairwise summation);
    # cumsum is sequential by definition, and its last element is the same float
    # C arrives at.
    zb = float(np.cumsum(dist)[-1]) if dist.size else 0.0

    # C tracks first appearance in new[] and permutes nr[] by it. Every cluster
    # is guaranteed non-empty (each medoid holds itself), so new[] always fills.
    seen: list[int] = []
    flagged = np.zeros(len(medoids), dtype=bool)
    for lab in labels:
        if not flagged[lab]:
            flagged[lab] = True
            seen.append(int(lab))
            if len(seen) == len(medoids):
                break
    return zb, medoids[np.array(seen, dtype=int)]


# ----------------------------------------------------------------------
# The driver (clara.c :: cl_clara)
# ----------------------------------------------------------------------


def _draw_sample(
    rng: _Randm, n: int, n_sam: int, kk: int, nrx: np.ndarray, kall: bool,
    jran: int, lrg_sam: bool,
) -> np.ndarray:
    """Draw one sub-sample's indices, sorted ascending (0-based).

    Reproduces the index-drawing block of ``cl_clara``. Once a valid sample has
    been seen, and unless we are sampling more than half the data, the running
    best medoids are seeded into the sample first and the remainder drawn around
    them — Kaufman & Rousseeuw's "each sub-dataset is forced to contain the
    medoids obtained from the best sub-dataset until then".

    ``n_sam`` is the number of indices to *draw*, which is the complement's size
    when ``lrg_sam``; the caller inverts.
    """
    def draw() -> int:
        # C: rand_k = 1 + (int)(rnn * randm(&nrun)), clamped; 1-based there.
        rand_k = int(n * rng())
        return min(rand_k, n - 1)

    nsel: list[int] = []

    # nunfs (the count of samples abandoned to missing data) is always 0 here --
    # we reject NaNs up front -- so C's `nunfs + 1 != jran` reduces to `jran != 1`.
    if kall and jran != 1 and not lrg_sam:
        nsel = sorted(int(v) for v in nrx[:kk])
    else:
        while True:
            rand_k = draw()
            if kall and rand_k in nrx[:kk]:
                continue
            break
        nsel.append(rand_k)
        if len(nsel) == n_sam:
            return np.array(nsel, dtype=int)

    # C runs this as a do-while, so it always adds at least one index.
    while True:
        while True:
            rand_k = draw()
            if kall and lrg_sam and rand_k in nrx[:kk]:
                continue
            # C walks nsel linearly for the first entry >= rand_k and inserts
            # there, keeping it sorted; a redraw on an exact hit. Same position.
            pos = bisect.bisect_left(nsel, rand_k)
            if pos < len(nsel) and nsel[pos] == rand_k:
                continue  # already sampled -- redraw
            nsel.insert(pos, rand_k)
            break
        if len(nsel) >= n_sam:
            return np.array(nsel, dtype=int)


def _cl_clara(x: np.ndarray, kk: int, nran: int, nsam: int) -> np.ndarray:
    """The resampling loop of ``cl_clara``; returns 0-based labels."""
    if not np.isfinite(x).all():
        # C's clara handles NA by tracking per-column missing codes and can bail
        # out with jstop; none of that is ported, so refuse rather than silently
        # disagree with R.
        raise ValueError(
            "clara() does not support missing or non-finite values; "
            "R's clara handles them by a separate code path that is not ported."
        )

    n = x.shape[0]
    nsamb = nsam * 2
    full_sample = n == nsam
    lrg_sam = n < nsamb  # sampling more than half -- draw the complement instead
    n_sam = n - nsam if lrg_sam else nsam

    rng = _Randm()
    kall = False
    zba = -1.0
    nrx = np.zeros(kk, dtype=int)

    for jran in range(1, nran + 1):
        if full_sample:
            nsel = np.arange(n)
        else:
            nsel = _draw_sample(rng, n, n_sam, kk, nrx, kall, jran, lrg_sam)
            if lrg_sam:
                # We hold the *unsampled* complement; invert it.
                mask = np.ones(n, dtype=bool)
                mask[nsel] = False
                nsel = np.flatnonzero(mask)

        dys = _pairwise_dys(x, nsel)
        # C maxes over dys[1..n_dys], skipping the permanently-zero dys[0]; with
        # non-negative distances that is just the matrix maximum.
        s = float(dys.max()) if dys.size else 0.0

        nrepr = _bswap2(kk, dys, s)
        zb, medoids = _selec(x, nsel[nrepr])

        if not kall or zba > zb:  # first proper sample, or a new best
            kall = True
            zba = zb
            nrx = medoids

        if full_sample:
            break  # further samples would be identical

    # resul(): assign the entire dataset to the winning medoids.
    labels, _ = _assign(x, nrx)
    return labels
