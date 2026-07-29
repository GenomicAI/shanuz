"""Declared bands for the numbers that move, and the handoff check they rest on.

Two numbers in the R-comparison suite are expected to differ from Seurat's, and
both were recorded as prose. `deseq2`'s overlap with Seurat's top 50 genes is a
*divergence measurement* — truecell tests pseudobulk samples where `DESeq2DETest`
tests cells — and JackStraw's PC cutoff moves with the seed, because R's
`JackRandom` seeds each replicate from its loop index and truecell seeds from its
`seed` argument. A sentence in a vignette does not fail, so a genuine regression
that landed anywhere inside the expected spread read as ordinary variation. This
module gives each of those a :class:`Band` with a stated reason, and gives the
numbers that should be *exact* a band as well, so the parity tables in
`de_vignette.md` and `dimreduc_vignette.md` are checked rather than described.

It also carries the precondition without which none of those bands mean
anything. Both `--report` paths align the R side onto the Python side by name —
`reindex` in the dim-reduction tutorial, an index intersection in the DE one —
and neither noticed when the R run predated the handoff it was answering. On
this working copy the DE report was comparing a 25 July Python run against a
19 July R run taken on a *different* cluster assignment, and printed a full
parity table: `pct.1` and `pct.2` are pure counts over the shared cells, and they
disagreed for 12,491 of 13,712 genes. The numbers looked like a `wilcox`
regression rather than a stale file. :class:`StaleReferenceError` is raised now
instead.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np

__all__ = [
    "Band", "Verdict", "StaleReferenceError",
    "check_bands", "render_verdicts",
    "PCT_ROUNDING_TOL", "check_shared_groups", "check_same_cells",
    "check_same_features",
]


class StaleReferenceError(RuntimeError):
    """The R reference on disk was not produced from the current handoff.

    Raised rather than warned: every parity number downstream of the handoff is
    meaningless once the two sides are describing different cells or genes, and
    a warning in the middle of a table of six-decimal figures does not read as
    "ignore all of these".
    """


# ---------------------------------------------------------------------------
# Bands
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Band:
    """An inclusive range a measurement is asserted to fall in, and why.

    ``why`` is not decoration. A band with no stated reason cannot be widened or
    tightened by the next reader without guessing whether the number moved
    because the port drifted or because the band was arbitrary to begin with.
    """

    low: float
    high: float
    why: str
    fmt: str = ".4g"

    def __post_init__(self) -> None:
        if self.high < self.low:
            raise ValueError(f"empty band: [{self.low}, {self.high}]")
        if not self.why.strip():
            raise ValueError("a band must say why it is where it is")

    @property
    def exact(self) -> bool:
        """True when the band admits a single value — a parity assertion."""
        return self.low == self.high

    def holds(self, value: float) -> bool:
        """Is ``value`` inside the band? ``NaN`` never is.

        A missing measurement is the case this whole module exists for: a report
        that prints ``nan`` and carries on is exactly how a stale reference went
        unnoticed, so absence fails rather than passing silently.
        """
        v = float(value)
        if math.isnan(v):
            return False
        return self.low <= v <= self.high

    def describe(self) -> str:
        if self.exact:
            return f"= {self.low:{self.fmt}}"
        if self.high == math.inf:
            return f">= {self.low:{self.fmt}}"
        if self.low == -math.inf:
            return f"<= {self.high:{self.fmt}}"
        return f"[{self.low:{self.fmt}}, {self.high:{self.fmt}}]"


@dataclass(frozen=True)
class Verdict:
    """One measurement judged against its band."""

    name: str
    value: float
    band: Band

    @property
    def ok(self) -> bool:
        return self.band.holds(self.value)


def check_bands(bands: dict[str, Band], measured: dict[str, float]) -> list[Verdict]:
    """Judge every declared band against ``measured``.

    A band with no measurement is reported as ``NaN`` and therefore fails, so
    deleting the line that computes a number cannot quietly retire its band.
    """
    return [Verdict(name, float(measured.get(name, float("nan"))), band)
            for name, band in bands.items()]


def render_verdicts(verdicts: Sequence[Verdict], title: str,
                    indent: str = "  ") -> bool:
    """Print the band table. Returns True when every band holds."""
    if not verdicts:
        return True
    width = max(len(v.name) for v in verdicts)
    print(f"\n{indent}{title}")
    print(f"{indent}{'-' * (width + 40)}")
    for v in verdicts:
        mark = "ok " if v.ok else "OUT"
        val = "nan" if math.isnan(v.value) else f"{v.value:{v.band.fmt}}"
        print(f"{indent}{mark}  {v.name:<{width}}  {val:>10}  "
              f"{v.band.describe()}")
    failed = [v for v in verdicts if not v.ok]
    if failed:
        print(f"\n{indent}{len(failed)} of {len(verdicts)} outside the declared band:")
        for v in failed:
            print(f"{indent}  {v.name}: {v.band.why}")
    return not failed


# ---------------------------------------------------------------------------
# The handoff check
# ---------------------------------------------------------------------------

# Seurat rounds pct.1/pct.2 to three decimals inside FindMarkers, so two runs
# over the same cells agree to half a unit in the last place and no closer.
PCT_ROUNDING_TOL = 5e-4


def check_shared_groups(py, r, *, source: str, tol: float = PCT_ROUNDING_TOL,
                        columns: Iterable[str] = ("pct.1", "pct.2")) -> float:
    """Verify an R marker table was computed on the current group assignment.

    ``pct.1`` and ``pct.2`` are the fraction of each group's cells in which the
    gene is detected: pure counting over the handoff in ``groups.csv``, with no
    statistics in the way. If the two sides agree about which cells are in which
    group they agree here to R's rounding, whatever they do about p-values —
    which makes this the one column that separates "the port changed" from "the
    R run is older than the handoff it is being compared against".

    Returns the largest absolute difference. Raises
    :class:`StaleReferenceError` when it exceeds ``tol``.
    """
    shared = py.index.intersection(r.index)
    if len(shared) == 0:
        raise StaleReferenceError(
            f"{source} shares no genes with the Python output — regenerate it "
            f"with `Rscript tutorials/pbmc3k_de_verify.R`.")
    worst = 0.0
    worst_col = ""
    for col in columns:
        if col not in py.columns or col not in r.columns:
            continue
        d = float((py.loc[shared, col] - r.loc[shared, col]).abs().max())
        if d > worst:
            worst, worst_col = d, col
    if worst > tol:
        raise StaleReferenceError(
            f"{source} was not computed on the current groups.csv: {worst_col} "
            f"differs by up to {worst:.4f}, and these are counts over the shared "
            f"cells, so they agree to {tol} when both sides test the same "
            f"cells.\nRe-run `python tutorials/pbmc3k_de_tutorial.py` and then "
            f"`Rscript tutorials/pbmc3k_de_verify.R` before comparing.")
    return worst


def check_same_cells(py_cells: Sequence[str], r_frame, *, source: str) -> None:
    """Verify an R table covers exactly the cells the Python run wrote out.

    ``report_concordance`` reindexes the R frame onto the Python cell order,
    which fills missing barcodes with ``NaN`` rather than complaining; the
    correlations downstream then come out ``NaN`` and print as ``nan``.
    """
    r_cells = list(r_frame.index)
    missing = [c for c in py_cells if c not in set(r_cells)]
    extra = [c for c in r_cells if c not in set(py_cells)]
    if missing or extra:
        raise StaleReferenceError(
            f"{source} covers a different cell set than figures_dimreduc/"
            f"cells.txt: {len(missing)} of the {len(py_cells)} Python cells are "
            f"absent and {len(extra)} are unexpected (e.g. "
            f"{(missing or extra)[:3]}).\nRe-run "
            f"`python tutorials/pbmc3k_dimreduc_tutorial.py` and then "
            f"`Rscript tutorials/pbmc3k_dimreduc_verify.R`.")


def check_same_features(py_features: Sequence[str], r_features: Sequence[str], *,
                        source: str, key=None) -> None:
    """Verify an R feature × dim table used the current HVG selection."""
    norm = key or (lambda x: x)
    py_set = {norm(f) for f in py_features}
    r_set = {norm(f) for f in r_features}
    if py_set != r_set:
        only_py = sorted(py_set - r_set)[:3]
        only_r = sorted(r_set - py_set)[:3]
        raise StaleReferenceError(
            f"{source} used a different feature set than "
            f"figures_dimreduc/hvg_features.txt: {len(py_set - r_set)} Python "
            f"features are absent (e.g. {only_py}) and {len(r_set - py_set)} are "
            f"unexpected (e.g. {only_r}).\nRe-run the Python tutorial and then "
            f"`Rscript tutorials/pbmc3k_dimreduc_verify.R`.")


def summarise_spread(values: Iterable[float], reference: Optional[float] = None) -> dict:
    """Min/max/median of a seed sweep, and the largest gap to ``reference``.

    Used to derive the JackStraw band from measurement rather than from the one
    run that happened to be in front of whoever wrote the vignette.
    """
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        raise ValueError("no values to summarise")
    out = {"n": int(arr.size), "min": float(arr.min()), "max": float(arr.max()),
           "median": float(np.median(arr))}
    if reference is not None:
        out["max_abs_gap"] = float(np.abs(arr - float(reference)).max())
    return out
