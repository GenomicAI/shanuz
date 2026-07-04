"""Composition / differential-abundance testing across conditions.

Answers "is cluster/cell-type X over-represented in condition A vs B?" — the
enrichment step every spatial/atlas Seurat analysis ends with. Reports a
*directional* log2 proportion ratio (unambiguous), a per-group Fisher exact
test (group vs rest), BH-adjusted q-values, and the overall chi-square.

Depends only on scipy (a core dependency); BH is computed inline so no
statsmodels requirement.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def composition_test(
    seurat,
    group_by: str,
    split_by: str,
    reference: Optional[str] = None,
) -> pd.DataFrame:
    """Directional abundance test of ``group_by`` categories across ``split_by``.

    Mirrors the enrichment table an analyst builds by hand: for a two-level
    ``split_by`` (e.g. condition), each ``group_by`` category (e.g. cluster) gets
    a ``log2(prop_test / prop_reference)`` and a Fisher exact p (that category vs
    all others). p-values are BH-adjusted. The overall chi-square p is stored in
    ``df.attrs['chisq_p']``.

    Parameters
    ----------
    group_by  : categorical metadata column tested for enrichment (rows).
    split_by  : metadata column with exactly two levels (the conditions).
    reference : which ``split_by`` level is the denominator (default: the first
                sorted level). log2 > 0 ⇒ enriched in the *other* level.

    Returns a DataFrame ordered by ``log2_ratio`` with columns:
    ``group, n_<ref>, n_<test>, prop_<ref>, prop_<test>, log2_ratio,
    odds_ratio, p, padj, sig, enriched_in``.
    """
    md = seurat.meta_data
    for col in (group_by, split_by):
        if col not in md.columns:
            raise KeyError(f"'{col}' not in meta_data.")
    tab = pd.crosstab(md[group_by].astype(str), md[split_by].astype(str))
    conds = list(tab.columns)
    if len(conds) != 2:
        raise ValueError(
            f"split_by='{split_by}' must have exactly 2 levels, found {conds}."
        )
    ref = reference if reference is not None else conds[0]
    if ref not in conds:
        raise ValueError(f"reference '{ref}' not a level of {split_by}: {conds}.")
    test = [c for c in conds if c != ref][0]

    n_ref, n_test = tab[ref].sum(), tab[test].sum()
    rows = []
    for grp in tab.index:
        a, b = tab.loc[grp, test], tab.loc[grp, ref]          # this group
        c, d = n_test - a, n_ref - b                          # all other groups
        odds, p = stats.fisher_exact([[a, b], [c, d]])
        prop_test = a / n_test if n_test else np.nan
        prop_ref = b / n_ref if n_ref else np.nan
        with np.errstate(divide="ignore"):
            log2 = np.log2(prop_test / prop_ref) if prop_ref else np.nan
        rows.append({
            "group": grp, f"n_{ref}": b, f"n_{test}": a,
            f"prop_{ref}": prop_ref, f"prop_{test}": prop_test,
            "log2_ratio": log2, "odds_ratio": odds, "p": p,
        })
    df = pd.DataFrame(rows)
    df["padj"] = _bh(df["p"].to_numpy())
    bins = [-np.inf, 0.001, 0.01, 0.05, np.inf]
    df["sig"] = pd.cut(df["padj"], bins, labels=["***", "**", "*", "ns"])
    df["enriched_in"] = np.where(df["log2_ratio"] > 0, test, ref)
    df = df.sort_values("log2_ratio", ascending=False).reset_index(drop=True)
    df.attrs["chisq_p"] = float(stats.chi2_contingency(tab.to_numpy())[1])
    df.attrs["reference"] = ref
    df.attrs["test"] = test
    return df
