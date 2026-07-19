"""Gene-set / signature scoring and cell-cycle scoring.

Mirrors Seurat's AddModuleScore() and CellCycleScoring().
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .lazy import is_lazy


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assay_data(seurat, assay: Optional[str], layer: str = "data"):
    """Return (matrix features×cells, feature_names) from an assay layer."""
    from .assay import Assay
    from .assay5 import Assay5
    from ._sparse import is_matrix_empty

    assay_obj = seurat.assays[assay or seurat.active_assay]
    if isinstance(assay_obj, Assay5):
        feats = assay_obj._all_feature_names
        if layer in assay_obj.layers:
            mat = assay_obj.layers[layer]
        elif "data" in assay_obj.layers:
            mat = assay_obj.layers["data"]
        else:
            mat = assay_obj.layers["counts"]
        return mat, feats
    else:
        feats = assay_obj._feature_names
        if layer == "counts":
            return assay_obj.counts, feats
        data = assay_obj.data
        return (data if not is_matrix_empty(data) else assay_obj.counts), feats


def _row(mat, idx) -> np.ndarray:
    row = mat[idx, :]
    if sp.issparse(row):
        return np.asarray(row.todense()).flatten()
    return np.asarray(row).flatten()


def _alnum(s: str) -> str:
    return "".join(ch for ch in s if ch.isalnum()).upper()


def _resolve_symbols(genes, feat_names, search: bool) -> list[str]:
    """Map requested gene symbols onto feature names actually in the assay.

    With ``search=False`` this is a plain membership filter. With ``search=True``
    it also matches case-insensitively and ignoring punctuation — a local
    stand-in for Seurat's ``UpdateSymbolList`` (no network lookup), enough to
    reconcile e.g. ``HLA-A`` vs ``HLA.A`` or lower/upper-case panels.
    """
    feat_set = set(feat_names)
    if not search:
        return [g for g in genes if g in feat_set]
    lut: dict[str, str] = {}
    for f in feat_names:
        lut.setdefault(f.upper(), f)
        lut.setdefault(_alnum(f), f)
    resolved = []
    for g in genes:
        if g in feat_set:
            resolved.append(g)
            continue
        cand = lut.get(g.upper()) or lut.get(_alnum(g))
        if cand is not None:
            resolved.append(cand)
    return resolved


# ---------------------------------------------------------------------------
# AddModuleScore  (mirrors R AddModuleScore())
# ---------------------------------------------------------------------------

def add_module_score(
    seurat,
    features: Union[Sequence[str], Sequence[Sequence[str]], dict],
    pool: Optional[list[str]] = None,
    nbin: int = 24,
    ctrl: int = 100,
    name: str = "Cluster",
    assay: Optional[str] = None,
    layer: str = "data",
    seed: int = 1,
    search: bool = False,
) -> "object":
    """Score one or more gene programs per cell.

    Mirrors R's AddModuleScore(): each program's score is the mean expression of
    its genes minus the mean expression of a control set drawn from the same
    average-expression bins (so highly/lowly expressed genes are controlled for).

    Parameters
    ----------
    features : a single gene list, a list of gene lists, or a name->list dict.
    pool     : genes to sample controls from (default: all features).
    nbin     : number of average-expression bins (Seurat default 24).
    ctrl     : control genes sampled per program gene (Seurat default 100).
    name     : metadata column prefix; programs become ``{name}1``, ``{name}2`` …
               (or the dict keys when ``features`` is a dict).
    seed     : RNG seed for control-gene sampling.
    search   : if True, resolve program genes not found verbatim by a
               case/punctuation-insensitive match (local ``UpdateSymbolList``).

    Stores one metadata column per program and returns ``seurat``.
    """
    rng = np.random.default_rng(seed)

    # Normalise `features` into (labels, list-of-lists).
    if isinstance(features, dict):
        labels = list(features.keys())
        programs = [list(v) for v in features.values()]
    elif len(features) > 0 and isinstance(features[0], (list, tuple, set)):
        programs = [list(p) for p in features]
        labels = [f"{name}{i + 1}" for i in range(len(programs))]
    else:
        programs = [list(features)]
        labels = [f"{name}1"]

    mat, feat_names = _assay_data(seurat, assay, layer)
    feat_idx = {f: i for i, f in enumerate(feat_names)}
    pool = list(pool) if pool is not None else list(feat_names)
    pool = [g for g in pool if g in feat_idx]

    # Average expression per pooled gene, then equal-frequency bins.
    pool_rows = [feat_idx[g] for g in pool]
    if sp.issparse(mat) or is_lazy(mat):
        data_avg = np.asarray(mat[pool_rows, :].mean(axis=1)).flatten()
    else:
        data_avg = np.asarray(mat)[pool_rows, :].mean(axis=1)
    # Tiny jitter breaks ties so qcut can form `nbin` equal-frequency bins.
    jitter = rng.standard_normal(len(data_avg)) / 1e30
    bins = pd.qcut(data_avg + jitter, q=min(nbin, len(pool)), labels=False, duplicates="drop")
    gene_to_bin = {g: int(b) for g, b in zip(pool, bins)}
    bin_to_genes: dict[int, list[str]] = {}
    for g, b in gene_to_bin.items():
        bin_to_genes.setdefault(b, []).append(g)

    n_cells = mat.shape[1]
    for label, genes in zip(labels, programs):
        used = _resolve_symbols(genes, feat_names, search)
        if not used:
            seurat.meta_data[label] = np.zeros(n_cells)
            continue

        # Control gene set: per program-gene, sample `ctrl` from its bin.
        ctrl_genes: set[str] = set()
        for g in used:
            b = gene_to_bin.get(g)
            if b is None:
                continue
            candidates = bin_to_genes.get(b, [])
            if not candidates:
                continue
            size = min(ctrl, len(candidates))
            picked = rng.choice(candidates, size=size, replace=False)
            ctrl_genes.update(picked.tolist())

        feat_scores = np.mean([_row(mat, feat_idx[g]) for g in used], axis=0)
        if ctrl_genes:
            ctrl_scores = np.mean([_row(mat, feat_idx[g]) for g in ctrl_genes], axis=0)
        else:
            ctrl_scores = np.zeros(n_cells)
        seurat.meta_data[label] = feat_scores - ctrl_scores

    return seurat


# ---------------------------------------------------------------------------
# Cell-cycle scoring  (mirrors R CellCycleScoring())
# ---------------------------------------------------------------------------

# Tirosh et al. 2016 cell-cycle gene sets (human), as shipped in Seurat's
# `cc.genes.updated.2019`.
CC_GENES = {
    "s_genes": [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2",
        "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
        "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
        "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
        "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
        "BRIP1", "E2F8",
    ],
    "g2m_genes": [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
        "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
        "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
        "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
        "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
        "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
        "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ],
}


def cell_cycle_scoring(
    seurat,
    s_features: Optional[list[str]] = None,
    g2m_features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    layer: str = "data",
    set_ident: bool = False,
    nbin: int = 24,
    ctrl: int = 100,
    seed: int = 1,
) -> "object":
    """Score S and G2/M phases and assign a discrete phase per cell.

    Mirrors R's CellCycleScoring(): runs AddModuleScore for the S and G2/M gene
    sets, writes ``S.Score`` / ``G2M.Score`` to metadata, and assigns ``Phase``:
    ``G1`` when both scores are ≤ 0, otherwise whichever of S / G2M is larger.
    Defaults to the Tirosh 2016 human gene sets (``CC_GENES``).

    If ``set_ident`` is True, the active identity is set to ``Phase``.
    """
    s_features = s_features if s_features is not None else CC_GENES["s_genes"]
    g2m_features = g2m_features if g2m_features is not None else CC_GENES["g2m_genes"]

    add_module_score(
        seurat,
        features={"S.Score": s_features, "G2M.Score": g2m_features},
        nbin=nbin, ctrl=ctrl, assay=assay, layer=layer, seed=seed,
    )

    s = seurat.meta_data["S.Score"].values.astype(float)
    g2m = seurat.meta_data["G2M.Score"].values.astype(float)

    phase = np.empty(len(s), dtype=object)
    for i in range(len(s)):
        if s[i] <= 0 and g2m[i] <= 0:
            phase[i] = "G1"
        elif s[i] > g2m[i]:
            phase[i] = "S"
        else:
            phase[i] = "G2M"
    seurat.meta_data["Phase"] = phase

    if set_ident:
        seurat.idents = list(phase)
    return seurat
