"""Differential expression / marker gene detection.

Mirrors Seurat's FindMarkers() and FindAllMarkers().
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import wilcoxon, ranksums


def find_markers(
    seurat,
    ident_1: Union[str, list[str]],
    ident_2: Optional[Union[str, list[str]]] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    test_use: str = "wilcox",
    only_pos: bool = False,
    min_pct: float = 0.1,
    min_pct_2: float = 0.0,
    logfc_threshold: float = 0.25,
    features: Optional[list[str]] = None,
    max_cells_per_ident: Optional[int] = None,
    random_seed: int = 1,
) -> pd.DataFrame:
    """Find differentially expressed marker genes.

    Mirrors R's FindMarkers(pbmc, ident.1 = 2).

    Parameters
    ----------
    ident_1         : cluster label(s) for group 1
    ident_2         : cluster label(s) for group 2 (None = all others)
    test_use        : statistical test: 'wilcox', 't', 'bimod', 'roc'
    only_pos        : only return positive markers
    min_pct         : minimum fraction cells expressing gene in either group
    logfc_threshold : minimum log2 fold-change filter
    features        : restrict to these genes (default: all)
    max_cells_per_ident : downsample each group to this many cells

    Returns
    -------
    DataFrame with columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj
    Sorted by p_val ascending.
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]
    cells = seurat.cell_names()
    idents = list(seurat.idents)

    # Resolve ident strings
    ident_1_set = {str(ident_1)} if isinstance(ident_1, str) else {str(i) for i in ident_1}
    if ident_2 is None:
        ident_2_set = {str(i) for i in set(idents) if str(i) not in ident_1_set}
    else:
        ident_2_set = {str(ident_2)} if isinstance(ident_2, str) else {str(i) for i in ident_2}

    # Cell indices for each group
    cells_1 = [c for c, i in zip(cells, idents) if str(i) in ident_1_set]
    cells_2 = [c for c, i in zip(cells, idents) if str(i) in ident_2_set]

    if not cells_1:
        raise ValueError(f"No cells found with ident {ident_1}.")
    if not cells_2:
        raise ValueError(f"No cells found for comparison group.")

    # Optional downsampling
    if max_cells_per_ident is not None:
        rng = np.random.default_rng(random_seed)
        if len(cells_1) > max_cells_per_ident:
            cells_1 = list(rng.choice(cells_1, max_cells_per_ident, replace=False))
        if len(cells_2) > max_cells_per_ident:
            cells_2 = list(rng.choice(cells_2, max_cells_per_ident, replace=False))

    # Get expression matrix for all genes (features × cells)
    data, feature_names = _get_expression_matrix(assay_obj, layer)

    cell_idx_map = {c: i for i, c in enumerate(cells)}
    idx_1 = [cell_idx_map[c] for c in cells_1]
    idx_2 = [cell_idx_map[c] for c in cells_2]

    if sp.issparse(data):
        mat1 = data[:, idx_1].toarray().astype(float)  # (features × n1)
        mat2 = data[:, idx_2].toarray().astype(float)  # (features × n2)
    else:
        mat1 = np.asarray(data)[:, idx_1].astype(float)
        mat2 = np.asarray(data)[:, idx_2].astype(float)

    # Restrict features
    if features is not None:
        feat_set = set(features)
        feat_mask = np.array([f in feat_set for f in feature_names])
    else:
        feat_mask = np.ones(len(feature_names), dtype=bool)

    # Percent cells expressing (> 0)
    pct1 = (mat1 > 0).mean(axis=1)
    pct2 = (mat2 > 0).mean(axis=1)

    # Pre-filter: gene must be expressed in at least min_pct of either group
    pct_mask = (pct1 >= min_pct) | (pct2 >= min_pct)
    combined_mask = feat_mask & pct_mask

    # Log2 fold change (log-normalized data is already log-scale, so exponentiate first)
    mean1 = mat1.mean(axis=1)
    mean2 = mat2.mean(axis=1)
    # Since data is log1p(CPM), exponentiate to get CPM then compute fold change
    fc_mask = np.zeros(len(feature_names))
    avg_log2fc = np.log2(np.expm1(mean1) + 1) - np.log2(np.expm1(mean2) + 1)

    # Pre-filter by logfc_threshold
    if logfc_threshold > 0:
        fc_mask_arr = np.abs(avg_log2fc) >= logfc_threshold
        combined_mask = combined_mask & fc_mask_arr

    test_indices = np.where(combined_mask)[0]
    if len(test_indices) == 0:
        return pd.DataFrame(
            columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
        )

    # Statistical tests
    p_vals = np.ones(len(test_indices))

    if test_use == "wilcox":
        for i, fi in enumerate(test_indices):
            x1 = mat1[fi, :]
            x2 = mat2[fi, :]
            if x1.sum() == 0 and x2.sum() == 0:
                p_vals[i] = 1.0
            else:
                _, p = ranksums(x1, x2)
                p_vals[i] = p if not np.isnan(p) else 1.0
    elif test_use == "t":
        from scipy.stats import ttest_ind
        for i, fi in enumerate(test_indices):
            x1 = mat1[fi, :]
            x2 = mat2[fi, :]
            _, p = ttest_ind(x1, x2, equal_var=False)
            p_vals[i] = p if not np.isnan(p) else 1.0
    else:
        raise ValueError(f"Unsupported test_use: {test_use!r}. Use 'wilcox' or 't'.")

    # Bonferroni correction (Seurat default: multiply by total gene count)
    n_total = len(feature_names)
    p_val_adj = np.minimum(p_vals * n_total, 1.0)

    results = pd.DataFrame(
        {
            "p_val": p_vals,
            "avg_log2FC": avg_log2fc[test_indices],
            "pct.1": pct1[test_indices],
            "pct.2": pct2[test_indices],
            "p_val_adj": p_val_adj,
        },
        index=[feature_names[i] for i in test_indices],
    )

    if only_pos:
        results = results[results["avg_log2FC"] > 0]

    return results.sort_values("p_val")


def find_all_markers(
    seurat,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    test_use: str = "wilcox",
    only_pos: bool = False,
    min_pct: float = 0.1,
    logfc_threshold: float = 0.25,
    max_cells_per_ident: Optional[int] = None,
    random_seed: int = 1,
) -> pd.DataFrame:
    """Find marker genes for each cluster vs all others.

    Mirrors R's FindAllMarkers(pbmc, only.pos = TRUE).

    Returns a single DataFrame with an extra 'cluster' column.
    """
    clusters = sorted(set(str(i) for i in seurat.idents))
    all_results = []

    for cluster in clusters:
        try:
            df = find_markers(
                seurat,
                ident_1=cluster,
                ident_2=None,
                assay=assay,
                layer=layer,
                test_use=test_use,
                only_pos=only_pos,
                min_pct=min_pct,
                logfc_threshold=logfc_threshold,
                max_cells_per_ident=max_cells_per_ident,
                random_seed=random_seed,
            )
            if len(df) > 0:
                df = df.copy()
                df["cluster"] = cluster
                df["gene"] = df.index
                all_results.append(df)
        except Exception as e:
            print(f"Warning: cluster {cluster} marker finding failed: {e}")

    if not all_results:
        return pd.DataFrame(
            columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]
        )

    combined = pd.concat(all_results, axis=0)
    combined = combined[["cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]]
    return combined.sort_values(["cluster", "p_val"])


# ------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------

def _get_expression_matrix(assay_obj, layer: Optional[str]):
    """Return (matrix features×cells, feature_names) using best available layer."""
    from .assay5 import Assay5
    from .assay import Assay

    if isinstance(assay_obj, Assay5):
        feature_names = assay_obj._all_feature_names
        if layer is not None and layer in assay_obj.layers:
            return assay_obj.layers[layer], feature_names
        for candidate in ("data", "counts"):
            if candidate in assay_obj.layers:
                return assay_obj.layers[candidate], feature_names
        raise ValueError("No expression data layer found in Assay5.")
    else:
        feature_names = assay_obj._feature_names
        from ._sparse import is_matrix_empty
        if layer == "counts":
            return assay_obj.counts, feature_names
        if layer == "scale_data" or layer == "scale.data":
            return assay_obj.scale_data, feature_names
        # Prefer log-normalized data
        if not is_matrix_empty(assay_obj.data):
            return assay_obj.data, feature_names
        return assay_obj.counts, feature_names
