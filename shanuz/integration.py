"""Batch correction & integration.

Mirrors Seurat's RunHarmony() (via the harmony R package) and the Seurat v5
IntegrateLayers() dispatch API.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd

from .dimreduc import DimReduc


def run_harmony(
    seurat,
    group_by: Union[str, list[str]],
    reduction: str = "pca",
    dims: Optional[Union[list[int], range]] = None,
    reduction_name: str = "harmony",
    reduction_key: str = "harmony_",
    theta: Optional[Union[float, list[float]]] = None,
    lambda_: Optional[Union[float, list[float]]] = None,
    sigma: float = 0.1,
    nclust: Optional[int] = None,
    max_iter_harmony: int = 10,
    assay: Optional[str] = None,
    seed: int = 0,
) -> None:
    """Run Harmony batch correction on an existing reduction.

    Mirrors R's ``RunHarmony(obj, group.by.vars = "batch")``. Takes the cell
    embeddings of ``reduction`` (PCA by default), removes batch effects with
    ``harmonypy``, and stores the corrected embeddings as a new DimReduc under
    ``reduction_name`` — same shape as the input, so it can be passed straight
    to ``find_neighbors(reduction="harmony")`` / ``run_umap(reduction="harmony")``.

    Parameters
    ----------
    group_by         : metadata column(s) identifying the batch(es) to correct
    reduction        : source reduction to correct (default 'pca')
    dims             : which dimensions of ``reduction`` to use (0-indexed;
                       default all available)
    reduction_name   : storage key for the corrected reduction
    theta            : diversity clustering penalty (harmonypy default when None)
    lambda_          : ridge regression penalty (harmonypy default when None)
    sigma            : soft-clustering width
    nclust           : number of Harmony clusters (harmonypy default when None)
    max_iter_harmony : maximum Harmony iterations
    seed             : random seed for reproducibility
    """
    try:
        import harmonypy
    except ImportError as exc:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "run_harmony requires 'harmonypy'. Install it with "
            "`pip install shanuz[integration]` or `pip install harmonypy`."
        ) from exc

    assay_name = assay or seurat.active_assay

    if reduction not in seurat.reductions:
        raise KeyError(
            f"Reduction '{reduction}' not found. Run run_pca() first."
        )
    dr = seurat.reductions[reduction]
    embeddings = dr.cell_embeddings  # (cells × dims)

    if dims is not None:
        embeddings = embeddings[:, list(dims)]

    group_vars = [group_by] if isinstance(group_by, str) else list(group_by)
    missing = [g for g in group_vars if g not in seurat.meta_data.columns]
    if missing:
        raise KeyError(
            f"group_by column(s) {missing} not found in meta_data."
        )
    meta = seurat.meta_data[group_vars]

    n_cells = embeddings.shape[0]

    np.random.seed(seed)
    harmony_obj = harmonypy.run_harmony(
        embeddings,
        meta,
        group_vars,
        theta=theta,
        lamb=lambda_,
        sigma=sigma,
        nclust=nclust,
        max_iter_harmony=max_iter_harmony,
        random_state=seed,
    )
    # harmonypy stores corrected embeddings as (dims × cells); orient robustly
    # to (cells × dims) by matching the known cell count.
    corrected = np.asarray(harmony_obj.Z_corr)
    if corrected.shape[0] != n_cells and corrected.shape[1] == n_cells:
        corrected = corrected.T

    cells = seurat.cell_names()
    dim_names = [f"{reduction_key}{i + 1}" for i in range(corrected.shape[1])]

    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=corrected,
        assay_used=assay_name,
        key=reduction_key,
        cell_names=cells,
        feature_names=dim_names,
    )


def integrate_layers(
    seurat,
    method: str = "harmony",
    orig_reduction: str = "pca",
    new_reduction: Optional[str] = None,
    group_by: Optional[Union[str, list[str]]] = None,
    assay: Optional[str] = None,
    **kwargs,
) -> None:
    """Integrate layers/batches (Seurat v5 ``IntegrateLayers`` dispatch API).

    Mirrors ``IntegrateLayers(obj, method = HarmonyIntegration,
    orig.reduction = "pca")``. A thin dispatcher over the individual
    integration routines.

    Parameters
    ----------
    method         : 'harmony', 'cca', or 'rpca'.
    orig_reduction : reduction to integrate (only used by 'harmony';
                     default 'pca')
    new_reduction  : storage key for the integrated reduction
                     (defaults to '{method}')
    group_by       : batch column identifying the layers/batches to integrate
                     (required for every method)
    """
    method = method.lower()
    new_reduction = new_reduction or method

    if method in ("harmony", "harmonyintegration"):
        if group_by is None:
            raise ValueError("method='harmony' requires group_by (batch column).")
        run_harmony(
            seurat,
            group_by=group_by,
            reduction=orig_reduction,
            reduction_name=new_reduction,
            assay=assay,
            **kwargs,
        )
    elif method in ("cca", "rpca", "ccaintegration", "rpcaintegration"):
        if group_by is None:
            raise ValueError(
                f"method={method!r} requires group_by (batch column)."
            )
        reduction = "rpca" if method.startswith("rpca") else "cca"
        _integrate_anchor_reduction(
            seurat,
            group_by=group_by,
            reduction=reduction,
            new_reduction=new_reduction,
            assay=assay,
            **kwargs,
        )
    else:
        raise ValueError(
            f"Unknown integration method {method!r}. "
            "Supported: 'harmony', 'cca', 'rpca'."
        )


def _integrate_anchor_reduction(
    seurat,
    group_by: str,
    reduction: str,
    new_reduction: str,
    assay: Optional[str] = None,
    n_pcs: int = 30,
    seed: int = 42,
    **kwargs,
) -> None:
    """CCA/RPCA layer integration → a corrected reduction (Seurat v5 path).

    Splits the object by ``group_by`` into per-batch datasets, finds anchors
    between them (:func:`shanuz.anchors.find_integration_anchors`), builds the
    batch-corrected ``"integrated"`` assay (:func:`shanuz.anchors.integrate_data`),
    then runs ``scale_data`` + ``run_pca`` on it and stores the resulting
    embedding under ``new_reduction`` (reindexed to the object's cell order).
    """
    from .anchors import find_integration_anchors, integrate_data
    from .preprocessing import scale_data
    from .reduction import run_pca

    if isinstance(group_by, (list, tuple)):
        if len(group_by) != 1:
            raise ValueError(
                "CCA/RPCA integration supports a single group_by column."
            )
        group_by = group_by[0]
    if group_by not in seurat.meta_data.columns:
        raise KeyError(f"group_by column {group_by!r} not found in meta_data.")

    k_weight = kwargs.pop("k_weight", 100)
    sd_weight = kwargs.pop("sd_weight", 1.0)

    groups = list(pd.unique(seurat.meta_data[group_by]))
    if len(groups) < 2:
        raise ValueError(
            f"group_by={group_by!r} has < 2 levels; nothing to integrate."
        )

    all_cells = seurat.cell_names()
    labels = seurat.meta_data[group_by]
    objects = [
        seurat.subset(cells=[c for c in all_cells if labels.loc[c] == g])
        for g in groups
    ]

    anchors = find_integration_anchors(
        objects, reduction=reduction, seed=seed, **kwargs
    )
    merged = integrate_data(anchors, k_weight=k_weight, sd_weight=sd_weight)

    scale_data(merged, assay="integrated")
    run_pca(
        merged,
        n_pcs=n_pcs,
        assay="integrated",
        reduction_name=new_reduction,
        seed=seed,
    )

    # Reindex the integrated embedding back to the original cell order.
    seurat.reductions[new_reduction] = merged.reductions[new_reduction].subset(
        cells=all_cells
    )
