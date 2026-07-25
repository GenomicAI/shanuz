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
    orig_reduction : reduction to integrate (default 'pca'). Every method
                     corrects this reduction and writes a new one of the same
                     shape — the anchor methods included, which is what makes
                     them interchangeable with Harmony here.
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
            orig_reduction=orig_reduction,
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
    orig_reduction: str = "pca",
    assay: Optional[str] = None,
    seed: int = 42,
    **kwargs,
) -> None:
    """CCA/RPCA layer integration → a corrected reduction (Seurat v5 path).

    Mirrors ``CCAIntegration`` / ``RPCAIntegration``: split the object by
    ``group_by``, find anchors between the batches, and hand them to
    :func:`shanuz.anchors.integrate_embeddings`, which corrects
    ``orig_reduction`` in place of the expression matrix.

    This is **not** the v4 ``IntegrateData`` workflow, and the difference is not
    cosmetic. v4 corrects expression and then re-scales and re-runs PCA, which
    lands you in a *new* basis; v5 corrects the existing embedding, so the
    output stays in ``orig_reduction``'s basis and keeps its loadings. Running
    the v4 workflow behind this API produced an embedding that agreed with
    Seurat's on 1 of 30 dimensions — not because the correction was wrong, but
    because it was a different object with the same shape.

    Two further details are carried over from the R:

    * ``k.filter`` is ``NA`` for **both** methods here, so the expression-space
      anchor filter never runs (v4 applies it to CCA);
    * ``RPCAIntegration`` runs ``ScaleData`` on each batch, while
      ``CCAIntegration`` slices the object's *existing* ``scale.data``. Scaling
      per batch centres each on its own mean, which is what reciprocal PCA
      needs and what CCA is deliberately not given.
    """
    from .anchors import find_integration_anchors, integrate_embeddings
    from .preprocessing import scale_data

    if isinstance(group_by, (list, tuple)):
        if len(group_by) != 1:
            raise ValueError(
                "CCA/RPCA integration supports a single group_by column."
            )
        group_by = group_by[0]
    if group_by not in seurat.meta_data.columns:
        raise KeyError(f"group_by column {group_by!r} not found in meta_data.")

    if orig_reduction not in seurat.reductions:
        raise KeyError(
            f"orig_reduction {orig_reduction!r} not found. "
            "IntegrateLayers corrects an existing reduction — run run_pca() first."
        )

    k_weight = kwargs.pop("k_weight", 100)
    sd_weight = kwargs.pop("sd_weight", 1.0)
    dims_to_integrate = kwargs.pop("dims_to_integrate", None)
    # CCAIntegration and RPCAIntegration both call FindIntegrationAnchors with
    # k.filter = NA. v4's default of 200 only applies to the object-list API.
    kwargs.setdefault("k_filter", None)

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

    # RPCAIntegration re-scales each batch (ScaleData per object → per-object
    # PCA), so every batch is centred on its own mean before the reciprocal
    # spaces are built; global scaling leaves the batch mean-shift in PC1 and
    # RPCA then under-finds mutual pairs. CCAIntegration does the opposite — it
    # assigns each batch a *slice* of the object's existing scale.data — and
    # subsetting already carried that slice through, so CCA is left alone.
    if reduction == "rpca":
        for obj in objects:
            scale_data(obj, assay=assay)

    # Seurat corrects the SMALLER dataset onto the larger one:
    # PairwiseIntegrateReference reverses the merge pair whenever the second
    # object has more cells, so the reference is whichever batch is biggest.
    # Taking the first batch instead is invisible on equal splits and pulls the
    # wrong way on real ones — ifnb is CTRL 6,548 vs STIM 7,451.
    reference = kwargs.pop(
        "reference", int(np.argmax([len(o.cell_names()) for o in objects]))
    )
    anchors = find_integration_anchors(
        objects, reduction=reduction, reference=reference, seed=seed, **kwargs
    )
    corrected = integrate_embeddings(
        anchors,
        seurat.reductions[orig_reduction],
        new_reduction=new_reduction,
        dims_to_integrate=dims_to_integrate,
        k_weight=k_weight,
        sd_weight=sd_weight,
    )

    # Reindex the corrected embedding back to the original cell order.
    seurat.reductions[new_reduction] = corrected.subset(cells=all_cells)
