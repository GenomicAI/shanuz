"""Pooled CRISPR-screen analysis — Seurat's ``CalcPerturbSig`` + ``RunMixscape``.

In a pooled CRISPR screen (Perturb-seq / CROP-seq / ECCITE-seq) every cell
receives one guide RNA targeting one gene, and the whole pool is sequenced
together. The promise is a knockout phenotype for hundreds of genes in a single
experiment; the catch is that **carrying a guide is not the same as being
perturbed**. A cell can pick up a guide and still escape the knockout — the edit
fails, one allele survives, the protein lingers — so the cells labelled
``"IFNGR2"`` are a *mixture* of true knockouts (KO) and non-perturbed escapers
(NP) that look just like controls. Averaging over that mixture dilutes, and can
entirely mask, the real phenotype. Mixscape (Papalexi, Mimitou et al. 2021) is
the method that separates the two, so downstream analysis runs on genuinely
perturbed cells. It comes in two steps, both ported here, plus a visualization.

**Step 1 — :func:`calc_perturb_sig` (Seurat's ``CalcPerturbSig``).** Guide
assignment is confounded by everything else that varies between cells: cell
cycle, sequencing depth, replicate, ambient RNA. To strip that away, each cell's
expression has subtracted from it the *average of its nearest non-targeting (NT)
control cells* — its ``num_neighbors`` (default 20) nearest neighbours **among the
NT cells**, in a reduction (default the first ``ndims`` PCs), optionally computed
within each ``split_by`` batch. What remains — the **local perturbation
signature** — is the deviation of this cell from the controls it most resembles,
with the shared technical variation cancelled out. The signatures are stored as a
new assay (default ``"PRTB"``) of the same genes × cells shape.

**Step 2 — :func:`run_mixscape` (Seurat's ``RunMixscape``).** Working on the
perturbation signature, each target gene is handled independently:

1. **Differential expression.** The gene's cells are tested against the NT cells
   (on the original ``de_assay``, default ``"RNA"``) with :func:`shanuz.find_markers`;
   the genes that pass (``|avg_log2FC| > logfc_threshold`` and
   ``p_val_adj < pval_cutoff``) are the axes along which this perturbation shows
   up. A gene with fewer than ``min_de_genes`` such axes has no detectable
   phenotype and all its cells are called NP.
2. **Iterative mixture classification.** Restricted to those DE genes, a
   *perturbation vector* is formed — the mean signature of the currently-KO cells
   minus the mean signature of the NT cells — and every cell is projected onto it
   to a single **perturbation score**. A two-component Gaussian mixture
   (``sklearn.mixture.GaussianMixture``) fit to those scores splits the low
   (non-perturbed, anchored by the NT cells) mode from the high (knockout) mode;
   a gene cell whose posterior for the high component exceeds 0.5 is called KO.
   Because the perturbation vector depends on which cells are currently KO, the
   fit is iterated — recompute the vector from the new KO set, re-project, re-fit
   — until the KO set stops changing or ``iter_num`` rounds elapse. This is the
   EM-style refinement that lets the true KO cells define their own axis.

Results are written to ``obj.meta_data`` under Seurat's names: ``mixscape_class``
(``"<gene> KO"`` / ``"<gene> NP"`` / ``"NT"``), which is also set as the active
identity; ``mixscape_class.global`` (``"KO"`` / ``"NP"`` / ``"NT"``); and a
posterior column ``mixscape_class_p_<type>`` (default ``mixscape_class_p_ko``).
Per-gene bookkeeping (DE-gene count, iterations, KO count) is stashed in
``obj.misc["mixscape"]``, along with a ``scores`` frame per gene holding that
gene's perturbation score for its own cells and the NT cells — R keeps the same
thing in the ``Tool(object, "RunMixscape")`` slot, and
:func:`~shanuz.plotting.plot_perturb_score` reads it back to draw the density.

**Visualization — :func:`mixscape_lda` (Seurat's ``MixscapeLDA``).** Where
:func:`run_mixscape` asks *which cells are perturbed*, this asks *how do the guide
populations differ from each other and from control*, and answers with one
supervised map on which each guide class forms its own cloud. It reads only the
perturbation signature and the raw guide labels — the KO/NP calls are not used —
so it can follow :func:`calc_perturb_sig` directly. See its docstring for the
per-guide PCA blocks and the discriminant fit over them.

Two deliberate, documented choices differ from a literal reading of R:

* **The mixture is fit on the pooled NT + gene-cell scores**, not the gene cells
  alone, so the non-perturbed component is anchored by the controls — which is the
  entire reason NT cells are carried into the fit. A gene whose cells are almost
  all genuinely perturbed still gets a well-defined NP mode from the NT scores.
* **The perturbation signature is read straight from the assay's ``data`` layer**;
  Seurat's tutorial first runs ``ScaleData`` (centre only) on the ``PRTB`` assay.
  Centring each gene shifts the perturbation vector's two group means by the *same* per-gene
  constant, leaving the vector unchanged, and shifts every projected score by one
  global constant — which the re-fit mixture absorbs. So the KO/NP calls are
  invariant to that centring and it is skipped.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ----------------------------------------------------------------------
# Public API — perturbation signature
# ----------------------------------------------------------------------


def calc_perturb_sig(
    seurat,
    assay: str = "RNA",
    features: Optional[Sequence[str]] = None,
    layer: str = "data",
    labels: str = "gene",
    nt_class: str = "NT",
    split_by: Optional[str] = None,
    num_neighbors: int = 20,
    reduction: str = "pca",
    ndims: int = 15,
    new_assay: str = "PRTB",
):
    """Compute each cell's local perturbation signature (Seurat's ``CalcPerturbSig``).

    Mirrors ``CalcPerturbSig(object, assay, gd.class, nt.cell.class, reduction,
    ndims, num.neighbors, new.assay.name = "PRTB")``. For every cell, the mean
    expression of its ``num_neighbors`` nearest non-targeting (NT) control cells
    — in the first ``ndims`` dimensions of ``reduction`` — is subtracted from its
    own expression. The residual (the deviation from the controls the cell most
    resembles) is stored as a new assay, ready for :func:`run_mixscape`.

    Parameters
    ----------
    seurat        : a :class:`~shanuz.Shanuz` object with a guide-assignment
                    metadata column and a computed ``reduction``.
    assay         : source expression assay (default ``"RNA"``).
    features      : genes to include (default: the assay's variable features, or
                    all features if none are set).
    layer         : expression layer to difference (default ``"data"``, the
                    log-normalized values).
    labels        : metadata column holding each cell's target-gene / guide class.
    nt_class      : the value in ``labels`` marking non-targeting control cells.
    split_by      : optional metadata column; neighbours are found only within the
                    same group (e.g. replicate), so batch is not mistaken for signal.
    num_neighbors : NT neighbours averaged per cell (Seurat default 20); capped at
                    the number of NT cells available in the group.
    reduction     : reduction whose embedding defines "nearest" (default ``"pca"``).
    ndims         : leading dimensions of ``reduction`` to use (default 15).
    new_assay     : name of the perturbation-signature assay to create (default
                    ``"PRTB"``).

    Returns
    -------
    Shanuz
        ``seurat``, with the perturbation-signature assay ``new_assay`` attached.
    """
    from .assay5 import create_assay5_object

    emb = seurat.embeddings(reduction)
    if ndims is not None:
        emb = emb[:, :ndims]
    emb = np.asarray(emb, dtype=float)

    data, feats, cells = _layer_matrix(seurat.assays[assay], layer)
    if features is not None:
        keep = [f for f in features if f in set(feats)]
        idx = [feats.index(f) for f in keep]
        data = data[idx, :]
        feats = keep
    if not feats:
        raise ValueError("No features selected for the perturbation signature.")

    labels_vec = _aligned_meta(seurat, labels, cells)
    nt_mask = labels_vec == nt_class
    if not nt_mask.any():
        raise ValueError(
            f"No non-targeting cells: column {labels!r} has no value {nt_class!r}."
        )

    groups = _split_groups(seurat, split_by, cells)

    signature = np.array(data, dtype=float, copy=True)
    for gidx in groups:
        nt_local = gidx[nt_mask[gidx]]
        if nt_local.size == 0:
            # No controls in this split — fall back to the whole-dataset NT mean.
            control = data[:, nt_mask].mean(axis=1, keepdims=True)
            signature[:, gidx] = data[:, gidx] - control
            continue
        k = int(min(num_neighbors, nt_local.size))
        neigh = _nt_neighbor_means(emb, gidx, nt_local, data, k)
        signature[:, gidx] = data[:, gidx] - neigh

    prtb = create_assay5_object(
        data=sp.csc_matrix(signature),
        feature_names=list(feats),
        cell_names=list(cells),
        key=f"{new_assay.lower()}_",
    )
    prtb.variable_features = list(feats)
    seurat.assays[new_assay] = prtb
    seurat.misc.setdefault("calc_perturb_sig", {})[new_assay] = {
        "assay": assay,
        "reduction": reduction,
        "ndims": ndims,
        "num_neighbors": num_neighbors,
        "n_features": len(feats),
    }
    return seurat


# ----------------------------------------------------------------------
# Public API — mixscape classification
# ----------------------------------------------------------------------


def run_mixscape(
    seurat,
    assay: str = "PRTB",
    labels: str = "gene",
    nt_class: str = "NT",
    de_assay: str = "RNA",
    layer: str = "data",
    min_de_genes: int = 5,
    min_cells: int = 5,
    logfc_threshold: float = 0.25,
    min_pct: float = 0.05,
    pval_cutoff: float = 0.05,
    iter_num: int = 10,
    prtb_type: str = "KO",
    new_class: str = "mixscape_class",
    de_test: str = "wilcox",
    seed: int = 0,
    verbose: bool = False,
):
    """Classify perturbed vs. escaping cells per guide (Seurat's ``RunMixscape``).

    Mirrors ``RunMixscape(object, assay = "PRTB", labels = "gene",
    nt.class.name = "NT", de.assay = "RNA", min.de.genes = 5, iter.num = 10,
    prtb.type = "KO")``. Operating on the perturbation signature from
    :func:`calc_perturb_sig`, each target gene's cells are split into knockout
    (``KO``) and non-perturbed (``NP``) by an iterative two-component Gaussian
    mixture over their projection onto the gene's perturbation vector (see the
    module docstring). NT cells stay ``NT``.

    The object is mutated in place: ``mixscape_class`` (also set as the active
    identity), ``mixscape_class.global``, and ``mixscape_class_p_<type>`` metadata
    columns are written, and per-gene bookkeeping is stashed in
    ``obj.misc["mixscape"]``.

    Parameters
    ----------
    seurat          : a :class:`~shanuz.Shanuz` object carrying the ``assay``
                      perturbation signature and a ``labels`` guide column.
    assay           : perturbation-signature assay (default ``"PRTB"``).
    labels          : metadata column of per-cell target-gene / guide class.
    nt_class        : value in ``labels`` marking non-targeting controls.
    de_assay        : assay used for the gene-vs-NT differential expression
                      (default ``"RNA"``).
    layer           : signature layer to project (default ``"data"``).
    min_de_genes    : a gene needs at least this many DE genes to be testable;
                      otherwise all its cells are NP (Seurat default 5).
    min_cells       : a gene needs at least this many cells; otherwise NP.
    logfc_threshold : ``|avg_log2FC|`` DE cutoff passed to ``find_markers`` (0.25).
    min_pct         : ``min.pct`` DE cutoff passed to ``find_markers`` (0.05).
    pval_cutoff     : adjusted-p cutoff a DE gene must clear (0.05).
    iter_num        : maximum mixture-refinement rounds per gene (Seurat 10).
    prtb_type       : label for the perturbed class (default ``"KO"``; use e.g.
                      ``"KD"`` for a knock-down screen). Also names the posterior
                      column ``mixscape_class_p_<type>``.
    new_class       : base name for the output columns / identity
                      (default ``"mixscape_class"``).
    de_test         : ``find_markers`` test for the gene-vs-NT DE (default
                      ``"wilcox"``).
    seed            : random state for the Gaussian mixture (determinism).
    verbose         : print each gene's DE-gene count and final KO count.

    Returns
    -------
    Shanuz
        ``seurat``, with the ``mixscape_class`` classification and identity.
    """
    sig, sig_feats, cells = _layer_matrix(seurat.assays[assay], layer)
    sig_feat_idx = {f: i for i, f in enumerate(sig_feats)}

    labels_vec = _aligned_meta(seurat, labels, cells)
    nt_idx = np.where(labels_vec == nt_class)[0]
    if nt_idx.size == 0:
        raise ValueError(
            f"No non-targeting cells: column {labels!r} has no value {nt_class!r}."
        )

    genes = sorted(
        {g for g in labels_vec.tolist() if isinstance(g, str) and g != nt_class}
    )
    if not genes:
        raise ValueError(f"No target genes to test in column {labels!r}.")

    n_cells = len(cells)
    mixscape_class = np.array(
        [str(x) for x in labels_vec], dtype=object
    )                                                # NT stays "NT"; genes filled below
    global_class = np.where(labels_vec == nt_class, "NT", "NP").astype(object)
    p_prtb = np.full(n_cells, np.nan, dtype=float)
    bookkeeping: dict[str, dict] = {}

    # find_markers reads the active identity, so drive the gene-vs-NT DE off the
    # guide labels — restored afterwards.
    saved_ident = pd.Categorical(list(seurat.idents))
    seurat.idents = list(labels_vec)
    try:
        for gene in genes:
            gene_local = np.where(labels_vec == gene)[0]
            info = {"n_cells": int(gene_local.size), "n_de": 0, "n_iter": 0, "n_ko": 0}

            if gene_local.size < min_cells:
                _assign(mixscape_class, global_class, gene_local, gene, "NP", prtb_type)
                info["scores"] = None
                bookkeeping[gene] = info
                if verbose:
                    print(f"[mixscape] {gene}: {gene_local.size} cells < min_cells → NP")
                continue

            de_genes = _de_genes(
                seurat, gene, nt_class, de_assay, de_test,
                logfc_threshold, min_pct, pval_cutoff, sig_feat_idx,
            )
            info["n_de"] = len(de_genes)
            if len(de_genes) < min_de_genes:
                _assign(mixscape_class, global_class, gene_local, gene, "NP", prtb_type)
                info["scores"] = None
                bookkeeping[gene] = info
                if verbose:
                    print(f"[mixscape] {gene}: {len(de_genes)} DE genes < min → NP")
                continue

            de_rows = [sig_feat_idx[g] for g in de_genes]
            ko_pos, post, n_iter, score = _mixscape_em(
                sig, de_rows, nt_idx, gene_local, iter_num, seed,
            )
            info["n_iter"] = n_iter
            info["n_ko"] = int(len(ko_pos))
            # R's gv data.frame: one `pvec` column plus the guide-label column,
            # indexed by cell — what PlotPerturbScore reads back out.
            if score is None:
                info["scores"] = None
            else:
                score_cells = [
                    cells[i] for i in np.concatenate([nt_idx, gene_local])
                ]
                info["scores"] = pd.DataFrame(
                    {
                        "pvec": score,
                        labels: [nt_class] * nt_idx.size + [gene] * gene_local.size,
                    },
                    index=score_cells,
                )

            ko_set = set(int(p) for p in ko_pos)
            for pos, cell in enumerate(gene_local):
                if pos in ko_set:
                    mixscape_class[cell] = f"{gene} {prtb_type}"
                    global_class[cell] = prtb_type
                else:
                    mixscape_class[cell] = f"{gene} NP"
                    global_class[cell] = "NP"
                p_prtb[cell] = post[pos]
            bookkeeping[gene] = info
            if verbose:
                print(
                    f"[mixscape] {gene}: {len(de_genes)} DE genes, "
                    f"{len(ko_set)}/{gene_local.size} {prtb_type} in {n_iter} iters"
                )
    finally:
        seurat.idents = saved_ident

    target = seurat.cell_names()

    def put(col, values):
        seurat.meta_data[col] = (
            pd.Series(list(values), index=cells).reindex(target).values
        )

    put(new_class, mixscape_class)
    put(f"{new_class}.global", global_class)
    put(f"{new_class}_p_{prtb_type.lower()}", p_prtb)
    seurat.idents = list(
        pd.Series(list(mixscape_class), index=cells).reindex(target).values
    )
    seurat.misc.setdefault("mixscape", {})[assay] = {
        "genes": bookkeeping,
        "nt_class": nt_class,
        "prtb_type": prtb_type,
    }
    return seurat


# ----------------------------------------------------------------------
# Public API — LDA visualization of perturbation classes
# ----------------------------------------------------------------------


def mixscape_lda(
    seurat,
    labels: str = "gene",
    nt_class: str = "NT",
    assay: str = "PRTB",
    de_assay: str = "RNA",
    layer: str = "data",
    npcs: int = 10,
    logfc_threshold: float = 0.25,
    min_pct: float = 0.1,
    pval_cutoff: float = 0.05,
    de_test: str = "wilcox",
    reduction_name: str = "lda",
    reduction_key: str = "LDA_",
    scale_max: float = 10.0,
    seed: int = 42,
    verbose: bool = False,
):
    """Linear-discriminant projection that separates the guide classes (Seurat's ``MixscapeLDA``).

    Mirrors ``MixscapeLDA(object, pc.assay = "PRTB", labels = "gene",
    nt.class.name = "NT", npcs = 10)``, which asks a complementary question to
    :func:`run_mixscape`: not *which cells are perturbed* but *how do the whole
    guide populations differ from one another and from control*. It builds a
    single supervised 2-D-ish map on which every guide class (and NT) forms its
    own cloud, the classic mixscape LDA plot. The only prerequisite is a
    perturbation-signature assay from :func:`calc_perturb_sig`; the mixscape KO/NP
    calls are **not** used — cells are grouped by their raw guide label.

    The construction (Seurat's ``PrepLDA`` → ``RunLDA``):

    1. **Per-guide feature blocks.** For each target gene, its cells are tested
       against NT (on ``de_assay``) to find that guide's response genes; a guide
       with fewer than ``npcs + 1`` such genes is dropped (it cannot support
       ``npcs`` components). Restricted to those genes on the perturbation-signature
       ``assay``, a PCA is fit on that guide's cells **plus** the NT cells, and then
       **every** cell in the dataset is projected onto that guide's ``npcs``-dim
       subspace. Each surviving guide thus contributes ``npcs`` columns describing
       where all cells fall along *its* perturbation axes.
    2. **One LDA over the concatenation.** The per-guide blocks are stacked side by
       side into one cell × (guides · ``npcs``) matrix and a single linear
       discriminant analysis (``sklearn`` ``LinearDiscriminantAnalysis``) is fit
       with the guide label as the class — finding the ``n_classes − 1`` directions
       that best separate the guide populations (including NT). The discriminant
       scores are stored as a reduction (default ``"lda"``, key ``"LDA_"``), and
       the per-cell class assignment and posteriors are written to metadata
       (``lda_assignments`` and ``LDAP_<class>``).

    Two choices differ from a literal reading of R, both documented:

    * **The per-guide subspace is read from the signature's ``data`` layer, scaled
      against the guide-plus-NT reference**, in place of Seurat's ``ScaleData`` →
      ``RunPCA`` → ``ProjectCellEmbeddings`` chain. The composition is the same map:
      centre/scale each response gene by the reference cells' mean and SD, project
      through the reference PCA loadings.
    * **The leave-one-out CV posterior** (MASS ``lda(..., CV = TRUE)``) that Seurat
      stashes in ``misc`` is not computed; only the resubstitution assignment and
      posterior are kept, which is all the plot and the metadata columns use.

    Parameters
    ----------
    seurat          : a :class:`~shanuz.Shanuz` object carrying the ``assay``
                      perturbation signature (from :func:`calc_perturb_sig`) and a
                      ``labels`` guide column.
    labels          : metadata column of per-cell target-gene / guide class.
    nt_class        : value in ``labels`` marking non-targeting controls.
    assay           : perturbation-signature assay projected for the LDA features
                      (default ``"PRTB"``).
    de_assay        : assay for the per-guide guide-vs-NT differential expression
                      (default ``"RNA"``).
    layer           : signature layer to project (default ``"data"``).
    npcs            : per-guide PCA components (Seurat default 10); a guide needs at
                      least ``npcs + 1`` DE genes to contribute.
    logfc_threshold : ``|avg_log2FC|`` DE cutoff passed to ``find_markers`` (0.25).
    min_pct         : ``min.pct`` DE cutoff passed to ``find_markers`` (0.1, as in
                      Seurat's ``TopDEGenesMixscape``).
    pval_cutoff     : adjusted-p cutoff a DE gene must clear (0.05).
    de_test         : ``find_markers`` test for the guide-vs-NT DE (``"wilcox"``).
    reduction_name  : key under which the LDA reduction is stored (default
                      ``"lda"``).
    reduction_key   : prefix for the discriminant dimension names (``"LDA_"``).
    scale_max       : clip scaled values to ``±scale_max`` before PCA, matching
                      Seurat's ``ScaleData`` (default 10; ``None`` to disable).
    seed            : random state for the per-guide PCA (determinism).
    verbose         : print each guide's DE-gene count and whether it contributed.

    Returns
    -------
    Shanuz
        ``seurat``, with the ``reduction_name`` LDA reduction and the
        ``lda_assignments`` / ``LDAP_<class>`` metadata columns.
    """
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    sig, sig_feats, cells = _layer_matrix(seurat.assays[assay], layer)
    sig_feat_idx = {f: i for i, f in enumerate(sig_feats)}

    labels_vec = _aligned_meta(seurat, labels, cells)
    nt_idx = np.where(labels_vec == nt_class)[0]
    if nt_idx.size == 0:
        raise ValueError(
            f"No non-targeting cells: column {labels!r} has no value {nt_class!r}."
        )

    genes = sorted(
        {g for g in labels_vec.tolist() if isinstance(g, str) and g != nt_class}
    )
    if not genes:
        raise ValueError(f"No target genes to test in column {labels!r}.")

    blocks: list[np.ndarray] = []
    feat_names: list[str] = []
    genes_used: list[str] = []

    # find_markers reads the active identity — drive the guide-vs-NT DE off the
    # guide labels, restored afterwards.
    saved_ident = pd.Categorical(list(seurat.idents))
    seurat.idents = list(labels_vec)
    try:
        for gene in genes:
            gene_local = np.where(labels_vec == gene)[0]
            de_genes = _de_genes(
                seurat, gene, nt_class, de_assay, de_test,
                logfc_threshold, min_pct, pval_cutoff, sig_feat_idx,
            )
            if len(de_genes) < npcs + 1:
                if verbose:
                    print(
                        f"[mixscape_lda] {gene}: {len(de_genes)} DE genes "
                        f"< npcs+1={npcs + 1} → skipped"
                    )
                continue
            de_rows = [sig_feat_idx[g] for g in de_genes]
            block = _lda_guide_block(
                sig, de_rows, nt_idx, gene_local, npcs, scale_max, seed,
            )
            blocks.append(block)
            feat_names.extend(f"{gene}_PC_{k + 1}" for k in range(block.shape[1]))
            genes_used.append(gene)
            if verbose:
                print(
                    f"[mixscape_lda] {gene}: {len(de_genes)} DE genes → "
                    f"{block.shape[1]} PCs"
                )
    finally:
        seurat.idents = saved_ident

    if not blocks:
        raise ValueError(
            f"No guide reached npcs+1={npcs + 1} DE genes; lower `npcs` or the DE "
            f"thresholds."
        )

    features = np.hstack(blocks)                     # cells × (guides · npcs)
    y = np.array([str(v) for v in labels_vec], dtype=object)
    n_comp = min(len(set(y.tolist())) - 1, features.shape[1])

    lda = LinearDiscriminantAnalysis(n_components=n_comp)
    lda.fit(features, y)
    embeddings = lda.transform(features)             # cells × n_comp
    loadings = np.asarray(lda.scalings_)[:, :n_comp]
    assignments = lda.predict(features)
    posterior = lda.predict_proba(features)          # cells × n_classes
    classes = [str(c) for c in lda.classes_]

    from .dimreduc import DimReduc

    dr = DimReduc(
        cell_embeddings=embeddings,
        cell_names=list(cells),
        feature_loadings=loadings,
        feature_names=list(feat_names),
        assay_used=assay,
        key=reduction_key,
        misc={
            "assignments": list(assignments),
            "posterior": posterior,
            "classes": classes,
            "genes_used": genes_used,
            "npcs": npcs,
        },
    )
    seurat.reductions[reduction_name] = dr

    target = seurat.cell_names()

    def put(col, values):
        seurat.meta_data[col] = (
            pd.Series(list(values), index=cells).reindex(target).values
        )

    put("lda_assignments", assignments)
    for j, cls in enumerate(classes):
        put(f"LDAP_{cls}", posterior[:, j])
    return seurat


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------


def _dense(mat) -> np.ndarray:
    if sp.issparse(mat):
        return np.asarray(mat.toarray(), dtype=float)
    return np.asarray(mat, dtype=float)


def _layer_matrix(assay_obj, layer: str):
    """Return ``(dense features × cells, feature_names, cell_names)`` for a layer."""
    from .assay5 import Assay5

    if isinstance(assay_obj, Assay5):
        if layer not in assay_obj.layers:
            raise KeyError(f"Layer {layer!r} not found on assay.")
        feats = list(assay_obj._layer_features.get(layer, assay_obj._all_feature_names))
        cells = list(assay_obj._layer_cells.get(layer, assay_obj._all_cell_names))
        mat = assay_obj.layers[layer]
    else:
        feats = list(assay_obj._feature_names)
        cells = list(assay_obj._cell_names)
        mat = getattr(assay_obj, layer.replace(".", "_"), None)
        if mat is None:
            mat = assay_obj.data
    return _dense(mat), feats, cells


def _aligned_meta(seurat, column: str, cells) -> np.ndarray:
    """The metadata ``column`` as an array aligned to the assay's ``cells``."""
    if column not in seurat.meta_data.columns:
        raise KeyError(f"Metadata column {column!r} not found.")
    return seurat.meta_data[column].reindex(cells).to_numpy()


def _split_groups(seurat, split_by: Optional[str], cells) -> list[np.ndarray]:
    """Index groups of cells to process independently (one group if no split)."""
    n = len(cells)
    if split_by is None:
        return [np.arange(n)]
    vals = _aligned_meta(seurat, split_by, cells)
    groups = []
    for level in pd.unique(vals):
        groups.append(np.where(vals == level)[0])
    return groups


def _nt_neighbor_means(emb, gidx, nt_local, data, k) -> np.ndarray:
    """Mean expression of each cell's ``k`` nearest NT neighbours (features × |gidx|)."""
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k).fit(emb[nt_local])
    _, idx = nn.kneighbors(emb[gidx])           # |gidx| × k, indices into nt_local
    nt_data = data[:, nt_local]                 # features × n_nt
    # For each query cell, average the k neighbour columns.
    means = nt_data[:, idx].mean(axis=2)        # features × |gidx|
    return means


def _de_genes(
    seurat, gene, nt_class, de_assay, de_test,
    logfc_threshold, min_pct, pval_cutoff, sig_feat_idx,
):
    """DE genes between ``gene`` cells and NT cells, restricted to signature genes."""
    from .markers import find_markers

    res = find_markers(
        seurat,
        ident_1=gene,
        ident_2=nt_class,
        assay=de_assay,
        test_use=de_test,
        logfc_threshold=logfc_threshold,
        min_pct=min_pct,
    )
    if res.empty or "p_val_adj" not in res.columns:
        return []
    passed = res.index[res["p_val_adj"] < pval_cutoff]
    return [g for g in passed if g in sig_feat_idx]


def _mixscape_em(sig, de_rows, nt_idx, gene_local, iter_num, seed):
    """Iterative 2-component mixture split of one gene's cells into KO / NP.

    Returns ``(ko_positions, posterior, n_iter, score)`` where ``ko_positions``
    index into ``gene_local`` (which cells are KO), ``posterior`` is the
    KO-component posterior for every gene cell (aligned to ``gene_local``), and
    ``score`` is the *first-iteration* perturbation score for every cell in
    ``concatenate([nt_idx, gene_local])`` — R's ``pvec``, kept for
    :func:`~shanuz.plotting.plot_perturb_score`. R stores it under
    ``if (n.iter == 0)``, i.e. from the round where the gene's cells are still
    taken whole, before any KO/NP split feeds back into the vector; later rounds
    would shift the axis and no longer match the published density plot.
    """
    from sklearn.mixture import GaussianMixture

    dat = sig[np.ix_(de_rows, np.concatenate([nt_idx, gene_local]))]
    n_nt = nt_idx.size
    n_gene = gene_local.size
    nt_cols = np.arange(n_nt)
    gene_cols = np.arange(n_nt, n_nt + n_gene)

    ko_pos = np.arange(n_gene)                  # start: all gene cells are KO
    post = np.full(n_gene, np.nan)
    nt_mean = dat[:, nt_cols].mean(axis=1)
    n_iter = 0
    score = None

    for it in range(1, iter_num + 1):
        n_iter = it
        if ko_pos.size == 0:
            break
        vec = dat[:, gene_cols[ko_pos]].mean(axis=1) - nt_mean
        proj = (vec @ dat).astype(float)        # perturbation score per local cell
        if it == 1:
            # R's ProjectVec is (v1 %*% v2) / (v2 %*% v2); the denominator is one
            # positive constant across cells, so the mixture split is unchanged by
            # it — but the score is a plotted axis, so scale it as R does.
            denom = float(vec @ vec)
            score = proj / denom if denom > 0 else proj.copy()
        if np.ptp(proj) <= 0:
            break
        try:
            gm = GaussianMixture(n_components=2, random_state=seed).fit(
                proj.reshape(-1, 1)
            )
        except Exception:
            break
        ko_comp = int(np.argmax(gm.means_.ravel()))
        prob_ko = gm.predict_proba(proj.reshape(-1, 1))[:, ko_comp]
        gene_prob = prob_ko[n_nt:]
        post = gene_prob
        new_ko = np.where(gene_prob > 0.5)[0]
        if new_ko.size == ko_pos.size and np.array_equal(np.sort(new_ko), np.sort(ko_pos)):
            ko_pos = new_ko
            break
        ko_pos = new_ko

    if np.all(np.isnan(post)):
        post = np.zeros(n_gene)
    return ko_pos, post, n_iter, score


def _lda_guide_block(sig, de_rows, nt_idx, gene_local, npcs, scale_max, seed):
    """One guide's LDA feature block: every cell projected onto its PCA subspace.

    The reference — this guide's cells plus the NT cells — sets both the per-gene
    centring/scaling and the PCA loadings; all cells in the dataset are then pushed
    through that same map, so the block says where each cell falls along *this*
    guide's perturbation axes. Returns ``n_all_cells × npcs``.
    """
    from sklearn.decomposition import PCA

    ref_cols = np.concatenate([nt_idx, gene_local])
    X = sig[de_rows, :]                              # n_de × n_all

    # ScaleData against the reference cells.
    ref = X[:, ref_cols]
    mu = ref.mean(axis=1, keepdims=True)
    sd = ref.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    scaled = (X - mu) / sd
    if scale_max is not None:
        scaled = np.clip(scaled, -scale_max, scale_max)

    pca = PCA(n_components=npcs, random_state=seed)
    pca.fit(scaled[:, ref_cols].T)                   # fit on guide + NT
    return pca.transform(scaled.T)                   # project every cell


def _assign(mixscape_class, global_class, gene_local, gene, kind, prtb_type):
    """Bulk-label a gene's cells (used for the all-NP shortcut paths)."""
    label = f"{gene} {prtb_type}" if kind == prtb_type else f"{gene} {kind}"
    for cell in gene_local:
        mixscape_class[cell] = label
        global_class[cell] = kind
