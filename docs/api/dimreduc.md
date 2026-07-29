# Dimensional reduction

Linear first, then the embeddings you look at. Every one of these writes a
[`DimReduc`](objects.md#truecell.dimreduc.DimReduc) into `obj.reductions` under a
key, and records the features it actually used — not the features it was asked
for, which are not always the same set.

`jack_straw` is the permutation test for how many PCs to keep. It is worth
reading its docstring before trusting the number: R's `JackRandom` seeds each
replicate from its loop index and is therefore deterministic, while this one
seeds from its `seed` argument and moves. Across 60 seeds on PBMC 3k it keeps
12–15 PCs, mode 13, which is R's answer. That spread is asserted as a band, not
described in prose — see [Fidelity](../fidelity.md#bands).

## Linear

::: truecell.reduction.run_pca

::: truecell.reduction.run_spca

::: truecell.reduction.run_ica

::: truecell.glmpca.glm_pca

## Non-linear embeddings

::: truecell.umap.run_umap

::: truecell.reduction.run_tsne

## How many components to keep

::: truecell.jackstraw.jack_straw

::: truecell.jackstraw.score_jackstraw
