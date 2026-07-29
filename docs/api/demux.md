# Demultiplexing and pooled screens

Two workflows that both start from a second assay of oligo counts and end in a
per-cell call written back into `meta_data`.

**Cell hashing** assigns pooled samples from hashtag counts. `hto_demux` is
Seurat's `HTODemux` — CLR, cluster into `k = n_hashtags + 1`, per-hashtag
negative-binomial background threshold, then singlet / doublet / negative.
`multiseq_demux` is the MULTI-seq alternative, a Gaussian-KDE quantile threshold
per barcode. On the cross-species ground truth they are 99.81 % call-concordant
with R.

**Mixscape** separates real CRISPR knockouts from escapers. `calc_perturb_sig`
subtracts each cell's nearest non-targeting controls; `run_mixscape` then fits
the two-component mixture per guide; `mixscape_lda` builds the supervised map on
which each guide population separates.

!!! note "The CLR margin defaults are deliberate"
    `hto_demux` and `multiseq_demux` normalize across **features** (margin 1),
    not across cells. That is what `HTODemux` does, and it is not the same
    default as the general-purpose CLR path. Changing it to 2 to "make them
    consistent" would break agreement with Seurat.

## Cell hashing

::: truecell.hto.hto_demux

::: truecell.multiseq.multiseq_demux

## Mixscape

::: truecell.mixscape.calc_perturb_sig

::: truecell.mixscape.run_mixscape

::: truecell.mixscape.mixscape_lda
