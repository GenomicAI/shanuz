# Loading data

## From 10x output

::: truecell.io.read_10x

## Bundled datasets

`truecell.datasets` stands in for R's `SeuratData`. Each loader downloads on first
call into `~/.truecell_data/` and returns a ready `Truecell` object; the whole set is
roughly 770 MB cached. These are the datasets the [tutorials](../tutorials/README.md)
run on, which is what makes each tutorial reproducible from a clean machine.

::: truecell.datasets.pbmc3k

::: truecell.datasets.pbmc8k

::: truecell.datasets.cbmc_citeseq

::: truecell.datasets.pbmc_hashing

::: truecell.datasets.ifnb

::: truecell.datasets.panc8

::: truecell.datasets.thp1_eccite

::: truecell.datasets.xenium_mouse_brain

::: truecell.datasets.visium_mouse_brain

## AnnData interoperability

::: truecell.compat.anndata.as_anndata

::: truecell.compat.anndata.from_anndata
