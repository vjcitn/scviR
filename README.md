# scviR

The scviR package provides an
experimental interface between R and [scvi-tools](https://docs.scvi-tools.org/en/stable/).

Our first release addresses the use of the [totalVI
model for CITE-seq data](https://docs.scvi-tools.org/en/stable/user_guide/models/totalvi.html).

- The scviR vignette works through a chunk of the colab tutorial
for [scvi-tools 0.19.0](https://colab.research.google.com/github/scverse/scvi-tutorials/blob/0.20.0/totalVI.ipynb); 0.20.0 employs muon, and this has not been addressed yet.
- scviR defines python infrastructure via the [basilisk](https://bioconductor.org/packages/basilisk)
discipline; the main python dependencies are declared in `R/basilisk.R`.
- We have collected a number of intermediate results so that the outputs of totalVI
can be explored without taking the time to fit the model.  The most important one
is the anndata instance with representations of the latent space, cluster
assignments, and UMAP projection:

```
> tot = get_totalVI_5k10k_adata() # retrieved on first call from Open Storage Network, cached
> tot
AnnData object with n_obs × n_vars = 10849 × 4000
    obs: 'n_genes', 'percent_mito', 'n_counts', 'batch', '_scvi_labels', '_scvi_batch', 'leiden_totalVI'
    var: 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches'
    uns: '_scvi_manager_uuid', '_scvi_uuid', 'hvg', 'leiden', 'log1p', 'neighbors', 'umap'
    obsm: 'X_totalVI', 'X_umap', 'denoised_protein', 'protein_expression', 'protein_foreground_prob'
    layers: 'counts', 'denoised_rna'
    obsp: 'connectivities', 'distances'
> table(tot$obs$batch)

 PBMC5k PBMC10k 
   3994    6855 
```
