# script to use scviR from Bioconductor to
# run most of the computations of https://colab.research.google.com/github/scverse/scvi-tutorials/blob/1.1.2/multimodal/totalVI.ipynb
#
# anndata will be acquired and cached by scvi.data
#
save_dir.name = "./total2024"
if (!dir.exists(save_dir.name)) {
  dir.create(save_dir.name)
}
library(scviR)
scvi = scviR()
adata = scvi$data$pbmcs_10x_cite_seq(save_path=save_dir.name)
adata

# these functions are defined in scviR ... this is not
# the way basilisk clients should work, because they
# pass python references around, but it will take time
# to properly encapsulate all the relevant tasks

ad = anndataR()
sc = scanpyR()

# use scanpy utilities to normalize

adata$layers["counts"] = adata$X # $copy()
sc$pp$normalize_total(adata)
sc$pp$log1p(adata)
adata$obs_names_make_unique()

# build up the protein component
protein_adata = ad$AnnData(adata$obsm["protein_expression"])
protein_adata$obs_names = adata$obs_names
adata$obsm["protein_expression"] = NULL
#del adata.obsm["protein_expression"]

# use mudata for multiple modalities
md = MuDataR()

mdata = md$MuData(dict("rna"= adata, "protein"= protein_adata))
mdata

# filter genes
sc$pp$highly_variable_genes(
    mdata$mod[["rna"]], #["rna"],
    n_top_genes=4000,
    flavor="seurat_v3",
    batch_key="batch",
    layer="counts",
)

# works: mdata$mod$rna[ ,mdata$mod$rna$var$highly_variable]

# add into mudata
mdata$mod$rna_subset = mdata$mod$rna[ ,mdata$mod$rna$var$highly_variable]$copy()

# crucial
mdata$update()

# perform 'setup' -- excluding protein_layer spec seems essential
scvi$model$TOTALVI$setup_mudata(
    mdata,
    rna_layer="counts",
#    protein_layer=NA,
    batch_key="batch",
    modalities=dict(
        "rna_layer"= "rna_subset",
        "protein_layer"= "protein",
        "batch_key"= "rna_subset"
    ),
)

# configure and train
model = scvi$model$TOTALVI(mdata)
model

model$train(max_epochs=50L)
