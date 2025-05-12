##!/usr/bin/env python
## coding: utf-8
#
## # CITE-seq analysis with totalVI
## 
## With totalVI, we can produce a joint latent representation of cells, denoised data for both protein and RNA, integrate datasets, and compute differential expression of RNA and protein. Here we demonstrate this functionality with an integrated analysis of PBMC10k and PBMC5k, datasets of peripheral blood mononuclear cells publicly available from 10X Genomics subset to the 14 shared proteins between them. The same pipeline would generally be used to analyze a single CITE-seq dataset.
## 
## If you use totalVI, please consider citing:
## 
## - Gayoso, A., Steier, Z., Lopez, R., Regier, J., Nazor, K. L., Streets, A., & Yosef, N. (2021). Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nature Methods, 18(3), 272-282.
#
## In[2]:
#
#
##get_ipython().system('pwd')
#
#
## ```{note}
## Running the following cell will install tutorial dependencies on Google Colab only. It will have no effect on environments other than Google Colab.
## ```
#
## In[3]:
#
#
##!pip install --quiet scvi-colab
##from scvi_colab import install
#
##install()
#
#
## In[4]:
#
#
#import os
#os.system('pip3 install mudata')
#
#
## In[5]:
#
#
#os.system('pip3 install muon')
#
#
## In[6]:
#
#
#os.system('pip3 install scvi-tools')
#
#
## In[7]:
#
#
#os.system('pip3 install packaging==24.2')
#
#
## In[8]:
#
#
#import tempfile
#
#import anndata as ad
#import matplotlib.pyplot as plt
#import mudata as md
#import muon
#import numpy as np
#import scanpy as sc
#import scvi
#import seaborn as sns
#import torch
#
#
## ## Imports and data loading
#
## In[9]:
#
#
#scvi.settings.seed = 0
#print("Last run with scvi-tools version:", scvi.__version__)
#
#
## ```{note}
## You can modify `save_dir` below to change where the data files for this tutorial are saved.
## ```
#
## In[10]:
#
#
#sc.set_figure_params(figsize=(6, 6), frameon=False)
#sns.set_theme()
#torch.set_float32_matmul_precision("high")
#save_dir = tempfile.TemporaryDirectory()
#print(save_dir)
##run_line_magic('config', 'InlineBackend.print_figure_kwargs={"facecolor": "w"}')
##run_line_magic('config', 'InlineBackend.figure_format="retina"')
#
#
## This data that is used in this tutorial is the same from the totalvi manuscript and is PBMCs from 10x Genomics profiled with RNA and protein that was already filtered as described in the totalVI manuscript (low quality cells, doublets, lowly expressed genes, etc.) and in this link from the [totalvi reproducibility repo](https://github.com/YosefLab/totalVI_reproducibility/blob/master/data/data_filtering_scripts/pbmc_10k/pbmc_10k.py). We do this for both dataests (5K and 10K). Note that this filtering script is old and perhaps not all functionalities work with recent updates, but the general purpose and ideas remain the same.
## We then use those 2 filtered datasets and merge them so that only shared proteins remains for each. This process can be seen in the _load_pbmcs_10x_cite_seq function [here](https://github.com/scverse/scvi-tools/blob/main/src/scvi/data/_built_in_data/_cite_seq.py)
## 
#
## In[11]:
#
#
#adata = scvi.data.pbmcs_10x_cite_seq(save_path=save_dir.name)
##adata = ad.read_h5ad("/Users/vincentcarey/BIOC_SOURCES/scviR/pbmc_10k_protein_v3.h5ad")
#
#
## We run the standard workflow for keeping count and normalized data together:
#
## In[12]:
#
#
#adata.layers["counts"] = adata.X.copy()
#sc.pp.normalize_total(adata)
#sc.pp.log1p(adata)
#adata.obs_names_make_unique()
#
#
## ```{important}
## In this tutorial we will show totalVI's compatibility with the [MuData](https://mudata.readthedocs.io/en/latest/api/generated/mudata.MuData.html#mudata.MuData) format, which is a container for multiple AnnData objects. MuData objects can be read from the outputs of CellRanger using `muon.read_10x_h5`.
## 
## Furthermore, AnnData alone can also be used by storing the protein count data in `.obsm`, which is how it already is. For the AnnData-only workflow, see the documentation for `setup_anndata` in `scvi.model.TOTALVI`.
## ```
#
## In[13]:
#
#
#protein_adata = ad.AnnData(adata.obsm["protein_expression"])
#protein_adata.obs_names = adata.obs_names
#del adata.obsm["protein_expression"]
#mdata = md.MuData({"rna": adata, "protein": protein_adata})
#mdata
#
#
## In[14]:
#
#
#os.system('pip3 install scikit-misc')
#
#
## In[15]:
#
#
#import skmisc
#
#
## In[16]:
#
#
#mdata
#
#
## In[17]:
#
#
#sc.pp.highly_variable_genes(
#    mdata.mod["rna"],
#    n_top_genes=4000,
#    flavor="seurat_v3",
#    #batch_key="batch",
#    layer="counts",
#)
## Place subsetted counts in a new modality
#mdata.mod["rna_subset"] = mdata.mod["rna"][:, mdata.mod["rna"].var["highly_variable"]].copy()
#
#
## In[ ]:
#
#
#
#
#
## In[18]:
#
#
#mdata.update()
#
#
## In[19]:
#
#
#mdata
#
#
## ### Using own data
#
## If you want to analyze your own data, there are multiple ways to upload the count matrix \(in h5 or mtx format\):
## * It can be upload directly to the colab session storage, use the `Upload to session sotrage` button on the upper left corner in the `Files` tab.
## * You can also upload the data to Google Drive, and mount it using the Mount Drive button in the Files tab (for Colab users only).
## * Download directly from public available data source using command line tools such as curl. This will work not only for Colab.
## 
## In this example we are downloading h5 format of single-cell multiomic data generated by Proteintech Genomics directly into the Google Colab session storage. The data is from human resting PBMCs stained with the MultiPro® Human Discovery Panel (HDP) followed by processing using 10x Genomics Flex chemistry with Feature Barcoding Technology.
#
## In[20]:
#
#
#os.system('curl -LO https://ptgngsdata.s3.us-west-2.amazonaws.com/counts/E44_1_restPBMC_DCpos_filtered_feature_bc_matrix.h5')
#
#
## We then run the same preprocess workflow we show before: 
#
## In[ ]:
#
#
#
#
#
## In[21]:
#
#
library(scviR)
mdata1 = muonR()$read_10x_h5("E44_1_restPBMC_DCpos_filtered_feature_bc_matrix.h5")
mdata1$mod["rna"]$var_names_make_unique()
#
#
## In[22]:
#
#
reticulate::py_run_string('r.mdata1.mod["rna"].layers["counts"] = r.mdata1.mod["rna"].X.copy()')
scr = scanpyR()
scr$pp$normalize_total(mdata1$mod["rna"])
scr$pp$log1p(mdata1$mod["rna"])
#
scr$pp$highly_variable_genes(
    mdata1$mod["rna"],
    n_top_genes=4000L,
    flavor="seurat_v3",
    layer="counts",
)

## Place subsetted counts in a new modality
#mdata1.mod["rna_subset"] = mdata1.mod["rna"][:, mdata1.mod["rna"].var["highly_variable"]].copy()
#
#
## Becuase of the filtering process we will re create the mdata here
#
## In[23]:
#
#
#mdata1 = md.MuData(mdata1.mod)
#
#
## In[24]:
#
#
## we need to work with dense and not sparse matrices:
#mdata1["prot"].X = mdata1["prot"].X.toarray()
#mdata1["rna_subset"].X = mdata1["rna_subset"].X.toarray()
#mdata1.mod["rna_subset"].layers["counts"] = mdata1.mod["rna_subset"].layers["counts"].toarray()
#
#
## In[25]:
#
#
#mdata1.update()
#
#
## In[26]:
#
#
#mdata1
#
#
## In[27]:
#
#
#mdata
#
#
## ### Setup mudata
#
## Now we run `setup_mudata`, which is the MuData analog to `setup_anndata`. The caveat of this workflow is that we need to provide this function which modality of the `mdata` object contains each piece of data. So for example, the batch information is in `mdata.mod["rna"].obs["batch"]`. Therefore, in the `modalities` argument below we specify that the `batch_key` can be found in the `"rna_subset"` modality of the MuData object.
## 
## Notably, we provide `protein_layer=None`. This means scvi-tools will pull information from `.X` from the modality specified in `modalities` (`"protein"` in this case). In the case of RNA, we want to use the counts, which we stored in `mdata.mod["rna"].layers["counts"]`.
#
## In[28]:
#
#
#scvi.model.TOTALVI.setup_mudata(
#    mdata,
#    rna_layer="counts",
#    protein_layer=None,
#    #batch_key="batch",
#    modalities={
#        "rna_layer": "rna_subset",
#        "protein_layer": "protein",
#        "batch_key": "rna_subset",
#    },
#)
#
#
## In[29]:
#
#
#scvi.model.TOTALVI.setup_mudata(
#    mdata1,
#    rna_layer="counts",
#    protein_layer=None,
#    modalities={
#        "rna_layer": "rna_subset",
#        "protein_layer": "prot",
#        "batch_key": "rna_subset",
#    },
#)
#
#
## ```{note}
## Specify the modality of each argument via the `modalities` dictionary, which maps layer/key arguments to MuData modalities.
## ```
#
## ## Prepare and run model
#
## For the rest of this tutorial we will use the model with the external data we loaded.
#
## In[30]:
#
#
#model = scvi.model.TOTALVI(mdata1)
#
#
## In[ ]:
#
#
#model.train()
#
#
## We will plot the loss curves for training and validation using auto alignment for the yaxis 
#
## In[ ]:
#
#
#last_val_valid = np.array(model.history["elbo_validation"])[-1]
#last_val_train = np.array(model.history["elbo_train"])[-1]
#global_min_loss = min(
#    np.min(model.history["elbo_train"]), np.min(model.history["elbo_validation"])
#)
#last_max_loss = max(last_val_train, last_val_valid)[0]
#global_max_loss = max(
#    np.max(model.history["elbo_train"]), np.max(model.history["elbo_validation"])
#)
#
## In[ ]:
#
#
## Compute the min and max of both train and validation losses
#min_loss = min(min(last_val_train, last_val_valid), global_min_loss)
#max_loss = max(max(last_val_train, last_val_valid), global_max_loss)
#ylim_min = 0.995 * min_loss  # 0.5% below the minimum
#ylim_max = min(
#    global_max_loss, ylim_min + (last_max_loss - ylim_min) * 4
#)  # keep it under the 25% part of figure
#
#
## In[ ]:
#
#
#fig, ax = plt.subplots(1, 1)
#model.history["elbo_train"].plot(ax=ax, label="train")
#model.history["elbo_validation"].plot(ax=ax, label="validation")
#ax.set(
#    title="Negative ELBO over training epochs", ylim=(ylim_min, ylim_max)
#)  # you can still plug in you numbers
#ax.legend()
#
#
## ## Analyze outputs
## 
## We use Scanpy and muon for clustering and visualization after running totalVI. It's also possible to save totalVI outputs for an R-based workflow.
#
## In[ ]:
#
#
#rna = mdata1.mod["rna_subset"]
#protein = mdata1.mod["prot"]
#
#
## In[ ]:
#
#
## arbitrarily store latent in rna modality
#TOTALVI_LATENT_KEY = "X_totalVI"
#rna.obsm[TOTALVI_LATENT_KEY] = model.get_latent_representation()
#
#
## In[ ]:
#
#
#rna_denoised, protein_denoised = model.get_normalized_expression(n_samples=25, return_mean=True)
#rna.layers["denoised_rna"] = rna_denoised
#protein.layers["denoised_protein"] = protein_denoised
#
#protein.layers["protein_foreground_prob"] = 100 * model.get_protein_foreground_probability(
#    n_samples=25, return_mean=True
#)
#parsed_protein_names = [p.split("_")[0] for p in protein.var_names]
#protein.var["clean_names"] = parsed_protein_names
#mdata1.update()
#
#
## Now we can compute clusters and visualize the latent space.
#
## In[ ]:
#
#
#TOTALVI_CLUSTERS_KEY = "leiden_totalVI"
#
#sc.pp.neighbors(rna, use_rep=TOTALVI_LATENT_KEY)
#sc.tl.umap(rna)
#sc.tl.leiden(rna, key_added=TOTALVI_CLUSTERS_KEY)
#
#
## In[ ]:
#
#
#mdata1.update()
#
#
## We can now use muon plotting functions which can pull data from either modality of the MuData object. We will show the umap of the model embeddings with leiden clusters (and batch inegration of the datasets if exists). Following that we will show the denoised protein values and the foreground probability of the 14 protein listed.
#
## In[ ]:
#
#
#muon.pl.embedding(
#    mdata1,
#    basis="rna_subset:X_umap",
#    color=[f"rna_subset:{TOTALVI_CLUSTERS_KEY}"],
#    frameon=False,
#    ncols=1,
#)
#
#
## ### Visualize denoised protein values
#
## In[ ]:
#
#
#max_prot_to_plot = 14
#protein.var_names[:max_prot_to_plot]
#
#
## In[ ]:
#
#
#muon.pl.embedding(
#    mdata1,
#    basis="rna_subset:X_umap",
#    color=protein.var_names[:max_prot_to_plot],
#    frameon=False,
#    ncols=3,
#    vmax="p99",
#    wspace=0.1,
#    layer="denoised_protein",
#)
#
#
## ### Visualize probability of foreground
## 
## Here we visualize the probability of foreground for the first 14 proteins in the list and cell (projected on UMAP).
## Some proteins are easier to disentangle than others. Some proteins end up being "all background".
## For example, CD15 does not appear to be captured well, when looking at the denoised values above we see little localization in the monocytes.
#
## ```{note}
## While the foreground probability could theoretically be used to identify cell populations, we recommend using the denoised protein expression, which accounts for the foreground/background probability, but preserves the dynamic range of the protein measurements. Consequently, the denoised values are on the same scale as the raw data and it may be desirable to take a transformation like log or square root.
## ```
#
## By viewing the foreground probability, we can get a feel for the types of cells in our dataset. For example, it's very easy to see a population of monocytes based on the CD14 foregroud probability.
#
## In[ ]:
#
#
#muon.pl.embedding(
#    mdata1,
#    basis="rna_subset:X_umap",
#    layer="protein_foreground_prob",
#    color=protein.var_names[:max_prot_to_plot],
#    frameon=False,
#    ncols=3,
#    vmax="p99",
#    wspace=0.1,
#    color_map="cividis",
#)
#
#
## ## Differential expression
#
## Here we do a one-vs-all DE test, where each cluster is tested against all cells not in that cluster. The results for each of the one-vs-all tests is concatenated into one DataFrame object. Inividual tests can be sliced using the "comparison" column. Genes and proteins are included in the same DataFrame.
#
## ```{important}
## We do not recommend using totalVI denoised values in other differential expression tools, as denoised values are a summary of a random quantity. The totalVI DE test takes into account the full uncertainty of the denoised quantities.
## ```
#
## In[ ]:
#
#
#de_df = model.differential_expression(
#    groupby="rna_subset:leiden_totalVI", delta=0.5, batch_correction=True
#)
#de_df.head(5)
#
#
## Now we filter the results such that we retain features above a certain Bayes factor (which here is on the natural log scale) and genes with greater than 10% non-zero entries in the cluster of interest.
#
## In[ ]:
#
#
#filtered_pro = {}
#filtered_rna = {}
#cats = rna.obs[TOTALVI_CLUSTERS_KEY].cat.categories
#for c in cats:
#    cid = f"{c} vs Rest"
#    cell_type_df = de_df.loc[de_df.comparison == cid]
#    cell_type_df = cell_type_df.sort_values("lfc_median", ascending=False)
#
#    cell_type_df = cell_type_df[cell_type_df.lfc_median > 0]
#
#    pro_rows = cell_type_df.index.str.contains("TotalSeqB")
#    data_pro = cell_type_df.iloc[pro_rows]
#    data_pro = data_pro[data_pro["bayes_factor"] > 0.7]
#
#    data_rna = cell_type_df.iloc[~pro_rows]
#    data_rna = data_rna[data_rna["bayes_factor"] > 3]
#    data_rna = data_rna[data_rna["non_zeros_proportion1"] > 0.1]
#
#    filtered_pro[c] = data_pro.index.tolist()[:3]
#    filtered_rna[c] = data_rna.index.tolist()[:2]
#
#
## We can also use general scanpy visualization functions
#
## In[ ]:
#
#
#sc.tl.dendrogram(rna, groupby=TOTALVI_CLUSTERS_KEY, use_rep=TOTALVI_LATENT_KEY)
## This is a bit of a hack to be able to use scanpy dendrogram with the protein data
#protein.obs[TOTALVI_CLUSTERS_KEY] = rna.obs[TOTALVI_CLUSTERS_KEY]
#protein.obsm[TOTALVI_LATENT_KEY] = rna.obsm[TOTALVI_LATENT_KEY]
#sc.tl.dendrogram(protein, groupby=TOTALVI_CLUSTERS_KEY, use_rep=TOTALVI_LATENT_KEY)
#
#
## Matrix plot displays totalVI denoised protein expression per leiden cluster.
#
## In[ ]:
#
#
#sc.pl.matrixplot(
#    protein,
#    protein.var["clean_names"],
#    groupby=TOTALVI_CLUSTERS_KEY,
#    gene_symbols="clean_names",
#    dendrogram=True,
#    swap_axes=True,
#    layer="denoised_protein",
#    cmap="Greens",
#    standard_scale="var",
#)
#
#
## This is a selection of some of the markers that turned up in the RNA DE test.
#
## In[ ]:
#
#
#sc.pl.umap(
#    rna,
#    color=[
#        TOTALVI_CLUSTERS_KEY,
#        "IGHD",
#        "FCER1A",
#        "SCT",
#        "GZMH",
#        "NOG",
#        "FOXP3",
#        "C1QA",
#        "SIGLEC1",
#        "XCL2",
#        "GZMK",
#    ],
#    legend_loc="on data",
#    frameon=False,
#    ncols=3,
#    layer="denoised_rna",
#    wspace=0.2,
#)
#
#
## In[ ]:
#
#
#
#
