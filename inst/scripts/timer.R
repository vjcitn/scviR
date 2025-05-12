## ----doe44--------------------------------------------------------------------
library(scviR)
HDP.h5 = cacheCiteseqHDPdata()
mdata1 = muonR()$read_10x_h5(HDP.h5)
mdata1$mod["rna"]$var_names_make_unique()
reticulate::py_run_string('r.mdata1.mod["rna"].layers["counts"] = r.mdata1.mod["rna"].X.copy()')
mdata1


## ----doscanpy-----------------------------------------------------------------
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


## ----domu---------------------------------------------------------------------
py_run_string('r.mdata1.mod["rna_subset"] = r.mdata1.mod["rna"][:, r.mdata1.mod["rna"].var["highly_variable"]].copy()')
mdata1 = MuDataR()$MuData(mdata1$mod)


## ----densify------------------------------------------------------------------
py_run_string('r.mdata1["prot"].X = r.mdata1["prot"].X.toarray()')
py_run_string('r.mdata1["rna_subset"].X = r.mdata1["rna_subset"].X.toarray()')
py_run_string('r.mdata1.mod["rna_subset"].layers["counts"] = r.mdata1.mod["rna_subset"].layers["counts"].toarray()')
mdata1$update()


## ----dosetup------------------------------------------------------------------
scviR()$model$TOTALVI$setup_mudata(
    mdata1,
    rna_layer="counts",
    protein_layer=reticulate::py_none(),
    modalities=list(
        "rna_layer"= "rna_subset",
        "protein_layer"= "prot",
        "batch_key"= "rna_subset"
    ), 
)


## ----what---------------------------------------------------------------------
model = scviR()$model$TOTALVI(mdata1)
model


## ----dotrain, message=FALSE, results="hide"-----------------------------------
n_epochs = 50L
acc = "cpu"
tchk = try(reticulate::import("torch"))
if (!inherits(tchk, "try-error") && tchk$backends$mps$is_available()) acc = "mps"
if (!inherits(tchk, "try-error") && tchk$backends$cuda$is_built()) acc = "gpu"
t1 = system.time(model$train(max_epochs=n_epochs, accelerator = acc))

## DO AGAIN

model = scviR()$model$TOTALVI(mdata1)
model


## ----dotrain, message=FALSE, results="hide"-----------------------------------
n_epochs = 50L
acc = "cpu"
tchk = try(reticulate::import("torch"))
if (!inherits(tchk, "try-error") && tchk$backends$mps$is_available()) acc = "mps"
if (!inherits(tchk, "try-error") && tchk$backends$cuda$is_built()) acc = "gpu"
t2 = system.time(model$train(max_epochs=n_epochs, accelerator = "cpu"))
