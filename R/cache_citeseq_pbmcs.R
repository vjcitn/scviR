#' grab scvi-tools-processed PBMC CITE-seq data in anndata format (gzipped) from Open Storage Network
#' @import BiocFileCache
#' @note Original h5ad files obtained using scvi-tools 0.18.0 scvi.data.pbmcs_10x_cite_seq,
#' then processed according to steps in the scviR vignette, which follow the
#' [scvi-tools tutorial](https://colab.research.google.com/github/scverse/scvi-tutorials/blob/0.18.0/totalVI.ipynb) by Gayoso et al.
#' @return invisibly, the path to the .h5ad file
#' @examples
#' h5path = cache_citeseq_5k10k_pbmcs()
#' cmeta = rhdf5::h5ls(h5path)
#' dim(cmeta)
#' head(cmeta, 17)
#' @export
cache_citeseq_5k10k_pbmcs = function() {
  ca = BiocFileCache()
  pa = bfcquery(ca, "demo2.h5ad")
# returns tibble
  if (nrow(pa)>1) stop("demo2.h5ad has multiple instances in cache, please inspect.")
  else if (nrow(pa)==1) return(pa$rpath)
# we need to retrieve if we get here
  gzdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/demo2.h5ad.gz"
  td = tempdir()
  targ = paste0(td, "/demo2.h5ad.gz")
  download.file(gzdat, targ)
  system(paste("gunzip", targ))  # bad?
  invisible(bfcrpath(ca, sub(".gz$", "", targ), action="move"))
}

#' grab scvi-tools VAE instance built on the PBMC datasets following the tutorial
#' @import BiocFileCache
#' @importFrom utils unzip
#' @note VAE construction followed tutorial at 
#' `https://docs.scvi-tools.org/en/stable/tutorials/notebooks/totalVI.html`.
#' @return invisibly, the path to the .zip file holding the fitted VAE and associated data
#' @examples
#' zpath = cache_citeseq_5k10k_tutvae()
#' td = tempdir()
#' utils::unzip(zpath, exdir=td)
#' vaedir = paste0(td, "/vae2_ov")
#' scvi = scviR()
#' adm = anndataR()
#' hpath = cache_citeseq_5k10k_pbmcs()
#' adata = adm$read(hpath)
#' mod = scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata, use_gpu=FALSE)
#' mod
#' @export
cache_citeseq_5k10k_tutvae = function() {
  ca = BiocFileCache()
  pa = bfcquery(ca, "vae2_ov.zip")
# returns tibble
  if (nrow(pa)>1) stop("vae2_ov.zip has multiple instances in cache, please inspect.")
  else if (nrow(pa)==1) return(pa$rpath)
  zdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/vae2_ov.zip"
  td = tempdir()
  targ = paste0(td, "/vae2_ov.zip")
  download.file(zdat, targ)
  invisible(bfcrpath(ca, targ, action="move"))
}
#
#
#
#
#Parameter validation failed:
#stvjc@stvjc-XPS-13-9300:~/YOSEF_Variational/scviR/R$ docker run -v /home/stvjc/OSN:/data -v /home/stvjc/RC2:/config/rclone -ti rclone/rclone:latest ls osn:/bir190004-bucket01/BiocScviR/
# 51842190 demo1.h5ad.gz
# 16477930 vae1_ov.zip
#st

#' helper to get the tutorial VAE for PBMCs from scvi-tools tutorial
#' @param use_gpu logical(1), defaulting to FALSE, passed to TOTALVI.load
#' @examples
#' get_citeseq_tutvae()
#' @export
get_citeseq_tutvae = function(use_gpu=FALSE) {
   zpath = cache_citeseq_5k10k_tutvae()
   td = tempdir()
   unzip(zpath, exdir=td)
   vaedir = paste0(td, "/vae2_ov")
   scvi = scviR()
   adm = anndataR()
   hpath = cache_citeseq_5k10k_pbmcs()
   adata = adm$read(hpath)
   mod = scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata, use_gpu=use_gpu)
   mod
}

#' helper to get the processed anndata for CITE-seq PBMCs from scvi-tools tutorial
#' @examples
#' get_citeseq_5k10k_pbmcs()
#' @export
get_citeseq_5k10k_pbmcs = function() {
   h5path = cache_citeseq_5k10k_pbmcs()
   anndataR()$read(h5path)
}

.osn_bucket_to_cache = function(entity, folder="BiocScviR",
    prefix="https://mghp.osn.xsede.org/bir190004-bucket01/",
    ca = BiocFileCache::BiocFileCache()) {
    pa = bfcquery(ca, entity)
    if (nrow(pa)>1) stop(sprintf("%s has multiple instances in cache, please inspect.",
               entity))
    else if (nrow(pa)==1) return(pa$rpath)
    target = paste0(prefix, folder, "/", entity)
    tf = tempfile()
    download.file(target, tf)
    bfcrpath(ca, target, action="move")
}

#' get an anndata reference to 5k10k protein after totalVI from tutorial
#' @examples
#' get_pro_5k10k_adata()
#' @export
get_pro_5k10k_adata = function() {
   ans = .osn_bucket_to_cache( "pbmc5k10k_pro_adata.h5ad" )
   anndataR()$read(ans)
}

#' get matrices of normalized quantifications from full totalVI 5k10k from tutorial
#' @examples
#' nmlist = get_totalVI_normalized_5k10k()
#' sapply(nmlist, dim)
#' @export
get_totalVI_normalized_5k10k = function() {
   ans = .osn_bucket_to_cache( "nmlzd_5k10k.rda" )
   load(ans, envir=.GlobalEnv)
   get("nmlzd_5k10k")
}

#' get anndata reference to full totalVI processing of 5k10k data
#' @examples
#' full = get_totalVI_5k10k_adata()
#' full
#' @export
get_totalVI_5k10k_adata = function() {
   ans = .osn_bucket_to_cache( "full_5k10k_totalVI.h5ad" )
   anndataR()$read(ans)
}


#
#get_totalVI_normalized_5k10k = function() {
#  ca = BiocFileCache()
#  pa = bfcquery(ca, "nmlzd_5k10k.rda")
#  if (nrow(pa)>1) stop("demo1.h5ad has multiple instances in cache, please inspect.")
## returns tibble
#  if (nrow(pa)>1) stop("demo1.h5ad has multiple instances in cache, please inspect.")
#  else if (nrow(pa)==1) return(pa$rpath)
## we need to retrieve if we get here
#  gzdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/demo1.h5ad.gz"
#  td = tempdir()
#  targ = paste0(td, "/demo1.h5ad.gz")
#  download.file(gzdat, targ)
#  system(paste("gunzip", targ))  # bad?
#  invisible(bfcrpath(ca, sub(".gz$", "", targ), action="move"))
#191170194 BiocScviR/nmlzd_5k10k.rda
#  3677260 BiocScviR/pbmc5k10k_pro_adata.h5ad
#
