#' grab scvi-tools-processed PBMC CITE-seq data in anndata format (gzipped) from Open Storage Network
#' @import BiocFileCache
#' @note Original h5ad files obtained using scvi-tools 0.18.0 scvi.data.pbmcs_10x_cite_seq,
#' then processed according to steps in the scviR vignette, which follow the
#' scvi-tools tutorial by Gayoso et al.
#' @return invisibly, the path to the .h5ad file
#' @examples
#' h5path = cache_citeseq_pbmcs()
#' cmeta = rhdf5::h5ls(h5path)
#' dim(cmeta)
#' head(cmeta, 17)
#' @export
cache_citeseq_pbmcs = function() {
  ca = BiocFileCache()
  pa = bfcquery(ca, "demo1.h5ad")
# returns tibble
  if (nrow(pa)>1) stop("demo1.h5ad has multiple instances in cache, please inspect.")
  else if (nrow(pa)==1) return(pa$rpath)
# we need to retrieve if we get here
  gzdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/demo1.h5ad.gz"
  td = tempdir()
  targ = paste0(td, "/demo1.h5ad.gz")
  download.file(gzdat, targ)
  system(paste("gunzip", targ))  # bad?
  invisible(bfcrpath(ca, sub(".gz$", "", targ), action="move"))
}

#' grab scvi-tools VAE instance built on the PBMC datasets following the tutorial
#' @import BiocFileCache
#' @note VAE construction followed tutorial at 
#' `https://docs.scvi-tools.org/en/stable/tutorials/notebooks/totalVI.html`.
#' @return invisibly, the path to the .zip file holding the fitted VAE and associated data
#' @examples
#' zpath = cache_citeseq_tutvae()
#' td = tempdir()
#' unzip(zpath, exdir=td)
#' vaedir = paste0(td, "/vae1_ov")
#' scvi = scviR()
#' adm = anndataR()
#' hpath = cache_citeseq_pbmcs()
#' adata = adm$read(hpath)
#' mod = scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata, use_gpu=FALSE)
#' mod
#' @export
cache_citeseq_tutvae = function() {
  ca = BiocFileCache()
  pa = bfcquery(ca, "vae1_ov.zip")
# returns tibble
  if (nrow(pa)>1) stop("vae1_ov.zip has multiple instances in cache, please inspect.")
  else if (nrow(pa)==1) return(pa$rpath)
  zdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/vae1_ov.zip"
  td = tempdir()
  targ = paste0(td, "/vae1_ov.zip")
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
#' @examples
#' get_citeseq_tutvae()
#' @export
get_citeseq_tutvae = function() {
   zpath = cache_citeseq_tutvae()
   td = tempdir()
   unzip(zpath, exdir=td)
   vaedir = paste0(td, "/vae1_ov")
   scvi = scviR()
   adm = anndataR()
   hpath = cache_citeseq_pbmcs()
   adata = adm$read(hpath)
   mod = scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata, use_gpu=FALSE)
   mod
}

#' helper to get the processed anndata for CITE-seq PBMCs from scvi-tools tutorial
#' @examples
#' get_citeseq_pbmcs()
#' @export
get_citeseq_pbmcs = function() {
   h5path = cache_citeseq_pbmcs()
   anndataR()$read(h5path)
}
