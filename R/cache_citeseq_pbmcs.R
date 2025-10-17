#' grab scvi-tools-processed PBMC CITE-seq data in anndata format (gzipped) from Open Storage Network
#' @import BiocFileCache
#' @importFrom utils download.file
#' @note Original h5ad files obtained using scvi-tools 0.18.0 scvi.data.pbmcs_10x_cite_seq,
#' then processed according to steps in the scviR vignette, which follow the
#' [scvi-tools tutorial](https://colab.research.google.com/github/scverse/scvi-tutorials/blob/0.18.0/totalVI.ipynb) by Gayoso et al.
#' @note It may be advantageous to set `options(timeout=3600)` or to allow an even greater
#' time for internet downloads, if working at a relatively slow network connection.
#' @return invisibly, the path to the .h5ad file
#' @examples
#' h5path <- cacheCiteseq5k10kPbmcs()
#' cmeta <- rhdf5::h5ls(h5path)
#' dim(cmeta)
#' head(cmeta, 17)
#' @export
cacheCiteseq5k10kPbmcs <- function() {
  ca <- BiocFileCache()
  pa <- bfcquery(ca, "pbmc_citeseq_5k10k.h5ad")
  # returns tibble
  if (nrow(pa) > 1) {
    stop("pbmc_citeseq_5k10k.h5ad has multiple instances in cache, please inspect.")
  } else if (nrow(pa) == 1) {
    return(pa$rpath)
  }
  # we need to retrieve if we get here
  gzdat <- "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/pbmc_citeseq_5k10k.h5ad.gz"
  td <- tempdir()
  targ <- paste0(td, "/pbmc_citeseq_5k10k.h5ad.gz")
  download.file(gzdat, targ)
  system(paste("gunzip", targ)) # bad?
  invisible(bfcrpath(ca, sub(".gz$", "", targ), action = "copy"))
}

#     new: character string: A suggestion for a replacement function.
#
# package: character string: The package to be used when suggesting
#          where the deprecated function might be listed.
#
#     msg: character string: A message to be printed, if missing a
#          default message is used.
#
#     old: character string specifying the function (default) or usage
#          which is being deprecated.


#' Deprecated: grab scvi-tools VAE instance built on the PBMC datasets following the tutorial
#' @import BiocFileCache
#' @importFrom utils unzip
#' @note the serialized model is obsolete
#' @note VAE construction followed tutorial at
#' `https://docs.scvi-tools.org/en/stable/tutorials/notebooks/totalVI.html`.
#' @note It may be advantageous to set `options(timeout=3600)` or to allow an even greater
#' time for internet downloads, if working at a relatively slow network connection.
#' @return invisibly, the path to the .zip file holding the fitted VAE and associated data
#' @examples
#' \dontrun{
#' zpath <- cacheCiteseq5k10kTutvae()
#' td <- tempdir()
#' utils::unzip(zpath, exdir = td)
#' vaedir <- paste0(td, "/vae2_ov")
#' scvi <- scviR()
#' adm <- anndataR()
#' hpath <- cacheCiteseq5k10kPbmcs()
#' adata <- adm$read_h5ad(hpath)
#' mod <- scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata) #, use_gpu = FALSE)
#' mod
#' }
#' @export
cacheCiteseq5k10kTutvae <- function() {
.Deprecated(new="cacheCiteseqHDPmodel", package="scviR",
   msg="the serialized model is obsolete", old="cacheCiteseq5k10kTutvae")
  ca <- BiocFileCache()
  pa <- bfcquery(ca, "vae2_ov.zip")
  # returns tibble
  if (nrow(pa) > 1) {
    stop("vae2_ov.zip has multiple instances in cache, please inspect.")
  } else if (nrow(pa) == 1) {
    return(pa$rpath)
  }
  zdat <- "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/vae2_ov.zip"
  td <- tempdir()
  targ <- paste0(td, "/vae2_ov.zip")
  download.file(zdat, targ)
  invisible(bfcrpath(ca, targ, action = "copy"))
}
#
#
#
#
# Parameter validation failed:
# stvjc@stvjc-XPS-13-9300:~/YOSEF_Variational/scviR/R$ docker run -v /home/stvjc/OSN:/data -v /home/stvjc/RC2:/config/rclone -ti rclone/rclone:latest ls osn:/bir190004-bucket01/BiocScviR/
# 51842190 demo1.h5ad.gz
# 16477930 vae1_ov.zip
# st

#' helper to get the tutorial VAE for PBMCs from scvi-tools tutorial
#' @param use_gpu logical(1), defaulting to FALSE, passed to TOTALVI.load
#' @return python reference to anndata
#' @note March 2024 use_gpu ignored
#' @examples
#' \dontrun{
#' getCiteseqTutvae()
#' }
#' @export
getCiteseqTutvae <- function(use_gpu = FALSE) {
  zpath <- cacheCiteseq5k10kTutvae()
  td <- tempdir()
  unzip(zpath, exdir = td)
  vaedir <- paste0(td, "/vae2_ov")
  scvi <- scviR()
  adm <- anndataR()
  hpath <- cacheCiteseq5k10kPbmcs()
  adata <- adm$read_h5ad(hpath)
  mod <- scvi$model$`_totalvi`$TOTALVI$load(vaedir, adata) #, use_gpu = use_gpu)
  mod
}

#' helper to get the processed anndata for CITE-seq PBMCs from scvi-tools tutorial
#' @return python reference to anndata
#' @note It may be advantageous to set `options(timeout=3600)` or to allow an even greater
#' time for internet downloads, if working at a relatively slow network connection.
#' @examples
#' getCiteseq5k10kPbmcs()
#' @export
getCiteseq5k10kPbmcs <- function() {
  h5path <- cacheCiteseq5k10kPbmcs()
  anndataR()$read_h5ad(h5path)
}

.osn_bucket_to_cache <- function(
    entity, folder = "BiocScviR",
    prefix = "https://mghp.osn.xsede.org/bir190004-bucket01/",
    ca = BiocFileCache::BiocFileCache()) {
  pa <- bfcquery(ca, entity)
  if (nrow(pa) > 1) {
    stop(sprintf(
      "%s has multiple instances in cache, please inspect.",
      entity
    ))
  } else if (nrow(pa) == 1) {
    return(pa$rpath)
  }
  target <- paste0(prefix, folder, "/", entity)
  tf <- tempfile(entity) # for metadata
  download.file(target, tf)
  bfcrpath(ca, tf, action = "copy")
}

#' get an anndata reference to 5k10k protein after totalVI from tutorial
#' @note It may be advantageous to set `options(timeout=3600)` or to allow an even greater
#' time for internet downloads, if working at a relatively slow network connection.
#' @return python reference to anndata
#' @examples
#' getPro5k10kAdata()
#' @export
getPro5k10kAdata <- function() {
  ans <- .osn_bucket_to_cache("pbmc5k10k_pro_adata.h5ad")
  anndataR()$read_h5ad(ans)
}

#' get matrices of normalized quantifications from full totalVI 5k10k from tutorial
#' @return list of matrices
#' @examples
#' nmlist <- getTotalVINormalized5k10k()
#' vapply(nmlist, dim, numeric(2))
#' @export
getTotalVINormalized5k10k <- function() {
  ans <- .osn_bucket_to_cache("nmlzd_5k10k.rda")
  load(ans, envir = .GlobalEnv)
  get("nmlzd_5k10k")
}

#' get anndata reference to full totalVI processing of 5k10k data
#' @return python reference to anndata
#' @examples
#' full <- getTotalVI5k10kAdata()
#' full
#' @export
getTotalVI5k10kAdata <- function() {
  ans <- .osn_bucket_to_cache("full_5k10k_totalVI.h5ad")
  anndataR()$read_h5ad(ans)
}

#' get SCE for 10k PBMC annotated as in OSCA book chapter 12
#' @param clear_cache logical(1) will delete relevant entries in available cache before continuing, defaults to FALSE
#' @note This is a SingleCellExperiment instance with data on 7472 cells from a 10x
#' CITE-seq experiment.  An altExp component includes
#' antibody-derived tag (ADT) counts on 17 proteins.  The data are acquired and
#' processed as described in ch 12 of the OSCA book, circa February 2023.
#' A metadata element (se.averaged) includes the result of averaging protein abundance
#' estimates within ADT-based clusters, as is done to give rise to Figure 12.8 of
#' the OSCA book.
#' @return SingleCellExperiment instance
#' @examples
#' ch12sce <- getCh12Sce()
#' ch12sce
#' @export
getCh12Sce <- function(clear_cache = FALSE) {
  if (clear_cache) {
    ca <- BiocFileCache::BiocFileCache()
    avail <- bfcquery(ca, "ch12sce.rda")
    if (nrow(avail) > 0) {
      to_kill <- avail$rid
      bfcremove(ca, to_kill)
    }
  }
  ans <- .osn_bucket_to_cache("ch12sce.rda")
  get(load(ans))
}

#' get list of cluster-specific SCE for 10k PBMC annotated as in OSCA book chapter 12
#' @note This is a list of SingleCellExperiment instances with data on a total of 7472 cells from a 10x
#' CITE-seq experiment.  An altExp component in each list element includes
#' antibody-derived tag (ADT) counts on 17 proteins.  The data are acquired and
#' processed as described in ch 12 of the OSCA book, circa February 2023.
#' List elements correspond to mRNA-based sub-clusters of ADT-based clusters.
#' @return SimpleList of SingleCellExperiment instances
#' @examples
#' ch12_allsce <- getCh12AllSce()
#' vapply(ch12_allsce, ncol, numeric(1))
#' @export
getCh12AllSce <- function() {
  ans <- .osn_bucket_to_cache("ch12_allsce.rda")
  get(load(ans))
}
