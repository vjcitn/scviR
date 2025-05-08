#' grab scvi-tools muon-oriented VAE instance built on 
#' the PBMC datasets following the tutorial
#' @import BiocFileCache
#' @importFrom utils unzip
#' @note VAE construction followed tutorial at
#' `https://docs.scvi-tools.org/en/stable/tutorials/notebooks/totalVI.html`.
#' @note We are using the scvi tutorial read early may 2025.  The
#' notebook uses "h5 format of single-cell multiomic data generated 
#' by Proteintech Genomics ...
#' The data is from human resting PBMCs 
#' stained with the MultiPro® Human Discovery 
#' Panel (HDP) followed by processing using 10x Genomics 
#' Flex chemistry with Feature Barcoding Technology."
#' @note It may be advantageous to set `options(timeout=3600)` or to allow an even greater
#' time for internet downloads, if working at a relatively slow network connection.
#' @return invisibly, the path to the .zip file holding the weights in pt format for the fitted VAE
#' @examples
#' zpath <- cacheCiteseqHDPmodel()
#' td <- tempdir()
#' utils::unzip(zpath, exdir = td)
#' vaedir <- paste0(td, "/vae3_pt")
#' dir(vaedir)
#' @export
cacheCiteseqHDPmodel <- function() {
  ca <- BiocFileCache()
  .osn_bucket_to_cache("vae3_pt.zip")
}

#' retrieve and cache a 349-protein CITE-seq dataset as employed in
#' scvi-tools tutorial
#' @export
cacheCiteseqHDPdata <- function() {
  ca <- BiocFileCache()
  .osn_bucket_to_cache("E44_1_restPBMC_DCpos_filtered_feature_bc_matrix.h5")
}
