#' grab scvi-tools-processed PBMC CITE-seq data in anndata format (gzipped) from Open Storage Network
#' @import BiocFileCache
#' @note Original h5ad files obtained using scvi-tools 0.18.0 scvi.data.pbmcs_10x_cite_seq,
#' then processed according to steps in the scviR vignette, which follow the
#' scvi-tools tutorial by Gayoso et al.
#' @examples
#' h5path = cache_citeseq_pbmcs()
#' rhdf5::h5ls(h5path)
#' @export
cache_citeseq_pbmcs = function() {
  ca = BiocFileCache()
  pa = bfcquery(ca, "BiocScviR/demo1.h5ad")
  if (nchar(pa)>0) return(pa)
  gzdat = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocScviR/demo1.h5ad.gz"
  td = tempdir()
  targ = paste0(td, "/demo1.h5ad.gz")
  download.file(gzdat, targ)
  system(paste("gunzip", targ))  # bad?
  invisible(bfcrpath(ca, sub(".gz$", "", targ), action="move"))
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
