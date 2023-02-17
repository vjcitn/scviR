#' produce a heatmap from a specialized CITE-seq SingleCellExperiment
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @param x SingleCellExperiment instance that has an `se.averaged` component
#' in its metadata
#' @param lb numeric(1) lower bound on 'breaks' sequence for ComplexHeatmap::pheatmap, defaults
#' to -3
#' @param ub numeric(1) upper bound on 'breaks' sequence for ComplexHeatmap::pheatmap,
#' defaults to 3
#' @param do_z logical(1) if TRUE, divide the residuals by their standard deviation across clusters,
#' defaults to false
#' @return ComplexHeatmap::pheatmap instance
#' @note See the OSCA book ch12.5.2 for the application.  
#' @examples
#' ch12sce = get_ch12sce()
#' adt_profiles(ch12sce)
#' adt_profiles(ch12sce, do_z = TRUE)
#' @export
adt_profiles = function(x, lb=-3, ub=3, do_z=FALSE) {
 stopifnot("se.averaged" %in% names(metadata(x)))
 avg = assay(metadata(x)$se.averaged)
 sds = rowSds(avg)
 quant = avg - rowMeans(avg)
 if (do_z) quant = t(t(quant)/sds)
 pheatmap::pheatmap(quant,
  breaks=seq(lb,ub, length.out=101))
}

