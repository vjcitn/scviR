
#' basic interface
#' @importFrom reticulate import
#' @examples
#' scvi = scviR()
#' scvi
#' scvi$model
#' @export
scviR = function() {
 proc = basilisk::basiliskStart(bsklenv)
 on.exit(basilisk::basiliskStop(proc))
 basilisk::basiliskRun(proc, function()
      reticulate::import("scvi"))
}


