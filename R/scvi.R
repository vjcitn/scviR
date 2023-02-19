
#' basic interface
#' @importFrom reticulate import
#' @return basiliskRun result with import from reticulate, typically a Module
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


