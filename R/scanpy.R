
#' basic interface
#' @importFrom reticulate import
#' @examples
#' sc = scanpyR()
#' sc
#' sc$pp
#' @export
scanpyR = function() {
  proc = basilisk::basiliskStart(bsklenv)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function()
    reticulate::import("scanpy"))
}


#' basic interface to anndata
#' @examples
#' ad = anndataR()
#' ad
#' ad$read
#' @export
anndataR = function() {
  proc = basilisk::basiliskStart(bsklenv)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function()
    reticulate::import("anndata"))
}
