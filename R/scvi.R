#' basic interface
#' @importFrom reticulate import
#' @return basiliskRun result with import from reticulate, typically a Module
#' @examples
#' scvi <- scviR()
#' scvi
#' scvi$model
#' @export
scviR <- function() {
  proc <- basilisk::basiliskStart(bsklenv, testload="anndata")
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function() {
   reticulate::py_run_string(" 
import sys
import pandas
sys.modules['pandas.core.indexes.numeric'] = pandas.core.indexes.base
pandas.core.indexes.base.Int64Index = pandas.core.indexes.base.Index")
  reticulate::import("scvi")
  })
}
