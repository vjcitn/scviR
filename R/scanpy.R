#' basic interface
#' @importFrom reticulate import
#' @return basiliskRun result with import from reticulate, typically a Module
#' @examples
#' sc <- scanpyR()
#' sc
#' sc$pp
#' @export
scanpyR <- function() {
#  proc <- basilisk::basiliskStart(bsklenv, testload="anndata")
#  on.exit(basilisk::basiliskStop(proc))
#  basilisk::basiliskRun(proc, function() {
#    reticulate::import("scanpy")
#  })
reticulate::py_require("scanpy")
reticulate::py_require("scikit-misc")
reticulate::import("skmisc") # but only implicitly used ...
reticulate::import("scanpy")
}


#' basic interface to anndata
#' @return basiliskRun result with import from reticulate, typically a Module
#' @examples
#' ad <- anndataR()
#' ad
#' ad$read_h5ad
#' @export
anndataR <- function() {
#  proc <- basilisk::basiliskStart(bsklenv, testload="anndata")
#  on.exit(basilisk::basiliskStop(proc))
#  basilisk::basiliskRun(proc, function() {
#    reticulate::import("anndata")
#  })
reticulate::py_require("anndata")
reticulate::import("anndata")
}

#' basic interface to MuData
#' @return basiliskRun result with import from reticulate, typically a Module
#' @examples
#' md <- MuDataR()
#' md
#' head(names(md))
#' @export
MuDataR <- function() {
#  proc <- basilisk::basiliskStart(bsklenv, testload="anndata")
#  on.exit(basilisk::basiliskStop(proc))
#  basilisk::basiliskRun(proc, function() {
#    reticulate::import("mudata")
#  })
reticulate::py_require("mudata")
reticulate::import("mudata")
}

#' basic interface to muon
#' @return basiliskRun result with import from reticulate, typically a Module
#' @examples
#' md <- muonR()
#' md
#' head(names(md))
#' @export
muonR <- function() {
#  proc <- basilisk::basiliskStart(bsklenv, testload="anndata")
#  on.exit(basilisk::basiliskStop(proc))
#  basilisk::basiliskRun(proc, function() {
#    reticulate::import("muon")
#  })
reticulate::py_require("muon")
reticulate::import("muon")
}
