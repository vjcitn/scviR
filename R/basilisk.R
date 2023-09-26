# necessary for python module control
#' python declarations
#' @import basilisk
bsklenv <- basilisk::BasiliskEnvironment(
  envname = "bsklenv",
  pkgname = "scviR",
  packages = c("numpy==1.24.4", "pandas==1.5.3"),
  pip = c(
    "scvi-tools==1.0.3", "scanpy==1.9.5",
    "matplotlib==3.7.3"
  )
)

# dropping scikit-misc on 18 Feb because it can't be used on M1 mac
