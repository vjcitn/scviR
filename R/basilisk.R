# necessary for python module control
#' python declarations
#' @import basilisk
bsklenv <- basilisk::BasiliskEnvironment(
  envname = "bsklenv",
  pkgname = "scviR",
  packages = c("python=3.10", "numpy=1.24.4", "pandas=1.5.3"),
  pip = c(
    "anndata==0.10.1",
    "docrep==0.3.2",
    "flax==0.8.1",
    "jax==0.4.23",
    "jaxlib==0.4.23",
    "optax==0.1.9",
    "scipy==1.12.0",
    "scikit-learn==1.4.0",
    "rich==13.7.1",
    "h5py==3.9.0",
    "torch==1.13.1",
    "lightning==2.1.4",
    "torchmetrics==0.11.0",
#    "pyro-ppl==1.6.0",
    "tqdm==4.66.2",
#    "numpyro==0.12.1",
    "ml-collections==0.1.1",
    "mudata==0.2.3",
    "scvi-tools==1.1.2", "scanpy==1.9.5",
    "matplotlib==3.7.3", "scikit-misc==0.3.1"
  )
)

# dropping scikit-misc on 18 Feb because it can't be used on M1 mac
