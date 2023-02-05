

# necessary for python module control
bsklenv <- basilisk::BasiliskEnvironment(envname="bsklenv",
    pkgname="scviR",
    packages=c("numpy==1.23.1", "pandas==1.4.4"),
      pip=c("scvi-tools==0.20.0", "scanpy==1.9.1", "scikit-misc==0.1.4",
       "matplotlib==3.6.3"))

