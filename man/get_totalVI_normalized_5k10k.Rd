% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_citeseq_pbmcs.R
\name{get_totalVI_normalized_5k10k}
\alias{get_totalVI_normalized_5k10k}
\title{get matrices of normalized quantifications from full totalVI 5k10k from tutorial}
\usage{
get_totalVI_normalized_5k10k()
}
\value{
list of matrices
}
\description{
get matrices of normalized quantifications from full totalVI 5k10k from tutorial
}
\examples{
nmlist = get_totalVI_normalized_5k10k()
vapply(nmlist, dim, numeric(2))
}
