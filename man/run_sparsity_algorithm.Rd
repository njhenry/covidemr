% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R
\name{run_sparsity_algorithm}
\alias{run_sparsity_algorithm}
\title{Run TMB precision matrix sparsity algorithm}
\usage{
run_sparsity_algorithm(adfun, verbose = FALSE)
}
\arguments{
\item{adfun}{The ADFunction to normalize}

\item{verbose}{[bool, default FALSE] return a message about normalization?}
}
\description{
Run "symbolic analysis", a set of algorithms that prune the
  precision matrix to increase sparsity and speed up optimization
}
