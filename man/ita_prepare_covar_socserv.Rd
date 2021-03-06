% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ita_covariate_prep_functions.R
\name{ita_prepare_covar_socserv}
\alias{ita_prepare_covar_socserv}
\title{Prepare social services covariate}
\usage{
ita_prepare_covar_socserv(covar_fp, model_years, location_table)
}
\arguments{
\item{covar_fp}{Filepath to the raw covariate}

\item{model_years}{Vector of years to consider for modeling}

\item{location_table}{Location code merge table for Italy}
}
\value{
data.table containing prepared covariate data and vector of indices
}
\description{
Covariate-specific prep function for social services coverage
  for families. This is a convenience function for specific use with datasets
  exported from IStat.
}
