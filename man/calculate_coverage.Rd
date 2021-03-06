% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictive_validity_functions.R
\name{calculate_coverage}
\alias{calculate_coverage}
\title{Calculate coverage across predictive draws}
\usage{
calculate_coverage(
  in_data,
  num_field,
  denom_field,
  draw_fields,
  coverage_levels = c(0.5, 0.8, 0.95),
  group_fields = NULL,
  pois_sim = TRUE
)
}
\arguments{
\item{in_data}{Input data.table}

\item{num_field}{Numerator field for the observed data}

\item{denom_field}{Denominator field for the observed data}

\item{draw_fields}{Character vector of fields containing predictive draws, in
rate space (e.g. mortality rates)}

\item{coverage_levels}{[optional, default c(.5, .8, .95)] Uncertainty
intervals to calculate from the posterior predictive draws}

\item{group_fields}{[optional, default NULL] If the predictive validity
metrics should be grouped, list the fields to group them by here. If NULL
(the default), the predictive validity metrics will be calculated across
the entire dataset}

\item{pois_sim}{[optional, default TRUE] Should the coverage be estimated
using realizations Poisson distribution centered at population * rate (as
opposed to the central value, population * p)?}
}
\value{
Data.table containing the following fields:
  - 'covg<X>': Empirical coverage for the X% uncertainty interval
  - Any grouping columns specified in the function arguments
}
\description{
Calculate how often observed data falls within X% uncertainty
  intervals of the posterior predictive draws. The observed data should come
 in the form of separate numerator and denominator columns, while the
 estimates should be a single field estimating a rate (i.e. already
 normalized by denominator)
}
