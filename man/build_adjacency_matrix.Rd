% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_data_functions.R
\name{build_adjacency_matrix}
\alias{build_adjacency_matrix}
\title{Build adjacency matrix}
\usage{
build_adjacency_matrix(
  poly_sp,
  style = "B",
  allow_zero_neighbors = TRUE,
  manually_add_links = list()
)
}
\arguments{
\item{poly_sp}{A SpatialPolygonsDataFrame object}

\item{style}{One of "B","W","C", or "S". For more information, see the
documentation in the \code{\link{spdep}} package:
https://www.rdocumentation.org/packages/spdep/versions/1.1-3/topics/nb2mat}

\item{allow_zero_neighbors}{[bool, default TRUE] Allow polygons with no
neighbors?}

\item{manually_add_links}{[list, default NULL] A list of links to add in
addition to those specified by polygon adjacency. This should be formatted
as a list where each item is a numeric vector of length 2, specifying the
indices of polygons that should be listed as adjacent. The numeric vectors
are not sensitive to ordering (so item `c(1, 3)` is the same as `c(3, 1)`)}
}
\value{
sparse dsCMatrix representing adjacency between polygons in the
  `poly_sp` object
}
\description{
Builds a matrix describing adjacency between polygons in the
  a spatial object. This function is a wrapper of two functions in the
  \code{\link{spdep}} package.
}
