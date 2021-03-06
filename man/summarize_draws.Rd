% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postestimation_functions.R
\name{summarize_draws}
\alias{summarize_draws}
\title{Summarize predictive draws}
\usage{
summarize_draws(draws)
}
\arguments{
\item{draws}{[matrix] matrix of dimensions (num obs) by (num draws)}
}
\value{
data.table with columns 'mean','median','upper','lower' (of 95% UI)
}
\description{
Summarize the mean and select quantiles of a matrix of posterior
  draws, where draws are stored in columns
}
