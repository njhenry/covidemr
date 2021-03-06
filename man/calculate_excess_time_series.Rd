% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postestimation_functions.R
\name{calculate_excess_time_series}
\alias{calculate_excess_time_series}
\title{Calculate excess mortality time series for a single population}
\usage{
calculate_excess_time_series(baseline_mat, starting_pop, obs_deaths_vec)
}
\arguments{
\item{baseline_mat}{[numeric] A matrix of the baseline mortality rate by
week, where each row corresponds to one week of the time series (in order)
and each column corresponds to a predictive draw of the baseline mortality
*rate* for that week.}

\item{starting_pop}{[numeric] A scalar value giving the population size at
the beginning of the excess mortality analysis.}

\item{obs_deaths_vec}{[numeric] A vector of observed deaths for each week of
the excess mortality analysis. This vector should be the same length as the
number of rows in the `baseline_mat`.}
}
\value{
A list with three items:
  - 'pop': Vector of estimated adjusted population for each week of the time
      series. The first week population will always be equal to the value of
      `starting_pop`.
  - 'baseline_deaths': A numeric matrix of size (num weeks) x (num draws)
      giving the estimated NUMBER of baseline deaths, a count, for each week
      and draw.
  - 'smrs': A numeric matrix of size (num weeks) x (num draws) giving the
      estimated standardized mortality ratio associated with each week and
      draw.
  - 'excess_deaths': A numeric matrix of size (num weeks) x (num draws)
      giving the estimated number of excess deaths, a count, for each week
      and draw.
}
\description{
Given draws of a baseline mortality rate, a starting population,
  and weekly deaths, calculate excess deaths, SMRs, and estimated population
  by draw and over the study time period. NOTE: This is a child function to
  estimate these quantities for a single sub-population, and should be called
  repeatedly to estimate excess mortality across different groups.
}
\details{
To calculate excess mortality in a time series, the effect of
  mortality changes on the population denominator must be considered. Without
  complete information about how many people are entering an age group by
  aging in or (for the youngest age group) through birth as opposed to the
  number of people aging out of an age group or dying each week, we make the
  assumption that at baseline mortality, the size of this population subgroup
  would have remained approximately stable. This is a reasonable assumption
  for most populations on the time scale of a few weeks. "Baseline mortality"
  here is estimated as the mean mortality *rate* calculated across all
  predictive draws for baseline mortality, multiplied by the population. The
  difference between the observed number of deaths and the expected baseline
  mortality will be subtracted from the population of future weeks. For
  example:

    Starting population: 1,000
    Mortality:
     week | mean_baseline_rate | observed_deaths
     -----+--------------------+----------------
      1   |               0.01 |            200
      2   |               0.02 |            100
      3   |               0.01 |             10

    Week 1:
     * Estimated population size = starting population, 1,000
     * Baseline deaths = 1,000 * 0.01 = 10
     * Excess deaths = 200 - 10 = 190
    Week 2:
     * Estimated population size = 1,000 - 190 = 810
     * Baseline deaths = 810 * 0.02 = 16.2
     * Excess deaths = 100 - 16.2 = 83.8
    Week 3:
     * Estimated population size = 810 - 83.8 = 726.2
     * Baseline deaths = 726.2 * 0.01 = 7.262
     * Excess deaths = 10 - 7.262 = 2.738

 Note that in the example, there are excess deaths measured in week 3. This
 would not be the case for a time series analysis where the population was
 not adjusted for previous excess deaths.

 An excess mortality analysis is first run comparing against the mean
 baseline, to estimate population over the time series, and then comparing to
 baseline mortality by draw to calculate SMRs and excess deaths in a way that
 preserves uncertainty.
}
