
#' Calculate RMSE and RSE for predictive estimates
#'
#' @description Given a dataset, calculate Root Mean Squared Error and Relative
#'  Squared Error (both weighted and unweighted) for predictions across a
#'  dataset. The observed data should come in the form of separate numerator and
#'  denominator columns, while the estimates should be a single field estimating
#'  a rate (i.e. already normalized by denominator)
#'
#' @param in_data Input data.table
#' @param num_field Numerator field for the observed data
#' @param denom_field Denominator field for the observed data
#' @param est_field Estimator field, presented as a rate (num/denom)
#' @param group_fields [optional, default NULL] If the predictive validity
#'   metrics should be grouped, list the fields to group them by here. If NULL
#'   (the default), the predictive validity metrics will be calculated across
#'   the entire dataset
#'
#' @return A data.table with the following fields:
#'   - 'rmse': Root mean squared error, unweighted
#'   - 'rmse_weighted': Root mean squared error, weighted by population size
#'   - 'rse': Relative squared error, unweighted
#'   - 'rse_weighted': Relative squared error, weighted by population size
#'   - Any grouping columns specified in the function arguments
#'
#' @import data.table
#' @export
calculate_rmse_rse <- function(
  in_data, num_field, denom_field, est_field, group_fields = NULL
){
  # Ensure that there are no missing columns
  reqd_cols <- c(num_field, denom_field, est_field, group_fields)
  missing_cols <- setdiff(reqd_cols, colnames(in_data))
  if(length(missing_cols) > 0){
    stop("Missing fields for RMSE calculation:" , paste(missing_cols, collapse=', '))
  }
  # If no grouping field is specified, make a dummy field
  oos_cols <- c('year', 'age_group_code', 'location_code', 'week')
  keep_cols <- intersect(colnames(in_data), c(reqd_cols, oos_cols))
  dt_for_error <- data.table::copy(in_data[, ..keep_cols])
  if(length(group_fields) == 0){
    dt_for_error[, agg_dummy := 1 ]
    group_by <- 'agg_dummy'
  } else {
    group_by <- group_fields
  }
  # Set all columns needed for error calculation
  dt_for_error[, obsfld := get(num_field) / get(denom_field) ]
  dt_for_error[, estfld := get(est_field) ]
  dt_for_error[, wgtfld := get(denom_field)/sum(get(denom_field)), by = group_by]
  # Get the mean across observations as a baseline for RSE
  dt_for_error[, baseline := sum(get(num_field))/sum(get(denom_field)), by = group_by]
  # Calculate an alternate baseline value, which is mean across observations
  #  by age group only
  if('age_group_code' %in% colnames(dt_for_error)){
    dt_for_error[, bl_by_age := sum(get(num_field))/sum(get(denom_field)), by = age_group_code]
  } else {
    dt_for_error[, bl_by_age := as.numeric(NA) ]
  }
  # Calculate an OOS baseline, which is the mean observation across all other
  #  years for a particular location, week, and age group
  compare_years <- na.omit(unique(dt_for_error$year))
  if(all(oos_cols %in% colnames(dt_for_error)) & (length(compare_years) > 0)){
    oos_grps <- setdiff(oos_cols, 'year')
    oos_compare_dt <- rbindlist(lapply(compare_years, function(oosyr){
      oosdt <- copy(dt_for_error[
        year != oosyr,
        bl_oos := sum(get(num_field))/sum(get(denom_field)),
        by = oos_grps
      ])
      oosdt[, year := oosyr ]
      return(oosdt)
    }))
    # Merge out-of-sample estimates back onto the original value
    dt_for_error[oos_compare_dt, bl_oos := i.bl_oos, on = oos_cols]
  } else {
    dt_for_error[, bl_oos := as.numeric(NA) ]
  }
  # Aggregate to get RMSE and RSE, weighted and unweighted
  dt_for_error[, sq_resid := (estfld - obsfld)**2]
  error_dt <- dt_for_error[, .(
      rmse = sqrt(sum(sq_resid)/.N),
      rmse_weighted = sqrt(sum(wgtfld * sq_resid) / sum(wgtfld)),
      rse = sum(sq_resid) / sum((baseline - obsfld)**2),
      rse_weighted = (sum(wgtfld*sq_resid) / sum(wgtfld*(baseline-obsfld)**2)),
      rse_by_age = sum(sq_resid) / sum((bl_by_age-obsfld)**2),
      rse_by_age_weighted = (sum(wgtfld*sq_resid) / sum(wgtfld*(bl_by_age-obsfld)**2)),
      rse_vs_oos = sum(sq_resid) / sum((bl_oos-obsfld)**2),
      rse_vs_oos_weighted = (sum(wgtfld*sq_resid) / sum(wgtfld*(bl_oos-obsfld)**2))
    ),
    by = group_by
  ]
  # Clean up and return
  if(length(group_fields) == 0) error_dt[, agg_dummy := NULL ]
  return(error_dt)
}


#' Calculate coverage across predictive draws
#'
#' @description Calculate how often observed data falls within X% uncertainty
#'   intervals of the posterior predictive draws. The observed data should come
#'  in the form of separate numerator and denominator columns, while the
#'  estimates should be a single field estimating a rate (i.e. already
#'  normalized by denominator)
#' @param in_data Input data.table
#' @param num_field Numerator field for the observed data
#' @param denom_field Denominator field for the observed data
#' @param draw_fields Character vector of fields containing predictive draws, in
#'   rate space (e.g. mortality rates)
#' @param coverage_levels [optional, default c(.5, .8, .95)] Uncertainty
#'   intervals to calculate from the posterior predictive draws
#' @param group_fields [optional, default NULL] If the predictive validity
#'   metrics should be grouped, list the fields to group them by here. If NULL
#'   (the default), the predictive validity metrics will be calculated across
#'   the entire dataset
#' @param pois_sim [optional, default TRUE] Should the coverage be estimated
#'   using realizations Poisson distribution centered at population * rate (as
#'   opposed to the central value, population * p)?
#'
#' @return Data.table containing the following fields:
#'   - 'covg<X>': Empirical coverage for the X% uncertainty interval
#'   - Any grouping columns specified in the function arguments
#'
#' @import data.table matrixStats
#' @export
calculate_coverage <- function(
  in_data, num_field, denom_field, draw_fields, coverage_levels = c(.5, .8, .95),
  group_fields = NULL, pois_sim = TRUE
){
  dt_for_covg <- data.table::copy(in_data)
  # Ensure that there are no missing columns
  missing_cols <- setdiff(
    c(num_field, denom_field, draw_fields, group_fields), colnames(in_data)
  )
  if(length(missing_cols) > 0){
    stop("Missing fields for coverage:" , paste(missing_cols, collapse=', '))
  }
  # If no grouping field is specified, make a dummy field
  if(length(group_fields) == 0){
    dt_for_covg[, agg_dummy := 1 ]
    group_by <- 'agg_dummy'
  } else {
    group_by <- group_fields
  }

  ## RUN DIFFERENT COVERAGE TESTS DEPENDING ON POIS_SIM
  if(pois_sim){
    ## CASE: Using Poisson simulation
    dt_for_covg[, obsfld := get(num_field) ]
    # Get all draws and denominators matrices with the same dimensions
    ndraws <- length(draw_fields)
    draws_mat <- as.matrix(dt_for_covg[, ..draw_fields])
    denoms_mat <- matrix(rep(dt_for_covg[[denom_field]], ndraws), ncol=ndraws)
    pois_samples <- matrix(
      rpois(length(draws_mat), denoms_mat * draws_mat),
      ncol=ndraws
    )
    dt_for_covg[, draw_fields] <- as.data.table(pois_samples)
  } else {
    ## CASE: Using probabilities directly
    # Get coverage by observation at all requested UI levels
    dt_for_covg[, obsfld := get(num_field) / get(denom_field) ]
  }

  # Compare observation field to draw fields
  covg_fields <- character(0)
  for(this_lev in coverage_levels){
    # Example coverage field name: 95% UI -> covg0950
    this_lev_field <- gsub('\\.', '', sprintf('covg%.3f', this_lev))
    covg_fields <- c(covg_fields, this_lev_field)
    # Example UI: 95% -> c(0.025, 0.975)
    alpha <- (1 - this_lev)/2
    dt_for_covg$c_lower <- matrixStats::rowQuantiles(
      as.matrix(dt_for_covg[, ..draw_fields]), probs = alpha
    )
    dt_for_covg$c_upper <- matrixStats::rowQuantiles(
      as.matrix(dt_for_covg[, ..draw_fields]), probs = 1 - alpha
    )
    # Check whether each observation is in the specified interval
    dt_for_covg[,
      (this_lev_field) := as.integer((obsfld >= c_lower) & (obsfld <= c_upper))
    ]
    dt_for_covg[, c('c_lower', 'c_upper') := NULL ]
  }

  # Aggregate across observations
  coverage_dt <- dt_for_covg[, lapply(.SD, mean), .SDcols=covg_fields, by=group_by]
  # Clean up and return
  if(length(group_fields) == 0) coverage_dt[, agg_dummy := NULL ]
  return(coverage_dt)
}
