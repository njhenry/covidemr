
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
  missing_cols <- setdiff(
    colnames(in_data), c(num_field, denom_field, est_field, group_fields)
  )
  if(length(missing_cols) > 0){
    stop("Missing fields for RMSE calculation:" , paste(missing_cols, collapse=', '))
  }
  # If no grouping field is specified, make a dummy field
  dt_for_error <- data.table::copy(in_data)
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
  # Aggregate to get RMSE and RSE, weighted and unweighted
  error_dt <- dt_for_error[, .(
      rmse = sqrt(sum((estfld - obsfld)**2)),
      rmse_weighted = sqrt(sum(wgtfld * (estfld - obsfld)**2) / sum(wgtfld)),
      rse = sum((estfld - obsfld)**2) / sum((baseline - obsfld)**2),
      rse_weighted = (
        sum(wgtfld * (estfld - obsfld)**2) / sum(wgtfld * (baseline - obsfld)**2)
      )
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
#'
#' @return Data.table containing the following fields:
#'   - 'covg<X>': Empirical coverage for the X% uncertainty interval
#'   - Any grouping columns specified in the function arguments
#'
#' @import data.table matrixStats
#' @export
calculate_coverage <- function(
  in_data, num_field, denom_field, draw_fields, coverage_levels = c(.5, .8, .95),
  group_fields = NULL
){
  # Ensure that there are no missing columns
  missing_cols <- setdiff(
    colnames(in_data), c(num_field, denom_field, draw_fields, group_fields)
  )
  if(length(missing_cols) > 0){
    stop("Missing fields for coverage:" , paste(missing_cols, collapse=', '))
  }
  # If no grouping field is specified, make a dummy field
  dt_for_covg <- data.table::copy(in_data)
  dt_for_covg[, obsfld := get(num_field) / get(denom_field) ]
  if(length(group_fields) == 0){
    dt_for_covg[, agg_dummy := 1 ]
    group_by <- 'agg_dummy'
  } else {
    group_by <- group_fields
  }
  # Get coverage by observation at all requested UI levels
  covg_fields <- character(0)
  for(this_lev in coverage_levels){
    # Example coverage field name: 95% UI -> covg0950
    this_lev_field <- gsub('\\.', '', sprintf('covg%.3f', this_lev))
    covg_fields <- c(covg_fields, this_lev_field)
    # Example UI: 95% -> c(0.025, 0.975)
    alpha <- (1 - this_lev)/2
    dt_for_covg$c_lower <- matrixStats::rowQuantiles(
      dt_for_covg[, ..draw_fields], probs = alpha
    )
    dt_for_covg$c_upper <- matrixStats::rowQuantiles(
      dt_for_covg[, ..draw_fields], probs = 1 - alpha
    )
    # Check whether each observation is in the specified interval
    dt_for_covg[,
      (this_lev_field) := as.integer((obsfld > c_lower) & (obsfld < c_upper))
    ]
    dt_for_covg[, c('c_lower', 'c_upper') := NULL ]
  }
  # Aggregate across observations
  coverage_dt <- dt_for_covg[, lapply(.SD, mean), .SDcols=covg_fields, by=group_by]
  # Clean up and return
  if(length(group_fields) == 0) coverage_dt[, agg_dummy := NULL ]
  return(coverage_dt)
}
