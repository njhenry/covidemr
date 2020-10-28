
#' Get Variance Inflation Factors for a set of covariates
#'
#' @description Get Variance Inflation Factors (VIFs), indicating
#'   multicollinearity associated with each covariate in a dataset. A VIF value
#'   over ~5 indicates high collinearity and might warrant exclusion from a
#'   dataset
#'
#' @param in_data Input dataset, in data.table format
#' @param covar_names [char] vector of covariate names, each of which should
#'   appear as columns in the dataset
#'
#' @return data.table with two fields: 'covar_name' (representing the covariate
#'   name) and 'vif' (giving the variance inflation factor for each covariate).
#'
#' @import car data.table stats
#' @export
get_covar_vif <- function(in_data, covar_names){
  # Construct a simple model where all covariates predict a dummy outcome
  missing_covars <- covar_names[ !covar_names %in% names(in_data) ]
  if(length(missing_covars) > 0){
    stop("Missing covars: ", paste(missing_covars, collapse=', '))
  }
  training_data <- copy(in_data)[, ..covar_names][, outcome := 1 ]
  dummy_model <- stats::lm(outcome ~ . , data = training_data)

  # Run VIF analysis
  vif_results <- car::vif(dummy_model)

  # Format as data.table and return
  return(data.table(covar_name = names(vif_results), vif = unname(vif_results)))
}
