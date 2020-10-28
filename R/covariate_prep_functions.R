
#' Extend a covariate time series forward
#'
#' @description Helper function to extend a covariate time series to the last
#'   year modeled by copying forward the most recent year of data available.
#'
#' @param covar_data data.table containing covariates with a 'year' field
#' @param model_year Vector of years to be included in modeling
#'
#' @return covariate data.table, with years extended to the most recent year of
#'   data
#'
#' @import data.table
#' @export
extend_covar_time_series <- function(covar_data, model_years){
  # Validate inputs: ensure that the year field exists
  if(!('year' %in% names(covar_data))) stop("'year' field missing from data")

  missing_years <- setdiff(model_years, covar_data$year)
  # If no years are missing, return the original data.table
  if(length(missing_years) == 0) return(covar_data)

  # Extend forward
  most_recent_year <- covar_data[year==max(year),]
  full_data <- rbindlist(c(
      list(covar_data),
      lapply(missing_years, function(yr) copy(most_recent_year)[, year := yr])
  ))[order(year)]

  message("    Extended forward for years ",paste(missing_years, collapse=', '))
  return(full_data)
}


#' Validate prepared covariate data
#'
#' @description Given prepared covariate data, validate that the correct fields
#'   are present and that no data is missing
#'
#' @param prepped_covar data.table containing the prepped covariate
#' @param covar_name Short name for the covariate. This will be used to send the
#'   covariate to a particular prep function.
#' @param covar_indices Vector of identifiers that should be used to merge onto
#'   the modeling dataset. The prepared covariate dataset should only include
#'   identifier columns and the covariate value, specified by the covariate name
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return NULL (fails if there are any issues in prepared covariate data)
#'
#' @import data.table
check_covar_validity <- function(
  prepped_covar, covar_name, covar_indices, model_years, location_table
){
  # Check that both exist
  if(is.null(prepped_covar)) stop("Missing prepared covariate")
  if(is.null(covar_indices)) stop("Missing prepared covariate indices")
  # Check for correct columns
  all_colnames <- c(covar_name, covar_indices)
  missing_cols <- setdiff(all_colnames, names(prepped_covar))
  extra_cols <- setdiff(names(prepped_covar), all_colnames)
  if(length(c(missing_cols, extra_cols)) > 0){
    if(missing_cols) message("Missing columns: ", paste(missing_cols, collapse=', '))
    if(extra_cols) message("Extra columns: ", paste(extra_cols, collapse=', '))
    stop("Resolve column issues before proceeding")
  }
  # If an index exists, check that all needed values are included
  # Extra values are fine but will be dropped when merging on the data
  # Helper function to check for missing values
  check_missing <- function(idx_col, all_vals){
    any_missing <- FALSE
    if(idx_col %in% covar_indices){
      missing_vals <- setdiff(all_vals, unique(prepped_covar[[idx_col]]))
      if(length(missing_vals) > 0){
        message("Missing ", idx_col, ": ", paste0(missing_vals, collapse=', '))
        any_missing <- TRUE
      }
    }
    # Return TRUE if there is an error, and FALSE otherwise
    return(any_missing)
  }
  any_missing <- c(
    check_missing('year', model_years),
    check_missing('sex', c('male','female')),
    check_missing('week', 1:52),
    check_missing('location_code', location_table$location_code)
  )
  if(any(any_missing)) stop("Resolve missingness errors to proceed.")

  # Validation passed!
  invisible(NULL)
}


#' Rescale observed covariates
#'
#' @description Rescale covariates so that each covariate has mean zero and
#'   variance 1 when measured at all observations in the dataset. Optionally,
#'   this function allows for rescaling based only on observations in the
#'   training subset of the data, as denoted by a field in the dataset
#'
#' @param input_data Data.table of the full dataset containing observations
#'   along with covariate data
#' @param covar_names [char] vector of covariate names. All of these should be
#'   included as fields in the `input_data`.
#' @param subset_field [char, optional] Fill this option to rescale based on
#'   only observations from the training data (not the test data). This
#'   argument gives the field to subset on
#' @param subset_field_values [optional] Fill this option to rescale based on
#'   only observations from the training data (not the test data). This argument
#'   gives the values of `subset_field` indicating that an observation was in
#'   the training dataset and should be used for rescaling
#'
#' @return Named list containing two data.tables:
#'     - 'data_rescaled': Prepared dataset with rescaled covariate values
#'     - 'covariate_scaling_factors': Mean and SD of the original dataset (use
#'          to back-transform the scaled covariates)
#'
#' @import data.table
#' @export
rescale_prepped_covariates <- function(
  input_data, covar_names, subset_field = NULL, subset_field_values = NULL
){
  # Check that all covariate names are in the dataset
  missing_covars <- covar_names[!covar_names %in% names(input_data)]
  if(length(missing_covars) > 0){
    stop("Missing covariates: ", paste(missing_covars, collapse=', '))
  }

  data_for_scale <- copy(input_data)
  data_rescaled <- copy(input_data)

  # Optionally take scaling factors from a subset of the data
  if(!is.null(subset_field) & (length(subset_field_values) > 0)){
    message("Subsetting data based on values of ", subset_field, "...")
    old_nrow <- nrow(data_for_scale)
    data_for_scale <- data_for_scale[
      data_for_scale[[subset_field]] %in% subset_field_values,
    ]
    new_nrow <- nrow(data_for_scale)
    message(paste("  - Subset from", old_nrow, "to", new_nrow, "rows."))
  }
  if(nrow(data_for_scale) == 0) stop("Scaling dataset has no rows.")

  # Create a list that will store scaling factors
  scaling_factors_list <- list()
  # Iteratively scale covariates
  for(covar_name in covar_names){
    if(data_for_scale[, uniqueN(get(covar_name)) ] < 2){
      stop("Not enough unique values to get rescaling factor for ", covar_name)
    }
    cov_mean <- mean(data_for_scale[[covar_name]], na.rm = TRUE)
    cov_sd <- sd(data_for_scale[[covar_name]], na.rm = TRUE)
    if(cov_sd < 1E-4){
      message("WARNING: scaling factor for ",covar_name," will be large (>1E4)")
    }
    # Rescale
    data_rescaled[, (covar_name) := (get(covar_name) - cov_mean) / cov_sd ]
    # Add rescaling factors to list
    scaling_factors_list[[covar_name]] <- list(
      cov_name = covar_name, cov_mean = cov_mean, cov_sd = cov_sd
    )
  }
  # Compile the full list of rescaling factors
  covariate_scaling_factors <- rbindlist(scaling_factors_list)

  # Return dataset with rescaled covariates and covariate scaling factors
  return(list(
    data_rescaled = data_rescaled,
    covariate_scaling_factors = covariate_scaling_factors
  ))
}
