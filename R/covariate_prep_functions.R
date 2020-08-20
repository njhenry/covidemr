
#' Prepare an Italy covariate for modeling
#'
#' @description Generic function to prepare a covariate pulled from IStat data.
#'   This function is a wrapper for sub-functions called for specific covariates
#'   and will fail if a non-initialized covariate has been passed. It returns
#'   both the prepared covariate and a vector of expected indices that can be
#'   used to merge back onto the original dataset.
#'
#' @param covar_data data.table containing the covariate
#' @param covar_name Short name for the covariate. This will be used to send the
#'   covariate to a particular prep function.
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table
#' @export
ita_prepare_covariate <- function(
  covar_data, covar_name, model_years, location_table
){
  # Check that an appropriate covariate function exists
  covar_func <- get_covar_prep_function(covar_name)
  # Apply the function to the data
  message(paste(" - Preparing",covar_name))
  covar_out_list <- do.call(
    covar_func,
    args = list(
      covar_data = covar_data,
      model_years = model_years,
      location_table = location_table
    )
  )
  # Check that output is valid
  check_covar_validity(
    prepped_covar = covar_out_list$prepped_covar,
    covar_name = covar_name,
    covar_indices = covar_out_list$covar_indices,
    model_years = model_years,
    location_table = location_table
  )
  # Return output
  return(covar_out_list)
}


#' Pull a specific prep function for an Italy covariate
#'
#' @description Helper function for \code{\link{prepare_covariate}} that checks
#'   whether a covariate-specific prep function exists yet by searching the
#'   package namespace. All covariate-specific prep functions should take the
#'   name format "ita_prep_covar_<COVAR NAME>". Errors if the function does not
#'   exist; returns the function name if the function does exist.
#'
#' @param covar_name Short name for the covariate.
#'
#' @return Name of the function to prepare this covariate
#'
get_covar_prep_function <- function(covar_name){
  # Search for prep functions
  f_pattern <- 'ita_prepare_covar_'
  potential_f_name <- paste0(f_pattern, covar_name)
  all_prep_functions <- ls('package:covidemr', pattern=sprintf("^%s",f_pattern))
  if(potential_f_name %in% all_prep_functions){
    return(potential_f_name)
  } else {
    stop(sprintf(
      "The function %s does not exist! Valid options include:\n - %s",
      potential_f_name,
      paste(all_prep_functions, collapse='\n - ')
    ))
  }
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

  message("Extended forward for years ",paste(missing_years, collapse=', '))
  return(full_data)
}


#' Prepare TFR covariate
#'
#' @description Covariate-specific prep function for TFR (total fertility rate).
#'   This is a convenience function for specific use with datasets exported from
#'   IStat.
#'
#' @param covar_data data.table containing the raw covariate data
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_tfr <- function(covar_data, model_years, location_table){

  setnames(covar_data, c('ITTER107','TIME','Value'), c('icode', 'year', 'tfr'))

  # Merge on location codes
  covar_data_merged <- merge(
    covar_data,
    location_table[, .(icode, location_code)],
    by='icode'
  )
  # Drop unnecessary columns
  covar_indices <- c('location_code', 'year')
  all_cols <- c(covar_indices, 'tfr')
  covar_data_merged <- covar_data_merged[, ..all_cols ]

  # Extend years
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare unemployment covariate
#'
#' @description Covariate-specific prep function for unemployment by sex and by
#'   quarter. This is a convenience function for specific use with datasets
#'   exported from IStat.
#'
#' @param covar_data data.table containing the raw covariate data
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_unemp <- function(covar_data, model_years, location_table){

  covar_indices <- c('location_code','sex','year')

  # Update names
  setnames(covar_data, c('ITTER107','TIME','Value'), c('icode', 'year', 'unemp'))

  # Keep only annual values and convert to integer
  suppressWarnings(covar_data[, year := as.integer(year) ])
  covar_data <- covar_data[ !is.na(year), ]

  # Format sex identifiers
  covar_data[, sex := gsub('s', '', Gender) ]

  # Fix Bolzano and Trento location codes, then merge on standard code table
  covar_data[ icode == 'ITD1', icode := 'ITD10' ] # Bolzano
  covar_data[ icode == 'ITD2', icode := 'ITD20' ] # Trento
  covar_data_merged <- merge(
    covar_data,
    location_table[, .(location_code, icode)],
    by='icode'
  )

  # Subset columns and extend to 2020 and return
  covar_data_merged <- covar_data_merged[, c(covar_indices, 'unemp'), with=FALSE]
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare social services covariate
#'
#' @description Covariate-specific prep function for social services coverage
#'   for families. This is a convenience function for specific use with datasets
#'   exported from IStat.
#'
#' @param covar_data data.table containing the raw covariate data
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_socserv <- function(covar_data, model_years, location_table){

  covar_indices <- c('year', 'location_code')

  # Update names
  setnames(covar_data, c('ITTER107', 'TIME', 'Value'), c('icode', 'year', 'socserv'))
  # Subset to just proportion of families/children receiving at home social services
  covar_data <- covar_data[(TIPUTENZA1=='FAM') & (TIPSERVSOC=='HOMECARE'), ]

  # Fix for Sud Sardegna, which is represented by its old defunct provinces
  covar_ssd <- covar_data[ icode %in% c('ITG2B', 'ITG2C'), .(icode, year, socserv)]
  covar_ssd[, icode := 'IT111' ]
  covar_ssd <- covar_ssd[, .(socserv=mean(socserv)), by=.(icode, year)]
  covar_data <- rbindlist(list(covar_data, covar_ssd), use.names=TRUE, fill=TRUE)

  # Merge on location codes
  covar_data_merged <- merge(
    covar_data,
    location_table[, .(icode, location_code)],
    by='icode'
  )

  covar_data_merged <- covar_data_merged[, c(covar_indices, 'socserv'), with=FALSE]
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}
