#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
#' @import matrixcalc
#' @import Matrix
#' @export
rmvnorm_prec <- function(mu, prec, n.sims) {
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Matrix::Cholesky(prec, super=TRUE)
  return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z)))
}


#' Get weekly harmonic frequency
#'
#' @description Get the constant that translates from the unit circle to a
#'   52-week harmonic cycle. Used for fourier analysis of seasonality.
#'
#' @return A float constant
#'
#' @export
weekly_harmonic_freq <- function(){
  return(2. * pi / 52.)
}


#' Generate weekly fluctuations in mortality based on fit harmonics
#'
#' @description This function takes fit coefficients for `sin(.)` and `cos(.)`
#'   from fourier analysis in TMB and returns the fitted estimates of
#'   seasonality for each week.
#'
#' @param fourier_coefs A vector of length 2 * N, where N is the level of
#'   harmonics used to fit seasonality. The first of each pair of terms is the
#'   coefficient on `sin(f(x))`, and the second is the coeffient on `cos(f(x))`,
#'   where `f(x) = (level) * (yearly harmonic frequency) * (week)`
#'
#' @return Vector of length 52 (corrected number of weeks in a modeled year)
#'   with fit seasonality estimate for each week
#'
#' @export
get_fourier_seasonality_fit <- function(fourier_coefs){
  # Input data validation
  if(length(fourier_coefs) == 0 | class(fourier_coefs) != "numeric"){
    stop("Fourier coefficients must be a numeric vector with length > 0")
  }
  if(length(fourier_coefs) %% 2 != 0){
    stop("Length of coefficients must be a multiple of 2")
  }

  # Set week IDs and output vector
  week_ids <- 1:52
  seasonal_fit <- rep(0, length(week_ids))
  harmonics_level <- length(fourier_coefs) / 2

  # Add each level of harmonics
  for(lev in 1:harmonics_level){
    seasonal_fit = (
      seasonal_fit +
      fourier_coefs[2 * lev - 1] * sin(lev * week_ids * weekly_harmonic_freq()) +
      fourier_coefs[2 * lev] * cos(lev * week_ids * weekly_harmonic_freq())
    )
  }

  return(seasonal_fit)
}


#' Create post-estimation predictive draws
#'
#' @description Given the outputs from a fitted TMB model object, create
#'   an object with posterior predictive draws for all groupings specified by a
#'   template data.table
#'
#' @param tmb_sdreport output of `TMB::sdreport()` on the fitted model object.
#'   Should include a joint precision matrix (by specifying
#'   `getJointPrecision = TRUE` in the call to `sdreport()`). This object will be
#'   parsed to check for fixed effects, random effects, and the Fourier time
#'   series terms.
#' @param num_draws [int] How many posterior predictive samples to take?
#' @param covariate_names [char] All covariate field names, including 'intercept'
#' @param template_dt [data.table] table containing at least the following fields:
#'   - idx_age: Zero-indexed age grouping
#'   - idx_week: Zero-indexed week grouping (typically week - 1)
#'   - idx_year: Zero-indexed year grouping (typically year - min(year))
#'   - idx_loc: Zero-indexed location grouping
#'   - idx_fourier: Indexer for fourier analysis
#' @param fourier_harmonics_level [int, default NULL] Number of levels used to
#'   fit seasonality in the model. This parameter will be ignored if there are
#'   no `Z_fourier` parameters in the fitted output
#'
#' @return A named list with two items:
#'   - XX TODO
#'
#' @import tictoc
#' @import data.table
#' @export
generate_stwa_draws <- function(
  tmb_sdreport, num_draws, covariate_names, template_dt, fourier_harmonics_level
){
  # Copy input data
  templ <- data.table::copy(template_dt)

  # Get parameter names
  mu <- c(tmb_sdreport$par.fixed, tmb_sdreport$par.random)
  parnames <- names(mu)

  ## Input data checks
  # Check that joint precision matrix was created
  if(!"jointPrecision" %in% names(tmb_sdreport)) stop("Missing joint precision matrix")
  # Check that all the covariate names are there
  missing_covs <- setdiff(covariate_names, names(templ))
  if(length(missing_covs) > 0){
    stop("Missing covariates: ", paste0(missing_covs, collapse=', '))
  }
  if(length(covariate_names) != sum(parnames == 'beta_covs')){
    stop("Wrong number of covariates in model fit")
  }
  # Check that template data.table has all required columns
  template_req_cols <- paste0('idx_', c('loc','year','week','age','fourier'))
  missing_templ_cols <- setdiff(template_req_cols, names(templ))
  if(length(missing_templ_cols) > 0){
    stop("Missing columns: ", paste0(missing_templ_cols, collapse=', '))
  }

  ## Get parameter draws
  message(sprintf(" - Generating %i parameter draws...", num_draws))
  tictoc::tic(" - Parameter draw generation")
  param_draws <- rmvnorm_prec(
    mu = mu,
    prec = tmb_sdreport$jointPrecision,
    n.sims = num_draws
  )
  rownames(param_draws) <- parnames
  tictoc::toc()

  ## Generate predictive draws from parameter draws
  # Sort by id fields beforehand to add random effects more easily
  templ[, row_id := .I ]
  templ <- templ[order(idx_age, idx_week, idx_year, idx_loc)]
  templ[, sorted_row_id := .I ]

  # Prediction = exp( Covar FEs + age FEs + REs + sometimes seasonality )
  # Covariate fixed effects
  cov_fes <- as.matrix(template_dt[, ..use_covs]) %*% param_draws[parnames=='beta_covs', ]
  # Age fixed effects
  age_fes <- rbind(0, param_draws[parnames=='beta_ages', ])[ templ$idx_age + 1 ]

  # Random effects -- structure varies
  res <- matrix(0., nrow=nrow(templ), ncol=num_draws)

  if(any(parnames=='Z_stwa')){
    message(" - Adding joint location-year-week-age random effect")
    if(nrow(res) != nrow(param_draws["Z_stwa", ])) stop("Z_stwa dims issue")
    res <- res + param_draws[parnames=="Z_stwa", ]
  }

  if(any(parnames=='Z_sta')){
    message(" - Adding joint location-year-age random effect")
    n_ages <- max(templ$idx_age) + 1
    n_weeks <- max(templ$idx_week) + 1
    n_years <- max(templ$idx_year) + 1
    n_locs <- max(templ$idx_loc) + 1
    if(sum(parnames == 'Z_sta') != (n_ages * n_years * n_locs)){
      stop("Z_sta dims issue")
    }
    # Repeat year-location random effects (# weeks) times within each age group
    # This replicates the age-week-year-location ordering of the template DT
    z_sta_idx <- (
      rep(1:(n_years * n_locs), times = n_ages * n_weeks) +
      rep(((1:n_ages)-1) * (n_years*n_locs), each = n_years * n_locs * n_weeks)
    )
    res <- res + param_draws[parnames=='Z_sta',][z_sta_idx,]
  }

  if(any(parnames=='Z_fourier')){
    message("  - Adding seasonality effect from fourier analysis")
    if(is.null(fourier_harmonics_level)) stop("Harmonics level cannot be NULL")
    f_ncol <- fourier_harmonics_level * 2
    if(sum(parnames=='Z_fourier') %% (f_ncol) != 0){
      stop("Number of fourier coefficients is not evenly divisible by harmonics level")
    }
    num_f_groups <- sum(parnames=='Z_fourier') / f_ncol
    # Create a matrix of size (num weeks * fourier groups) by (num draws)
    z_fourier <- Reduce("rbind", lapply(1:num_f_groups, function(f_group){
      terms_idx <- 0:(f_ncol - 1) * num_f_groups + f_group
      # Return (52) fit weekly values by draw
      draws_list <- lapply(
        1:num_draws, function(dd){
          get_fourier_seasonality_fit(
            param_draws[parnames=='Z_fourier',][terms_idx, dd]
          )
        })
      return(Reduce('cbind', draws_list))
    }))
    # Create the associated lookup table
    fourier_lookup <- CJ(idx_fourier = (1:num_f_groups) - 1, idx_week = 0:51)
    fourier_lookup[, f_row := .I ]
    templ[fourier_lookup, on=c('idx_fourier', 'idx_week'), f_row := f_row ]
    # Add seasonality random effect
    res <- res + z_fourier[templ$f_row, ]
  }

  preds <- exp(cov_fes + age_fes + res)
  # Reorder by the original row ordering
  templ <- templ[order(row_id)]
  preds <- preds[templ$sorted_row_id, ]
  # Return parameter draws and predictive draws
  return(list(param_draws = param_draws, predictive_draws = preds))
}


#' Summarize predictive draws
#'
#' @description Summarize the mean and select quantiles of a matrix of posterior
#'   draws, where draws are stored in columns
#'
#' @param draws [matrix] matrix of dimensions (num obs) by (num draws)
#'
#' @return data.table with columns 'mean','median','upper','lower' (of 95% UI)
#'
#' @import matrixStats
#' @export
summarize_draws <- function(draws){
  summs <- cbind(
    rowMeans(draws), matrixStats::rowQuantiles(draws, probs=c(0.5, 0.025, 0.975))
  )
  colnames(summs) <- c('mean','median','lower','upper')
  return(as.data.table(summs))
}