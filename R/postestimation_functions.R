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
#' @param keep_params [char] Vector of parameter names to keep when generating
#'   draws from the precision matrix. Keep only parameters needed for prediction
#' @param num_draws [int] How many posterior predictive samples to take?
#' @param covariate_names [char] All covariate field names, including 'intercept'
#' @param template_dt [data.table] table containing at least the following fields:
#'    - idx_age: Zero-indexed age grouping
#'    - idx_week: Zero-indexed week grouping (typically week - 1)
#'    - idx_year: Zero-indexed year grouping (typically year - min(year))
#'    - idx_loc: Zero-indexed location grouping
#'    - idx_fourier: Indexer for fourier analysis
#' @param rescale_covars [bool, default TRUE] Should covariates in the template
#'   data.table be rescaled based on a set of normalizing factors?
#' @param covar_scaling_factors [default NULL] If `rescale_covars` is TRUE, this
#'   argument should be set to a data.table with three fields:
#'    - 'cov_name': The name of each covariate in `covariate_names`
#'    - 'cov_mean': The mean that will be subtracted for normalization
#'    - 'cov_sd': The standard deviation that will be divided for normalization
#' @param fourier_harmonics_level [int, default NULL] Number of levels used to
#'   fit seasonality in the model. This parameter will be ignored if there are
#'   no `Z_fourier` parameters in the fitted output
#'
#' @return A named list with three items:
#'    - 'param_names': Vector of parameter names in the order they have been
#'         extracted
#'    - 'param_draws': Matrix of parameter draws
#'    - 'pred_draws': Matrix of mortality predictive draws, taken at the
#'         observation points specified in the `template_dt`
#'
#' @import data.table
#' @export
generate_stwa_draws <- function(
  tmb_sdreport, keep_params, num_draws, covariate_names, template_dt,
  rescale_covars = FALSE, covar_scaling_factors = NULL,
  fourier_harmonics_level = NULL
){
  # Copy input data
  templ <- data.table::copy(template_dt)

  # Get parameter names
  mu <- c(tmb_sdreport$par.fixed, tmb_sdreport$par.random)
  param_subset <- which(names(mu) %in% keep_params)
  mu <- mu[ param_subset ]
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
  prec_mat <- tmb_sdreport$jointPrecision
  prec_subset <- which(colnames(prec_mat) %in% keep_params)
  prec_mat <- prec_mat[ prec_subset, prec_subset ]
  if(any(colnames(prec_mat) != parnames )) stop("Issue with parameter ordering")
  param_draws <- rmvnorm_prec(
    mu = mu,
    prec = prec_mat,
    n.sims = num_draws
  )
  rownames(param_draws) <- parnames

  ## Generate predictive draws from parameter draws
  # Sort by id fields beforehand to add random effects more easily
  templ[, row_id := .I ]
  templ <- templ[order(idx_age, idx_week, idx_year, idx_loc)]
  templ[, sorted_row_id := .I ]

  # If covariates need to be rescaled, do so now
  if(rescale_covars){
    message("     - Rescaling covariates in template data.table...")
    for(this_cov in setdiff(covariate_names, 'intercept')){
      cov_mean <- covar_scaling_factors[cov_name == this_cov, cov_mean]
      cov_sd <- covar_scaling_factors[cov_name == this_cov, cov_sd]
      templ[, (this_cov) := (get(this_cov) - cov_mean) / cov_sd ]
    }
  }

  # Prediction = inverse.logit( Covar FEs + age FEs + REs + sometimes seasonality )
  # Covariate fixed effects
  cov_fes <- (
    as.matrix(templ[, ..covariate_names]) %*% param_draws[parnames=='beta_covs', ]
  )
  # Age fixed effects
  age_fes <- rbind(0, param_draws[parnames=='beta_ages', ])[ templ$idx_age + 1 ]

  # Random effects -- structure varies
  res <- matrix(0., nrow=nrow(templ), ncol=num_draws)

  if(any(parnames=='Z_sta')){
    message("     - Adding joint location-year-age random effect")
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
    # Truncate if weeks are missing from the final year
    z_sta_idx <- z_sta_idx[1:nrow(templ)]
    res <- res + param_draws[parnames=='Z_sta',][z_sta_idx,]
  }

  if(any(parnames=='Z_fourier')){
    message("     - Adding seasonality effect from fourier analysis")
    if(is.null(fourier_harmonics_level)) stop("Harmonics level cannot be NULL")
    f_ncol <- fourier_harmonics_level * 2
    if(sum(parnames=='Z_fourier') %% (f_ncol) != 0){
      stop("Number of fourier coefficients is not evenly divisible by harmonics level")
    }
    num_f_groups <- sum(parnames=='Z_fourier') / f_ncol
    # Create a matrix of size (num weeks * fourier groups) by (num draws)
    z_fourier <- Reduce(
      'rbind',
      lapply(1:num_f_groups, function(f_group){
        terms_idx <- 0:(f_ncol - 1) * num_f_groups + f_group
        # Return (52) fit weekly values by draw
        draws_list <- lapply(
          1:num_draws, function(dd){
            get_fourier_seasonality_fit(
              param_draws[parnames=='Z_fourier',][terms_idx, dd]
            )
          })
        return(matrix(unlist(draws_list), ncol=num_draws))
      })
    )
    # Create the associated lookup table
    fourier_lookup <- CJ(idx_fourier = (1:num_f_groups) - 1, idx_week = 0:51)
    fourier_lookup[, f_row := .I ]
    templ[fourier_lookup, on=c('idx_fourier', 'idx_week'), f_row := f_row ]
    # Add seasonality random effect
    res <- res + z_fourier[templ$f_row, ]
  }

  log_preds <- cov_fes + age_fes + res
  preds <- exp(log_preds)
  # Reorder by the original row ordering
  templ <- templ[order(row_id)]
  preds <- preds[templ$sorted_row_id, ]
  # Return parameter draws and predictive draws
  return(list(
    param_names = parnames,
    param_draws = param_draws,
    predictive_draws = preds
  ))
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
#' @import data.table matrixStats
#' @export
summarize_draws <- function(draws){
  if('data.table' %in% class(draws)) draws <- as.matrix(draws)
  summs <- cbind(
    rowMeans(draws), matrixStats::rowQuantiles(draws, probs=c(0.5, 0.025, 0.975))
  )
  colnames(summs) <- c('mean','median','lower','upper')
  return(as.data.table(summs))
}

#' Calculate excess mortality time series for a single population
#'
#' @description Given draws of a baseline mortality rate, a starting population,
#'   and weekly deaths, calculate excess deaths, SMRs, and estimated population
#'   by draw and over the study time period. NOTE: This is a child function to
#'   estimate these quantities for a single sub-population, and should be called
#'   repeatedly to estimate excess mortality across different groups.
#'
#' @details To calculate excess mortality in a time series, the effect of
#'   mortality changes on the population denominator must be considered. Without
#'   complete information about how many people are entering an age group by
#'   aging in or (for the youngest age group) through birth as opposed to the
#'   number of people aging out of an age group or dying each week, we make the
#'   assumption that at baseline mortality, the size of this population subgroup
#'   would have remained approximately stable. This is a reasonable assumption
#'   for most populations on the time scale of a few weeks. "Baseline mortality"
#'   here is estimated as the mean mortality *rate* calculated across all
#'   predictive draws for baseline mortality, multiplied by the population. The
#'   difference between the observed number of deaths and the expected baseline
#'   mortality will be subtracted from the population of future weeks. For
#'   example:
#'
#'     Starting population: 1,000
#'     Mortality:
#'      week | mean_baseline_rate | observed_deaths
#'      -----+--------------------+----------------
#'       1   |               0.01 |            200
#'       2   |               0.02 |            100
#'       3   |               0.01 |             10
#'
#'     Week 1:
#'      * Estimated population size = starting population, 1,000
#'      * Baseline deaths = 1,000 * 0.01 = 10
#'      * Excess deaths = 200 - 10 = 190
#'     Week 2:
#'      * Estimated population size = 1,000 - 190 = 810
#'      * Baseline deaths = 810 * 0.02 = 16.2
#'      * Excess deaths = 100 - 16.2 = 83.8
#'     Week 3:
#'      * Estimated population size = 810 - 83.8 = 726.2
#'      * Baseline deaths = 726.2 * 0.01 = 7.262
#'      * Excess deaths = 10 - 7.262 = 2.738
#'
#'  Note that in the example, there are excess deaths measured in week 3. This
#'  would not be the case for a time series analysis where the population was
#'  not adjusted for previous excess deaths.
#'
#'  An excess mortality analysis is first run comparing against the mean
#'  baseline, to estimate population over the time series, and then comparing to
#'  baseline mortality by draw to calculate SMRs and excess deaths in a way that
#'  preserves uncertainty.
#'
#' @param baseline_mat [numeric] A matrix of the baseline mortality rate by
#'   week, where each row corresponds to one week of the time series (in order)
#'   and each column corresponds to a predictive draw of the baseline mortality
#'   *rate* for that week.
#' @param starting_pop [numeric] A scalar value giving the population size at
#'   the beginning of the excess mortality analysis.
#' @param obs_deaths_vec [numeric] A vector of observed deaths for each week of
#'   the excess mortality analysis. This vector should be the same length as the
#'   number of rows in the `baseline_mat`.
#'
#' @return A list with three items:
#'   - 'pop': Vector of estimated adjusted population for each week of the time
#'       series. The first week population will always be equal to the value of
#'       `starting_pop`.
#'   - 'baseline_deaths': A numeric matrix of size (num weeks) x (num draws)
#'       giving the estimated NUMBER of baseline deaths, a count, for each week
#'       and draw.
#'   - 'smrs': A numeric matrix of size (num weeks) x (num draws) giving the
#'       estimated standardized mortality ratio associated with each week and
#'       draw.
#'   - 'excess_deaths': A numeric matrix of size (num weeks) x (num draws)
#'       giving the estimated number of excess deaths, a count, for each week
#'       and draw.
#'
#' @export
calculate_excess_time_series <- function(
  baseline_mat, starting_pop, obs_deaths_vec
){
  ## Check that input data is formatted correctly
  errmsg <- function(msg) stop("ERROR IN TIME SERIES CALCULATION: ", msg)
  if(!is.matrix(baseline_mat)) errmsg("Mortality rates should be a matrix.")
  if(nrow(baseline_mat) != length(obs_deaths_vec)){
    errmsg("Number of weeks differ across observed deaths and baseline mortality")
  }
  if(length(starting_pop) != 1) errmsg("Starting population should have length 1")

  # Get data dimensions
  ndraw <- ncol(baseline_mat)
  nweek <- nrow(baseline_mat)

  ## Adjust population for excess mortality
  # - Calculate mean expected baseline mortality across all draws
  bl_mean_vec <- rowMeans(baseline_mat)
  # - Set up default population vector
  adj_pop_vec <- rep(starting_pop, nweek)
  # - Adjust all populations after the first week
  if(nweek > 1){
    for(this_wk in 2:nweek){
      last_wk <- this_wk - 1
      excess_last_week <- (
        obs_deaths_vec[last_wk] - (bl_mean_vec[last_wk] * adj_pop_vec[last_wk])
      )
      adj_pop_vec[this_wk] <- adj_pop_vec[last_wk] - excess_last_week
    }
  }

  ## Calculate baseline deaths, SMRs, and excess deaths by week
  baseline_deaths <- baseline_mat * adj_pop_vec
  smrs <- obs_deaths_vec / baseline_deaths
  excess_deaths <- obs_deaths_vec - baseline_deaths

  ## Return as list
  return(list(
    pop = adj_pop_vec,
    baseline_deaths = baseline_deaths,
    smrs = smrs,
    excess_deaths = excess_deaths
  ))
}


#' Calculate time series of excess mortality across multiple subpopulations
#'
#' @description Given draws of baseline mortality, observed starting population,
#'   and observed deaths by week, calculate excess mortality across many
#'   sub-population groupings. Sub-populations are identified by unique sets
#'   of ID columns. This is a wrapper for `calculate_excess_time_series()`,
#'   which calculates excess-adjusted population, excess deaths, and SMRs for a
#'   single population subgroup.
#'
#' @param experiment_draw_dt data.table containing the experimental data ONLY
#'   for the time period when excess mortality is to be calculated. Should
#'   contain all of the fields specificied by the other parameters.
#' @param baseline_draw_cols [character] columns in `experiment_draw_dt` that
#'   contain the draws for predicted baseline mortality rate.
#' @param group_cols [character] columns in `experiment_draw_dt` with unique
#'   grouping identifiers. There should be no NULLs in any of these columns.
#' @param week_col [character] name of the column containing week ID
#' @param obs_death_col [character] name of the column containing observed
#'   deaths by week
#' @param pop_col [character] name of the column containing estimated population
#'   by week, not adjusted (yet) to account for excess mortality.
#'
#' @return Named list of three items:
#'   - 'baseline_deaths': Data.table of baseline death draws with identifiers
#'   - 'smrs': Data.table of SMR draws with identifiers
#'   - 'excess_deaths': Data.table of excess death draws with identifiers
#'
#' @import data.table
#' @export
calculate_excess_time_series_by_group <- function(
  experiment_draw_dt, baseline_draw_cols, group_cols, week_col='week',
  obs_death_col='deaths', pop_col='pop'
){
  # Helper function specifying error reporting in this function
  errmsg <- function(msg) stop("Error in grouped excess function: ", msg)

  ## Validate input data
  if(!data.table::is.data.table(experiment_draw_dt)) errmsg(
    "Dataset must be a data.table."
  )
  reqd_cols <- c(baseline_draw_cols, group_cols, week_col, obs_death_col, pop_col)
  missing_cols <- setdiff(reqd_cols, colnames(experiment_draw_dt))
  if(length(missing_cols) > 0) errmsg(
    paste0("Missing required cols: ", missing_cols, collapse=', ')
  )
  if(any(is.na(experiment_draw_dt[, ..group_cols]))) errmsg(
    "Some grouping columns have NAs."
  )
  # The grouping columns should not include any of the other columns
  dupe_cols <- intersect(
    group_cols,
    c(baseline_draw_cols, week_col, obs_death_col, pop_col)
  )
  if(length(dupe_cols) > 0) errmsg("Group cols and measure cols shouldn't overlap.")

  ## Create data.table idenfying unique subpopulations
  subpop_id_dt <- unique(experiment_draw_dt[, ..group_cols])
  subpop_id_dt[, group_id := .I ]
  num_groups <- nrow(subpop_id_dt)

  ## Calculate excess mortality and SMRs for each subpopulation
  out_data_names <- c('baseline_deaths', 'smrs', 'excess_deaths')
  lists_by_subpop <- lapply(
    out_data_names, function(x) vector('list', length=num_groups)
  )
  names(lists_by_subpop) <- out_data_names

  for(group_idx in 1:num_groups){
    # Subset
    idx_dt <- copy(subpop_id_dt[ group_id == group_idx, ])
    data_sub_dt <- merge(x = idx_dt, y = experiment_draw_dt, by = group_cols)
    data_sub_dt <- data_sub_dt[order(get(week_col))]
    # Check that week order is correct and not duplicated
    if(any(diff(data_sub_dt[[week_col]]) != 1)){
      message(idx_dt)
      errmsg("Issue with week ordering in the subgroup listed above.")
    }
    # Run the sub-population excess mortality calculation
    this_group_excess_list <- calculate_excess_time_series(
      baseline_mat = as.matrix(data_sub_dt[, ..baseline_draw_cols]),
      starting_pop = max(data_sub_dt[[pop_col]]),
      obs_deaths_vec = data_sub_dt[[obs_death_col]]
    )
    # Create the data.tables for each data type, with identifiers added back to
    #  the dataset
    for(out_data_name in out_data_names){
      lists_by_subpop[[out_data_name]][[group_idx]] <- cbind(
        idx_dt,
        week = data_sub_dt[[week_col]],
        deaths = data_sub_dt[[obs_death_col]],
        pop = this_group_excess_list$pop,
        this_group_excess_list[[out_data_name]]
      )
    }
  }

  ## Combine SMR and excess deaths lists into unified datasets and return
  out_list <- lapply(
    out_data_names,
    function(out_data_name){
      out_dt <- data.table::rbindlist(lists_by_subpop[[out_data_name]])
      data.table::setnames(
        out_dt, c('week', 'deaths', 'pop'), c(week_col, obs_death_col, pop_col)
      )
      return(out_dt)
    }
  )
  names(out_list) <- out_data_names
  return(out_list)
}


#' Aggregate and summarize excess mortality estimates
#'
#' @description Given draws of baseline death COUNTS (not rates) and observed
#'   death counts, calculate summary metrics of excess mortality, optionally
#'   grouping by fields in the data
#'
#' @param baseline_deaths_dt data.table containing draws of baseline death
#'   COUNTS (not rates), observed death counts, and optional identifier columns
#'   for the period when excess mortality is being calculated
#' @param aggregate [bool, default FALSE] should the data be aggregated across
#'   grouping columns?
#' @param group_cols [char, default NULL] If `aggregate` is TRUE, group by these
#'   columns when aggregating baseline and observed deaths. If aggregation is
#'   specified and this argument is left as NULL, the aggregation will occur
#'   across all rows
#' @param baseline_draw_cols [character, default c('draw1', ... 'draw1000')]
#'   columns in `baseline_deaths_dt` that contain the draws for predicted
#'   baseline death counts
#' @param obs_death_col [character, default 'deaths'] name of the column in
#'   `baseline_deaths_dt` containing observed deaths by week
#' @param alpha [numeric, default 0.05] Cutoff for two-tailed uncertainty
#'    interval. For example, the default of 0.05 will use the 95% uncertainty
#'    interval for each outcome
#'
#' @return data.table containing at least the following fields:
#'    - Any grouping columns, if specified
#'    - `obs_death_col` as specified in the function argument
#'    - 'bl_<mean/lower/upper>': Mean and UI bounds for baseline death counts
#'    - 'smr_<mean/lower/upper>': Mean and UI bounds for Standardized Mortality
#'        Ratios, calculated as observed / baseline
#'    - 'ex_d_<mean/lower/upper>': Mean and UI bounds for Excess Deaths,
#'        calculated as observed - baseline
#'    - 'sig_under1': Did the observed deaths fall BELOW the baseline deaths UI?
#'    - 'sig_over1': Did the observed deaths fall ABOVE the baseline deaths UI?
#'
#' @import data.table matrixStats
#' @export
aggregate_summarize_excess <- function(
  baseline_deaths_dt, aggregate = FALSE, group_cols = NULL,
  baseline_draw_cols = paste0('draw',1:1000), obs_death_col = 'deaths',
  alpha = 0.05
){
  ## Check that inputs are valid:
  errmsg <- function(msg) stop("ERROR in excess mortality summary: ", msg)
  # Dataset must be in data.table format
  if(!is.data.table(baseline_deaths_dt)) errmsg("Dataset not in data.table format")
  # All required fields must be present
  reqd_cols <- c(group_cols, baseline_draw_cols, obs_death_col)
  missing_cols <- setdiff(reqd_cols, colnames(baseline_deaths_dt))
  if(length(missing_cols) > 0) errmsg(
    paste0("Missing fields: ", missing_cols, collapse=', ')
  )
  # No NAs allowed in grouping columns
  for(gcol in group_cols){
    if(any(is.na(baseline_deaths_dt[[gcol]]))) errmsg(
      paste('Grouping column', gcol, 'has NA values.')
    )
  }
  # Alpha must be between zero and 1
  if(!is.numeric(alpha) | alpha < 0 | alpha > 1) errmsg('Alpha must be on [0,1].')

  ## Optionally aggregate data
  summ_dt <- data.table::copy(baseline_deaths_dt)
  setnames(summ_dt, obs_death_col, 'deaths')
  if(aggregate == TRUE){
    summ_dt <- summ_dt[,
      lapply(.SD, sum),
      .SDcols = c(baseline_draw_cols, 'deaths'),
      by = group_cols
    ]
  }

  ## Calculate summary metrics
  # Baseline death summaries
  ui_percentiles <- c(alpha/2, 1-alpha/2)
  summ_dt$bl_mean <- rowMeans(summ_dt[, ..baseline_draw_cols])
  summ_dt$bl_lower <- matrixStats::rowQuantiles(
    as.matrix(summ_dt[, ..baseline_draw_cols]), probs=ui_percentiles[1], na.rm=T
  )
  summ_dt$bl_upper <- matrixStats::rowQuantiles(
    as.matrix(summ_dt[, ..baseline_draw_cols]), probs=ui_percentiles[2], na.rm=T
  )
  # SMRs and excess deaths
  summ_dt[, `:=` (
    smr_mean = deaths / bl_mean, smr_lower = deaths / bl_upper,
    smr_upper = deaths / bl_lower, ex_d_mean = deaths - bl_mean,
    ex_d_lower = deaths - bl_upper, ex_d_upper = deaths - bl_lower,
    sig_under1 = as.integer(deaths < bl_lower),
    sig_over1 = as.integer(deaths > bl_upper)
  )]

  ## Clean up and return
  # Drop draw columns
  summ_dt[, (baseline_draw_cols) := NULL ]
  setnames(summ_dt, 'deaths', obs_death_col)

  return(summ_dt)
}
