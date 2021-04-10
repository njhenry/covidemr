
#' Generate an ICAR precision matrix based on an adjacency matrix
#'
#' @description Generate a precision matrix for the intrinsic correlated autoregressive
#'  (ICAR) model specification, a special case of the correlated autoregressive (CAR)
#'  class of Markov random field models. This precision matrix is usually denoted as "Q".
#'
#' @details The precision matrix is fully specified by the adjacency weights, matrix W,
#'   defined as W = {w_ij} where w_ij is 1 if i and j are neighbors, and 0 otherwise. The
#'   precision matrix Q is defined as Q = D_w - W, where D_w is a diagonal matrix with
#'   each diagonal term d_ii equal to the sum of row i in W.
#'
#'   Note that the ICAR model is improper, in that the conditional distributions
#'   specified by the precision matrix do not determine a full joint distribution that
#'   integrates to 1; in other words, the precision matrix Q is not invertible. The ICAR
#'   precision matrix can still be used as a prior in a hierarchical model.
#'
#'   This function includes optional argument `scale_variance`. If set to `TRUE` (the
#'   default), the function will rescale the precision matrix to have a generalized
#'   variance of 1, which may aid in prior specifications that are comparable across
#'   areal spatial models with different geometries.
#'
#'   For more details, see:
#'   Banerjee, Carlin, and Gelfand (2015). Hierarchical Modeling and Analysis for Spatial
#'     Data, 2nd Edition. Section 6.4.3.3: CAR models and their difficulties.
#'   Riebler et al. (2016). An intuitive Bayesian sptial model for disease mapping that
#'     accounts for scaling. Statistical methods in medical research, 25(4):1145-65.
#'
#' @param W Adjacency matrix, with w_ij = w_ji = 1 if areal units i and j are neighbors,
#'   and zero otherwise. See function details for more information
#' @param scale_variance [default TRUE] Should the precision matrix be rescaled so that
#'  the generalized variance is equal to 1? Setting to TRUE may help with prior
#'  specification.
#'
#' @return Sparse ICAR precision matrix Q. See function details for more information.
#'
#' @import Matrix INLA
#' @export
icar_precision_from_adjacency <- function(W, scale_variance = TRUE){
  # Generate and return sparse precision matrix
  Q <- Matrix::Diagonal(n = nrow(W), x = Matrix::rowSums(W)) - W
  if(scale_variance){
    # Scale model to have generalized variance of 1
    constraint_matrix <- matrix(1, nrow = 1, ncol = ncol(Q))
    Q <- INLA::inla.scale.model(Q, constr = list(A = constraint_matrix, e = 0))
  }
  return(Q)
}


#' Assign seasonality grouping IDs
#'
#' @description Create a vector that assigns a (zero-indexed) seasonality
#'   grouping ID based on other categories in the data. If no categories are
#'   assigned, all rows get the same seasonality grouping ID (0). If categories
#'   are based on other ID fields, then a group is assigned to every unique
#'   combination of values in those fields. For example, if grouings were
#'   assigned based on the "age" and "location" fields, and the data contained
#'   four age groupings and 10 unique locations, then 20 unique seasonality
#'   groupings would be assigned.
#'
#' @param input_data Data.table of the full input dataset. Must contain fields
#'   for any grouping categories that are assigned
#' @param grouping_fields [optional, default NULL] character vector listing any
#'   fields in the dataset that should be used to create unique seasonality
#'   groups
#'
#' @return integer vector with the same number of rows as `input_data` assigning
#'   holdout IDs.
#'
#' @import data.table
#' @export
assign_seasonality_ids <- function(input_data, grouping_fields = NULL){

  # Check for missing fields
  missing <- setdiff(grouping_fields, names(input_data))
  if(length(missing) > 0) stop("Missing fields: ", paste(missing, collapse=', '))

  # Case of no groupings: Assign all zeroes
  if(length(grouping_fields) == 0){
    return(rep(0, nrow(input_data)))
  }

  # Case of groupings: create a grouping data.table and merge it into the data
  to_cj <- lapply(grouping_fields, function(fld) sort(unique(input_data[[fld]])))
  names(to_cj) <- grouping_fields
  holdout_id_dt <- do.call('CJ', to_cj)
  # Add holdout ID, zero-indexed
  holdout_id_dt[, idx_season := .I - 1 ]
  # Merge back onto dataset
  templ_dt <- copy(input_data)[, ..grouping_fields]
  templ_dt[holdout_id_dt, idx_season := idx_season, on = grouping_fields]

  # Check that the field is the correct length and return
  if(length(templ_dt$idx_season) != nrow(input_data)) stop("Holdout merge issue")
  return(templ_dt$idx_season)
}


#' Get Maximum a Priori (MAP) parameter estimates
#'
#' @description Find the Maximum a Priori (MAP) estimates for fixed effect
#'   parameter values, including age-specific fixed effects, in a simplified
#'   version of a GLM. Starting the full TMB model with fixed effect starting
#'   values set at the MAP has been shown to improve performance and run time.
#'
#' @param in_data Input data.table, including only the data used to train the
#'    model. Fields must include a numerator, a denominator, fields named for
#'    all covariates, and (optionally) an age ID field
#' @param events_field Field name for the number of events (eg deaths)
#' @param exposure_field Field name for the relative exposure (eg. pop-time)
#' @param covar_names Names of all covariates, including 'intercept' if desired
#' @param distribution_family Name of the distribution to to fit the GLM with,
#'   which also determines the link function to be used. Should be a valid
#'   argument to pass to `stats::glm(family = distribution_family)`
#' @param grouping_field [optional] Field that groups observations by ID
#'
#' @return Named list with two items:
#'     - "glm_fit": Full model fit for the GLM
#'     - "fixed_effects_map": Maximum a priori estimates for covariate fixed
#'          effects, organized as a named numeric vector
#'     - "fixed_effects_grouping": Maximum a priori estimates for grouped fixed
#'          effects, organized as a vector of length(num groups). This list item
#'          is NULL if the argument `grouping_field` was not specified
#'
#' @import data.table glue stats
#' @export
find_glm_map_parameter_estimates <- function(
  in_data, events_field, exposure_field, covar_names, distribution_family,
  grouping_field = NULL
){
  # Ensure that all columns are available in input data
  reqd <- c(events_field, exposure_field, covar_names, grouping_field)
  missing_fields <- setdiff(reqd, colnames(in_data))
  if(length(missing_fields) > 0){
    stop("MAP input data missing fields: ", paste(missing_fields, collapse=', '))
  }
  in_data <- na.omit(in_data, cols = reqd)

  # Add group-based fixed effects, if specified
  if(!is.null(grouping_field) & in_data[, uniqueN(get(grouping_field)) ] > 1){
    grp_vals <- sort(unique(in_data[[grouping_field]]))
    # The first group is set as 'default' - others vary with a fixed effect
    for(grp_val in grp_vals[2:length(grp_vals)]){
      in_data[, paste0('grp',grp_val) := 0 ]
      in_data[ get(grouping_field) == grp_val, paste0('grp',grp_val) := 1 ]
    }
    grp_cols <- paste0('grp', grp_vals[2:length(grp_vals)])
  } else {
    grp_cols <- c()
  }

  # Get a field representing successes as a proportion of exposure
  in_data[, rate_success := get(events_field) / get(exposure_field) ]

  # Set up formula with an offset based on link function, then run the GLM
  formula_char <- glue::glue(
    "rate_success ~ 0 + {paste(c(covar_names, grp_cols), collapse = ' + ')}"
  )
  .env <- environment()
  formula_parsed <- as.formula(formula_char, env = .env)
  glm_fit <- stats::glm(
    formula_parsed, data = in_data, family = distribution_family,
    weights = in_data[[exposure_field]]
  )

  # Return the full GLM fit, the covariate fixed effects, and (optionally) the
  #  grouped fixed effects
  covs_map <- glm_fit$coefficients[ covar_names ]
  if(length(grp_cols) > 0){
    grps_map <- c(0, glm_fit$coefficients[ grp_cols ])
    names(grps_map) <- paste0('grp',grp_vals)
  } else {
    grps_map <- NULL
  }
  return(list(
    glm_fit = glm_fit,
    fixed_effects_map = covs_map,
    fixed_effects_grouping = grps_map
  ))
}


#' Run TMB normalization
#'
#' @description Get the TMB normalization constant from random effects
#'
#' @param adfun The ADFunction to normalize
#' @param flag [char] Flagging variable to indicate whether random effects only have been
#'   incorporated into the joint negative log-likelihood
#' @param value [int, default 0] Value that the flagging variable takes when the data
#'   should NOT be incorporated into the jnll
#' @param verbose [bool, default FALSE] return a message about normalization?
#'
#' @return Normalized ADFunction
#'
#' @import TMB
#' @import tictoc
#' @export
normalize_adfun <- function(adfun, flag, value=0, verbose=FALSE){
  if(verbose) message(" - Running normalization")
  if(verbose) tictoc::tic("    Normalization")
  normalized <- TMB::normalize(adfun, flag=flag, value=value)
  if(verbose) tictoc::toc()
  return(normalized)
}


#' Run TMB precision matrix sparsity algorithm
#'
#' @description Run "symbolic analysis", a set of algorithms that prune the
#'   precision matrix to increase sparsity and speed up optimization
#'
#' @param adfun The ADFunction to normalize
#' @param verbose [bool, default FALSE] return a message about normalization?
#'
#' @import TMB
#' @import tictoc
#' @export
run_sparsity_algorithm <- function(adfun, verbose=FALSE){
  if(verbose) message(" - Running symbolic analysis to reduce run time")
  tictoc::tic("    Symbolic analysis")
  suppressMessages(suppressWarnings(TMB::runSymbolicAnalysis(adfun)))
  tictoc::toc()
  invisible()
}

#' Set up and run TMB
#'
#' @description Generic TMB model run handler. Sets up the ADFun object, applies
#'   model speedups and fixes as specified, and optimizes using `nlminb`. This
#'   is meant to be a helper function run by more specific model functions.
#'
#' @param tmb_data_stack List containing all data inputs expected in the TMB
#'   CPP file
#' @param params_list List containing all parameters expected in the TMB CPP file
#' @param tmb_random Character vector containing all random effects that will be
#'   optimized in the inner optimizer
#' @param tmb_map Named list containing parameters that will be treated in a
#'   particular way by the optimizer
#' @param tmb_outer_maxsteps Max number of steps taken by the outer optimizer
#' @param tmb_inner_maxsteps Max number of steps taken by the inner optimizer
#'   in a single outer optimizer step
#' @param normalize [boolean, default FALSE] Run TMB's automatic process normalization
#'   function? Adds two additional items to the TMB data stack: `auto_normalize` (1) and
#'   `early_return` (0). Both are used by the \code{\link{normalize_adfun}} function to
#'   get the normalizing constant prior to optimization. For more details, see
#'   \link{\code{TMB::normalize}}.
#' @param run_symbolic_analysis [boolean, default FALSE] run symbolic analysis
#'   to speed up model run time?
#' @param set_limits [boolean, default FALSE] Set limits for fixed effects?
#' @param limit_max [numeric, default 10] If limits are set in the model,
#'   maximum value for any fixed effect
#' @param limit_min [numeric, default -limit_max] If limits are set in the model,
#'   minimum value of any fixed effect
#' @param model_name [char, default "model"] name of the model
#' @param verbose [boolean, default FALSE] Should this function return logging
#'   information about the stage of model fitting, including the outer optizimer
#'   sampling? This will be somewhat long (>100 lines)
#' @param inner_verbose [boolean, default FALSE] Should this function return
#'   logging information about inner optimizer sampling? This can be useful for
#'   debugging, but will show very verbose (>10k lines) output.
#'
#' @return list of two objects: obj (ADFunction object), and opt (optimized
#'   nlminb object)
#'
#' @useDynLib covidemr
#' @import TMB glue tictoc optimx
#' @export
setup_run_tmb <- function(
  tmb_data_stack, params_list, tmb_random, tmb_map, DLL, tmb_outer_maxsteps,
  tmb_inner_maxsteps, normalize=FALSE, run_symbolic_analysis=FALSE,
  set_limits=FALSE, limit_max=10, limit_min=-limit_max,
  optimization_methods = c('nlminb','L-BFGS-B','Rcgmin','spg','bobyqa','CG','Nelder-Mead'),
  model_name="model", verbose=FALSE, inner_verbose=FALSE
){
  # Helper function to send a message only if verbose
  vbmsg <- function(x) if(verbose) message(x)
  # Setup
  vbmsg(paste0(c("\n",rep("*",nchar(model_name)+14)),collapse=''))
  vbmsg(glue::glue("***  {model_name} RUN  ***"))
  vbmsg(paste0(c(rep("*",nchar(model_name)+14),"\n"),collapse=''))

  # Set up openmp threads
  threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
  if(threads != '') {
    vbmsg(sprintf('Detected %s threads in OMP environmental variable.',threads))
    openmp(as.numeric(threads))
  } else {
    vbmsg("Did not detect environmental OMP variable, defaulting to 2 cores. \n
           You can set this using OMP_NUM_THREADS.")
    openmp(2)
  }

  # Add flags to the data input stack indicating whether automatic process normalization
  #  should be run
  if(normalize){
    tmb_data_stack$auto_normalize <- 1L
    tmb_data_stack$early_return <- 0L
  }

  # Try optimizing using a variety of algorithms (all fit in optimx)
  for(this_method in optimization_methods){
    # Make Autodiff function
    vbmsg("Constructing ADFunction...")
    tictoc::tic("  Making Model ADFun")
    obj <- TMB::MakeADFun(
      data = tmb_data_stack,
      parameters = params_list,
      random = tmb_random,
      map = tmb_map,
      DLL = 'covidemr',
      silent = inner_verbose
    )
    obj$env$tracemgc <- as.integer(verbose)
    obj$env$inner.control$trace <- as.integer(inner_verbose)
    tictoc::toc()
    # Optionally run a normalization fix for models with large random effect sets
    if(normalize) obj <- normalize_adfun(
      adfun=obj, flag='early_return', value=1, verbose=verbose
    )
    # Optionally run optimization algorithms to improve model run time
    if(run_symbolic_analysis) run_sparsity_algorithm(adfun=obj, verbose=verbose)
    # Optionally set upper and lower limits for fixed effects
    fe_names <- names(obj$par)
    if(set_limits==TRUE){
      fe_lower_vec = rep(limit_min, times=length(fe_names))
      fe_upper_vec = rep(limit_max, times=length(fe_names))
      names(fe_lower_vec) <- names(fe_upper_vec) <- fe_names
      vbmsg(glue::glue("Fixed effects limited to the range [{limit_min},{limit_max}]."))
    } else {
      fe_lower_vec <- -Inf
      fe_upper_vec <- Inf
    }
    # Optimize using nlminb
    tictoc::tic("  Optimization")
    message(glue("\n** OPTIMIZING USING METHOD {this_method} **"))
    opt <- optimx(
      par = obj$par, fn = obj$fn, gr = obj$gr,
      lower = fe_lower_vec, upper = fe_upper_vec,
      method = this_method,
      itnmax = tmb_outer_maxsteps,
      control = list(
        rel.tol = 1E-10,
        trace = as.integer(verbose),
        follow.on = FALSE,
        dowarn = as.integer(verbose),
        maxit = tmb_inner_maxsteps,
        starttests = FALSE
      )
    )
    if(opt$convcode == 0){
      message(glue("Optimization converged using method {this_method}!"))
      break()
    } else {
      message(glue("Optimization failed using method {this_method} (code {opt$convcode})\n"))
    }
  }
  conv_code <- opt$convcode
  vbmsg(glue::glue(
    "{model_name} optimization finished with convergence code {conv_code}.\n"
  ))
  tictoc::toc()
  vbmsg(glue::glue("*** {model_name} RUN COMPLETE **************\n\n"))
  return(list(obj=obj, opt=opt))
}
