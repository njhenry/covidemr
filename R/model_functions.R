
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
#' @param DLL Full filepath to the TMB C++ file to use for optimization
#' @param tmb_outer_maxsteps Max number of steps taken by the outer optimizer
#' @param tmb_inner_maxsteps Max number of steps taken by the inner optimizer
#'   in a single outer optimizer step
#' @param run_symbolic_analysis [boolean, default FALSE] run symbolic analysis
#'   to speed up model run time?
#' @param normalize [boolean, default FALSE] Run normalize fix for models with
#'   large random effects sets? Only run this option if you know what you are
#'   doing and have explicitly coded this fix into the TMB code.
#' @param set_limits [boolean, default FALSE] Set limits for fixed effects?
#' @param limit_max [numeric, default 10] If limits are set in the model,
#'   maximum value for any fixed effect
#' @param limit_min [numeric, default -limit_max] If limits are set in the model,
#'   minimum value of any fixed effect
#' @param model_name [char, default "model"] name of the model
#' @param verbose [boolean, default FALSE] Should this function return
#'
#' @return list of two objects: obj (ADFunction object), and opt (optimized
#'   nlminb object)
#'
#' @import TMB glue tictoc
#' @export
setup_run_tmb <- function(
  tmb_data_stack, params_list, tmb_random, tmb_map, DLL, tmb_outer_maxsteps,
  tmb_inner_maxsteps, run_symbolic_analysis=FALSE, normalize=FALSE,
  set_limits=FALSE, limit_max=10, limit_min=-limit_max, model_name="model",
  verbose=FALSE
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
    vbmsg("Did not detect environmental OMP variable, defaulting to 4 cores. \n
           You can set this using OMP_NUM_TREADS or when launching singularity image.")
    openmp(4)
  }

  # Add a flag to the data input stack if normalize is specified
  if(normalize) tmb_data_stack$flag <- 1
  # Make Autodiff function
  vbmsg("Constructing ADFunction...")
  TMB::config(tape.parallel=0, DLL=DLL)
  tictoc::tic("  Making Model ADFun")
  obj <- TMB::MakeADFun(
    data=tmb_data_stack,
    parameters=params_list,
    random=tmb_random,
    map=tmb_map,
    DLL=DLL,
    silent=TRUE
  )
  obj$env$tracemgc <- as.integer(verbose)
  obj$env$inner.control$trace <- 0
  tictoc::toc()
  # Optionally run a normalization fix for models with large random effect sets
  if(normalize){
    vbmsg("Running normalization fix:")
    tictoc::tic("  Normalization")
    obj <- TMB::normalize(obj, flag='flag')
    tictoc::toc()
  }
  # Optionally run optimization algorithms to improve model run time
  if(run_symbolic_analysis){
    vbmsg("Running symbolic analysis to reduce run time:")
    tictoc::tic("  Running symbolic analysis")
    symbolic_output <- suppressMessages(suppressWarnings(utils::capture.output(
      TMB::runSymbolicAnalysis(obj)
    )))
    tictoc::toc()
  }
  # Optionally set upper and lower limits for fixed effects
  fe_names <- names(obj$par)
  if(set_limits==TRUE){
    fe_lower_vec = rep(limit_min, times=length(fe_names))
    fe_upper_vec = rep(limit_max, times=length(fe_names))
    names(fe_lower_vec) <- names(fe_upper_vec) <- fe_names
  } else {
    fe_lower_vec <- -Inf
    fe_upper_vec <- Inf
  }
  vbmsg(glue::glue("Fixed effects limited to the range [{limit_min},{limit_max}]."))
  # Optimize using nlminb
  vbmsg("Optimizing using nlminb:")
  tictoc::tic("  Optimization")
  opt <- do.call("nlminb", list(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = fe_lower_vec,
    upper = fe_upper_vec,
    control = list(
      eval.max=tmb_outer_maxsteps,
      iter.max=tmb_inner_maxsteps,
      trace=as.integer(verbose)
    )
  ))
  conv_code <- opt$convergence
  vbmsg(glue::glue(
    "{model_name} optimization finished with convergence code {conv_code}.\n"
  ))
  tictoc::toc()
  vbmsg(glue::glue("*** {model_name} RUN COMPLETE **************\n\n"))
  return(list(obj=obj, opt=opt))
}


#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
#' @import matrixcalc
#' @export
rmvnorm_prec <- function(mu, prec, n.sims, find_near_pd=FALSE) {
  if(find_near_pd){
    if( !matrixcalc::is.positive.definite(as.matrix(prec))){
      message("WARNING: Joint Precision matrix is not positive definite.")
      message("Finding closest positive definite matrix...\n")
      prec <- Matrix(
        Matrix::nearPD(prec)$mat,
        sparse=TRUE
      )
    }
  }
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Cholesky(prec, super=TRUE)
  return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as(L_inv, "Matrix")), z)))
}
