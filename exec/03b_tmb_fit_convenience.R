tmb_outer_maxsteps = 3000; tmb_inner_maxsteps = 3000; normalize=FALSE; run_symbolic_analysis=FALSE; set_limits=FALSE; limit_max=10; limit_min=-limit_max; optimization_methods = 'nlminb'; model_name="model"; verbose=TRUE; inner_verbose=TRUE; this_method = 'nlminb'
vbmsg <- function(x) if(verbose) message(x)
# Setup
vbmsg(paste0(c("\n",rep("*",nchar(model_name)+14)),collapse=''))
vbmsg(glue::glue("***  {model_name} RUN  ***"))
vbmsg(paste0(c(rep("*",nchar(model_name)+14),"\n"),collapse=''))
vbmsg("Constructing ADFunction...")
tictoc::tic("  Making Model ADFun")
if(normalize){
  tmb_data_stack$auto_normalize <- 1L
  tmb_data_stack$early_return <- 0L
}
obj <- TMB::MakeADFun(
  data = tmb_data_stack,
  parameters = params_list,
  random = tmb_random,
  random.start = expression(rep(0, length(random))),
  map = tmb_map,
  DLL = 'covidemr',
  silent = inner_verbose
)
tictoc::toc()
obj$env$tracemgc <- as.integer(verbose)
obj$env$inner.control$trace <- as.integer(inner_verbose)
# Optionally run a normalization fix for models with large random effect sets
if(normalize) obj <- normalize_adfun(
  adfun=obj, flag='early_return', value=1, verbose=verbose
)

# Optimize using nlminb
tictoc::tic("  Optimization")
message(glue("\n** OPTIMIZING USING METHOD {this_method} **"))
opt <- optimx(
  par = obj$par, fn = obj$fn, gr = obj$gr,
  lower = -Inf, upper = Inf, method = this_method, itnmax = tmb_outer_maxsteps,
  hessian = FALSE,
  control = list(
    trace = as.integer(verbose), follow.on = FALSE, kkt = FALSE,
    dowarn = as.integer(verbose), maxit = tmb_inner_maxsteps, starttests = FALSE
  )
)
conv_code <- opt$convcode
vbmsg(glue::glue(
  "{model_name} optimization finished with convergence code {conv_code}.\n"
))
tictoc::toc()
vbmsg(glue::glue("*** {model_name} RUN COMPLETE **************\n\n"))
model_fit <- list(obj = obj, opt = opt)
