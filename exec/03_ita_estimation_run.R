## -----------------------------------------------------------------------------
##
## 03: Space-time excess mortality modeling code, using prepared data as input
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

message("================= Space-time mortality fitting script ================")
message("Start time: ", Sys.time())
message("Script start")
library(argparse)
library(data.table)
library(parallel)
library(sp)
library(sf)
library(tictoc)

dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Load run-specific settings from command line
ap <- argparse::ArgumentParser(
  description='COVID Excess Mortality: Model fitting script',
  allow_abbrev=FALSE
)
ap$add_argument('--run-sex', type='character', help='Sex-specific model to run')
ap$add_argument('--data-version', type='character', help='Prepped data version date')
ap$add_argument('--model-version', type='character', help='Model run version date')
ap$add_argument('--holdout', type='integer', default=0, help='Holdout number')
ap$add_argument(
  '--use-covs', type='character', nargs='+', help='Covariates to use (include intercept)'
)
ap$add_argument('--use-Z-stwa', help='Use space-time-week-age RE?', action='store_true')
ap$add_argument('--use-Z-sta', help='Use space-time-age RE?', action='store_true')
ap$add_argument('--use-Z-fourier', help='Use periodic-fit RE?', action='store_true')
ap$add_argument('--use-nugget', help='Add a nugget?', action='store_true')
ap$add_argument(
  '--fourier-levels', help='Number of levels for Fourier analysis',
  type='integer', default=2
)
ap$add_argument(
  '--fourier-groups', type='character', nargs='*', default=NULL,
  help='Grouping fields for seasonality (default one group for all data)'
)
args <- ap$parse_args(commandArgs(TRUE))
# args <- list(
#   run_sex = 'male', data_version = '20201026', model_version = '20201103f2',
#   holdout = 0, use_covs = c(
#     'intercept', 'tfr', 'unemp', 'socserv', 'tax_brackets', 'hc_access',
#     'elevation', 'temperature'
#   ), use_Z_stwa = FALSE, use_Z_sta = TRUE, use_Z_fourier = TRUE,
#   use_nugget = FALSE, fourier_levels = 2,
#   fourier_groups = NULL
# )
message(str(args))
use_covs <- args$use_covs # Shorten for convenience


## Load and prepare data -------------------------------------------------------

# Helper functions for loading prepared data
prep_dir <- file.path(config$paths$prepped_data, args$data_version)
get_prep_fp <- function(ff) file.path(prep_dir, config$prepped_data_files[[ff]])

prepped_data <- data.table::fread(get_prep_fp('full_data_rescaled'))
template_dt <- data.table::fread(get_prep_fp('template'))
covar_scaling_factors <- data.table::fread(get_prep_fp('covar_scaling_factors'))
location_table <- data.table::fread(get_prep_fp('location_table'), na.strings='')
adjmat <- readRDS(get_prep_fp('adjacency_matrix'))

## Subset data to sex being modeled; merge indices on prepared data
template_dt <- template_dt[sex==args$run_sex, ]
template_dt$idx_fourier <- assign_seasonality_ids(
  input_data = template_dt, grouping_fields = args$fourier_groups
)

prepped_data <- prepped_data[sex==args$run_sex, ]
# Set holdout IDs
prepped_data$idx_fourier <- assign_seasonality_ids(
  input_data = prepped_data, grouping_fields = args$fourier_groups
)

# Subset to only input data for this model
in_data_final <- prepped_data[(deaths<pop) & (pop>0) & (in_baseline==1), ]


## Define input data stack (the data used to fit the model) --------------------

# Define data stack
data_stack <- list(
  holdout = args$holdout,
  y_i = in_data_final$deaths,
  n_i = in_data_final$pop,
  days_exp_i = in_data_final$observed_days,
  X_ij = as.matrix(in_data_final[, ..use_covs]),
  idx_loc = in_data_final$idx_loc,
  idx_year = in_data_final$idx_year,
  idx_week = in_data_final$idx_week,
  idx_age = in_data_final$idx_age,
  idx_fourier = in_data_final$idx_fourier,
  idx_holdout = in_data_final$idx_holdout,
  loc_adj_mat = as(adjmat, 'dgTMatrix'),
  use_Z_stwa = as.integer(args$use_Z_stwa),
  use_Z_sta = as.integer(args$use_Z_sta),
  use_Z_fourier = as.integer(args$use_Z_fourier),
  use_nugget = as.integer(args$use_nugget),
  harmonics_level = as.integer(args$fourier_levels)
)


## Define parameters (what is fit in the model) and their constraints ----------

# Find MAP approximation of all model fixed effects. This will be used to set
#  the starting values for the model
max_a_priori_list <- covidemr::find_binomial_map_parameter_estimates(
  in_data = copy(in_data_final),
  numerator_field = 'deaths',
  denominator_field = 'pop',
  covar_names = use_covs,
  grouping_field = 'idx_age'
)
map_glm_fit <- max_a_priori_list$glm_fit

params_list <- list(
  # Fixed effects
  beta_covs = unname(max_a_priori_list$fixed_effects_map),
  beta_ages = unname(max_a_priori_list$fixed_effects_grouping),
  # Rho parameters
  rho_loc_trans = 0.0, rho_year_trans = 0.0, rho_week_trans = 0.0,
  rho_age_trans = 0.0,
  # Variance parameters
  log_sigma_loc = 0, log_sigma_year = 0, log_sigma_week = 0,
  log_sigma_age = 0, log_sigma_nugget = 0,
  # Structured random effect
  Z_stwa = array(0.0,
    dim = c(
      nrow(location_table), # Number of locations
      length(config$model_years), # Number of unique modeled years
      diff(config$model_week_range) + 1, # Number of weeks in each year
      length(config$age_cutoffs) # Number of age groups
    )
  ),
  Z_sta = array(0.0,
    dim = c(
      nrow(location_table), # Number of locations
      length(config$model_years), # Number of unique modeled years
      length(config$age_cutoffs) # Number of age groups
    )
  ),
  Z_fourier = array(0.0,
    dim = c(
      max(in_data_final$idx_fourier) + 1,
      2 * args$fourier_levels
    )
  ),
  # Unstructured random effect
  nugget = rep(0.0, length(data_stack$n_i))
)

# Fix particular parameter values using the TMB map
tmb_map <- list()
if(length(params_list$beta_ages) > 1){
  tmb_map$beta_ages <- as.factor(c(NA, 2:length(params_list$beta_ages)))
}
if(!args$use_Z_stwa){
  tmb_map$Z_stwa <- rep(as.factor(NA), length(params_list$Z_stwa))
  tmb_map$rho_week_trans <- as.factor(NA)
  tmb_map$log_sigma_week <- as.factor(NA)
}
if(!args$use_Z_sta) tmb_map$Z_sta <- rep(as.factor(NA), length(params_list$Z_sta))
if(!args$use_Z_fourier) tmb_map$Z_fourier <- rep(as.factor(NA), length(params_list$Z_fourier))
if(!args$use_nugget){
  tmb_map$nugget <- rep(as.factor(NA), length(params_list$nugget))
  tmb_map$log_sigma_nugget <- as.factor(NA)
}

# Set random effects
tmb_random <- character(0)
if(args$use_nugget) tmb_random <- c(tmb_random, 'nugget')
if(args$use_Z_stwa) tmb_random <- c(tmb_random, 'Z_stwa')
if(args$use_Z_sta) tmb_random <- c(tmb_random, 'Z_sta')
if(args$use_Z_fourier) tmb_random <- c(tmb_random, 'Z_fourier')


## Run modeling ----------------------------------------------------------------

tictoc::tic("Full TMB model fitting")
model_fit <- covidemr::setup_run_tmb(
  tmb_data_stack=data_stack,
  params_list=params_list,
  tmb_random=tmb_random,
  tmb_map=tmb_map,
  normalize = TRUE, run_symbolic_analysis = TRUE,
  tmb_outer_maxsteps=3000, tmb_inner_maxsteps=3000,
  model_name="ITA deaths model",
  verbose=FALSE, inner_verbose=FALSE,
  optimization_methods = c('nlminb', 'L-BFGS-B')
)

message("Getting sdreport and joint precision matrix...")
sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)
tictoc::toc()


## Create post-estimation predictive objects -----------------------------------

# Draws of parameters and baseline modeled deaths (assuming no COVID)
tictoc::tic("Full model postestimation")
postest_list <- vector('list', length = ceiling(config$num_draws / 50))
for(ii in 1:length(postest_list)){
  postest_list[[ii]] <- covidemr::generate_stwa_draws(
    tmb_sdreport = sdrep,
    keep_params = c("beta_covs", "beta_ages", "Z_stwa", "Z_sta", "Z_fourier"),
    num_draws = min(50, config$num_draws - (ii - 1) * 50),
    covariate_names = use_covs,
    template_dt = template_dt,
    rescale_covars = TRUE,
    covar_scaling_factors = covar_scaling_factors,
    fourier_harmonics_level = args$fourier_levels
  )
}
cbindlist <- function(a_list) setDT(unlist(a_list, recursive=FALSE))
param_draws <- cbindlist(c(
  list(data.table(parameter = postest_list[[1]]$param_names)),
  lapply(postest_list, function(sl) as.data.table(sl$param_draws))
))
colnames(param_draws) <- c('parameter', paste0('V',1:config$num_draws))
pred_draws <- cbindlist(lapply(postest_list, function(sl) as.data.table(sl$predictive_draws)))
colnames(pred_draws) <- paste0('V',1:config$num_draws)
rm(postest_list)
# Summarize draws
pred_summary <- cbind(template_dt, summarize_draws(pred_draws))

if(args$holdout == 0){
  # For in-sample runs, generate additional predictive estimates:
  # Draws of excess mortality (true deaths - baseline)
  excess_draws_list <- covidemr::get_excess_death_draws(
    template_dt = template_dt,
    baseline_draws = pred_draws,
    death_data = prepped_data[in_baseline==0,]
  )
  excess_index <- excess_draws_list$obs_deaths
  excess_draws <- excess_draws_list$excess_draws
  proportion_draws <- excess_draws_list$proportion_draws
  # Summarize draws
  prop_summ <- summarize_draws(proportion_draws)
  names(prop_summ) <- paste0('prop_',names(prop_summ))
  excess_summary <- cbind(excess_index, summarize_draws(excess_draws), prop_summ)
} else {
  # FOR OUT-OF-SAMPLE RUNS ONLY: Keep only the OOS data subset
  pred_sub_idx <- which(pred_summary$idx_holdout == args$holdout)
  pred_draws <- pred_draws[ pred_sub_idx , ]
}
tictoc::toc()


## Save outputs ----------------------------------------------------------------

out_file_stub <- sprintf("%s_holdout%i", args$run_sex, args$holdout)
out_dir <- file.path(config$paths$model_results, args$model_version)
dir.create(out_dir, showWarnings = FALSE)

## Save all model output files
model_results_dir <- config$path$model_results
model_run_version <- args$model_version
run_sex <- args$run_sex
holdout <- args$holdout

# Always write out args and config
yaml::write_yaml(args, file = glue::glue(
  '{model_results_dir}/{model_run_version}/{run_sex}_{holdout}_model_args.yaml'
))
yaml::write_yaml(config, file = glue::glue(
  '{model_results_dir}/{model_run_version}/{run_sex}_{holdout}_model_config.yaml'
))

if(holdout != 0){
  save_objs <- c('model_fit','pred_draws')
} else {
  save_objs <- names(config$results_files)
}

for(obj_str in save_objs){
  # Parse output file
  out_fp <- glue::glue(config$results_files[[obj_str]])
  message(glue::glue("Saving {obj_str} to {out_fp}"))
  # Save to file differently depending on format
  if(endsWith(tolower(out_fp), 'rds')){
    saveRDS(get(obj_str), file = out_fp)
  } else if(endsWith(tolower(out_fp), 'csv')){
    fwrite(get(obj_str), file = out_fp, row.names = FALSE)
  } else {
    stop(sprintf("Save function not enabled for file type of %s", out_fp))
  }
}

message("Model fitting script COMPLETE")
message("End time: ", Sys.time())
