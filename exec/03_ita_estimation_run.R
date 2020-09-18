## -----------------------------------------------------------------------------
## 
## 03: Space-time excess mortality modeling code, using prepared data as input
## 
## For more details, see README at https://github.com/njhenry/covidemr/
## 
## -----------------------------------------------------------------------------


library(data.table)
library(sp)
library(sf)
library(parallel)

dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Settings 
## TODO: Convert to command line
run_sex <- 'male'
prepped_data_version <- '20200917'
model_run_version <- '20200917'
holdout <- 0
use_covs <- c('intercept') # , 'tfr','unemp','socserv'
use_Z_stwa <- FALSE
use_Z_sta <- !use_Z_stwa
use_Z_fourier <- !use_Z_stwa
fourier_levels <- 2


## Load and prepare data -------------------------------------------------------

# Load data
prep_file <- function(key){
  prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
  return(file.path(prep_dir, config$prepped_data_files[[key]]))
}
prepped_data <- fread(prep_file('full_data'))
template_dt <- fread(prep_file('template'))
location_table <- fread(prep_file('location_table'))
adjmat <- readRDS(prep_file('adjacency_matrix'))

# Prepare the following required inputs for modeling:
# Input data stack:
#  - Response y_i and population n_i, as vectors
#  - Days of exposure exp_i, for normalization
#  - Fixed effects (with intercept) as a matrix
#  - Zero-indexed random effects: year, location code (province), age
#  - Adjacency matrix
# Parameters list:
#  - Vector of covariate fixed effects
#  - Vector of age fixed effects
#  - Array of space (x) time (x) age-structured random effects
#  - Vector of iid random effects
#  - Random effect hyperparameters
# Additional data:
#  - Table associating province/age/time values with random effect indices

## Subset data to sex being modeled; merge indices on prepared data
template_dt <- template_dt[sex==run_sex, ]
template_dt[, idx_fourier := 0 ]

prepped_data <- prepped_data[sex==run_sex, ]
# Set holdout IDs
prepped_data[, idx_fourier := 0 ] # Alternate option: group by age
prepped_data$idx_holdout <- 1

# Subset to only input data for this model
in_data_final <- prepped_data[(deaths<pop) & (pop>0) & (in_baselin=e=1), ]


data_stack <- list(
  holdout = holdout,
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
  use_Z_stwa = as.integer(use_Z_stwa),
  use_Z_sta = as.integer(use_Z_sta),
  use_Z_fourier = as.integer(use_Z_fourier),
  harmonics_level = as.integer(fourier_levels)
)

params_list <- list(
  # Fixed effects
  beta_covs = rep(0.0, length(use_covs)),
  beta_ages = rep(0.0, length(unique(template_dt$idx_age))),
  # Structured random effect
  Z_stwa = array(
    0.0, 
    dim = c(
      nrow(location_table), # Number of locations
      length(config$model_years), # Number of unique modeled years
      diff(config$model_week_range) + 1, # Number of weeks in each year
      length(config$age_cutoffs) # Number of age groups
    )
  ),
  Z_sta = array(
    0.0,
    dim = c(
      nrow(location_table), # Number of locations
      length(config$model_years), # Number of unique modeled years
      length(config$age_cutoffs) # Number of age groups
    )
  ),
  Z_fourier = array(
    0.0,
    dim = c(
      max(in_data_final$idx_fourier) + 1,
      2 * fourier_levels
    )
  ),
  # Unstructured random effect
  nugget = rep(0.0, length(data_stack$n_i)),
  # Rho parameters
  rho_loc_trans = 0.0, rho_year_trans = 0.0, rho_week_trans = 0.0,
  rho_age_trans = 0.0,
  # Variance parameters
  log_sigma_loc = 0, log_sigma_year = 0, log_sigma_week = 0,
  log_sigma_age = 0, log_sigma_nugget = 0
)


## Hold first age fixed effect constant
tmb_map <- list()
if(length(params_list$beta_ages) > 1){
  tmb_map$beta_ages <- as.factor(c(NA, 2:length(params_list$beta_ages)))
}
if( !use_Z_stwa ){
  tmb_map$Z_stwa <- rep(as.factor(NA), prod(dim(params_list$Z_stwa)))
}
if( !use_Z_sta ){
  tmb_map$Z_sta <- rep(as.factor(NA), prod(dim(params_list$Z_sta)))
}
if( !use_Z_fourier ){
  tmb_map$Z_fourier <- rep(as.factor(NA), prod(dim(params_list$Z_fourier)))
}
tmb_map$nugget <- rep(as.factor(NA), length(params_list$nugget))
tmb_map$rho_week_trans <- as.factor(NA)
tmb_map$log_sigma_week <- as.factor(NA)
tmb_map$log_sigma_nugget <- as.factor(NA)

tmb_random <- c('nugget')
if(use_Z_stwa) tmb_random <- c(tmb_random, 'Z_stwa')
if(use_Z_sta) tmb_random <- c(tmb_random, 'Z_sta')
# if(use_Z_fourier) tmb_random <- c(tmb_random, "Z_fourier")

## Run modeling!
model_fit <- setup_run_tmb(
  tmb_data_stack=data_stack,
  params_list=params_list,
  tmb_random=tmb_random,
  tmb_map=tmb_map,
  normalize = TRUE, run_symbolic_analysis = TRUE,
  tmb_outer_maxsteps=3000, tmb_inner_maxsteps=3000, 
  model_name="ITA deaths model", verbose=FALSE,
  optimization_methods = c('nlminb', 'L-BFGS-B')
)

sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)



## Create post-estimation predictive objects -----------------------------------

# Draws of parameters and baseline modeled deaths (assuming no COVID)
postest_draws_list <- covidemr::generate_stwa_draws(
  tmb_sdreport = sdrep,
  num_draws = config$num_draws,
  covariate_names = use_covs,
  template_dt = template_dt,
  fourier_harmonics_level = fourier_levels
)
param_draws <- postest_draws_list$param_draws
pred_draws <- postest_draws_list$predictive_draws
# Summarize draws
pred_summary <- cbind(template_dt, summarize_draws(pred_draws))

# Generate draws of excess mortality (true deaths - baseline)
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



## Save outputs ----------------------------------------------------------------

out_file_stub <- sprintf("%s_holdout%i", run_sex, holdout)
out_dir <- file.path(config$paths$model_results, model_run_version)
dir.create(out_dir, showWarnings = FALSE)

## Save all model output files
model_results_dir <- config$path$model_results

for(obj_str in names(config$results_files)){
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

message("*** FIN ***")
