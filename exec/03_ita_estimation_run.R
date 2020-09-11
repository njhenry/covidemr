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
run_sex <- 'female'
prepped_data_version <- '20200909'
model_run_version <- '20200909'
holdout <- 0
use_covs <- c('intercept','tfr','unemp','socserv')

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

prepped_data <- prepped_data[sex==run_sex, ]
# Set holdout IDs
prepped_data$idx_holdout <- 1

# Subset to only input data for this model
in_data_final <- prepped_data[(deaths<pop) & (pop>0) & (in_baseline==1), ]


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
  idx_holdout = in_data_final$idx_holdout,
  loc_adj_mat = as(adjmat, 'dgTMatrix')
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
if(length(params_list$beta_ages) == 1){
  tmb_map <- list()
} else {
  tmb_map = list(
    beta_ages = as.factor(c(NA, 2:length(params_list$beta_ages)))
  )
}

## Run modeling!
model_fit <- setup_run_tmb(
  tmb_data_stack=data_stack,
  params_list=params_list,
  tmb_random=c('Z_stwa','nugget'),
  tmb_map=tmb_map,
  normalize = TRUE, run_symbolic_analysis = TRUE,
  tmb_outer_maxsteps=1000, tmb_inner_maxsteps=1000, 
  model_name="ITA deaths model", verbose=TRUE
)

sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)


## Save outputs ----------------------------------------------------------------

out_file_stub <- sprintf("%s_holdout%i", run_sex, holdout)
out_dir <- file.path(config$paths$model_results, model_run_version)
dir.create(out_dir, showWarnings = FALSE)

## Create post-estimation predictive objects -----------------------------------

mu <- c(sdrep$par.fixed, sdrep$par.random)
parnames <- names(mu)
keep_fields <- which(parnames %in% c('beta_covs', 'beta_ages', 'Z_stwa'))
parnames_sub <- parnames[keep_fields]

if(sum(parnames=='Z_stwa') != nrow(template_dt)) stop("Issue with draw dimensions")

tictoc::tic(sprintf("%i draws", config$num_draws))
pryr::mem_change(
  draws <- covidemr::rmvnorm_prec(
    mu = mu[keep_fields],
    prec = sdrep$jointPrecision[keep_fields, keep_fields],
    n.sims = config$num_draws
  )
)
tictoc::toc()


## Generate predictions by location-year-week-age

# Reorder template matrix by age-week-year-location
# This is the reverse ordering of the random effect dimensions, accounting for 
# how arrays are translated into vectors
template_dt <- template_dt[order(idx_age, idx_week, idx_year, idx_loc)]

# Draws == logit^-1( Covariate effects + age fixed effect + stwa RE + nugget )
# Covariate effects
cov_fes <- as.matrix(template_dt[, ..use_covs]) %*% draws[parnames_sub=='beta_covs',]
# Age fixed effect
age_fes <- rbind(0, draws[parnames_sub=='beta_ages',])[ template_dt$idx_age + 1, ]
# Add it all together and take the inverse logit
preds <- plogis(cov_fes + age_fes + draws[parnames_sub=='Z_stwa', ])

# Get mean, lower, upper
summs <- cbind(
  rowMeans(preds), matrixStats::rowQuantiles(preds, probs=c(0.5, 0.025, 0.975))
)
colnames(summs) <- c('mean','median','lower','upper') 
summs <- cbind(template_dt, summs)


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
