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

dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Settings 
## TODO: Convert to command line
run_sex <- 'male'
prepped_data_version <- '20200818'
holdout <- 0
use_covs <- c('intercept','tfr','unemp','socserv')

## Load and prepare data -------------------------------------------------------

# Load data
prep_file <- function(key){
  prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
  return(file.path(prep_dir, config$prepped_data_files[[key]]))
}
prepped_data <- fread(prep_file('full_data'))
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

## Set up random effects indices
age_groups <- create_age_groups(config$age_cutoffs)
age_groups[, idx_age := age_group_code - min(age_group_code) ]

year_dt <- data.table(year=sort(config$model_years))
year_dt[, idx_year := .I - 1 ]

week_dt <- data.table(week = min(config$model_week_range):max(config$model_week_range))
week_dt[, idx_week := .I - 1 ]

location_table[, idx_loc := .I - 1 ]

# Merge on data
prepped_data <- merge(
  prepped_data, location_table[, .(location_code, idx_loc)], by='location_code'
)
prepped_data <- merge(prepped_data, year_dt, by='year')
prepped_data <- merge(prepped_data, week_dt, by='week')
prepped_data <- merge(
  prepped_data, age_groups[,.(age_group_code, idx_age)], by='age_group_code'
)

# Set holdout IDs
prepped_data$idx_holdout <- 1

prepped_data[, intercept := 1 ]

# Subset to only input data for this model
in_data_final <- prepped_data[(deaths<pop) & (sex==run_sex) & (in_baseline==1), ]


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
  beta_ages = rep(0.0, nrow(age_groups)),
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
  rho_loc_trans = 1, rho_year_trans = 1, rho_week_trans = 1, rho_age_trans = 1,
  # Variance parameters
  log_sigma_Z = 0, log_sigma_nugget = 0 
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

