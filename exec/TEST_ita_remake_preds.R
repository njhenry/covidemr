## -----------------------------------------------------------------------------
##
## 03: Space-time excess mortality modeling code, using prepared data as input
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

rm(list=ls())

message("================= Space-time mortality fitting script ================")
message("Start time: ", Sys.time())
message("Script start")
library(data.table)
library(glue)
library(tictoc)

dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Load run-specific settings from command line
args <- list(
  run_sex = 'male', data_version = '20201026', model_version = '20201103f3fageloc',
  holdout = 0, use_covs = c(
    'intercept', 'tfr', 'unemp', 'socserv', 'tax_brackets', 'hc_access',
    'elevation', 'temperature'
  ), use_Z_stwa = FALSE, use_Z_sta = TRUE, use_Z_fourier = TRUE,
  use_nugget = FALSE, fourier_levels = 3,
  fourier_groups = c('age_group_code', 'location_code')
)
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

out_file_stub <- sprintf("%s_holdout%i", args$run_sex, args$holdout)
out_dir <- file.path(config$paths$model_results, args$model_version)
dir.create(out_dir, showWarnings = FALSE)

## Save all model output files
model_results_dir <- config$path$model_results
model_run_version <- args$model_version
run_sex <- args$run_sex
holdout <- args$holdout

model_fit <- readRDS(glue(config$results_files$model_fit))
sdrep <- readRDS(glue(config$results_files$sdrep))


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

save_objs <- c(
  'param_draws', 'pred_draws', 'pred_summary', 'excess_draws',
  'proportion_draws', 'excess_summary'
)

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

message("End time: ", Sys.time())
