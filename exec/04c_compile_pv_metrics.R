## #############################################################################
##
## 04c: Compile predictive validity metrics from multiple runs
##
## Purpose: Combine tables of predictive validity metrics across many model runs
##
## #############################################################################

library(data.table)
library(glue)

dev_fp <- '~/repos/covidemr/'
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

# TODO: convert to CLI
compare_runs <- c(
  '20201212f3fal', '20201212f2fal', '20201212f1fal', '20201212f3fa',
  '20201212f2fa', '20201212f1fa', '20201212f3fl', '20201212f2fl',
  '20201212f1fl', '20201212f3', '20201212f2', '20201212f1', '20201212nof'
)
is_oos <- 'is'

results_list <- list(
  full_noagg = vector('list', length=length(compare_runs)),
  provmonth = vector('list', length=length(compare_runs))
)

for(ii in 1:length(compare_runs)){
  run <- compare_runs[ii]
  pv_dir <- file.path(config$paths$model_results, run, 'pred_metrics')
  for(aggtype in names(results_list)){
    results_list[[aggtype]][[ii]] <- rbindlist(list(
      cbind(
        data.table(run_date = run, sex = 'male'),
        data.table::fread(glue::glue('{pv_dir}/male_{is_oos}_error_{aggtype}.csv')),
        data.table::fread(glue::glue('{pv_dir}/male_{is_oos}_coverage_{aggtype}.csv'))
      ),
      cbind(
        data.table(run_date = run, sex = 'female'),
        data.table::fread(glue::glue('{pv_dir}/female_{is_oos}_error_{aggtype}.csv')),
        data.table::fread(glue::glue('{pv_dir}/female_{is_oos}_coverage_{aggtype}.csv'))
      )
    ))
    results_list[[aggtype]][[ii]]$agg_type <- aggtype
  }
}

final_results <- rbindlist(lapply(results_list, rbindlist))

out_fp <- glue::glue('{config$paths$model_results}/compare_model_pv_{Sys.Date()}.csv')
fwrite(final_results, file=out_fp)
