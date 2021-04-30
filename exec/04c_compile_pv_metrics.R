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

compare_runs <- c(
  '20210421_f2fal', '20210421_f2fal_nn', '20210421_f2fal_ns',
  '20210421_f2fl_ns', '20210421_f3fal_ns' # '20210421_f2fal_ns_nn'
)
is_oos <- c('is','oos')

for(isval in is_oos){
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
          data.table::fread(glue::glue('{pv_dir}/male_{isval}_error_{aggtype}.csv')),
          data.table::fread(glue::glue('{pv_dir}/male_{isval}_coverage_{aggtype}.csv'))
        ),
        cbind(
          data.table(run_date = run, sex = 'female'),
          data.table::fread(glue::glue('{pv_dir}/female_{isval}_error_{aggtype}.csv')),
          data.table::fread(glue::glue('{pv_dir}/female_{isval}_coverage_{aggtype}.csv'))
        )
      ))
      results_list[[aggtype]][[ii]]$agg_type <- aggtype
    }
  }

  final_results <- rbindlist(lapply(results_list, rbindlist))

  out_fp <- glue::glue('{config$paths$model_results}/compare_model_{isval}_pv_{Sys.Date()}.csv')
  fwrite(final_results, file=out_fp)
}
