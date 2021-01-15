## #############################################################################
##
## 04b: GENERATE PREDICTIVE VALIDITY METRICS FOR A MODEL RUN
##
## Estimate root mean squared error, relative squared error, and coverage for
##   a model run (either in-sample or out-of-sample). Various aggregations of
##   model results before estimating error:
##    - No aggregation (most detailed)
##    - Results aggregated by age group, most detailed otherwise
##    - Results aggregated by "month" (4-week grouping), most detailed otherwise
##    - Results aggregated by region, most detailed otherwise
##    - Results aggregated by region and age group, most detailed otherwise
##
## Report predictive validity metrics by grouping for most detailed run:
##    - By week
##    - By "month" (4-week grouping)
##    - By age group
##    - By region
##    - By region and age group
##
## #############################################################################

library(argparse)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Load run-specific settings from command line
ap <- argparse::ArgumentParser(
  description='COVID Excess Mortality: Generate predictive validity metrics'
)
ap$add_argument('--run-sex', type='character', help='Sex-specific model to run')
ap$add_argument('--data-version', type='character', help='Prepped data version date')
ap$add_argument('--model-version', type='character', help='Model run version date')
ap$add_argument(
  '--oos',
  help='Create out-of-sample PV metrics? Default in-sample',
  action='store_true'
)
args <- ap$parse_args(commandArgs(TRUE))
# args <- list(
#   run_sex = 'male', data_version ='20201203', model_version = '20201211f2fageloc',
#   oos = TRUE
# )
# Set some variables as globals for use with `glue()` later in script
message(str(args))
run_sex <- args$run_sex
model_run_version <- args$model_version

## Create output directory
model_results_dir <- config$paths$model_results
model_vers_dir <- file.path(model_results_dir, model_run_version)

if(!dir.exists(model_vers_dir)){
  stop("Model directory ", model_vers_dir, " does not exist.")
}
pv_dir <- file.path(model_vers_dir, 'pred_metrics')
viz_dir <- file.path(model_vers_dir, 'viz')
if(!dir.exists(pv_dir)) dir.create(pv_dir)
if(!dir.exists(viz_dir)) dir.create(viz_dir)


## Load model results of draws -------------------------------------------------

# Load template dataset and subset by sex
template_dt <- fread(file.path(
  config$paths$prepped_data, args$data_version, config$prepped_data_files$template
))[ sex == run_sex, ]
template_dt[, rowidx := .I ]
check_index <- template_dt$rowidx

# Add on baseline population and deaths counts
num_denom_data <- fread(file.path(
  config$paths$prepped_data, args$data_version, config$prepped_data_files$full_data
))[ (sex == run_sex) & (in_baseline == 1), ]
data_merge_cols <- c(
  'location_code', 'year', 'week', 'age_group_code', 'region_code'
)
template_dt[
  num_denom_data,
  `:=` (deaths = i.deaths, pop = i.pop),
  on = data_merge_cols
]
if(nrow(template_dt) != length(check_index) | any(check_index != template_dt$rowidx)){
  stop("Adding numerator and denominator data changed row IDs")
}

# Load draw-level data and column bind with identifying data
if(args$oos){
  holdouts <- 1:5
  fp_prefix <- glue::glue('{run_sex}_oos')
} else {
  holdouts <- 0
  fp_prefix <- glue::glue('{run_sex}_is')
}

draws_list <- vector('list', length=length(holdouts))
for(ii in 1:length(holdouts)){
  holdout <- holdouts[ii]
  message("Loading draws for holdout ", holdout, "...")
  # Subset the template to only this holdout (most useful for OOS)
  templ_sub <- copy(template_dt)
  if(holdout != 0) templ_sub <- templ_sub[idx_holdout == holdout, ]
  # Load draws matrix and rename fields
  draws_mat <- fread(glue(config$results_files$pred_draws))
  draws_mat[draws_mat > 1] <- 1
  draw_col_names <- paste0('draw', 1:ncol(draws_mat))
  colnames(draws_mat) <- draw_col_names
  # Check that the template matches the number of draws
  if(nrow(templ_sub) != nrow(draws_mat)){
    stop("Mismatch in template and draw row counts")
  }
  # Column bind and assign to list
  draws_list[[ii]] <- cbind(templ_sub, draws_mat)
}

# Row-bind across samples
pred_data_full <- na.omit(
  data.table::rbindlist(draws_list), cols=c('deaths','pop',draw_col_names)
)
rm(list = c('template_dt','num_denom_data','templ_sub','draws_mat','draws_list'))

plot_years <- sort(unique(template_dt$year))
plot_colors <- RColorBrewer::brewer.pal(name = "Set2", n=length(plot_years))
names(plot_colors) <- as.character(plot_years)


## Generate predictive validity metrics: Most detailed input data --------------

pred_data_full$pred_mean <- rowMeans(pred_data_full[, ..draw_col_names])
pred_data_full[, obs := deaths / pop ]

error_full <- calculate_rmse_rse(
  in_data = pred_data_full, num_field = 'deaths', denom_field = 'pop',
  est_field = 'pred_mean'
)
coverage_full <- calculate_coverage(
  in_data = pred_data_full, num_field = 'deaths', denom_field = 'pop',
  draw_fields = draw_col_names, coverage_levels = c(.5, .8, .9, .95, .99),
  pois_sim = TRUE
)
## Save to file
fwrite(error_full, file = glue::glue('{pv_dir}/{fp_prefix}_error_full_noagg.csv'))
fwrite(coverage_full, file = glue::glue('{pv_dir}/{fp_prefix}_coverage_full_noagg.csv'))


## Plot observed vs expected
pred_data_full[, `:=` (obs_plot = obs * 1E5, pred_mean_plot = pred_mean * 1E5)]
# Age groups
age_groups <- unique(pred_data_full$age_group_name)
age_cols <- RColorBrewer::brewer.pal(name = "Set1", n=length(age_groups))
names(age_cols) <- age_groups


obs_scatter_all <- ggplot(data=pred_data_full) +
  geom_point(
    aes(x=obs_plot, y=pred_mean_plot, size=sqrt(pop)),
    alpha = .2
  ) +
  facet_wrap('age_group_name', ncol=ceiling(sqrt(uniqueN(pred_data_full$age_group_name)))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = '#AAAAAA') +
  scale_color_manual(values=age_cols) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  scale_size_continuous(
    breaks = c(5E3, 1E4, 5E4, 1E5), labels = c('5k', '10k', '50k', '100k')
  ) +
  labs(
    title = paste0('Weekly mortality rate among ',run_sex,'s'),
    subtitle = 'Reported as deaths per 100,000 person-weeks',
    x = 'Observed', y = 'Predicted (mean)', color = 'Age Group',
    size = 'Population\nsize'
  ) +
  theme_bw()
png(
  filename = glue::glue('{viz_dir}/{fp_prefix}_scatter_full_noagg.png'),
  height = 6, width = 6, units = 'in', res = 100
)
print(obs_scatter_all)
dev.off()


## Aggregate across all observations across age groups and "month", then get
##   predictive validity metrics
pred_data_full[, wkgrp := ceiling(week/4)]

data_by_prov_month <- aggregate_data_and_draws(
  in_data = pred_data_full, num_field = 'deaths', denom_field = 'pop',
  draw_fields = c(draw_col_names, 'obs'),
  group_fields = c('year', 'wkgrp', 'location_code', 'location_name')
)
# Merge population back on
pop_prepped <- fread(file.path(
  config$paths$prepped_data, args$data_version, config$prepped_data_files$population
))
pop_agg <- pop_prepped[, .(pop=sum(pop)), by=c('year', 'location_code')]
data_by_prov_month[pop_agg, pop := i.pop, on=c('year', 'location_code')]

# Get mean and predicted weekly mortality rates
data_by_prov_month$pred_mean <- rowMeans(data_by_prov_month[, ..draw_col_names])
data_by_prov_month[, deaths := obs * pop ]

error_provmonth <- calculate_rmse_rse(
  in_data = data_by_prov_month, num_field = 'deaths', denom_field = 'pop',
  est_field = 'pred_mean'
)
coverage_provmonth <- calculate_coverage(
  in_data = data_by_prov_month, num_field = 'deaths', denom_field = 'pop',
  draw_fields = draw_col_names, coverage_levels = c(.5, .8, .9, .95, .99),
  pois_sim = TRUE
)
# Save to file
fwrite(error_provmonth, file = glue::glue('{pv_dir}/{fp_prefix}_error_provmonth.csv'))
fwrite(coverage_provmonth, file = glue::glue('{pv_dir}/{fp_prefix}_coverage_provmonth.csv'))


data_by_prov_month[, `:=` (obs_plot = obs * 1E5, pred_mean_plot = pred_mean * 1E5)]
xybreaks <- seq(10, 50, by=10)

obs_scatter_prov_month <- ggplot(data=data_by_prov_month[year < 2020 | wkgrp < 3,]) +
  geom_point(
    aes(x=obs_plot, y=pred_mean_plot, size=sqrt(pop), color=as.character(year)),
    alpha = .2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = '#AAAAAA') +
  scale_color_manual(values=plot_colors) +
  scale_x_continuous(limits = range(xybreaks), breaks = xybreaks, labels = xybreaks, trans='sqrt') +
  scale_y_continuous(limits = range(xybreaks), breaks = xybreaks, labels = xybreaks, trans='sqrt') +
  scale_size_continuous(guide = 'none') +
  labs(
    title = paste0('Mortality rate among ',run_sex,'s'),
    subtitle = 'Aggregated by month and province',
    x = 'Observed', y = 'Predicted (mean)', color = 'Year',
    size = 'Population\nsize'
  ) +
  theme_dark()
png(
  filename = glue::glue('{viz_dir}/{fp_prefix}_scatter_provmonth.png'),
  height = 6, width = 6, units = 'in', res = 100
)
print(obs_scatter_prov_month)
dev.off()
