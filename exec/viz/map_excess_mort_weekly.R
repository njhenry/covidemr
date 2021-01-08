## -----------------------------------------------------------------------------
##
## VISUALIZATION: Map excess mortality by week
##
## -----------------------------------------------------------------------------

## Load required packages and inputs

# library(covidemr)
library(argparse)
library(data.table)
library(RColorBrewer)
library(sp)
library(sf)
library(tictoc)
# DEVELOPMENT: rebuild library
dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Run settings
parser <- argparse::ArgumentParser(
  description='COVID Excess Mortality: viz script'
)
parser$add_argument('--data-version', type='character', help='Prepped data version')
parser$add_argument('--model-version', type='character', help='Model run date')
args <- parser$parse_args(commandArgs(TRUE))
prepped_data_version <- args$data_version
model_run_version <- args$model_version

# Helper function to create a filepath for a particular prepped data object
prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
viz_dir <- file.path(config$paths$model_results, model_run_version, 'viz')
dir.create(viz_dir, showWarnings = FALSE)
get_prep_fp <- function(ff) file.path(prep_dir, config$prepped_data_files[[ff]])

# Load input data
location_table <- data.table::fread(get_prep_fp('location_table'), na.strings="")
template_dt <- fread(get_prep_fp('template'))
data_full <- fread(get_prep_fp('full_data'))
pop_dt <- fread(get_prep_fp('population'))
polys_sf <- readRDS(get_prep_fp('shapefile_sf'))
covid_mort <- fread(get_prep_fp('covid_deaths'))

# Merge extra location identifiers onto the spatial polygons
merge_on_cols <- c('location_code', 'location_name', 'region_code', 'region_name')
polys_sf <- merge(polys_sf, location_table[, ..merge_on_cols], by='location_code')
# Create a regional sf object by dissolving provinces
regions_sf <- aggregate(polys_sf, by=list(polys_sf$region_code), FUN=first)


# Set week info
start_week <- covidemr::week_id_from_date(config$first_death_date)
end_week <- covidemr::week_id_from_date(config$final_obs_date)
ndraw <- config$num_draws
draw_col_names <- paste0('draw', 1:ndraw)

data_full[ year < 2020 | week < start_week, in_baseline := 1 ]
data_full[ year == 2020 & week >= start_week & is.na(in_baseline), in_baseline := 0 ]

# Prepare a table of week start dates for plotting purposes
wk_starts <- data.table(
  year = 2020, week = start_week:end_week
)
wk_starts[, day_id := (week - 1) * 7
  ][, date := as.Date(day_id, origin = '2020-01-01')
  ][, date_plot := strftime(date, format = '%d %b %Y')]


## Load draw-level data --------------------------------------------------------

# Start timer for this section
tictoc::tic('Loading data')

draws_list <- list(male = NULL, female = NULL)
subpop_cols <- c('location_code', 'age_group_code', 'year', 'sex')

for(run_sex in c('male', 'female')){
  templ_sub <- template_dt[sex == run_sex, ]
  keep_idx <- which(templ_sub[,year == 2020 & week >= start_week & week <= end_week])
  # Combine relevant index and draws
  id_cols <- c(subpop_cols, 'week')
  draws_list[[run_sex]] <- cbind(
    templ_sub[keep_idx, ..id_cols],
    fread(glue::glue(
      config$results_files$pred_draws,
      .envir = list(
        model_results_dir = config$paths$model_results,
        model_run_version = model_run_version,
        run_sex = run_sex,
        holdout = 0
      )
    ))[keep_idx, ]
  )
  setnames(draws_list[[run_sex]], paste0('V', 1:ndraw), draw_col_names)
  # Merge on population and deaths
  draws_list[[run_sex]][
    data_full[in_baseline == 0, ],
    `:=` (deaths = i.deaths, pop = i.pop, observed_days = i.observed_days),
    on = id_cols
  ]
}

exp_draws <- rbindlist(draws_list)
rm(draws_list)

# Adjust mortality draws for weeks that were not completely observed
exp_draws[,
  (draw_col_names) := lapply(.SD, function(col) col / 7 * observed_days),
  .SDcols = draw_col_names
]

# End timer for this section
tictoc::toc()


## Calculate most-detailed SMRs and excess deaths ------------------------------

tictoc::tic("Calculating excess mortality")

most_detailed_list <- calculate_excess_time_series_by_group(
  experiment_draw_dt = exp_draws,
  baseline_draw_cols = draw_col_names,
  group_cols = subpop_cols,
  week_col = 'week',
  obs_death_col = 'deaths',
  pop_col = 'pop'
)
smrs_dt <- most_detailed_list$smrs
excess_deaths_dt <- most_detailed_list$excess_deaths
rm(most_detailed_list)

tictoc::toc()


## Aggregate to the province-week level and estimate excess --------------------

agg_draws <- aggregate_data_and_draws(
  in_data = exp_draws[ !is.na(pop), ], num_field = 'deaths', denom_field = 'pop',
  draw_fields = draw_col_names, group_fields = c('location_code', 'year', 'week'),
  summarize = TRUE
)
agg_draws[
  pop_dt[, .(pop = sum(pop)), by = .(location_code, year)],
  pop := i.pop,
  on = c('location_code', 'year')
]
agg_summ <- agg_draws[, -draw_col_names, with=FALSE]

# Estimate excess
agg_summ[, obs := deaths / pop
  ][, `:=` (
      ratio_mean=obs/pred_mean, ratio_lower=obs/pred_upper, ratio_upper=obs/pred_lower,
      ex_mean=deaths-pred_mean*pop, ex_lower=deaths-pred_upper*pop,
      ex_upper=deaths-pred_lower*pop
    )
  ][, `:=` (sig_under1 = as.integer(ratio_upper < 1), sig_over1 = as.integer(ratio_lower > 1))
  ][
    location_table,
    `:=` (location_name = i.location_name, region_name = i.region_name),
    on = 'location_code'
  ]


## Map excess mortality rate by week -------------------------------------------

fill_grps <- c(
  'Signif. lower', 'Not signif.', '+ 0% to 20%', '+ 20 to 50%', '+ 50% to 100%',
  '+ 100 to 200%', '+ >200%'
)
fill_colors <- c('#b3ffec','#b3b3b3','#FFFFB2','#FECC5C','#FD8D3C','#F03B20','#BD0026')
names(fill_colors) <- fill_grps

for(week_id in start_week:end_week){
  plot_title <- paste("Week of", wk_starts[week == week_id, date_plot ])
  plot_dt <- agg_summ[week == week_id, ]
  # Assign fill colors
  plot_dt[, fill_grp := fill_grps[2] ]
  plot_dt[sig_under1==1, fill_grp:=fill_grps[1] ]
  plot_dt[(sig_over1==1) & (ratio_mean<1.2), fill_grp:=fill_grps[3] ]
  plot_dt[(sig_over1==1) & (ratio_mean>=1.2) & (ratio_mean<1.5), fill_grp:=fill_grps[4] ]
  plot_dt[(sig_over1==1) & (ratio_mean>=1.5) & (ratio_mean<2.0), fill_grp:=fill_grps[5] ]
  plot_dt[(sig_over1==1) & (ratio_mean>=2.0) & (ratio_mean<3.0), fill_grp:=fill_grps[6] ]
  plot_dt[(sig_over1==1) & (ratio_mean>=3.0), fill_grp:=fill_grps[7] ]
  # Plot it
  covidemr::map_ita_choropleth(
    province_sf = polys_sf, region_sf = regions_sf,
    in_data = plot_dt,
    map_field = 'fill_grp',
    fill_list = list(
      values = fill_colors, limits = rev(fill_grps)
    ),
    titles_list = list(title = plot_title, fill = 'Excess Mortality'),
    fill_type = 'manual',
    save_fp = sprintf("%s/map_excess_rate_week%02d.png", viz_dir, week_id)
  )
}

## Aggregate to the national level and plot excess -----------------------------

agg_draws_week <- aggregate_data_and_draws(
  in_data = exp_draws[ !is.na(pop) , ], num_field = 'deaths', denom_field = 'pop',
  draw_fields = draw_col_names, group_fields = c('year', 'week'),
  summarize = TRUE
)
agg_draws_week[pop_dt[, .(pop=sum(pop)), by=year], pop := i.pop, on='year']
# Convert from rates to counts
for(dcol in c(draw_col_names, 'pred_mean', 'pred_lower', 'pred_upper')){
  agg_draws_week[[dcol]] <- agg_draws_week[[dcol]] * agg_draws_week$pop
}
agg_summ_week <- agg_draws_week[, -draw_col_names, with = FALSE ]

## Plot it
to_plot_natl_weekly <- merge(x=agg_summ_week, y=wk_starts[, .(week, date)], by='week')

mort_plot_colors <- c(
  'COVID-19' = '#336600',
  'Predicted' = '#993300',
  'Observed' = '#000000',
  'Excess (all-cause)' = '#F03B20'
)

date_breaks <- as.Date(sprintf('2020-%02d-01',3:8))
date_labs <- strftime(date_breaks, '%b %d')
y_breaks <- seq(0, 2.5E4, by=5E3)
y_labs <- c('0', '5,000', '10,000', '15,000', '20,000', '25,000')


natl_weekly_plot <- ggplot(data=to_plot_natl_weekly, aes(x=date)) +
  geom_ribbon(
    aes(ymin=pred_lower, ymax=pred_upper),
    color = NA, fill = mort_plot_colors['Predicted'], alpha = .5
  ) +
  geom_line(aes(y=pred_mean, color='Predicted'), lwd = .75) +
  geom_line(aes(y=deaths, color='Observed'), lwd = .75) +
  geom_vline(xintercept=wk_starts[week==9,date], color='#888888', linetype=2) +
  geom_vline(xintercept=wk_starts[week==22,date], color='#888888', linetype=2) +
  scale_x_date(limits = range(to_plot_natl_weekly$date), breaks=date_breaks, labels=date_labs) +
  scale_y_continuous(limits = range(y_breaks), breaks = y_breaks, labels = y_labs) +
  scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors)) +
  labs(
    x = 'Week starting', y = 'Deaths',
    title = 'Italy: all-cause mortality',
    subtitle = 'March through August 2020',
    color = '') +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position='bottom')
png(
  file.path(viz_dir, 'national_all_cause_march_to_august.png'),
  height=5, width=8, units='in', res=200
)
plot(natl_weekly_plot)
dev.off()


## Plot COVID deaths and excess mortality by week
natl_covid_and_excess <- merge(
  x = to_plot_natl_weekly[ week %in% 9:22, ],
  y = covid_mort[, .(covid_deaths = sum(covid_deaths)), by=week],
  by = 'week'
)[, `:=` (ex_mean=deaths-pred_mean, ex_lower=deaths-pred_upper, ex_upper=deaths-pred_lower)]

y_breaks <- seq(0, 1.2E4, by=2E3)
y_labs <- c('0 (baseline)', sprintf('%d,000', setdiff(y_breaks, 0)/1E3))

excess_natl_plot <- ggplot(data=natl_covid_and_excess, aes(x=date)) +
  geom_hline(yintercept = 0, color='black', lwd=.5) +
  geom_ribbon(
    aes(ymin=ex_lower, ymax=ex_upper),
    color = NA, fill = mort_plot_colors['Excess (all-cause)'], alpha = .5
  ) +
  geom_line(aes(y=covid_deaths, color='COVID-19'), lwd = .75) +
  geom_line(aes(y=ex_mean, color='Excess (all-cause)'), lwd = .75) +
  scale_x_date(limits = range(natl_covid_and_excess$date), breaks=date_breaks, labels=date_labs) +
  scale_y_continuous(
    limits = c(min(natl_covid_and_excess$ex_lower), max(natl_covid_and_excess$ex_upper)),
    breaks = y_breaks, labels = y_labs
  ) +
  scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors)) +
  labs(
    x='Week starting', y='Deaths',
    title='Italy: excess mortality and COVID-19 deaths',
    subtitle = 'March through May 2020',
    color='') +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position='bottom')
png(
  file.path(viz_dir, 'national_excess_march_to_may.png'),
  height=5, width=8, units='in', res=200
)
plot(excess_natl_plot)
dev.off()

# Plot with a weekly intercept
for(week_id in 9:21){
  week_start_date <- wk_starts[week == week_id, date]
  plot_subtitle <- paste("Week of", wk_starts[week == week_id, date_plot])
  excess_plot_weekly <- suppressWarnings(
    excess_natl_plot +
      geom_vline(xintercept = week_start_date, color='#666666', linetype=2) +
      labs(subtitle = plot_subtitle)
  )
  png(
    sprintf('%s/national_excess_week%02d.png', viz_dir, week_id),
    height=5, width=8, units='in', res=200
  )
  plot(excess_plot_weekly)
  dev.off()
}


## Convert to draws of counts for aggregation across weeks ---------------------

agg_overall <- copy(agg_draws_week)[ week %in% 9:21, ]
agg_overall <- agg_overall[, lapply(.SD, sum), .SDcols=c(draw_col_names, 'deaths')]
agg_overall$pred_mean <- rowMeans(agg_overall[, ..draw_col_names])
agg_overall$pred_lower <- rowQuantiles(as.matrix(agg_overall[, ..draw_col_names]), probs=.025)
agg_overall$pred_upper <- rowQuantiles(as.matrix(agg_overall[, ..draw_col_names]), probs=.975)
agg_overall[, `:=` (
  excess_mean = deaths - pred_mean,
  excess_lower = deaths - pred_upper,
  excess_upper = deaths - pred_lower
)]
agg_summ_overall <- agg_overall[, -draw_col_names, with = FALSE ]

# Get cumulative deaths across experimental weeks
agg_wk_cumul <- copy(agg_draws_week)
for(dcol in c(draw_col_names, 'deaths')){
  agg_wk_cumul[[dcol]] <- cumsum(agg_wk_cumul[[dcol]])
}
agg_wk_cumul$pred_mean <- rowMeans(agg_wk_cumul[, ..draw_col_names])
agg_wk_cumul$pred_lower <- rowQuantiles(as.matrix(agg_wk_cumul[, ..draw_col_names]), probs=.025)
agg_wk_cumul$pred_upper <- rowQuantiles(as.matrix(agg_wk_cumul[, ..draw_col_names]), probs=.975)
agg_wk_cumul[, `:=` (
  excess_mean = deaths - pred_mean,
  excess_lower = deaths - pred_upper,
  excess_upper = deaths - pred_lower
)]
agg_wk_cumul_summ <- agg_wk_cumul[, -draw_col_names, with = FALSE ]


## Aggregate deaths by province across the 9-21 week range ---------------------

agg_counts <- copy(agg_draws)[ year == 2020 & week %in% 9:21, ]
for(dcol in draw_col_names) agg_counts[, (dcol) := get(dcol) * pop ]
agg_counts$ct_mean <- rowMeans(agg_counts[, ..draw_col_names])
agg_counts$ct_lower <- rowQuantiles(as.matrix(agg_counts[, ..draw_col_names]), probs=.025)
agg_counts$ct_upper <- rowQuantiles(as.matrix(agg_counts[, ..draw_col_names]), probs=.975)

# Aggregate by province
agg_by_prov <- agg_counts[,
  lapply(.SD, sum),
  .SDcols = c('deaths', 'ct_mean', 'ct_lower', 'ct_upper', draw_col_names),
  by = location_code
][pop_dt[year==2020, .(pop=sum(pop)), by=location_code], pop:=i.pop, on='location_code']
# Get ranking by excess counts
agg_by_prov_top90 <- agg_by_prov[, -draw_col_names, with=FALSE
  ][, `:=` (ex_mean=deaths-ct_mean, ex_lower=deaths-ct_upper, ex_upper=deaths-ct_lower)
  ][order(-ex_mean)
  ][, cumul_ex := cumsum(ex_mean) / sum(ex_mean) ]
agg_by_prov_top90[, rank := .I ]

# The top 26 provinces account for >75% of all excess mortality
# The top 9 provinces account for >50% of all excess mortality
# The top 3 provinces account for nearly 30% of all excess mortality
agg_by_prov_top90[, top_grouping := 'Other'
  ][rank <= 26, top_grouping := 'Top 26 provinces'
  ][rank <= 9, top_grouping := 'Top 9 provinces'
  ][rank <= 3, top_grouping := 'Top 3 provinces' ]

top_prov_dt <- agg_by_prov_top90[, .(location_code, ex_mean, rank, top_grouping, pop)]
top_prov_dt <- merge(location_table, top_prov_dt, by='location_code')[order(rank)]
top_prov_dt[, cumul_prop_pop := cumsum(pop)/sum(pop) ]
top_prov_dt[, cumul_prop_ex := cumsum(ex_mean)/sum(ex_mean) ]
fwrite(top_prov_dt, file=file.path(viz_dir, 'top_provinces_by_excess.csv'))


## Map the top provinces by excess
prov_grps <- c('Top 3 provinces', 'Top 9 provinces', 'Top 26 provinces', 'Other')
prov_colors <- c('#120d32', '#7d2482', '#fea873', '#b3b3b3')
names(prov_colors) <- prov_grps

covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = top_prov_dt,
  map_field = 'top_grouping',
  fill_list = list(
    values = prov_colors, limits = rev(prov_grps)
  ),
  titles_list = list(title = '', fill = 'Ranking by\nExcess Deaths'),
  fill_type = 'manual',
  save_fp = file.path(viz_dir, "province_top_burden_map.png")
)


## Plot the excess curve of the top provinces
weekly_by_province_data <- merge(
  agg_counts[, .(location_code, year, week, ct_mean, deaths)],
  top_prov_dt[, .(location_code, top_grouping)],
  by='location_code'
)[, .(ex_mean = sum(deaths) - sum(ct_mean)), by = .(top_grouping, week)]
weekly_by_province_data <- merge(
  x = weekly_by_province_data,
  y = data.table(top_grouping = as.character(prov_grps), grp_order = 1:4),
  by = 'top_grouping'
)[ order(week, grp_order)
  ][, ex_to_plot := cumsum(ex_mean), by=week
  ][ wk_starts[, .(week, date)], date := i.date, on='week'
  ][ order(grp_order) ][
  ][, top_grouping := factor(top_grouping, levels = rev(prov_grps)) ]


excess_by_prov_grp_plot <- ggplot(data = weekly_by_province_data, aes(x=date)) +
  geom_hline(yintercept = 0, color='black', lwd=.5) +
  geom_ribbon(aes(fill=top_grouping, ymax=ex_to_plot), ymin=0, color=NA) +
  scale_x_date(limits = range(natl_covid_and_excess$date), breaks=date_breaks, labels=date_labs) +
  scale_y_continuous(
    limits = c(min(natl_covid_and_excess$ex_lower), max(natl_covid_and_excess$ex_upper)),
    breaks = y_breaks, labels = y_labs
  ) +
  scale_fill_manual(values=prov_colors, breaks=names(prov_colors)) +
  labs(
    x='Week starting', y='Deaths',
    title='Italy: excess mortality by province grouping',
    subtitle = 'March through May 2020',
    fill = ''
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position='bottom')
png(
  file.path(viz_dir, 'province_top_burden_lines.png'),
  height=5, width=8, units='in', res=200
)
plot(excess_by_prov_grp_plot)
dev.off()


## Get aggregates by age -------------------------------------------------------

fill_grps <- c(
  'Signif. lower', 'Not signif.', '+ 0% to 20%', '+ 20 to 50%', '+ 50% to 100%',
  '+ 100 to 200%', '+ >200%'
)
fill_colors <- c('#b3ffec','#b3b3b3','#FFFFB2','#FECC5C','#FD8D3C','#F03B20','#BD0026')
names(fill_colors) <- fill_grps

age_groups_dt <- covidemr::create_age_groups(config$age_cutoffs)

## Get deviation from baseline mortality rate by province and age group in March
ag_subset_weeks <- 9:21

agg_by_prov_age_march <- aggregate_data_and_draws(
  in_data = exp_draws[ !is.na(pop) & week %in% ag_subset_weeks, ],
  num_field = 'deaths', denom_field = 'pop',
  draw_fields = draw_col_names,
  group_fields = c('location_code', 'year', 'age_group_code'),
  summarize = TRUE
)
agg_by_prov_age_march <- agg_by_prov_age_march[
    pop_dt[year == 2020, ], pop := i.pop, on = c('location_code', 'age_group_code')
  ][, obs := (deaths / length(ag_subset_weeks)) / pop # Get average across 5 weeks
  ][, `:=` (
      ratio_mean = obs / pred_mean,
      ratio_upper = obs / pred_lower,
      ratio_lower = obs / pred_upper
    )
  ][, `:=` (sig_over1 = as.numeric(ratio_lower>1), sig_under1 = as.numeric(ratio_upper<1))
  ][, fill_grp := fill_grps[2]
  ][ sig_under1==1, fill_grp:=fill_grps[1]
  ][ (sig_over1==1) & (ratio_mean<1.2), fill_grp:=fill_grps[3]
  ][ (sig_over1==1) & (ratio_mean>=1.2) & (ratio_mean<1.5), fill_grp:=fill_grps[4]
  ][ (sig_over1==1) & (ratio_mean>=1.5) & (ratio_mean<2.0), fill_grp:=fill_grps[5]
  ][ (sig_over1==1) & (ratio_mean>=2.0) & (ratio_mean<3.0), fill_grp:=fill_grps[6]
  ][ (sig_over1==1) & (ratio_mean>=3.0), fill_grp:=fill_grps[7]
  ][, -draw_col_names, with=FALSE
  ][
    age_groups_dt,
    age_group_name := i.age_group_name,
    on = 'age_group_code'
  ]

## Create four plots by age group
for(this_ag in unique(agg_by_prov_age_march$age_group_name)){
  covidemr::map_ita_choropleth(
    province_sf = polys_sf, region_sf = regions_sf,
    in_data = agg_by_prov_age_march[age_group_name == this_ag, ],
    map_field = 'fill_grp',
    fill_list = list(
      values = fill_colors, limits = rev(fill_grps)
    ),
    titles_list = list(title = this_ag, fill = 'Excess Mortality'),
    fill_type = 'manual',
    save_fp = sprintf("%s/ag_excess_rate_map_%s.png", viz_dir, gsub(" ","_",this_ag))
  )
}


## Aggregate by age group and week, then create line plot

age_week_dt <- copy(exp_draws[ !is.na(pop) & week %in% 9:22 , ])
for(dcol in draw_col_names) age_week_dt[, (dcol) := get(dcol) * pop ]
age_week_dt <- age_week_dt[,
    lapply(.SD, sum), .SDcols = c('deaths', draw_col_names),
    by=c('age_group_code', 'week')
  ][age_groups_dt, age_group_name := i.age_group_name, on='age_group_code']

age_week_dt$pred_mean <- rowMeans(age_week_dt[, ..draw_col_names])
age_week_dt$pred_lower <- rowQuantiles(as.matrix(age_week_dt[, ..draw_col_names]), probs=.025)
age_week_dt$pred_upper <- rowQuantiles(as.matrix(age_week_dt[, ..draw_col_names]), probs=.975)
age_week_dt <- age_week_dt[,
    `:=` (ex_mean=deaths-pred_mean, ex_lower=deaths-pred_upper, ex_upper=deaths-pred_lower)
  ][, -draw_col_names, with = FALSE
  ][wk_starts, date := i.date, on='week' ]

excess_age_week_plot <- ggplot(data=age_week_dt, aes(x=date)) +
  facet_wrap('age_group_name', ncol=2, scales = 'fixed') +
  geom_hline(yintercept = 0, color='black', lwd=.5) +
  geom_ribbon(
    aes(ymin=ex_lower, ymax=ex_upper),
    color = NA, fill = mort_plot_colors['Excess (all-cause)'], alpha = .5
  ) +
  geom_line(aes(y=ex_mean, color='Excess (all-cause)'), lwd = .75) +
  scale_x_date(limits = range(age_week_dt$date), breaks=date_breaks, labels=date_labs) +
  scale_y_continuous(
    limits = c(min(age_week_dt$ex_lower), max(age_week_dt$ex_upper)),
    breaks = c(0, 5E3, 1E4), labels = c('0', '5,000', '10,000')
  ) +
  scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors), guide=FALSE) +
  labs(
    x='Week starting', y='Deaths',
    title='Italy: excess mortality by age group',
    subtitle = 'March through May 2020',
    color='') +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position='bottom')
png(
  file.path(viz_dir, 'ag_excess_lines_plot.png'),
  height=6, width=10, units='in', res=300
)
plot(excess_age_week_plot)
dev.off()


## Create functions to generate excess and SMRs, then save to file
agg_ex_smrs_save <- function(dt, group_cols, out_fp, printout = TRUE){
  dt_agg <- dt[,
    lapply(.SD, sum), .SDcols = c('deaths', draw_col_names), by = group_cols
  ]
  dt_agg$cf_mean <- rowMeans(dt_agg[, ..draw_col_names])
  dt_agg$cf_lower <- rowQuantiles(as.matrix(dt_agg[, ..draw_col_names]), probs=.025)
  dt_agg$cf_upper <- rowQuantiles(as.matrix(dt_agg[, ..draw_col_names]), probs=.975)
  dt_agg[, `:=` (
    ex_mean=deaths-cf_mean, ex_lower=deaths-cf_upper, ex_upper=deaths-cf_lower,
    smr_mean=deaths/cf_mean, smr_lower=deaths/cf_upper, smr_upper=deaths/cf_lower
  )]
  keep_cols <- c(
    group_cols, 'ex_mean','ex_lower','ex_upper','smr_mean','smr_lower','smr_upper'
  )
  fwrite(dt_agg[, ..keep_cols], file = out_fp)
  if(printout){
    knitr::kable(dt_agg[, ..keep_cols])
  } else {
    invisible()
  }
}

# Get overall excess with UIs by age group and sex *nationally*
age_sex_dt <- copy(exp_draws[ !is.na(pop) & week %in% 9:21 , ])
for(dcol in draw_col_names) age_sex_dt[, (dcol) := get(dcol) * pop * observed_days / 7 ]
age_sex_dt[
    location_table,
    `:=` (region_code = i.region_code, region_name = i.region_name, location_name = i.location_name),
    on = 'location_code'
  ][age_groups_dt, age_group_name := i.age_group_name, on='age_group_code'
  ][, dummy := 1 ]

agg_ex_smrs_save(age_sex_dt, group_cols=c('dummy'), out_fp=file.path(viz_dir,'natl_excess_all.csv'))
agg_ex_smrs_save(age_sex_dt, group_cols=c('age_group_name','sex'), out_fp=file.path(viz_dir,'natl_excess_by_age_sex.csv'))
agg_ex_smrs_save(age_sex_dt, group_cols=c('sex'), out_fp=file.path(viz_dir,'natl_excess_by_sex.csv'))
agg_ex_smrs_save(age_sex_dt, group_cols=c('age_group_name'), out_fp=file.path(viz_dir,'natl_excess_by_age.csv'))

agg_ex_smrs_save(
  age_sex_dt, group_cols=c('region_code','region_name','age_group_name','sex'),
  out_fp=file.path(viz_dir,'regional_excess_by_age_sex.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('region_code','region_name','sex'),
  out_fp=file.path(viz_dir,'regional_excess_by_sex.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('region_code','region_name','age_group_name'),
  out_fp=file.path(viz_dir,'regional_excess_by_age.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('region_code','region_name'),
  out_fp=file.path(viz_dir,'regional_excess.csv'), printout = TRUE
)

agg_ex_smrs_save(
  age_sex_dt, group_cols=c('location_code','location_name','age_group_name','sex'),
  out_fp=file.path(viz_dir,'provincial_excess_by_age_sex.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('location_code','location_name','sex'),
  out_fp=file.path(viz_dir,'provincial_excess_by_sex.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('location_code','location_name','age_group_name'),
  out_fp=file.path(viz_dir,'provincial_excess_by_age.csv'), printout = FALSE
)
agg_ex_smrs_save(
  age_sex_dt, group_cols=c('location_code','location_name'),
  out_fp=file.path(viz_dir,'provincial_excess.csv'), printout = TRUE
)

# Get excess and SMRs only for the north of Italy (regions 1-9)
northern_subset <- age_sex_dt[region_code %in% 1:9, ]
agg_ex_smrs_save(northern_subset, group_cols='dummy', out_fp=file.path(viz_dir,'northern_excess_overall.csv'))
agg_ex_smrs_save(northern_subset, group_cols=c('age_group_name','sex'), out_fp=file.path(viz_dir,'northern_excess_by_age_sex.csv'))
agg_ex_smrs_save(northern_subset, group_cols=c('sex'), out_fp=file.path(viz_dir,'northern_excess_by_sex.csv'))
agg_ex_smrs_save(northern_subset, group_cols=c('age_group_name'), out_fp=file.path(viz_dir,'northern_excess_by_age.csv'))


## Make regional map showing COVID deaths per 100k population prior to Aug 31 --

# Get population by region
pop_reg <- merge(
  x = location_table[, .(region_code, location_code)],
  y = pop_dt[year==2020, .(location_code, pop)],
  by = 'location_code'
)[, .(pop = sum(pop)), by=region_code]

# Get COVID deaths per 100k population prior to August 31
covid_deaths_per_100k <- covid_mort[
    week <= 35, .(covid_deaths = sum(covid_deaths)), by=region_code
  ][ pop_reg, pop := i.pop, on = 'region_code'
  ][, deaths_per_100k := covid_deaths / pop * 1E5 ]

message(covid_mort[week <= 35, sum(covid_deaths)], ' COVID deaths prior to August 31')
message('Of those, ',covid_mort[week < 22, sum(covid_deaths)],' COVID deaths prior to May 27')
message(
  "During those weeks, registered COVID deaths made up ",
  round(covid_mort[week < 22, sum(covid_deaths)] / 6E7 * 1E5),
  " deaths per 100k population"
)
message("Mortality was highest in Lombardy, with nearly 17k deaths")
message("This is equivalent to 167 deaths per 100k population")

covidemr::map_ita_choropleth_region(
  region_sf = regions_sf, in_data = covid_deaths_per_100k,
  map_field='deaths_per_100k',
  fill_list = list(
    limits = c(0, 167),
    breaks = c(10, 20, 50, 100, 150),
    labels = c('< 10', '20', '50', '100', '> 150'),
    colors = RColorBrewer::brewer.pal('Greens', n=9)[c(1, 4:9)]
  ),
  fill_type = 'gradientn',
  titles_list = list(
    title = 'Italy: COVID mortality by region prior to 31 August',
    fill = 'Deaths\nper 100,000\npopulation'),
  save_fp = file.path(viz_dir, 'COVID_deaths_per_100k.png')
)


## Plot excess mortality lines by week for all provinces -----------------------

for(loc_code in unique(location_table$location_code)){
  this_loc_name <- location_table[location_code==loc_code, location_name]
  this_reg_name <- location_table[location_code==loc_code, region_name]
  plot_title <- sprintf('%s (%s)', this_loc_name, this_reg_name)
  plot_fp <- glue::glue(
    '{viz_dir}/em_lines_prov_{sprintf("%03d",loc_code)}_',
    '{gsub(" ","_",this_loc_name)}_{gsub(" ","_",this_reg_name)}.png'
  )

  to_plot_dt <- agg_summ[location_code == loc_code & week %in% 9:22,
    ][wk_starts, date := i.date, on='week']

  prov_line_plot <- ggplot(data=to_plot_dt, aes(x=date)) +
    geom_hline(yintercept = 0, color='black', lwd=.5) +
    geom_ribbon(
      aes(ymin=ex_lower, ymax=ex_upper),
      color = NA, fill = mort_plot_colors['Excess (all-cause)'], alpha = .5
    ) +
    geom_line(aes(y=ex_mean, color='Excess (all-cause)'), lwd = .75) +
    scale_x_date(limits = range(to_plot_dt$date), breaks=date_breaks, labels=date_labs) +
    scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors), guide=FALSE) +
    labs(
      x='Week starting', y='Deaths',
      title=plot_title,
      color='') +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  png(plot_fp, height=4, width=6, units='in', res=300)
  plot(prov_line_plot)
  dev.off()
}


## Map the week of peak excess across northern provinces -----------------------

peak_excess_dt <- agg_summ[ week %in% 9:21,
  ][, week_rank := frank(-ex_mean), by = location_code
  ][ week_rank == 1,
  ][ wk_starts, date := i.date, on='week'
  ][location_table, region_code := i.region_code, on = 'location_code']

wk_labs <- c(
  '11 Mar or prior', '18 Mar', '25 Mar', '01 Apr', '08 Apr', '15 Apr or later')
wk_colors <- RColorBrewer::brewer.pal('PuOr', n=11)[c(2,3,5,7,9,10)]
names(wk_colors) <- wk_labs

peak_excess_dt[
    week <= 11, week_label := wk_labs[1]
  ][ week == 12, week_label := wk_labs[2]
  ][ week == 13, week_label := wk_labs[3]
  ][ week == 14, week_label := wk_labs[4]
  ][ week == 15, week_label := wk_labs[5]
  ][ week >= 16, week_label := wk_labs[6]]

covidemr::map_ita_choropleth(
  province_sf = polys_sf[polys_sf$region_code %in% 1:9, ],
  region_sf = regions_sf,
  in_data = peak_excess_dt[region_code %in% 1:9, ],
  map_field = 'week_label',
  fill_list = list(
    values = wk_colors, limits = rev(wk_labs)
  ),
  titles_list = list(
    title = 'Week of peak excess mortality',
    subtitle = 'March through May 2020',
    fill = 'Week beginning'),
  fill_type = 'manual',
  save_fp = file.path(viz_dir, 'peak_excess_week_map.png')
)


## Map the maximum excess across northern provinces ----------------------------

fill_grps <- c(
  '+ 0% to 50%', '+ 50% to 100%', '+ 100% to 200%', '+ 200% to 300%',
  '+ 300% to 400%', '+ 400% to 500%', '+ >500%'
)
fill_colors <- c(RColorBrewer::brewer.pal('OrRd', n=9)[c(3,5:9)], '#2e051a')
names(fill_colors) <- fill_grps
fill_brks <- c(1., 1.5, 2., 3., 4., 5., 6., 100)

peak_excess_dt[, fill_grp := fill_grps[cut(ratio_mean, fill_brks, labels=FALSE)]]


covidemr::map_ita_choropleth(
  province_sf = polys_sf[polys_sf$region_code %in% 1:9, ],
  region_sf = regions_sf,
  in_data = peak_excess_dt[region_code %in% 1:9, ],
  map_field = 'fill_grp',
  fill_list = list(
    values = fill_colors, limits = rev(fill_grps)
  ),
  titles_list = list(
    title = 'Peak excess mortality over baseline',
    subtitle = 'March through May 2020',
    fill = 'Excess Mortality'),
  fill_type = 'manual',
  save_fp = file.path(viz_dir, 'peak_excess_amount_map.png')
)
