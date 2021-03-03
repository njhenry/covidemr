## -----------------------------------------------------------------------------
##
## VISUALIZATION: Map excess mortality by week
##
## -----------------------------------------------------------------------------

## Load required packages and inputs

# library(covidemr)
library(argparse)
library(data.table)
library(gridExtra)
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
top_prov_labs <- fread(get_prep_fp('province_labels'))
top_reg_labs <- fread(get_prep_fp('region_labels'))
top_reg_labs[, label_name := gsub('\\\\n','\n', label_name)]


# Merge extra location identifiers onto the spatial polygons
merge_on_cols <- c('location_code', 'location_name', 'region_code', 'region_name')
polys_sf <- merge(polys_sf, location_table[, ..merge_on_cols], by='location_code')
# Create a regional sf object by dissolving provinces
regions_sf <- aggregate(polys_sf, by=list(polys_sf$region_code), FUN=first)


# Set week info
start_week <- covidemr::week_id_from_date(config$first_death_date)
end_week <- covidemr::week_id_from_date(config$final_obs_date)
focus_end_week <- 21 # Final week for the period highlighted in the paper

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

# Create age group metadata
age_groups_dt <- covidemr::create_age_groups(config$age_cutoffs)


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
  .SDcols = c(draw_col_names)
]

# End timer for this section
tictoc::toc()


## Calculate most-detailed SMRs and excess deaths ------------------------------

tictoc::tic("Calculating most-detailed excess mortality")

most_detailed_list <- calculate_excess_time_series_by_group(
  experiment_draw_dt = exp_draws,
  baseline_draw_cols = draw_col_names,
  group_cols = subpop_cols,
  week_col = 'week',
  obs_death_col = 'deaths',
  pop_col = 'pop'
)
baseline_deaths_dt <- most_detailed_list$baseline_deaths
smrs_dt <- most_detailed_list$smrs
excess_deaths_dt <- most_detailed_list$excess_deaths
rm(most_detailed_list)

tictoc::toc()


## Define fill colors and breaks -----------------------------------------------

# Breaks used for most choropleths
# fill_grps <- c(
#   'Signif. lower', 'Not signif.', '+ 0% to 20%', '+ 20 to 50%', '+ 50% to 100%',
#   '+ 100 to 200%', '+ >200%'
# )
# fill_colors <- c('#b3ffec','#b3b3b3','#FFFFB2','#FECC5C','#FD8D3C','#F03B20','#BD0026')
fill_grps <- c(
  'Not signif.\nhigher', '+ 0% to 20%', '+ 20 to 50%', '+ 50% to 100%', '+ 100 to 200%', '+ >200%'
)
fill_colors <- c('#b3b3b3','#FFFFB2','#FECC5C','#FD8D3C','#F03B20','#BD0026')
names(fill_colors) <- fill_grps

# Breaks used when highlighting northern provinces
fill_grps_highlight <- c(
  '+ 0% to 50%', '+ 50% to 100%', '+ 100% to 200%', '+ 200% to 300%',
  '+ 300% to 400%', '+ 400% to 500%', '+ >500%'
)
fill_colors_highlight <- c(RColorBrewer::brewer.pal('OrRd', n=9)[c(3,5:9)], '#2e051a')
names(fill_colors_highlight) <- fill_grps_highlight
fill_brks_highlight <- c(1., 1.5, 2., 3., 4., 5., 6., 1E8)


## Aggregate across various dimensions and estimate excess ---------------------

tictoc::tic("Aggregating excess mortality for plotting")

# Add additional identifiers onto the baseline deaths data.table
baseline_deaths_dt <- merge(
  x = baseline_deaths_dt,
  y = location_table[, .(location_code, location_name, region_code, region_name)],
  by = 'location_code',
  all.x = TRUE
)
baseline_deaths_dt <- merge(
  x = baseline_deaths_dt,
  y = age_groups_dt,
  by = 'age_group_code',
  all.x = TRUE
)

# Grouped by province and week
summ_by_prov_week_dt <- aggregate_summarize_excess(
  baseline_deaths_dt, aggregate = TRUE,
  group_cols = c('location_code', 'location_name', 'region_name', 'week')
)
# Grouped by week only
summ_by_week_dt <- aggregate_summarize_excess(
  baseline_deaths_dt, aggregate=TRUE, group_cols=c('week')
)[wk_starts, date := i.date, on = 'week']
# Update so that the final (incompletely-observed week) trend of expected and
#  observed deaths is correct
final_week_obs_days <- as.integer(as.Date('2020-08-31') - max(wk_starts$date)) + 1
count_cols <- c(
  'deaths', 'bl_mean', 'bl_lower', 'bl_upper', 'ex_d_mean', 'ex_d_lower',
  'ex_d_upper'
)
for(ccol in count_cols){
  suppressWarnings(
    summ_by_week_dt[week == max(week), (ccol) := get(ccol)*7/final_week_obs_days ]
  )
}

# Grouped by province, across the duration of the "focal period" (March to May)
summ_by_prov_focus_dt <- aggregate_summarize_excess(
  baseline_deaths_dt[ week %in% (start_week:focus_end_week), ],
  aggregate = TRUE, group_cols = c('location_code')
)[
  pop_dt[year==2020, .(pop=sum(pop)), by=location_code],
  pop := i.pop,
  on = 'location_code'
]

# Grouped by province and age group, across the duration of the focal period
summ_by_prov_age_focus <- aggregate_summarize_excess(
  baseline_deaths_dt[ week %in% (start_week:focus_end_week), ],
  aggregate = TRUE, group_cols = c('location_code', 'age_group_code', 'age_group_name')
)[, fill_grp := fill_grps[1]
  ][ (sig_over1==1) & (smr_mean<1.2), fill_grp:=fill_grps[2]
  ][ (sig_over1==1) & (smr_mean>=1.2) & (smr_mean<1.5), fill_grp:=fill_grps[3]
  ][ (sig_over1==1) & (smr_mean>=1.5) & (smr_mean<2.0), fill_grp:=fill_grps[4]
  ][ (sig_over1==1) & (smr_mean>=2.0) & (smr_mean<3.0), fill_grp:=fill_grps[5]
  ][ (sig_over1==1) & (smr_mean>=3.0), fill_grp:=fill_grps[6] ]

# Grouped by age group and week across the focal period
summ_by_age_week_focus <- aggregate_summarize_excess(
  baseline_deaths_dt[ week %in% (start_week:(focus_end_week + 1)), ],
  aggregate = TRUE, group_cols = c('age_group_code', 'age_group_name', 'week')
)[ order(week, age_group_code)
  ][ ex_d_mean < 0, ex_d_mean := 0
  ][, ex_d_plot := cumsum(ex_d_mean), by = week
  ][wk_starts, date := i.date, on = 'week'
  ][ order(age_group_code, week) ]


tictoc::toc()

## Map excess mortality rate by province and week ------------------------------

tictoc::tic("Making plots")

for(week_id in start_week:end_week){
  plot_title <- paste("Week of", wk_starts[week == week_id, date_plot ])
  plot_dt <- data.table::copy(summ_by_prov_week_dt[week == week_id, ])
  # Assign fill colors
  plot_dt[, fill_grp := fill_grps[1] ]
  plot_dt[(sig_over1==1) & (smr_mean<1.2), fill_grp:=fill_grps[2] ]
  plot_dt[(sig_over1==1) & (smr_mean>=1.2) & (smr_mean<1.5), fill_grp:=fill_grps[3] ]
  plot_dt[(sig_over1==1) & (smr_mean>=1.5) & (smr_mean<2.0), fill_grp:=fill_grps[4] ]
  plot_dt[(sig_over1==1) & (smr_mean>=2.0) & (smr_mean<3.0), fill_grp:=fill_grps[5] ]
  plot_dt[(sig_over1==1) & (smr_mean>=3.0), fill_grp:=fill_grps[6] ]
  # Plot it
  covidemr::map_ita_choropleth(
    province_sf = polys_sf, region_sf = regions_sf,
    in_data = plot_dt,
    map_field = 'fill_grp',
    fill_list = list(
      values = fill_colors, limits = rev(fill_grps)
    ),
    titles_list = list(title=plot_title, fill='Relative change\nfrom baseline', x='', y=''),
    fill_type = 'manual',
    save_fp = sprintf("%s/map_excess_rate_week%02d.png", viz_dir, week_id)
  )
}


## Plot excess deaths against COVID deaths at the national level ---------------

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

focus_date_lims <- as.Date(c('2020-02-26', '2020-05-28'))
focus_date_breaks <- as.Date(c(
  '2020-02-26', '2020-03-01', '2020-04-01', '2020-05-01', '2020-05-27'
))
focus_date_labs <- c('', 'Mar 01', 'Apr 01', 'May 01', 'May 27')

natl_weekly_plot <- ggplot(data=summ_by_week_dt, aes(x=date)) +
  geom_ribbon(
    aes(ymin=bl_lower, ymax=bl_upper),
    color = NA, fill = mort_plot_colors['Predicted'], alpha = .5
  ) +
  geom_line(aes(y=bl_mean, color='Predicted'), lwd = .75) +
  geom_line(aes(y=deaths, color='Observed'), lwd = .75) +
  # Add vertical dotted lines for the focus period
  geom_vline(xintercept=wk_starts[week==start_week,date], color='#888888', linetype=2) +
  geom_vline(xintercept=wk_starts[week==(focus_end_week+1),date], color='#888888', linetype=2) +
  scale_x_date(limits = range(summ_by_week_dt$date), breaks=date_breaks, labels=date_labs) +
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
  x = summ_by_week_dt[ week %in% start_week:(focus_end_week+1), ],
  y = covid_mort[, .(covid_deaths = sum(covid_deaths)), by=week],
  by = 'week'
)

y_breaks <- seq(0, 1.2E4, by=2E3)
y_labs <- c('0 (baseline)', sprintf('%d,000', setdiff(y_breaks, 0)/1E3))

excess_natl_plot <- ggplot(data=natl_covid_and_excess, aes(x=date)) +
  geom_hline(yintercept = 0, color='black', lwd=.5) +
  geom_ribbon(
    aes(ymin=ex_d_lower, ymax=ex_d_upper),
    color = NA, fill = mort_plot_colors['Excess (all-cause)'], alpha = .5
  ) +
  geom_line(aes(y=covid_deaths, color='COVID-19'), lwd = .75) +
  geom_line(aes(y=ex_d_mean, color='Excess (all-cause)'), lwd = .75) +
  scale_x_date(
    limits = focus_date_lims,
    breaks = focus_date_breaks,
    labels = focus_date_labs
  ) +
  scale_y_continuous(
    limits = range(
      c(natl_covid_and_excess$ex_d_lower, natl_covid_and_excess$ex_d_upper)
    ),
    breaks = y_breaks, labels = y_labs
  ) +
  scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors)) +
  labs(
    x='Date', y='Deaths',
    title='COVID-19 deaths and excess mortality',
    subtitle = 'Italy, March through May 2020',
    color = ''
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = c(.94, .98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill='white', color='#CCCCCC', size=.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(.55, 'inches'),
    legend.key.height = unit(.25, 'inches'),
    legend.margin = margin(2.5, 4, 4, 4)
  )
png(
  file.path(viz_dir, 'national_excess_march_to_may.png'),
  height=5, width=8, units='in', res=200
)
plot(excess_natl_plot)
dev.off()

# Plot with a weekly intercept
for(week_id in start_week:focus_end_week){
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


## Plot provinces with the most excess deaths during the focus time period -----

top_excess_provinces_dt <- data.table::copy(summ_by_prov_focus_dt)
top_excess_provinces_dt <- top_excess_provinces_dt[
    order(-ex_d_mean)
  ][ ex_d_mean < 0, ex_d_mean := 0
  ][, cumul_ex := cumsum(ex_d_mean) / sum(ex_d_mean)
  ][, rank := .I ]

# Get top 3 provinces, provinces accounting for >50% of excess deaths, and
#  provinces accounting for more than >75% of excess deaths
num_under_75pct <- nrow(top_excess_provinces_dt[cumul_ex <= .75 ])
num_under_50pct <- nrow(top_excess_provinces_dt[cumul_ex <= .50 ])

prov_grps <- c(
  'Top 3 provinces',
  paste('Top', num_under_50pct, 'provinces'),
  paste('Top', num_under_75pct, 'provinces'),
  'Other'
)
prov_colors <- c('#120d32', '#7d2482', '#fea873', '#b3b3b3')
names(prov_colors) <- prov_grps

# Assign grouping names
top_excess_provinces_dt[, top_grouping := prov_grps[4]
  ][cumul_ex <= .75, top_grouping := prov_grps[3]
  ][cumul_ex <= .5, top_grouping := prov_grps[2]
  ][rank <= 3, top_grouping := prov_grps[1] ]

top_prov_dt <- top_excess_provinces_dt[, .(location_code, ex_d_mean, rank, top_grouping, pop)]
top_prov_dt <- merge(location_table, top_prov_dt, by='location_code')[order(rank)]
top_prov_dt[, cumul_prop_pop := cumsum(pop)/sum(pop) ]
top_prov_dt[, cumul_prop_ex := cumsum(ex_d_mean)/sum(ex_d_mean) ]
fwrite(top_prov_dt, file=file.path(viz_dir, 'top_provinces_by_excess.csv'))


## Map the top provinces by excess
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = top_prov_dt,
  map_field = 'top_grouping',
  fill_list = list(
    values = prov_colors, limits = rev(prov_grps)
  ),
  titles_list = list(title = '', fill = 'Ranking by\nExcess Deaths', x='', y=''),
  fill_type = 'manual',
  labels_dt = top_prov_labs,
  save_fp = file.path(viz_dir, "province_top_burden_map.png")
)


## Plot the excess curve of the top provinces
weekly_by_province_data <- merge(
  summ_by_prov_week_dt[, .(location_code, week, ex_d_mean)],
  top_prov_dt[, .(location_code, top_grouping)],
  by = 'location_code'
)[, .(ex_d_mean = sum(ex_d_mean)), by = .(week, top_grouping) ]
weekly_by_province_data[ ex_d_mean < 0, ex_d_mean := 0 ]
weekly_by_province_data <- merge(
  x = weekly_by_province_data,
  y = data.table(top_grouping = as.character(prov_grps), grp_order = 1:4),
  by = 'top_grouping'
)[ order(week, grp_order)
  ][, ex_to_plot := cumsum(ex_d_mean), by=week
  ][ wk_starts[, .(week, date)], date := i.date, on='week'
  ][ order(grp_order) ][
  ][, top_grouping := factor(top_grouping, levels = rev(prov_grps)) ]


excess_by_prov_grp_plot <- ggplot(data = weekly_by_province_data, aes(x=date)) +
  geom_hline(yintercept = 0, color='black', lwd=.25, alpha=.5) +
  geom_ribbon(aes(fill=top_grouping, ymax=ex_to_plot), ymin=0, color=NA) +
  scale_x_date(
    limits = focus_date_lims,
    breaks = focus_date_breaks,
    labels = focus_date_labs
  ) +
  scale_y_continuous(
    limits = range(
      c(natl_covid_and_excess$ex_d_lower, natl_covid_and_excess$ex_d_upper)
    ),
    breaks = y_breaks, labels = y_labs
  ) +
  scale_fill_manual(values=prov_colors, breaks=names(prov_colors)) +
  labs(
    x='Date', y='Deaths',
    title='Excess deaths by province grouping',
    subtitle = 'Italy, March through May 2020',
    fill = ''
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = c(.94, .98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill='white', color='#CCCCCC', size=.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(.55, 'inches'),
    legend.key.height = unit(.25, 'inches'),
    legend.margin = margin(2.5, 4, 4, 4)
  )
png(
  file.path(viz_dir, 'province_top_burden_lines.png'),
  height=5, width=8, units='in', res=200
)
plot(excess_by_prov_grp_plot)
dev.off()


## Get aggregates by age -------------------------------------------------------

## Create SMR choropleth maps by age group
for(this_ag in unique(summ_by_prov_age_focus$age_group_name)){
  covidemr::map_ita_choropleth(
    province_sf = polys_sf, region_sf = regions_sf,
    in_data = summ_by_prov_age_focus[age_group_name == this_ag, ],
    map_field = 'fill_grp',
    fill_list = list(
      values = fill_colors, limits = rev(fill_grps)
    ),
    titles_list = list(title = this_ag, fill = 'Relative change\nfrom baseline'),
    fill_type = 'manual',
    save_fp = sprintf("%s/ag_excess_rate_map_%s.png", viz_dir, gsub(" ","_",this_ag))
  )
}

age_plot_colors <- RColorBrewer::brewer.pal('YlGn', n=nrow(age_groups_dt) + 1)[-1]
names(age_plot_colors) <- age_groups_dt$age_group_name

age_groups_dt$age_group_name_plot <- factor(
  age_groups_dt$age_group_name, levels = rev(age_groups_dt$age_group_name)
)
summ_by_age_week_focus[
  age_groups_dt, age_group_name_plot := i.age_group_name_plot, on='age_group_name'
]

## Aggregate by age group and week, then create line plot
excess_age_week_plot <- ggplot(data=summ_by_age_week_focus, aes(x=date)) +
  geom_hline(yintercept = 0, color='black', lwd=.25, alpha=.5) +
  geom_ribbon(
    aes(fill = age_group_name_plot, ymax = ex_d_plot), ymin = 0, color = NA
  ) +
  scale_x_date(
    limits = focus_date_lims,
    breaks = focus_date_breaks,
    labels = focus_date_labs
  ) +
  scale_y_continuous(
    limits = range(
      c(natl_covid_and_excess$ex_d_lower, natl_covid_and_excess$ex_d_upper)
    ),
    breaks = y_breaks, labels = y_labs
  ) +
  scale_fill_manual(
    values=age_plot_colors, breaks=names(age_plot_colors)
  ) +
  labs(
    x='Date', y='Deaths',
    title='Excess deaths by age group',
    subtitle = 'Italy, March through May 2020',
    fill='') +
  theme_light() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = c(.94, .98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill='white', color='#CCCCCC', size=.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(.55, 'inches'),
    legend.key.height = unit(.25, 'inches'),
    legend.margin = margin(2.5, 4, 4, 4)
  )
png(
  file.path(viz_dir, 'ag_excess_lines_plot.png'),
  height=6, width=10, units='in', res=300
)
plot(excess_age_week_plot)
dev.off()


# Add plot that includes the population distribution and 2019 death distribution
#  by age group

ag_dummy <- CJ(age_group_code = sort(unique(pop_dt$age_group_code)), dummy=1:2)

pop_distrib_dt <- pop_dt[ year==2020, .(pop=sum(pop)), by = age_group_code
  ][ order(age_group_code)
  ][, rel_pop := pop / sum(pop)
  ][, to_plot := cumsum(rel_pop)
  ][ age_groups_dt, age_group_name_plot:=i.age_group_name_plot, on='age_group_code']
fwrite(pop_distrib_dt, file=file.path(viz_dir, 'pop_dist_jan_2020.csv'))
pop_distrib_dt <- merge(pop_distrib_dt, ag_dummy, by = 'age_group_code')

d_distrib_dt <- data_full[year < 2020, .(deaths=sum(deaths)), by=age_group_code
  ][ order(age_group_code)
  ][, rel_deaths := deaths / sum(deaths)
  ][, to_plot := cumsum(rel_deaths)
  ][ age_groups_dt, age_group_name_plot:=i.age_group_name_plot, on='age_group_code']
fwrite(d_distrib_dt, file=file.path(viz_dir, 'deaths_dist_2015_to_2019.csv'))
d_distrib_dt <- merge(d_distrib_dt, ag_dummy, by = 'age_group_code')

# Helper function to make similar plots for death and population distributions
make_distrib_fig <- function(dt, title){
  out_fig <- ggplot(data=dt, aes(x = dummy)) +
    geom_ribbon(ymin=0, aes(ymax = to_plot, fill=age_group_name_plot)) +
    scale_x_continuous(limits=c(1,2), breaks=c(1,2), labels=c('May 27','')) +
    scale_y_continuous(
      limits = c(-0.23,1.15), breaks = c(0, .25, .5, .75, 1.),
      labels = c('0%', '25%', '50%', '75%', '100%')
    ) +
    scale_fill_manual(
      values=age_plot_colors, breaks=names(age_plot_colors)
    ) +
    labs(title=title, x="", y="", color="") +
    theme_light() +
    theme(
      legend.position='none',
      axis.text.x = element_text(colour='white', angle=45, hjust=1),
      axis.text.y = element_text(size=8),
      plot.title = element_text(size=8, hjust=0.5),
      plot.margin = unit(c(0.2,0.5,0.2,0), 'cm'),
      axis.ticks.x = element_blank()
    )
  return(out_fig)
}

pop_distrib_fig <- make_distrib_fig(pop_distrib_dt, '\n\nPopulation\nJan 2020')
d_distrib_fig <- make_distrib_fig(d_distrib_dt, '\n\nDeaths\n2015-19')

## Plot it
png(
  file.path(viz_dir, 'ag_excess_lines_plot_with_refs.png'),
  height=6, width=12, units='in', res=300
)
gridExtra::grid.arrange(
  excess_age_week_plot, pop_distrib_fig, d_distrib_fig,
  layout_matrix = matrix(c(rep(1,8), 2, 3), nrow=1),
  padding = unit(0, 'line')
)
dev.off()


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
    fill = 'Deaths\nper 100,000\npopulation', x='', y=''),
  labels_dt = top_reg_labs,
  save_fp = file.path(viz_dir, 'COVID_deaths_per_100k.png')
)


## Plot excess mortality lines by week for each province -----------------------

tictoc::tic("Plotting mortality lines by province")

for(loc_code in unique(location_table$location_code)){
  this_loc_name <- location_table[location_code==loc_code, location_name]
  this_reg_name <- location_table[location_code==loc_code, region_name]
  plot_title <- sprintf('%s (%s)', this_loc_name, this_reg_name)
  plot_fp <- glue::glue(
    '{viz_dir}/em_lines_prov_{sprintf("%03d",loc_code)}_',
    '{gsub(" ","_",this_loc_name)}_{gsub(" ","_",this_reg_name)}.png'
  )

  to_plot_dt <- summ_by_prov_week_dt[location_code == loc_code & week %in% 9:22,
    ][wk_starts, date := i.date, on='week']

  prov_line_plot <- ggplot(data=to_plot_dt, aes(x=date)) +
    geom_hline(yintercept = 0, color='black', lwd=.5) +
    geom_ribbon(
      aes(ymin=ex_d_lower, ymax=ex_d_upper),
      color = NA, fill = mort_plot_colors['Excess (all-cause)'], alpha = .5
    ) +
    geom_line(aes(y=ex_d_mean, color='Excess (all-cause)'), lwd = .75) +
    scale_x_date(
      limits = focus_date_lims,
      breaks = focus_date_breaks,
      labels = focus_date_labs
    ) +
    scale_color_manual(values=mort_plot_colors, breaks=names(mort_plot_colors), guide=FALSE) +
    labs(
      x='Date', y='Deaths',
      title=plot_title,
      color='') +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  png(plot_fp, height=4, width=6, units='in', res=300)
  plot(prov_line_plot)
  dev.off()
}

tictoc::toc() # END plotting mortality lines by province


## Map the week of peak excess across northern provinces -----------------------

peak_excess_dt <- summ_by_prov_week_dt[ week %in% 9:21,
  ][, week_rank := frank(-ex_d_mean), by = location_code
  ][ week_rank == 1,
  ][ wk_starts, date := i.date, on='week'
  ][location_table, region_code := i.region_code, on = 'location_code']

wk_labs <- c(
  '11 Mar or prior', '18 Mar', '25 Mar', '01 Apr', '08 Apr', '15 Apr or later')
wk_colors <- RColorBrewer::brewer.pal('PuOr', n=11)[c(2,3,5,7,9,10)]
wk_breaks <- c(.99, 12:16 - .5, 52.01)
names(wk_colors) <- wk_labs

peak_excess_dt[, week_label := wk_labs[cut(week, wk_breaks, labels=FALSE)]]

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
    fill = 'Week beginning', x='', y=''),
  fill_type = 'manual',
  save_fp = file.path(viz_dir, 'peak_excess_week_map.png')
)


## Map the maximum excess across northern provinces ----------------------------


peak_excess_dt$fill_grp <- fill_grps_highlight[
  cut(peak_excess_dt$smr_mean, fill_brks_highlight, labels=FALSE)
]

covidemr::map_ita_choropleth(
  province_sf = polys_sf[polys_sf$region_code %in% 1:9, ],
  region_sf = regions_sf,
  in_data = peak_excess_dt[region_code %in% 1:9, ],
  map_field = 'fill_grp',
  fill_list = list(
    values = fill_colors_highlight, limits = rev(fill_grps_highlight)
  ),
  titles_list = list(
    title = 'Peak excess mortality over baseline',
    subtitle = 'March through May 2020',
    fill = 'Excess Mortality', x='', y=''),
  fill_type = 'manual',
  save_fp = file.path(viz_dir, 'peak_excess_amount_map.png')
)

tictoc::toc() # END all maps and plots


## Aggregate across various dimensions and save to CSV files -------------------

tictoc::tic("Creating excess mortality aggregate files")

# Helper function to aggregate data and then save to file
agg_then_save <- function(dt, group_cols, out_fp, printout = TRUE){
  dt_agg <- aggregate_summarize_excess(
    dt, aggregate = TRUE, group_cols = group_cols
  )
  fwrite(dt_agg, file = out_fp)
  if(printout){
    print(knitr::kable(dt_agg))
    message("")
  } else {
    invisible()
  }
}

# Run aggregations for focus weeks (all provinces)
bl_focus <- baseline_deaths_dt[week %in% (start_week:focus_end_week), ]
vizfp <- function(basefile) file.path(viz_dir, basefile)

agg_then_save(bl_focus, group_cols=NULL, out_fp=vizfp('natl_excess_all.csv'))
agg_then_save(bl_focus, group_cols=c('age_group_name','sex'), out_fp=vizfp('natl_excess_by_age_sex.csv'))
agg_then_save(bl_focus, group_cols='sex', out_fp=vizfp('natl_excess_by_sex.csv'))
agg_then_save(bl_focus, group_cols='age_group_name', out_fp=vizfp('natl_excess_by_age.csv'))
agg_then_save(
  bl_focus, group_cols=c('region_code','region_name','age_group_name','sex'),
  out_fp=vizfp('regional_excess_by_age_sex.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('region_code','region_name','sex'),
  out_fp=vizfp('regional_excess_by_sex.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('region_code','region_name','age_group_name'),
  out_fp=vizfp('regional_excess_by_age.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('region_code','region_name'),
  out_fp=vizfp('regional_excess.csv'), printout = TRUE
)
agg_then_save(
  bl_focus, group_cols=c('location_code','location_name','age_group_name','sex'),
  out_fp=vizfp('provincial_excess_by_age_sex.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('location_code','location_name','sex'),
  out_fp=vizfp('provincial_excess_by_sex.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('location_code','location_name','age_group_name'),
  out_fp=vizfp('provincial_excess_by_age.csv'), printout = FALSE
)
agg_then_save(
  bl_focus, group_cols=c('location_code','location_name'),
  out_fp=vizfp('provincial_excess.csv'), printout = TRUE
)

# Get excess and SMRs only for the north of Italy (regions 1-9)
northern_subset <- bl_focus[region_code %in% 1:9, ]
agg_then_save(northern_subset, group_cols=NULL, out_fp=vizfp('northern_excess_overall.csv'))
agg_then_save(northern_subset, group_cols=c('age_group_name','sex'), out_fp=vizfp('northern_excess_by_age_sex.csv'))
agg_then_save(northern_subset, group_cols=c('sex'), out_fp=vizfp('northern_excess_by_sex.csv'))
agg_then_save(northern_subset, group_cols=c('age_group_name'), out_fp=vizfp('northern_excess_by_age.csv'))

tictoc::toc() # END creating excess mortality aggregate files
