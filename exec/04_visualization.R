## -----------------------------------------------------------------------------
##
## 04: Data and model visualization script
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(glue)
library(sf)
library(RColorBrewer)

dev_fp <- '~/Documents/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Settings
## TODO: Convert to command line
prepped_data_version <- '20201019'
model_run_version <- '20201019'
holdout <- 0


## Load all data ---------------------------------------------------------------

# Load excess mortality and baseline mortality
model_results_dir <- config$paths$model_results
pred_list <- excess_list <- list(male=NULL, female=NULL)
for(run_sex in names(pred_list)){
  pred_list[[run_sex]] <- fread(glue(config$results_files$pred_summary))
  excess_list[[run_sex]] <- fread(glue(config$results_files$excess_summary))
}
preds_dt <- rbindlist(pred_list)[location_code < 109, ]
excess_dt <- rbindlist(excess_list)[location_code < 109, ]
# Some basic renaming
summ_names <- c('mean','median','lower','upper')
setnames(excess_dt, summ_names, paste0('ex_',summ_names))

# Load provinces spatial object
shp_sf <- readRDS(file.path(
  config$paths$prepped_data, prepped_data_version,
  config$prepped_data_files$shapefile_sf
))


## Make maps showing greatest weekly excess across all weeks -------------------

viz_dir <- file.path(config$paths$prepped_data, prepped_data_version, 'viz')
dir.create(viz_dir, showWarnings = FALSE)

ita_viz <- function(
  in_data, in_sf, merge_col, viz_col, title, color_breaks, color_labels,
  color_title, colors
){
  # Subset data down to breaks limits
  to_merge <- copy(in_data[, c(merge_col, viz_col), with=F])
  setnames(to_merge, viz_col, 'to_plot')
  to_merge[ to_plot >= max(color_breaks) - 1E-8, to_plot := max(color_breaks) - 1E-8 ]
  to_merge[ to_plot <= min(color_breaks) + 1E-8, to_plot := min(color_breaks) + 1E-8 ]
  # Merge on data
  viz_sf <- merge(
    x = in_sf[, c(merge_col, 'geometry')],
    y = to_merge,
    by = merge_col,
    all.x = TRUE
  )
  # Visualize
  fig <- ggplot(data=viz_sf, aes(fill=to_plot)) +
    geom_sf(lwd=.5, color='#444444') +
    scale_fill_gradientn(
      colors=colors, breaks=color_breaks, labels=color_labels, na.value='#BBBBBB'
    ) +
    labs(title=title, x=NULL, y=NULL, fill=color_title) +
    theme_minimal()
  return(fig)
}

pdf(file.path(viz_dir, 'max_excess_proportion.pdf'), height=11, width=8.5)
# As a PROPORTION of baseline mortality
for(this_sex in unique(excess_dt$sex)){
  for(this_age_group in unique(excess_dt$age_group_name)){
    dt_sub <- copy(excess_dt[sex==this_sex & age_group_name==this_age_group, ])
    dt_sub <- dt_sub[, .(max_prop = max(prop_mean)), by=location_code ]
    color_breaks <- seq(1, 4, by=.5)
    color_labels <- paste0(color_breaks * 100, '%')
    print(ita_viz(
      in_data = dt_sub, in_sf = shp_sf, merge_col = 'location_code',
      viz_col = 'max_prop',
      title = glue('Highest proportion above baseline: {this_sex}, ages {this_age_group}'),
      color_breaks = color_breaks, color_labels = color_labels,
      color_title = 'Max ratio\nto baseline',
      colors = rev(RColorBrewer::brewer.pal('Spectral', n=7))
    ))
  }
}
dev.off()


pdf(file.path(viz_dir, 'max_excess_count.pdf'), height=11, width=8.5)
# Total COUNT of excess mortality
for(this_sex in unique(excess_dt$sex)){
  for(this_age_group in unique(excess_dt$age_group_name)){
    dt_sub <- copy(excess_dt[sex==this_sex & age_group_name==this_age_group, ])
    dt_sub <- dt_sub[, .(excess = sum(ex_mean), pop=mean(pop)), by=location_code]
    dt_sub[, excess_per_100k := excess / pop * 1E5 ]
    color_labels <- color_breaks <- seq(0, 3000, by=1000)
    print(ita_viz(
      in_data = dt_sub, in_sf = shp_sf, merge_col = 'location_code',
      viz_col = 'excess_per_100k',
      title = glue('Excess mortality per 100k population: {this_sex}, ages {this_age_group}'),
      color_breaks = color_breaks, color_labels = color_labels,
      color_title = 'Excess mortality',
      colors = rev(RColorBrewer::brewer.pal('Spectral', n=7))
    ))
  }
}
dev.off()


pdf(file.path(viz_dir, 'excess_proportion_lineplots.pdf'), height=16, width=22)
for(this_sex in unique(excess_dt$sex)){
  for(this_age_group in unique(excess_dt$age_group_name)){
    dt_sub <- copy(excess_dt[sex==this_sex & age_group_name==this_age_group, ])
    dt_sub[ prop_lower > 5, prop_lower := 4.9997 ]
    dt_sub[ prop_median > 5, prop_median := 4.9998 ]
    dt_sub[ prop_upper > 5, prop_upper := 4.9999 ]
    print(
      ggplot(data=dt_sub, aes(x=week)) +
        facet_wrap('location_name', ncol=13) +
        geom_hline(yintercept = 1, linetype = 2, lwd=.5) +
        geom_ribbon(aes(ymin=prop_lower, ymax=prop_upper), color=NA, fill='#AA2288', alpha=.5) +
        geom_line(aes(y=prop_median)) +
        lims(y=c(0, 5)) +
        labs(
          y='Mortality\nproportional to\nbaseline', x='Week of 2020',
          title=glue('Excess mortality line plots by province: {this_sex}, {this_age_group}')
        ) +
        theme_minimal()
    )
  }
}
dev.off()
