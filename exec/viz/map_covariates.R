## -----------------------------------------------------------------------------
##
## VISUALIZATION: Make maps of covariates
##
## -----------------------------------------------------------------------------

## Load required packages and inputs

# library(covidemr)
library(data.table)
library(RColorBrewer)
library(sp)
library(sf)
# DEVELOPMENT: rebuild library
dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Run settings
prepped_data_version <- '20201104'

# Helper function to create a filepath for a particular prepped data object
prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
viz_dir <- file.path(prep_dir, 'viz')
dir.create(viz_dir, showWarnings = FALSE)
get_prep_fp <- function(ff) file.path(prep_dir, config$prepped_data_files[[ff]])


## Load and pre-format data ----------------------------------------------------

# Load data
location_table <- data.table::fread(get_prep_fp('location_table'), na.strings="")
polys_sf <- readRDS(get_prep_fp('shapefile_sf'))
covs_list <- readRDS(get_prep_fp('covars_list'))

# Merge extra location identifiers onto the spatial polygons
merge_on_cols <- c('location_code', 'location_name', 'region_code', 'region_name')
polys_sf <- merge(polys_sf, location_table[, ..merge_on_cols], by='location_code')

# Create a regional sf object by dissolving provinces
regions_sf <- aggregate(polys_sf, by=list(polys_sf$region_code), FUN=first)


## Make plots ------------------------------------------------------------------

## TFR
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$tfr$prepped_covar[year==2020],
  map_field = 'tfr',
  fill_lims = c(.90, 1.6),
  fill_list = list(
    limits = c(.9, 1.6),
    breaks = c(1.0, 1.1, 1.3, 1.45),
    labels = c('< 1.0', '1.0 - 1.2', '1.2 - 1.4', '> 1.4'),
    colors = c("#EDF8E9","#EDF8E9","#BAE4B3","#74C476","#31A354","#006D2C","#006D2C")
  ),
  titles_list = list(fill = 'Total\nfertility\nrate'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_tfr.png')
)

## Unemployment among men
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$unemp$prepped_covar[sex=='male' & year==2020],
  map_field = 'unemp',
  fill_lims = c(0, 22),
  fill_list = list(
    breaks = c(0, 5, 10, 15, 20),
    labels = c('0%', '5%', '10%', '15%', '> 20%'),
    colors = c("#FEEDDE","#FDD0A2","#FDAE6B","#FD8D3C","#E6550D","#A63603","#A63603")
  ),
  titles_list = list(fill = 'Unemployment\nrate (men)'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_unemployment_men.png')
)

## Unemployment among women
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$unemp$prepped_covar[sex=='female' & year==2020],
  map_field = 'unemp',
  fill_lims = c(0, 22),
  fill_list = list(
    breaks = c(0, 5, 10, 15, 20),
    labels = c('0%', '5%', '10%', '15%', '> 20%'),
    colors = c("#FEEDDE","#FDD0A2","#FDAE6B","#FD8D3C","#E6550D","#A63603","#A63603")
  ),
  titles_list = list(fill = 'Unemployment\nrate (women)'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_unemployment_women.png')
)

## Social services
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$socserv$prepped_covar[year==2020, ],
  map_field = 'socserv',
  fill_lims = c(0, .35),
  fill_list = list(
    breaks = c(0, .1, .2, .3),
    labels = c('0%', '10%', '20%', '> 30%'),
    colors = c("#F0F9E8","#CCEBC5","#A8DDB5","#7BCCC4","#43A2CA","#0868AC")
  ),
  titles_list = list(fill = 'At-home social\nservices: proportion\nreached'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_social_services.png')
)

## Proportion making under 10k Euros
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$tax_brackets$prepped_covar[year==2020, ],
  map_field = 'tax_brackets',
  fill_lims = c(.2, .5),
  fill_list = list(
    breaks = c(.2, .3, .4, .5),
    labels = c('< 20%', '30%', '40%', '50%'),
    colors = c("#F2F0F7","#DADAEB","#BCBDDC","#9E9AC8","#756BB1","#54278F")
  ),
  titles_list = list(fill = 'Households\nwith taxable\nincome under\n10k Euros/year'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_tax_brackets.png')
)


## Travel time to healthcare
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$hc_access$prepped_covar,
  map_field = 'hc_access',
  fill_lims = c(1, 10),
  fill_list = list(
    breaks = c(1, 2, 5, 10),
    labels = c('1 min', '2 min', '5 min', '10 min'),
    colors = RColorBrewer::brewer.pal("Greys", n = 6)
  ),
  titles_list = list(fill = 'Transit time\nto nearest\nhealth facility'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_hc_travel_time.png')
)

## Elevation
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$elevation$prepped_covar,
  map_field = 'elevation',
  fill_lims = c(0, 500),
  fill_list = list(
    breaks = seq(0, 500, by=100),
    labels = c('0m', '100m', '200m', '300m', '400m', '> 500m'),
    colors = RColorBrewer::brewer.pal("YlOrBr", n = 6)
  ),
  titles_list = list(fill = 'Mean\nresidential\nelevation'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_elevation.png')
)

## Temperature on a particular date
covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$temperature$prepped_covar[year == 2020 & week == 1, ],
  map_field = 'temperature',
  fill_lims = c(0, 30),
  fill_list = list(
    breaks = seq(0, 30, by=10),
    labels = c('0 C', '10 C', '20 C', '30 C'),
    colors = rev(RColorBrewer::brewer.pal("RdYlBu", n = 9))
  ),
  titles_list = list(fill = 'Temperature,\nweek of\n1 Jan 2020'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_temperature_jan1.png')
)

covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf,
  in_data = covs_list$temperature$prepped_covar[year == 2020 & week == 26, ],
  map_field = 'temperature',
  fill_lims = c(0, 30),
  fill_list = list(
    breaks = seq(0, 30, by=10),
    labels = c('0 C', '10 C', '20 C', '30 C'),
    colors = rev(RColorBrewer::brewer.pal("RdYlBu", n = 9))
  ),
  titles_list = list(fill = 'Temperature,\nweek of\n30 June 2020'),
  fill_type = 'gradientn',
  save_fp = glue::glue('{viz_dir}/covar_temperature_june30.png')
)
