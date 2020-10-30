## -----------------------------------------------------------------------------
##
## 02: Data Preparation Code for Italy ISTAT deaths, population, and covariates
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

## Load required packages and inputs

# library(covidemr)
library(data.table)
library(sp)
library(sf)
# DEVELOPMENT: rebuild library
dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## TODO: Load using command line argument
parser <- argparse::ArgumentParser(
  description='COVID Excess Mortality: Pre-modeling data prep script',
  allow_abbrev=FALSE
)
parser$add_argument('--data-version', type='character', help='Prepped data version')
args <- parser$parse_args(commandArgs(TRUE))
prepped_data_version <- args$data_version

# Helper function to create a filepath for a particular prepped data object
prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
get_prep_fp <- function(ff) file.path(prep_dir, config$prepped_data_files[[ff]])

# Load prepared location table
location_table <- data.table::fread(get_prep_fp('location_table'))
# Get age groups
age_groups <- create_age_groups(config$age_cutoffs)


## Load and format deaths data -------------------------------------------------

deaths_raw <- covidemr::ita_read_istat_csv(config$paths$raw_deaths)
deaths_prepped <- covidemr::ita_prepare_deaths_all_years(
  deaths_raw = deaths_raw,
  age_cutoffs = config$age_cutoffs,
  model_years = config$model_years,
  first_covid_death_date = as.Date(config$first_death_date)
)

# Save prepared deaths data to file
data.table::fwrite(deaths_prepped, file = get_prep_fp('deaths'))
rm(deaths_raw)


## Load and format population data ---------------------------------------------

pop_raw <- covidemr::ita_read_istat_csv(config$path$raw_pop)
pop_prepped <- covidemr::ita_prepare_pop(pop_raw, age_cutoffs=config$age_cutoffs)

# Save to file
data.table::fwrite(pop_prepped, file = get_prep_fp('population'))


## Load spatial data: polygons, adjacency matrix, population raster ------------

# Spatial polygons
polys_list <- covidemr::load_format_spatial_polys(
  shp_in_fp = config$paths$shp_generalized,
  abbrev_field = 'SIGLA',
  location_table = location_table
)
# Adjacency matrix
# Add manual adjacencies for islands:
#   - Messina, Sicilia (83) <-> Reggio Calabria, mainland (80) - proximity
#   - Cagliari, Sardinia (92) <-> Roma, mainland (58) - most common flights/ferries
adjmat <- covidemr::build_adjacency_matrix(
  poly_sp = polys_list$shp_sp,
  allow_zero_neighbors = FALSE,
  manually_add_links = list(c(83, 80), c(92, 58))
)
# Population raster
pop_raster <- covidemr::load_format_pop_raster(
  pop_fp_format = config$paths$pop_raster_layers,
  model_years = config$model_years,
  country_polys = polys_list$shp_sf,
  projection = config$projection
)

# Save to file
saveRDS(polys_list$shp_sf, file = get_prep_fp('shapefile_sf'))
saveRDS(polys_list$shp_sp, file = get_prep_fp('shapefile_sp'))
saveRDS(adjmat, file = get_prep_fp('adjacency_matrix'))
raster::writeRaster(pop_raster, file = get_prep_fp('pop_raster'), overwrite=TRUE)


## Load and format covariates --------------------------------------------------

covar_fps <- config$paths$raw_covars
covar_fps$pop_density <- get_prep_fp('population')
covar_names <- names(covar_fps)

# Prepare each covariate separately using a standard function interface
covars_prepped <- lapply(covar_names, function(covar_name){
  covidemr::ita_prepare_covariate(
    covar_name = covar_name,
    covar_fp = covar_fps[[covar_name]],
    model_years = config$model_years,
    location_table = data.table::copy(location_table),
    pop_raster = pop_raster,
    polys_sf = polys_list$shp_sf,
    projection = config$projection
  )
})
names(covars_prepped) <- covar_names

# Save to file
saveRDS(covars_prepped, file = get_prep_fp('covars_list'))

## Create template dataset containing IDs for all categories to be modeled -----

template_dt <- create_template_dt(
  model_years = config$model_years,
  model_week_range = config$model_week_range,
  location_table = location_table,
  age_groups = age_groups
)

# Merge on covariates
template_dt[, intercept := 1 ]
for(covar_name in covar_names){
  template_dt <- merge(
    x = template_dt,
    y = covars_prepped[[covar_name]]$prepped_covar,
    by = covars_prepped[[covar_name]]$covar_indices,
    all.x = TRUE
  )
  if(any(is.na(template_dt[(year < 2020) | (week <= 26)][[covar_name]]))){
    stop("Missing values for covar ",covar_name," in template dataset")
  }
}

# Save out
fwrite(template_dt, file = get_prep_fp('template'))


## Merge deaths and population onto the template to get full input dataset -----

in_data <- merge(
  x = template_dt, y = pop_prepped,
  by = c('location_code', 'year', 'sex', 'age_group_code'),
  all.x = TRUE
)
if(any(is.na(in_data$pop))) stop("Some populations missing from input data")

in_data <- merge(
  x = in_data, y = deaths_prepped,
  by = c('location_code', 'year', 'sex', 'week', 'age_group_code'),
  all.x = TRUE
)
if(nrow(in_data[is.na(deaths)])/nrow(in_data) > .2){
  stop("Deaths missing from more than 20% of observations - check data!")
}
in_data[is.na(deaths), `:=` (deaths=0, observed_days=7)]

# Save to file
fwrite(in_data, file = get_prep_fp('full_data'))


## Rescale all covariates to N(0, 1) from the training dataset -----------------

rescale_list <- rescale_prepped_covariates(
  input_data = in_data,
  covar_names = covar_names,
  subset_field = 'in_baseline',
  subset_field_values = 1
)

fwrite(rescale_list$data_rescaled, get_prep_fp('full_data_rescaled'))
fwrite(rescale_list$covariate_scaling_factors, get_prep_fp('covar_scaling_factors'))

message("** Data prep complete. **")
