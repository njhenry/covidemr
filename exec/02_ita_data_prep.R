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
prepped_data_version <- '20201026'

# Helper function to create a filepath for a particular prepped data object
get_prep_fp <- function(file_type){
  return(file.path(
    config$paths$prepped_data,
    prepped_data_version,
    config$prepped_data_files[[file_type]]
  ))
}

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

data.table::fwrite(pop_prepped, file = get_prep_fp('population'))


## Load spatial data: polygons, adjacency matrix, population raster ------------

polys_list <- covidemr::load_format_spatial_polys(
  shp_in_fp = config$paths$shp_generalized,
  abbrev_field = 'SIGLA',
  location_table = location_table
)
# Create adjacency matrix
# Add manual adjacencies for islands:
#   - Messina, Sicilia (83) <-> Reggio Calabria, mainland (80) - proximity
#   - Cagliari, Sardinia (92) <-> Roma, mainland (58) - most common flights/ferries
adjmat <- covidemr::build_adjacency_matrix(
  poly_sp = polys_list$shp_sp,
  allow_zero_neighbors = FALSE,
  manually_add_links = list(c(83, 80), c(92, 58))
)

# Load and format population raster
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


## Create template dataset containing IDs for all categories to be modeled -----

sex_index <- data.table(sex = c('male','female'), c_j = 1)[, idx_sex := .I - 1 ]

years_index <- data.table(
  year = sort(config$model_years),
  c_j = 1
)[, idx_year := .I - 1]

week_dt <- data.table(
  week = min(config$model_week_range):max(config$model_week_range),
  c_j = 1
)[, idx_week := .I - 1]

# Full template dataset contains sex, location, year, week, age
# idx_<index> is a zero-indexed code for each ID column
template_dt <- location_table[order(location_code)
  ][, `:=` (c_j = 1, idx_loc = .I - 1)
  ][ sex_index, on='c_j', allow.cartesian=TRUE
  ][ age_groups[, `:=` (c_j=1, idx_age=.I-1)], on='c_j', allow.cartesian=TRUE
  ][ years_index, on='c_j', allow.cartesian=TRUE
  ][ week_dt, on='c_j', allow.cartesian=TRUE
  ][, c_j := NULL ]

# Merge on covariates
template_dt[, intercept := 1 ]
for(covar_name in covar_names){
  template_dt <- merge(
    x = template_dt,
    y = covars_prepped[[covar_name]]$prepped_covar,
    by = covars_prepped[[covar_name]]$covar_indices,
    all.x = TRUE
  )
}

# Quick validation
for(covar_name in covar_names){
  if(any(is.na(template_dt[[covar_name]]))) stop(sprintf(
    "Missing values for covar %s in template dataset", covar_name
  ))
}

# Save out
write.csv(template_dt, file = get_prep_fp('template'), row.names = FALSE)


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
write.csv(in_data, file = get_prep_fp('full_data'), row.names = FALSE)

message("** Data prep complete. **")
