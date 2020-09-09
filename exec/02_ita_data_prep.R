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
dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## TODO: Load using command line argument
prepped_data_version <- '20200909'

# Load prepared location table
prepped_data_dir <- file.path(config$paths$prepped_data, prepped_data_version)
location_table <- data.table::fread(file.path(
  prepped_data_dir, config$prepped_data_files$location_table
))

# Get age groups
age_groups <- create_age_groups(config$age_cutoffs)

# Define a convenience function for loading CSVs
ita_fread <- function(fp){
  fread(fp, encoding=config$encoding, na.strings=c("","NA","n.d."))
}


## Load and format deaths data -------------------------------------------------

deaths_raw <- ita_fread(config$paths$raw_deaths)

deaths_prepped <- ita_prepare_deaths(
  deaths_raw,
  age_cutoffs = config$age_cutoffs,
  model_years = config$model_years,
  first_covid_death_date = as.Date(config$first_death_date)
)
# Save prepared deaths data to file
write.csv(
  deaths_prepped,
  file = file.path(prepped_data_dir, config$prepped_data_files$deaths),
  row.names = FALSE
)


## Load and format population data ---------------------------------------------

pop_raw <- ita_fread(config$path$raw_pop)

pop_prepped <- ita_prepare_pop(pop_raw, age_cutoffs=config$age_cutoffs)

write.csv(
  pop_prepped,
  file = file.path(prepped_data_dir, config$prepped_data_files$population),
  row.names = FALSE
)


## Load and format covariates --------------------------------------------------

covars_raw <- lapply(config$paths$raw_covars, ita_fread)
covar_names <- names(config$paths$raw_covars)
names(covars_raw) <- covar_names

# Prepare
covars_prepped <- lapply(covar_names, function(covar_name){
  ita_prepare_covariate(
    copy(covars_raw[[covar_name]]),
    covar_name,
    model_years = config$model_years,
    location_table = location_table
  )
})
names(covars_prepped) <- covar_names



## Merge deaths, population, and covariates to get prepped input dataset -------

in_data <- merge(
  x = deaths_prepped,
  y = pop_prepped, 
  by = c('location_code', 'year', 'sex', 'age_group_code'),
  all = TRUE
)
for(covar_name in covar_names){
  in_data <- merge(
    x = in_data,
    y = covars_prepped[[covar_name]]$prepped_covar,
    by = covars_prepped[[covar_name]]$covar_indices,
    all = TRUE
  )
}

# Quick validation
if(any(is.na(in_data$deaths))) stop("Missing deaths in final dataset")
if(any(is.na(in_data$pop))) stop("Missing population in final dataset")
for(covar_name in covar_names){
  if(any(is.na(in_data[[covar_name]]))) stop(sprintf(
    "Missing values for covar %s in final dataset", covar_name
  ))
}

# Save to file
write.csv(
  in_data,
  file = file.path(prepped_data_dir, config$prepped_data_files$full_data),
  row.names = FALSE
)


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

# Save out
write.csv(
  template_dt,
  file = file.path(prepped_data_dir, config$prepped_data_files$template),
  row.names = FALSE
)


## Load shapefile, create adjacency matrix, and cache both ---------------------

shp_sf <- sf::st_read(config$paths$shp_generalized)
# Order by location code
shp_sf$location_code <- sapply(
  shp_sf$SIGLA, 
  function(SIG) location_table[abbrev==SIG, location_code]
)
shp_sf <- shp_sf[order(shp_sf$location_code), ]

# Create SP version
shp_sp <- as(shp_sf, "Spatial")
# Reset polygon IDs to reduce confusion
for(ii in 1:nrow(shp_sp@data)){
  shp_sp@polygons[[ii]]@ID <- as.character(shp_sp@data[[ii, 'location_code']])
}

# Create adjacency matrix
adjmat <- build_adjacency_matrix(shp_sp)

# Add manual adjacencies for islands:
#   - Messina, Sicilia (83) <-> Reggio Calabria, mainland (80) - proximity
#   - Cagliari, Sardinia (92) <-> Roma, mainland (58) - most common flights/ferries
add_links <- list(c(83, 80), c(92, 58))
for(add_link in add_links){
  adjmat[add_link[1], add_link[2]] <- 1
  adjmat[add_link[2], add_link[1]] <- 1
}

# Save to file
saveRDS(shp_sf, file=file.path(prepped_data_dir, config$prepped_data_files$shapefile_sf))
saveRDS(shp_sp, file=file.path(prepped_data_dir, config$prepped_data_files$shapefile_sp))
saveRDS(adjmat, file=file.path(prepped_data_dir, config$prepped_data_files$adjacency_matrix))

message("** Data prep complete. **")
