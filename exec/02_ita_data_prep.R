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
# DEVELOPMENT: rebuild library
dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

# Load prepared location table
prepped_data_dir <- file.path(config$paths$prepped_data, config$prepped_data_version)
location_table <- data.table::fread(file.path(
  prepped_data_dir, config$prepped_data_files$location_table
))


## Load and format deaths data -------------------------------------------------

deaths_raw <- data.table::fread(
  config$paths$raw_deaths, encoding=config$encoding, na.strings=c("","NA","n.d.")
)
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

pop_raw <- fread(config$path$raw_pop, encoding=config$encoding, na.strings=c('','NA','n.d.'))

pop_prepped <- ita_prepare_deaths(pop_raw, age_cutoffs=config$age_cutoffs)

write.csv(
  pop_prepped,
  file = file.path(prepped_data_dir, config$prepped_data_files$population),
  row.names = FALSE
)


## Load and format covariates --------------------------------------------------

