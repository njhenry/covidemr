## -----------------------------------------------------------------------------
##
## 02a: Data diagnostics for input deaths, population, and covariates
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

## Load required packages and inputs

# library(covidemr)
library(data.table)
library(knitr)
library(sf)
# DEVELOPMENT: rebuild library
dev_fp <- '~/repos/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## TODO: Load using command line argument
prepped_data_version <- '20201026'

# Helper functions to create a filepath for a particular prepped data object or
#  diagnostic filepath
prep_dir <- file.path(config$paths$prepped_data, prepped_data_version)
diag_dir <- file.path(prep_dir, 'diagonstics')
dir.create(diag_dir, showWarnings=FALSE)
get_prep_fp <- function(ff) file.path(prep_dir, config$prepped_data_files[[ff]])
get_diag_fp <- function(ff) file.path(diag_dif, ff)

## Load input data
location_table <- data.table::fread(get_prep_fp('location_table'))
age_groups <- create_age_groups(config$age_cutoffs)
polys_sf <- readRDS(get_prep_fp('shapefile_sf'))
covars_prepped <- readRDS(get_prep_fp('covars_list'))
covar_names <- names(covars_prepped)
data_rescaled <- data.table::fread(get_prep_fp('full_data_rescaled'))


## Run variance inflation factor analysis on data ------------------------------

message("\n\n\n*** Running VIF analysis ***")
vif_results <- get_covar_vif(
  in_data = data_rescaled[ in_baseline == 1, ],
  covar_names = covar_names
)
message("VIF results:", appendLF=FALSE)
print(knitr::kable(vif_results))


## Univariate covariate visualizations -----------------------------------------




## Visualize covariates against each other -------------------------------------



## Visualize covariates against mortality --------------------------------------