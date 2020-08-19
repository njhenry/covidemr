## -----------------------------------------------------------------------------
## 
## 00: Create standard list of Provinces
## 
## For more details, see README at https://github.com/njhenry/covidemr/
## 
## -----------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------

library(data.table)

# DEVELOPMENT: load package functions and config
dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))


## Create standard location list -----------------------------------------------

