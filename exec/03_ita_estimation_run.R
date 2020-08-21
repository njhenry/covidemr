## -----------------------------------------------------------------------------
## 
## 03: Space-time excess mortality modeling code, using prepared data as input
## 
## For more details, see README at https://github.com/njhenry/covidemr/
## 
## -----------------------------------------------------------------------------


library(data.table)
library(sp)
library(sf)

dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
devtools::load_all(dev_fp)
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))


## Load and prepare data -------------------------------------------------------

