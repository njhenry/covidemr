## -----------------------------------------------------------------------------
## 
## 00: Create standard list of Provinces
## 
## For more details, see README at https://github.com/njhenry/covidemr/
## 
## -----------------------------------------------------------------------------

## Load packages

library(data.table)

# DEVELOPMENT: load package functions and config
dev_fp <- '/ihme/code/covid-19/user/nathenry/covidemr/'
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## DEVELOPMENT: set globals
## TODO: Convert to command line argument
prepped_data_version <- '20200823'

## Create standard location list

prov_table <- foreign::read.dbf(gsub('.shp$', '.dbf', config$paths$shp_generalized))
prov_table <- prov_table[, c('COD_REG','COD_PROV','DEN_UTS','SIGLA')]
names(prov_table) <- c('region_code','location_code','location_name','abbrev')
prov_table <- as.data.table(prov_table)

## Add on ITTER107 code
pop_raw <- fread(
  config$paths$raw_pop, encoding=config$encoding, na.strings=c("","NA","n.d.")
)
setnames(pop_raw, c('Territory','ITTER107'), c('location_name','icode'))
pop_sub <- unique(pop_raw[, .(location_name, icode)])
pop_sub[, location_name := iconv(
  location_name, from=tolower(gsub('-','',config$encoding)), to='ASCII//TRANSLIT'
)]
# Add necessary renames
renames <- list(
  c("Valle d'Aosta / VallA(C)e d'Aoste", "Aosta"),
  c("Bolzano / Bozen", "Bolzano"),
  c("ForlA!-Cesena", "Forli'-Cesena"),
  c("Massa-Carrara", "Massa Carrara")
)
for(rename in renames) pop_sub[ location_name == rename[1], location_name := rename[2]]
old_provs <- c('ITG2C','ITG2B','ITG2A','ITG29')
pop_sub <- pop_sub[!(icode %in% old_provs), ]

# Merge tables to add code
prov_table_full <- merge(prov_table, pop_sub, by='location_name', all=TRUE)
if(nrow(prov_table_full) != 107)) stop("ISSUE: Wrong number of provinces!")
prov_table_full <- prov_table_full[order(location_code)]

## Create output directory and save

out_dir <- file.path(config$paths$prepped_data, prepped_data_version)
dir.create(out_dir, showWarnings=FALSE)

out_fp <- file.path(out_dir, config$prepped_data_files$location_table)
write.csv(prov_table_full, file=out_fp, row.names=FALSE)
