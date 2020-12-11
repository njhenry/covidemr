## -----------------------------------------------------------------------------
##
## 00: Create standard list of Provinces
##
## For more details, see README at https://github.com/njhenry/covidemr/
##
## -----------------------------------------------------------------------------

# Load packages
library(data.table)

# DEVELOPMENT: load package functions and config
dev_fp <- '~/repos/covidemr/'
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))

## Load input data version from command line
ap <- argparse::ArgumentParser(
  description='COVID Excess Mortality: Create input locations list',
  allow_abbrev=FALSE
)
ap$add_argument('--data-version', type='character', help='Prepped data version date')
args <- ap$parse_args(commandArgs(TRUE))
prepped_data_version <- args$data_version


## Create standard location list -----------------------------------------------

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
  c("Valle d'Aosta / Vall~A(c)e d'Aoste", "Aosta"),
  c("Bolzano / Bozen", "Bolzano"),
  c("ForlA!-Cesena", "Forli'-Cesena"),
  c("Forl~Anot-Cesena", "Forli'-Cesena"),
  c("Massa-Carrara", "Massa Carrara")
)
for(rename in renames) pop_sub[ location_name == rename[1], location_name := rename[2]]
old_provs <- c('ITG2C','ITG2B','ITG2A','ITG29')
pop_sub <- pop_sub[!(icode %in% old_provs), ]

# Merge tables to add code
prov_table_with_icode <- merge(prov_table, pop_sub, by='location_name', all=TRUE)
if(nrow(prov_table_with_icode) != 107) stop("ISSUE: Wrong number of provinces!")

## Prepare and merge on region and macroregion identifiers
region_codes <- unique(data.table::fread(
  config$paths$raw_location_codes, na.strings = c("", "n.d.")
)[, c(9, 10, 1, 11, 15)])
colnames(region_codes) <- c(
  'macroregion_code', 'macroregion_name', 'region_code', 'region_name', 'abbrev'
)
# Replace some region names with special characters
region_codes[region_name=="Valle d'Aosta/Vall\xe9e d'Aoste", region_name := "Valle d'Aosta"]
region_codes[region_name=="Trentino-Alto Adige/S\xfcdtirol", region_name := "Trentino-Alto Adige"]
# Merge onto province-level data
prov_table_full <- merge(
  x = region_codes, y = prov_table_with_icode, by = c('region_code', 'abbrev'), all=TRUE
)
if(nrow(prov_table_full) != 107){
  stop("ISSUE: Wrong number of provinces after region merge!")
}
prov_table_full <- prov_table_full[order(location_code)]


## Create output directory and save --------------------------------------------

out_dir <- file.path(config$paths$prepped_data, prepped_data_version)
dir.create(out_dir, showWarnings=FALSE)

out_fp <- file.path(out_dir, config$prepped_data_files$location_table)
data.table::fwrite(prov_table_full, file=out_fp)
