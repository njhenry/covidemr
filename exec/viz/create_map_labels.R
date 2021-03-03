## ---------------------------------------------------------------------------------------
##
## VISUALIZATION: Create the tables of label names and positions to plot in maps of Italy
##
## ---------------------------------------------------------------------------------------

## Setup -------------------------------------------------------------------------------->

library(data.table)
library(dplyr)
library(ggplot2)
library(sp)
library(sf)

## Which provinces should be labeled?
top_provs <- c('Milano', 'Bergamo', 'Brescia')

## Filepaths
dev_fp <- '~/repos/covidemr/'
prepped_data_version <- '20210113'
devtools::load_all(dev_fp)


## Data loading
config <- yaml::read_yaml(file.path(dev_fp, 'inst/extdata/config.yaml'))
config$paths$prepped_data <- '~/temp_data/ita_emr/prepped_data/'
get_prep_fp <- function(ff) file.path(
  config$paths$prepped_data, prepped_data_version,
  config$prepped_data_files[[ff]]
)

location_table <- data.table::fread(get_prep_fp('location_table'), na.strings="")
polys_sf <- readRDS(get_prep_fp('shapefile_sf'))

# Merge extra location identifiers onto the spatial polygons
merge_on_cols <- c('location_code', 'location_name', 'region_code', 'region_name')
polys_sf <- merge(polys_sf, location_table[, ..merge_on_cols], by='location_code')
# Create a regional sf object by dissolving provinces
regions_sf <- aggregate(polys_sf, by=list(polys_sf$region_code), FUN=first)


## PREP PART 1: Default labels located at centroids of locations - fields for starting and
##  ending line segments

# Regions
region_points <- suppressWarnings(
  sf::st_coordinates(sf::st_centroid(regions_sf))
)
region_labels_default <- data.table::data.table(
  label_name = regions_sf$region_name,
  x = region_points[,'X'], y = region_points[,'Y'],
  x1_lseg = as.numeric(NA), y1_lseg = as.numeric(NA),
  x2_lseg = as.numeric(NA), y2_lseg = as.numeric(NA),
  color = '#000000'
)
data.table::fwrite(
  region_labels_default,
  file = gsub('.csv', '_default.csv', get_prep_fp('region_labels'))
)

# Provinces
province_points <- suppressWarnings(
  sf::st_coordinates(sf::st_centroid(polys_sf))
)
province_labels_default <- data.table::data.table(
  label_name = polys_sf$location_name,
  x = province_points[,'X'], y = province_points[,'Y'],
  x1_lseg = as.numeric(NA), y1_lseg = as.numeric(NA),
  x2_lseg = as.numeric(NA), y2_lseg = as.numeric(NA),
  color = '#000000'
)
data.table::fwrite(
  province_labels_default,
  file = gsub('.csv', '_default.csv', get_prep_fp('province_labels'))
)


## PREP PART 2 -------------------------------------------------------------------------->
## Test mapping for updated labels

# Create blank plots
polys_sf$dummy <- 'a'
regions_sf$dummy <- 'a'

prov_fig <- covidemr::map_ita_choropleth(
  province_sf = polys_sf, region_sf = regions_sf, in_data = as.data.table(polys_sf),
  map_field = 'dummy', fill_list = list(values = '#FFFFFF', limits = 'a'),
  titles_list = list(title = "Testing"), fill_type = 'manual'
)
reg_fig <- covidemr::map_ita_choropleth_region(
  region_sf = regions_sf, in_data = as.data.table(regions_sf), map_field='dummy',
  fill_list = list(values = '#FFFFFF', limits = 'a'),
  titles_list = list(title = "Testing"), fill_type = 'manual'
)


# Plot tests!
province_test_fp <- '~/Desktop/test_province_labels.png'
province_labels <- data.table::fread(get_prep_fp('province_labels'))
prov_with_labels <- prov_fig +
  geom_text(
    data = province_labels,
    aes(x=x, y=y, label=label_name), hjust='middle', vjust='center', size=3
  ) +
  geom_segment(
    data = province_labels,
    aes(x=x1_lseg, y=y1_lseg, xend=x2_lseg, yend=y2_lseg)
  )
png(province_test_fp, height = 7.25, width = 6, units='in', res=300)
print(prov_with_labels)
dev.off()

region_test_fp <- '~/Desktop/test_region_labels.png'
region_labels <- data.table::fread(get_prep_fp('region_labels'))
region_labels[, label_name := gsub('\\\\n','\n', label_name)]
reg_with_labels <- reg_fig +
  geom_text(
    data = region_labels,
    aes(x=x, y=y, label=label_name), hjust='middle', vjust='center', size=3, lineheight=.75
  ) +
  geom_segment(
    data = region_labels,
    aes(x=x1_lseg, y=y1_lseg, xend=x2_lseg, yend=y2_lseg)
  )
png(region_test_fp, height = 7.25, width = 6, units='in', res=300)
print(reg_with_labels)
dev.off()
