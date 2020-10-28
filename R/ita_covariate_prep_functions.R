
#' Prepare an Italy covariate for modeling
#'
#' @description Generic function to prepare a covariate for modeling. This
#'   function is a wrapper for sub-functions called for specific covariates and
#'   will fail if a non-initialized covariate has been passed. It returns both
#'   the prepared covariate and a vector of expected indices that can be used
#'   to merge back onto the original dataset.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param covar_name Short name for the covariate. This will be used to send the
#'   covariate to a particular prep function.
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#' @param pop_raster [optional] Population rasterBrick object, with one layer
#'   per modeling year. Only used to prepare raster covariates, default NULL.
#' @param polys_sf [optional] Polygon boundaries
#' @param projection [optional] Character vector giving the proj4 CRS
#'   definition. Example: "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
#'
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table
#' @export
ita_prepare_covariate <- function(
  covar_name, covar_fp, model_years, location_table, pop_raster = NULL,
  polys_sf = NULL, projection = NULL
){
  # Check that an appropriate covariate function exists
  covar_func <- get_covar_prep_function(covar_name)
  # Apply the function to the data
  message(" - Preparing ",covar_name)
  covar_out_list <- do.call(
    covar_func,
    args = list(
      covar_fp = covar_fp,
      model_years = model_years,
      location_table = location_table,
      pop_raster = pop_raster,
      polys_sf = polys_sf,
      projection = projection
    )[names(formals(covar_func))]
  )
  # Check that output is valid
  check_covar_validity(
    prepped_covar = covar_out_list$prepped_covar,
    covar_name = covar_name,
    covar_indices = covar_out_list$covar_indices,
    model_years = model_years,
    location_table = location_table
  )
  # Return output
  return(covar_out_list)
}


#' Pull a specific prep function for an Italy covariate
#'
#' @description Helper function for \code{\link{prepare_covariate}} that checks
#'   whether a covariate-specific prep function exists yet by searching the
#'   package namespace. All covariate-specific prep functions should take the
#'   name format "ita_prep_covar_<COVAR NAME>". Errors if the function does not
#'   exist; returns the function name if the function does exist.
#'
#' @param covar_name Short name for the covariate.
#'
#' @return Name of the function to prepare this covariate
#'
get_covar_prep_function <- function(covar_name){
  # Search for prep functions
  f_pattern <- 'ita_prepare_covar_'
  potential_f_name <- paste0(f_pattern, covar_name)
  all_prep_functions <- ls('package:covidemr', pattern=sprintf("^%s",f_pattern))
  if(potential_f_name %in% all_prep_functions){
    return(potential_f_name)
  } else {
    stop(sprintf(
      "The function %s does not exist! Valid options include:\n - %s",
      potential_f_name,
      paste(all_prep_functions, collapse='\n - ')
    ))
  }
}


#' Validate prepared covariate data
#'
#' @description Given prepared covariate data, validate that the correct fields
#'   are present and that no data is missing
#'
#' @param prepped_covar data.table containing the prepped covariate
#' @param covar_name Short name for the covariate. This will be used to send the
#'   covariate to a particular prep function.
#' @param covar_indices Vector of identifiers that should be used to merge onto
#'   the modeling dataset. The prepared covariate dataset should only include
#'   identifier columns and the covariate value, specified by the covariate name
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return NULL (fails if there are any issues in prepared covariate data)
#'
#' @import data.table
check_covar_validity <- function(
  prepped_covar, covar_name, covar_indices, model_years, location_table
){
  # Check that both exist
  if(is.null(prepped_covar)) stop("Missing prepared covariate")
  if(is.null(covar_indices)) stop("Missing prepared covariate indices")
  # Check for correct columns
  all_colnames <- c(covar_name, covar_indices)
  missing_cols <- setdiff(all_colnames, names(prepped_covar))
  extra_cols <- setdiff(names(prepped_covar), all_colnames)
  if(length(c(missing_cols, extra_cols)) > 0){
    if(missing_cols) message("Missing columns: ", paste(missing_cols, collapse=', '))
    if(extra_cols) message("Extra columns: ", paste(extra_cols, collapse=', '))
    stop("Resolve column issues before proceeding")
  }

  # If an index exists, check that all needed values are included
  # Extra values are fine but will be dropped when merging on the data
  # Helper function to check for missing values
  check_missing <- function(idx_col, all_vals){
    any_missing <- FALSE
    if(idx_col %in% covar_indices){
      missing_vals <- setdiff(all_vals, unique(prepped_covar[[idx_col]]))
      if(length(missing_vals) > 0){
        message("Missing ", idx_col, ": ", paste0(missing_vals, collapse=', '))
        any_missing <- TRUE
      }
    }
    # Return TRUE if there is an error, and FALSE otherwise
    return(any_missing)
  }
  any_missing <- c(
    check_missing('year', model_years),
    check_missing('sex', c('male','female')),
    check_missing('week', 1:52),
    check_missing('location_code', location_table$location_code)
  )
  if(any(any_missing)) stop("Resolve missingness errors to proceed.")

  # Validation passed!
  invisible(NULL)
}


#' Extend a covariate time series forward
#'
#' @description Helper function to extend a covariate time series to the last
#'   year modeled by copying forward the most recent year of data available.
#'
#' @param covar_data data.table containing covariates with a 'year' field
#' @param model_year Vector of years to be included in modeling
#'
#' @return covariate data.table, with years extended to the most recent year of
#'   data
#'
#' @import data.table
#' @export
extend_covar_time_series <- function(covar_data, model_years){
  # Validate inputs: ensure that the year field exists
  if(!('year' %in% names(covar_data))) stop("'year' field missing from data")

  missing_years <- setdiff(model_years, covar_data$year)
  # If no years are missing, return the original data.table
  if(length(missing_years) == 0) return(covar_data)

  # Extend forward
  most_recent_year <- covar_data[year==max(year),]
  full_data <- rbindlist(c(
      list(covar_data),
      lapply(missing_years, function(yr) copy(most_recent_year)[, year := yr])
  ))[order(year)]

  message("    Extended forward for years ",paste(missing_years, collapse=', '))
  return(full_data)
}


#' Prepare TFR covariate
#'
#' @description Covariate-specific prep function for TFR (total fertility rate).
#'   This is a convenience function for specific use with datasets exported from
#'   IStat.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_tfr <- function(covar_fp, model_years, location_table){

  covar_data <- ita_read_istat_csv(covar_fp)
  data.table::setnames(
    covar_data, c('ITTER107','TIME','Value'), c('icode', 'year', 'tfr')
  )

  # Merge on location codes
  covar_data_merged <- data.table::merge.data.table(
    covar_data,
    location_table[, .(icode, location_code)],
    by='icode'
  )
  # Apply Sud Sardegna fix
  covar_data_merged <- backfill_input_data(
    input_data = covar_data_merged,
    index_field = 'icode',
    check_vals = 'IT111'
  )
  # Drop unnecessary columns
  covar_indices <- c('location_code', 'year')
  all_cols <- c(covar_indices, 'tfr')
  covar_data_merged <- covar_data_merged[, ..all_cols ]

  # Extend years
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare unemployment covariate
#'
#' @description Covariate-specific prep function for unemployment by sex and by
#'   quarter. This is a convenience function for specific use with datasets
#'   exported from IStat.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_unemp <- function(covar_fp, model_years, location_table){

  covar_data <- ita_read_istat_csv(covar_fp)
  covar_indices <- c('location_code','sex','year')

  # Update names
  setnames(covar_data, c('ITTER107','TIME','Value'), c('icode', 'year', 'unemp'))

  # Keep only annual values and convert to integer
  suppressWarnings(covar_data[, year := as.integer(year) ])
  covar_data <- covar_data[ !is.na(year), ]

  # Format sex identifiers
  covar_data[, sex := gsub('s', '', Gender) ]

  # Fix Bolzano and Trento location codes, then merge on standard code table
  covar_data[ icode == 'ITD1', icode := 'ITD10' ] # Bolzano
  covar_data[ icode == 'ITD2', icode := 'ITD20' ] # Trento
  covar_data_merged <- data.table::merge.data.table(
    covar_data,
    location_table[, .(location_code, icode)],
    by='icode'
  )
  # Apply Sud Sardegna fix
  covar_data_merged <- backfill_input_data(
    input_data = covar_data_merged,
    index_field = 'icode',
    check_vals = 'IT111'
  )
  covar_data_merged <- covar_data_merged[, c(covar_indices, 'unemp'), with=FALSE]

  # Subset columns and extend to 2020 and return
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare social services covariate
#'
#' @description Covariate-specific prep function for social services coverage
#'   for families. This is a convenience function for specific use with datasets
#'   exported from IStat.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_socserv <- function(covar_fp, model_years, location_table){

  covar_data <- ita_read_istat_csv(covar_fp)
  covar_indices <- c('year', 'location_code')

  # Update names
  setnames(covar_data, c('ITTER107', 'TIME', 'Value'), c('icode', 'year', 'socserv'))
  # Subset to just proportion of families/children receiving at home social services
  covar_data <- covar_data[(TIPUTENZA1=='FAM') & (TIPSERVSOC=='HOMECARE'), ]

  # Fix for Sud Sardegna, which is represented by its old defunct provinces.
  # Note that this fix differs from others because Sud Sardegna is not
  # represented in any years -- instead, the average is taken between Medio
  # Campidano and Carbonia-Iglesias
  covar_ssd <- covar_data[ icode %in% c('ITG2B', 'ITG2C'), .(icode, year, socserv)]
  covar_ssd[, icode := 'IT111' ]
  covar_ssd <- covar_ssd[, .(socserv=mean(socserv)), by=.(icode, year)]
  covar_data <- rbindlist(list(covar_data, covar_ssd), use.names=TRUE, fill=TRUE)

  # Merge on location codes
  covar_data_merged <- data.table::merge.data.table(
    covar_data,
    location_table[, .(icode, location_code)],
    by='icode'
  )

  covar_data_merged <- covar_data_merged[, c(covar_indices, 'socserv'), with=FALSE]
  prepped_covar <- extend_covar_time_series(
    covar_data = covar_data_merged,
    model_years = model_years
  )

  # Return
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare covariate: proportion of households with income under E10k
#'
#' @description Covariate-specific prep function for proportion of households
#'   with taxable income under 10,000 Euros per year. A related presentation of
#'   this data is mean taxable income per household, prepared in another
#'   function. This function is meant specifically for use with IStat data.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_tax_brackets <- function(
  covar_fp, model_years, location_table
){

  covar_indices <- c('year', 'location_code')
  covar_data <- ita_read_istat_csv(covar_fp)
  setnames(covar_data, 'TIME', 'year')

  # Set the province code based on the comuna code
  covar_data[, ITTER107 := as.character(ITTER107)]
  covar_data[, loc_nc := nchar(ITTER107) ]
  covar_data[, location_code := as.integer(substr(ITTER107, 1, loc_nc - 3))]

  # Some provinces in Sardinia were reorganized in 2018: try to reassign the
  #  provinces for those municipalities based on 2018 data
  sard_provinces <- location_table[region_code == 20, location_code ]
  sard_merge_table <- unique(covar_data[
    (year == 2018) & (location_code %in% sard_provinces),
    .(Territory, location_code)
  ])

  old_provs <- 104:107
  existing_provinces <- covar_data[ !(location_code %in% old_provs), ]
  sard_pre_2018 <- covar_data[location_code %in% old_provs,]
  sard_pre_2018[
    sard_merge_table, on='Territory', location_code := i.location_code
  ]
  # Reconstitute the full dataset with updated province codes
  covar_data <- rbindlist(list(existing_provinces, sard_pre_2018), use.names=TRUE)

  # Get the number of households with income less than 15k Euros
  covar_data[, low_bk := IMPORTOEURO %in% c('E_UN0', 'E0-10000')]
  prepped_covar <- covar_data[,
    .(tax_brackets = .SD[low_bk==1, sum(Value, na.rm=T)] / sum(Value, na.rm=T)),
    by = covar_indices
  ]
  # Extend for 2019-2020 and return
  prepped_covar <- extend_covar_time_series(
    covar_data = prepped_covar,
    model_years = model_years
  )
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}

#' Prepare covariate: mean taxable income across all households
#'
#' @description Covariate-specific prep function for average taxable income
#'   across all households in a province. This covariate is derived from two
#'   IStat datasets, one listing the number of households per tax bracket and
#'   another listing total taxed income by bracket. A related presentation of
#'   this data is proportion of households with less than 10k Euros in taxable
#'   income, prepared in the function `ita_prepare_covar_tax_brackets()`.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#'
#' @return data.table containing prepared covariate data and vector of indices
#'
#' @import data.table
#' @export
ita_prepare_covar_tax_income <- function(
  covar_fp, model_years, location_table
){
  # Load data
  covar_indices <- c('year', 'location_code')
  covar_data_long <- rbindlist(
    lapply(covar_fp, ita_read_istat_csv)
  )[, .(ITTER107, Territory, TIPO_DATO_MEF, TIME, Value)]
  # Cast so that total income and number of households are in different fields
  covar_data <- data.table::dcast(
    covar_data_long,
    ITTER107 + Territory + TIME ~ TIPO_DATO_MEF,
    fun.aggregate = function(x) sum(x, na.rm=T),
    value.var = 'Value'
  )
  setnames(covar_data, c('TIME','AGGINCF','AGGINCR'), c('year','hh','income'))

  # Set the province code based on the comuna code
  covar_data[, ITTER107 := as.character(ITTER107)]
  covar_data[, loc_nc := nchar(ITTER107) ]
  covar_data[, location_code := as.integer(substr(ITTER107, 1, loc_nc - 3))]

  # Some provinces in Sardinia were reorganized in 2018: try to reassign the
  #  provinces for those municipalities based on 2018 data
  sard_provinces <- location_table[region_code == 20, location_code ]
  sard_merge_table <- unique(covar_data[
    (year == 2018) & (location_code %in% sard_provinces),
    .(Territory, location_code)
  ])

  old_provs <- 104:107
  existing_provinces <- covar_data[ !(location_code %in% old_provs), ]
  sard_pre_2018 <- covar_data[location_code %in% old_provs,]
  sard_pre_2018[
    sard_merge_table, on='Territory', location_code := i.location_code
  ]
  # Reconstitute the full dataset with updated province codes
  covar_data <- rbindlist(list(existing_provinces, sard_pre_2018), use.names=TRUE)

  # Get the number of households with income less than 15k Euros
  prepped_covar <- covar_data[,
    .(tax_income = sum(income) / sum(hh)),
    by = covar_indices
  ]
  # Extend for 2019-2020 and return
  prepped_covar <- extend_covar_time_series(
    covar_data = prepped_covar,
    model_years = model_years
  )
  return(list(prepped_covar = prepped_covar, covar_indices = covar_indices))
}


#' Prepare covariate: population density
#'
#' @description Covariate-specific prep function for all-ages population density
#'   by Italian province
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#' @param polys_sf [optional] Polygon boundaries
#'
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table sf
#' @export
ita_prepare_covar_pop_density <- function(
  covar_fp, model_years, location_table, polys_sf
){
  # Load annual population data
  pop_by_age <- data.table::fread(covar_fp)
  pop_agg <- pop_by_age[, .(pop = sum(pop)), by=.(location_code, year)]
  pop_agg <- pop_agg[order(location_code, year)]

  # Generate areas for each location based on the projected polygons
  loc_areas <- data.table::data.table(
    location_code = polys_sf$location_code,
    area_sqkm = as.vector(sf::st_area(polys_sf) / 1E6) # Coming from units in m
  )
  pop_agg[ loc_areas, area_sqkm := as.numeric(i.area_sqkm), on = 'location_code'
    ][, pop_density := pop / area_sqkm # Units: people / km2
    ][, c('pop','area_sqkm') := NULL ]

  return(list(
    prepped_covar = pop_agg,
    covar_indices = c('location_code','year')
  ))
}


#' Prepare covariate: Elevation
#'
#' @description Covariate-specific prep function for elevation by location.
#'   Specifically, this aggregated covariate represents the population-weighted
#'   mean elevation inhabited by residents of each province.
#'
#' @param covar_fp Filepath to the raw covariate
#' @param pop_raster [optional] Population rasterBrick object, with one layer
#'   per modeling year. Only used to prepare raster covariates, default NULL.
#' @param polys_sf [optional] Polygon boundaries
#' @param projection [optional] Character vector giving the proj4 CRS
#'   definition. Example: "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
#'
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table raster fasterize
#' @export
ita_prepare_covar_elevation <- function(
  covar_fp, pop_raster, polys_sf, projection
){
  # Load raw data
  elev_raw <- raster::raster(covar_fp)
  # Get population in final year (mean elevation should be stable across years)
  pop <- pop_raster[[dim(pop_raster)[3]]]

  # Align with population raster and mask
  elev_projected <- raster::projectRaster(from = elev_raw, crs = projection)
  elev_sub <- raster::crop(x = elev_projected, y = pop)
  elev_sub <- raster::resample(x = elev_sub, y = pop, method = 'bilinear')
  elev_sub <- raster::mask(x = elev_sub, mask = pop)

  locs_raster <- fasterize::fasterize(
    sf = polys_sf, raster = pop, field = 'location_code', fun = 'last'
  )
  # If all dimensions align, run the population-weighted aggregation
  if(any(dim(elev_sub) != dim(pop))) stop("Access raster not aligned")
  if(any(dim(pop) != dim(locs_raster))) stop("Location raster not aligned")
  prepped_covar <- na.omit(data.table::data.table(
    elevation = as.vector(elev_sub),
    population = as.vector(pop),
    location_code = as.vector(locs_raster)
  ))[, .(elevation = weighted.mean(elevation, w=population)), by=location_code]
  prepped_covar <- prepped_covar[order(location_code)]

  return(list(prepped_covar = prepped_covar, covar_indices = 'location_code'))
}


#' Prepare covariate: access to healthcare facility
#'
#' @description Covariate-specific prep function for the average person's travel
#'   time to the nearest healthcare facility (by motor vehicle).
#'
#' @param covar_fp Filepath to the raw covariate
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#' @param pop_raster [optional] Population rasterBrick object, with one layer
#'   per modeling year. Only used to prepare raster covariates, default NULL.
#' @param polys_sf [optional] Polygon boundaries
#' @param projection [optional] Character vector giving the proj4 CRS
#'   definition. Example: "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
#'
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table raster fasterize
#' @export
ita_prepare_covar_hc_access <- function(
  covar_fp, model_years, location_table, pop_raster, polys_sf, projection
){

  rast_raw <- raster::raster(covar_fp)
  # Only use the final year of population data, which aligns with the access
  #  dataset
  pop <- pop_raster[[length(model_years)]]

  # Crop in unprojected space first - this will make the projection much faster
  unproj_extent <- raster::extent(
    raster::projectRaster(pop, crs = crs(rast_raw))
  )
  buffered_extent <- raster::extent(
    unproj_extent@xmin - .1, unproj_extent@xmax + .1,
    unproj_extent@ymin - .1, unproj_extent@ymax + .1
  )
  rast_sub <- raster::crop(rast_raw, y = buffered_extent)

  # Project, crop, resample and mask to population raster
  rast_sub <- raster::projectRaster(rast_sub, crs = projection)
  rast_sub <- raster::crop(x = rast_sub, y = pop)
  rast_resampled <- raster::resample(x = rast_sub, y = pop, method = 'ngb')
  rast_resampled <- raster::mask(x = rast_resampled, mask = pop)

  # Rasterize polygons to get location code for each pixel, and get "average
  #  access" using a population-weighted aggregation
  locs_raster <- fasterize::fasterize(
    sf = polys_sf, raster = pop, field = 'location_code', fun = 'last'
  )
  # If all dimensions align, run the population-weighted aggregation
  if(any(dim(rast_resampled) != dim(pop))) stop("Access raster not aligned")
  if(any(dim(pop) != dim(locs_raster))) stop("Location raster not aligned")
  prepped_covar <- na.omit(data.table::data.table(
    hc_access = as.vector(rast_resampled),
    population = as.vector(pop),
    location_code = as.vector(locs_raster)
  ))[, .(hc_access = weighted.mean(hc_access, w=population)), by=location_code]
  prepped_covar <- prepped_covar[order(location_code)]

  return(list(prepped_covar = prepped_covar, covar_indices = 'location_code'))
}


#' Load and cache daily temperature data from API
#'
#' @description Pull daily temperature data from the Meteostat API. This
#'   function uses cacheing: if the `cache_file` already exists, it pulls the
#'   results from file rather than re-querying the API
#'
#' @param model_years Years of data to pull
#' @param api_key API key for querying the Meteostat API
#' @param prov_latlong Table indexing locations to latitute and longitude. Must
#'   contain at least the following three fields:
#'    - 'location_code': Numeric code identifying location
#'    - 'lat': Latitude
#'    - 'lon': Longitude
#' @param cache_file File where the pulled data will be cached
#' @param refresh_cache [optional, default FALSE] Should the data be re-pulled
#'   even if a cache file exists?
#'
#' @return Data.table containing full Meteostat daily API data as well as
#'  location code for all available years
#'
#' @import data.table glue httr
#'
ita_temperature_query_api <- function(
  model_years, api_key, prov_latlong, cache_file, refresh_cache = FALSE
){
  # Load and return cached file if it exists
  if(file.exists(cache_file) & !refresh_cache){
    return(fread(cache_file))
  }

  # Prepare a list to store individual API query results
  temp_data_list <- vector('list', length=length(model_years))

  # Iterate through years and locations to save on memory
  for(year_idx in 1:length(model_years)){
    model_year <- model_years[year_idx]
    message("  Querying all province data for year ", model_year)
    # For each province, query the API for daily weather patterns across the
    #  year
    temp_data_list[[year_idx]] <- vector('list', length=nrow(prov_latlong))
    for(pnt_idx in 1:nrow(prov_latlong)){
      cat('.'); flush.console();
      this_prov_code <- prov_latlong[pnt_idx, location_code]
      # Build query
      meteo_query <- glue::glue(
        "https://api.meteostat.net/v2/point/daily?lat={prov_latlong[pnt_idx, lat]}",
        "&lon={prov_latlong[pnt_idx, lon]}&start={model_year}-01-01",
        "&end={model_year}-12-31"
      )
      # Send query with header
      meteo_response <- httr::GET(meteo_query, add_headers('x-api-key'=api_key))
      # Wait between requests
      Sys.sleep(0.55)
      response_code <- httr::status_code(meteo_response)
      if(response_code == 200){
        # The query succeeded
        response_dataset <- suppressWarnings(data.table::rbindlist(
          httr::content(meteo_response)$data, use.names = TRUE
        ))
        response_dataset[, location_code := this_prov_code ]
        response_dataset[, point_id := pnt_idx ]
        response_dataset$lat <- prov_latlong[pnt_idx, lat]
        response_dataset$lon <- prov_latlong[pnt_idx, lon]

        if(nrow(response_dataset)==1){
          message("    W: No temperature data for province ",this_prov_code," in ",model_year)
        } else {
          temp_data_list[[year_idx]][[pnt_idx]] <- response_dataset
        }

      } else {
        # The query failed
        stop(
          "API request failed with code ", response_code, " on year ",
          model_year, " and location code ", this_prov_code
        )
      }
    }
  }

  # Compile full dataset, save to file, and return
  raw_data_full <- data.table::rbindlist(
    lapply(temp_data_list, rbindlist, fill=TRUE)
  )
  fwrite(raw_data_full, file = cache_file)
  return(raw_data_full)
}


#' Prepare covariate: weekly temperature
#'
#' @description Covariate-specific prep function for average experienced weekly
#'   temperature by province. Note that this function uses cacheing, so if the
#'   raw data has been pulled from the API already it will not be re-pulled
#'
#' @param covar_fp Filepath to an API key for Meteostat, which will be used to
#'   pull weekly templerature data
#' @param model_years Vector of years to consider for modeling
#' @param location_table Location code merge table for Italy
#' @param pop_raster [optional] Population rasterBrick object, with one layer
#'   per modeling year. Only used to prepare raster covariates, default NULL.
#' @param polys_sf [optional] Polygon boundaries
#''
#' @return A list of two items:
#'   - "prepped_covar": data.table containing the prepped covariate
#'   - "covar_indices": Vector of identifiers that should be used to merge onto
#'     the modeling dataset. The prepared covariate dataset should only include
#'     identifier columns and the covariate value, specified by the covariate
#'     name
#'
#' @import data.table raster
#' @export
ita_prepare_covar_temperature <- function(
  covar_fp, model_years, location_table, pop_raster, polys_sf
){
  ## Find the pixel with the highest population in each province
  ## The API will be queried at these points
  prov_indexing <- index_populated_grid_cells(
    pop = pop_raster[[length(model_years)]],
    polys_sf = polys_sf,
    loc_field = 'location_code'
  )
  prov_latlong <- prov_indexing[, .SD[frank(-pop) <= 3,], by=location_code
    ][, pop := NULL
    ][, `:=` (lat = round(lat, 3), lon = round(lon, 3))
    ][order(location_code)]
  if(nrow(prov_latlong) != nrow(location_table) * 3) stop("Missing some provinces")

  # Query the Meteostat API for daily temperatures in selected locations
  api_key <- readLines(covar_fp)
  meteostat_folder <- file.path(dirname(api_key), 'meteostat')
  dir.create(meteostat_folder, showWarnings = FALSE)
  raw_temp_data <- ita_temperature_query_api(
    model_years = model_years,
    api_key = api_key,
    prov_latlong = prov_latlong,
    cache_file = file.path(meteostat_folder, 'api_data_compiled.csv'),
    refresh_cache = FALSE
  )

  # Aggregate temperature from days to weeks
  temp_weekly <- raw_temp_data[, date := as.Date(date, format = '%Y-%m-%d')
    ][, year := as.integer(strftime(date, format='%Y'))
    ][, week := as.integer(ceiling(as.numeric(strftime(date, format='%j'))/7))
    ][ week > 52, week := 52
    ][, .(tavg = mean(tavg, na.rm=T)), by=.(location_code, point_id, year, week)]

  # Merge on a template to ensure that all points, weeks, and years are included
  template <- data.table::merge.data.table(
    x = data.table::data.table(
      location_code = rep(location_table$location_code, each = 3),
      point_id = 1:(3 * nrow(location_table)),
      dummy = 1
    ),
    y = data.table::CJ(year = model_years, week = 1:52, dummy = 1),
    by = 'dummy',
    allow.cartesian = TRUE
  )[, dummy := NULL ]
  temp_weekly <- data.table::merge.data.table(
    x = template, y = temp_weekly, by = names(template), all.x=TRUE
  )[order(location_code, point_id, year, week)]
  # Subset to only included weeks (up through 2020 week 26)
  temp_weekly <- temp_weekly[(year < 2020) | (week <= 26), ]

  # For small gaps (less than 4 weeks), linearly interpolate available data
  num_na_start <- sum(is.na(temp_weekly$tavg))

  temp_weekly[, allna := all(is.na(tavg)), by = point_id]
  temp_weekly[
    allna == FALSE,
    temp_interp := stats::approx(x = .I, y=tavg, xout=.I, method='linear')$y,
    by=point_id
  ]
  temp_weekly[, na_grouping := data.table::rleid(is.na(tavg)), by = point_id
    ][, num_missing := .N, by = .(na_grouping, point_id)
    ][ is.na(tavg) & (num_missing < 4), tavg := temp_interp
    ][, c('num_missing', 'temp_interp', 'na_grouping') := NULL ]
  num_na_end <- sum(is.na(temp_weekly$tavg))

  # Aggregate to the province level
  temp_by_province <- temp_weekly[,
    .(tavg = mean(tavg, na.rm=T)),
    by=.(location_code, year, week)
  ]

  # When a weekly temperature is not available for a given province, pull it
  #  from a neighboring province
  reassign_dt <- data.table::rbindlist(lapply(list(
    c(1,2), c(4,9), c(5,18), c(6,9), c(10,9), c(11,45), c(14,13), c(21,25),
    c(22,25), c(29,28), c(33,19), c(34,19), c(35,20), c(36,20), c(37,39),
    c(38,39), c(42,41), c(43,41), c(44,41), c(48,50), c(51,54), c(52,54),
    c(55,54), c(57,54), c(64,62), c(65,63), c(66,69), c(67,69), c(68,69),
    c(71,77), c(72,74), c(76,77), c(78,77), c(79,77), c(82,83), c(84,81),
    c(85,81), c(88,87), c(91,90), c(95,111), c(98,18), c(100,47), c(101,80),
    c(102,77), c(109,41), c(110,74)
  ), as.list))
  names(reassign_dt) <- c('target', 'source')
  temp_by_province[ reassign_dt, source := i.source, on = c(location_code = 'target')
    ][
      temp_by_province,
      tavg_neighbor := i.tavg,
      on = c(source = 'location_code', week = 'week', year = 'year')
    ][, from_neighbor := as.integer(is.na(tavg))
    ][ from_neighbor == 1, tavg := tavg_neighbor
    ][, c('source','tavg_neighbor','from_neighbor') := NULL ]

  # Rename value column
  data.table::setnames(temp_by_province, 'tavg', 'temperature')

  # Return
  return(list(
    prepped_covar = temp_by_province,
    covar_indices = c('location_code','year','week')
  ))
}

