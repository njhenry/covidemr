#' Load and format spatial polygons for modeling
#'
#' @description Load a spatial polygons object in SF format, assign location
#'   codes based on each province abbreviation, and return polygons objects in
#'   both SF and SP formats
#'
#' @param shp_in_fp Input filepath to a spatial object (usually a shapefile)
#' @param abbrev_field Field listing two-character province abbreviations
#' @param location_table Data.frame listing province information. Should contain
#'   at least the following fields:
#'     - 'abbrev': Two-character province abbreviations (eg 'TO', 'VC')
#'     - 'location_code': Province codes that will be assigned to the spatial
#'          polygons object
#'
#' @return List with two items:
#'   - 'shp_sf': Formatted spatial polygons, in SF format
#'   - 'shp_sp': Formatted spatial polygons, in SP format
#'
#' @import sp sf
#' @export
load_format_spatial_polys <- function(shp_in_fp, abbrev_field, location_table){
  shp_sf <- sf::st_read(shp_in_fp)
  # Order by location code
  shp_sf$location_code <- unlist(sapply(
    shp_sf[[abbrev_field]],
    function(this_abbrev){
      location_table[location_table$abbrev==this_abbrev, 'location_code']
    }
  ))
  shp_sf <- shp_sf[order(shp_sf$location_code), ]

  # Create SP version
  shp_sp <- as(shp_sf, "Spatial")
  # Reset polygon IDs to reduce confusion
  for(ii in 1:nrow(shp_sp@data)){
    shp_sp@polygons[[ii]]@ID <- as.character(shp_sp@data[[ii, 'location_code']])
  }

  return(list(shp_sf = shp_sf, shp_sp = shp_sp))
}


#' Build adjacency matrix
#'
#' @description Builds a matrix describing adjacency between polygons in the
#'   a spatial object. This function is a wrapper of two functions in the
#'   \code{\link{spdep}} package.
#'
#' @param poly_sp A SpatialPolygonsDataFrame object
#' @param style One of "B","W","C", or "S". For more information, see the
#'   documentation in the \code{\link{spdep}} package:
#'   https://www.rdocumentation.org/packages/spdep/versions/1.1-3/topics/nb2mat
#' @param allow_zero_neighbors [bool, default TRUE] Allow polygons with no
#'   neighbors?
#' @param manually_add_links [list, default NULL] A list of links to add in
#'   addition to those specified by polygon adjacency. This should be formatted
#'   as a list where each item is a numeric vector of length 2, specifying the
#'   indices of polygons that should be listed as adjacent. The numeric vectors
#'   are not sensitive to ordering (so item `c(1, 3)` is the same as `c(3, 1)`)
#'
#' @return sparse dsCMatrix representing adjacency between polygons in the
#'   `poly_sp` object
#'
#' @import spdep
#' @export
build_adjacency_matrix <- function(
    poly_sp, style='B', allow_zero_neighbors=TRUE, manually_add_links=list()
){
  ## Generate adjacency matrix for the polygon
  adjmat <- spdep::nb2mat(
    neighbours = spdep::poly2nb(poly_sp),
    style = style,
    zero.policy = allow_zero_neighbors
  )
  adjrows <- nrow(adjmat)
  adjcols <- ncol(adjmat)
  message("Created adjacency matrix with dimensions ",adjrows,"x",adjcols)

  ## Manually add links when specified
  if(length(manually_add_links) > 0){
    for(add_link in manually_add_links){
      # Check that the link specification is valid
      if(length(add_link) != 2) stop("Manual link has length != 2")
      for(idx in c(1, 2)){
        if((add_link[idx] < 0) | (add_link[idx] > adjrows)) stop(
          "Invalid manual link ID (must be positive integer less than ",adjrows,
          "): ", add_link[idx]
        )
      }
      adjmat[add_link[1], add_link[2]] <- 1
      adjmat[add_link[2], add_link[1]] <- 1
    }
  }

  return(adjmat)
}


#' Load and format population raster
#'
#' @description Load a population raster from file, clip it to the given
#'   shapefile, and return as a raster brick for the years specified
#'
#' @param pop_fp_format Each of the worldpop rasters are downloaded separately
#'   by year. This format specifies the location of each raster, with the
#'   '{year}' standing in for a particular year. For example:
#'   `pop_filepath_format = '/path/to/pop/raster/my_raster_{year}.tif'`
#' @param model_years The years needed for the model. The function will load
#'   separate raster layers for each model year based on the filepath template
#' @param country_polys The spatial polygons object (in sf format) that the
#'  raster will be clipped to
#' @param projection Projection used for spatial data throughout the model
#'
#' @return rasterBrick clipped to the polygons object
#'
#' @import raster glue
#' @export
load_format_pop_raster <- function(
  pop_fp_format, model_years, country_polys, projection
){
  # Load raster data
  raster_list <- vector('list', length=length(model_years))
  for(year_idx in 1:length(model_years)){
    year <- model_years[year_idx]
    this_rast_layer <- raster::raster(glue::glue(pop_fp_format))
    if(class(this_rast_layer) != 'RasterLayer') stop("Bad pop raster for ", year)
    raster_list[[year_idx]] <- this_rast_layer
  }

  # Convert to raster brick
  pop_brick <- raster::brick(raster_list)
  # Switch to modeling projection (also used by the spatial object)
  pop_brick <- raster::projectRaster(pop_brick, crs = projection)
  # Clip and mask raster brick to polygon outlines
  pop_brick <- raster::crop(x = pop_brick, y = country_polys)
  pop_brick <- raster::mask(x = pop_brick, mask = country_polys)

  # Return the formatted population raster brick
  return(pop_brick)
}


#' Get populated grid cells by location
#'
#' @description Pull the latitude, longitude, and population for each grid cell
#'   associated with each location in a shapefile
#'
#' @param pop Population rasterLayer (representing a single year)
#' @param polys_sf sf object representing polygons, with a field that has a
#'    unique identifier for each location
#' @param loc_field Location identifier field in the polygons sf object
#'
#' @return A data.table containing four columns, where each row represents a
#'    grid cell:
#'     - <loc_field>: The location identifier
#'     - pop: Population
#'     - lat: Latitude
#'     - lon: Longitude
#'
#' @import data.table raster fasterize
#' @export
index_populated_grid_cells <- function(pop, polys_sf, loc_field){

  # Reproject population and province IDs to lat-longs
  dd_crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  prov_rast <- fasterize::fasterize(
    sf = polys_sf, raster = pop, field = loc_field
  )
  pop_reproj <- raster::projectRaster(from = pop, crs = dd_crs)
  prov_rast_reproj <- raster::projectRaster(
    from = prov_rast, crs = dd_crs, method = 'ngb'
  )

  # Get the latitude and longitude associated with each populated pixel in each
  #  province
  prov_coords <- raster::rasterToPoints(x = prov_rast_reproj)
  prov_indexing <- data.table::data.table(
    loc_dummy = as.vector(prov_rast_reproj),
    pop = as.vector(pop_reproj)
  )[!is.na(loc_dummy)]
  if(nrow(prov_indexing) != nrow(prov_coords)) stop("Issue with dimensions")
  prov_indexing$lon <- prov_coords[, 'x']
  prov_indexing$lat <- prov_coords[, 'y']
  prov_indexing <- na.omit(prov_indexing)
  setnames(prov_indexing, 'loc_dummy', loc_field)

  return(prov_indexing)
}
