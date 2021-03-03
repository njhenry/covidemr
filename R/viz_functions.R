#' Get largest Italian cities
#'
#' @description Returns a data.table with all Italian cities with an estimated
#'   population of over 500,000
#'
#' @return Data.table with three fields: 'name', 'lon', 'lat'
#' @import data.table
#' @export
italian_cities_dt <- function(){
  return(data.table(
    name = c('Rome', 'Milan', 'Naples', 'Turin', 'Palermo', 'Genoa'),
    lon =  c( 12.50,    9.18,    14.25,    7.67,     13.37,    8.92),
    lat =  c( 41.88,   45.49,    40.83,   45.07,     38.12,   44.40)
  ))
}


#' Make a choropleth map of Italian provinces
#'
#' @description Create a choropleth map for Italian provinces based on
#'   continuous input data by province
#'
#' @param province_sf [sf] Spatial object of provinces. Should have a
#'   `location_code` field that will be matched to data
#' @param region_sf [sf] Spatial object of regions. Will be added as an outline
#'   to map
#' @param in_data [data.table] Dataset to be plotted. Should only include one
#'   row of data per province
#' @param map_field [char] Field in `in_data` to be mapped
#' @param fill_list list of fill parameters that should be used, in the format
#'   of arguments to the `scale_fill_<fill_type>`
#' @param fill_lims [numeric, default NULL] Vector of upper and lower values for
#'   color scale.
#' @param titles_list [default empty list] list of any titles that should be
#'   included, in the format of arguments to the ggplot `labs` command
#' @param fill_type [optional, default 'continuous'] Scale fill type eg.
#'   'discrete' or 'continuous'
#' @param show_legend Should the legend be featured on the plot?
#' @param labels_dt [optional] Data.table used to label locations on the graph. Should be
#'   a data.table containing at least the fields `label_name`, `x`, `y`, `x1_lseg`,
#'   `x2_lseg`, `y1_lseg`, `y2_lseg`
#' @param save_fp [optional, default NULL] if not NULL (the default), saves to
#'   file rather than returning an object. If filled, should be a PNG extension
#' @param file_height [optional, default 8] map height in inches
#' @param file_width [optional, default 6] map width in inches
#'
#' @return Either returns a ggplot object or, if `save_fp` is not NULL, saves to
#'   file and returns NULL
#'
#' @import data.table ggplot2 glue
#' @export
map_ita_choropleth <- function(
  province_sf, region_sf, in_data, map_field, fill_list, fill_lims = NULL,
  titles_list = list(), fill_type = 'continuous', show_legend = TRUE, labels_dt = NULL,
  save_fp = NULL, file_height = 7.25, file_width = 6
){
  loc_codes <- unique(province_sf$location_code)
  # Ensure the correct data dimensions
  if(in_data[, .N, by=location_code][, max(N)] > 1) stop("Location duplicates in map")
  if(length(setdiff(loc_codes, in_data$location_code)) > 0){
    message("Missing values.. will fill with grey")
  }

  # Merge mapping field onto province polygons
  map_data <- copy(in_data[, c('location_code', map_field), with = FALSE])
  setnames(map_data, map_field, 'to_map')
  sf_merged <- merge(province_sf, map_data, by='location_code', all.x=TRUE)

  if(!is.null(fill_lims)){
    sf_merged[ sf_merged$to_map > max(fill_lims), 'to_map'] <- max(fill_lims) - 1E-7
    sf_merged[ sf_merged$to_map < min(fill_lims), 'to_map'] <- min(fill_lims) + 1E-7
    fill_list$limits <- fill_lims
  }

  line_color <- '#444444'
  na_color <- '#AAAAAA'

  # Construct map
  fig <- ggplot(data = sf_merged) +
    geom_sf(data = sf_merged, aes(fill=to_map), color=line_color, lwd=.2) +
    geom_sf(data = region_sf, color=line_color, fill=NA, lwd=.5) +
    do.call(glue::glue('scale_fill_{fill_type}'), fill_list) +
    do.call('labs', titles_list) +
    lims(x=c(3.5E5, 1.27E6), y=c(4.1E6, 5.2E6)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  if(!is.null(labels_dt)){
    fig <- fig +
      geom_text(
        data = labels_dt, aes(x=x, y=y, label=label_name),
        hjust='middle', vjust='center', size=3, lineheight=.75
      ) +
      geom_segment(
        data = labels_dt, aes(x=x1_lseg, y=y1_lseg, xend=x2_lseg, yend=y2_lseg)
      )
  }
  if(show_legend & (fill_type == 'manual') & !is.null(fill_list$values)){
    if(!is.null(names(fill_list$values))){
      dummy_dt <- data.table(fill_vals = names(fill_list$values))
      fig <- fig + geom_blank(data = dummy_dt, aes(fill = fill_vals))
    }
  }
  if(show_legend){
    fig <- fig + guides(fill = guide_legend(
        label.position="left", label.hjust=1, reverse=TRUE, title.hjust = .5
    )) + theme(
      legend.position = c(.99, .99),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill='white', color=line_color, size=.5),
      legend.title = element_text(size = 12, hjust = .5),
      legend.text = element_text(size = 10),
      legend.key.width = unit(.55, 'inches'),
      legend.key.height = unit(.25, 'inches'),
      legend.margin = margin(12, 12, 12, 12)
    )
  } else {
    fig <- fig + guides(fill = FALSE)
  }

  # Either save to file or return
  if(!is.null(save_fp)){
    png(save_fp, height = file_height, width = file_width, units='in', res=300)
    print(fig)
    dev.off()
    invisible()
  } else {
    return(fig)
  }
}



#' Make a choropleth map of Italian REGIONS
#'
#' @description Create a choropleth map for Italian regions based on
#'   continuous input data by region
#'
#' @param region_sf [sf] Spatial object of regions. Will be added as an outline
#'   to map
#' @param in_data [data.table] Dataset to be plotted. Should only include one
#'   row of data per region
#' @param map_field [char] Field in `in_data` to be mapped
#' @param fill_list list of fill parameters that should be used, in the format
#'   of arguments to the `scale_fill_<fill_type>`
#' @param fill_lims [numeric, default NULL] Vector of upper and lower values for
#'   color scale.
#' @param titles_list [default empty list] list of any titles that should be
#'   included, in the format of arguments to the ggplot `labs` command
#' @param fill_type [optional, default 'continuous'] Scale fill type eg.
#'   'discrete' or 'continuous'
#' @param show_legend Should the legend be featured on the plot?
#' @param labels_dt [optional] Data.table used to label locations on the graph. Should be
#'   a data.table containing at least the fields `label_name`, `x`, `y`, `x1_lseg`,
#'   `x2_lseg`, `y1_lseg`, `y2_lseg`
#' @param save_fp [optional, default NULL] if not NULL (the default), saves to
#'   file rather than returning an object. If filled, should be a PNG extension
#' @param file_height [optional, default 8] map height in inches
#' @param file_width [optional, default 6] map width in inches
#'
#' @return Either returns a ggplot object or, if `save_fp` is not NULL, saves to
#'   file and returns NULL
#'
#' @import data.table ggplot2 glue
#' @export
map_ita_choropleth_region <- function(
  region_sf, in_data, map_field, fill_list, fill_lims = NULL,
  titles_list = list(), fill_type = 'continuous', show_legend = TRUE, labels_dt = NULL,
  save_fp = NULL, file_height = 7.25, file_width = 6
){
  loc_codes <- unique(region_sf$region)
  # Ensure the correct data dimensions
  if(in_data[, .N, by=region_code][, max(N)] > 1) stop("Location duplicates in map")
  if(length(setdiff(loc_codes, in_data$region_code)) > 0){
    message("Missing values.. will fill with grey")
  }

  # Merge mapping field onto region polygons
  map_data <- copy(in_data[, c('region_code', map_field), with = FALSE])
  setnames(map_data, map_field, 'to_map')
  sf_merged <- merge(region_sf, map_data, by='region_code', all.x=TRUE)

  if(!is.null(fill_lims)){
    sf_merged[ sf_merged$to_map > max(fill_lims), 'to_map'] <- max(fill_lims) - 1E-7
    sf_merged[ sf_merged$to_map < min(fill_lims), 'to_map'] <- min(fill_lims) + 1E-7
    fill_list$limits <- fill_lims
  }

  line_color <- '#444444'
  na_color <- '#AAAAAA'

  # Construct map
  fig <- ggplot(data = sf_merged) +
    geom_sf(data = sf_merged, aes(fill=to_map), color=line_color, lwd=.5) +
    do.call(glue::glue('scale_fill_{fill_type}'), fill_list) +
    do.call('labs', titles_list) +
    lims(x=c(3.5E5, 1.27E6), y=c(4.1E6, 5.2E6)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  if(!is.null(labels_dt)){
    fig <- fig +
      geom_text(
        data = labels_dt, aes(x=x, y=y, label=label_name),
        hjust='middle', vjust='center', size=3, lineheight=.75
      ) +
      geom_segment(
        data = labels_dt, aes(x=x1_lseg, y=y1_lseg, xend=x2_lseg, yend=y2_lseg)
      )
  }
  if(show_legend & (fill_type == 'manual') & !is.null(fill_list$values)){
    if(!is.null(names(fill_list$values))){
      dummy_dt <- data.table(fill_vals = names(fill_list$values))
      fig <- fig + geom_blank(data = dummy_dt, aes(fill = fill_vals))
    }
  }
  if(show_legend){
    fig <- fig + guides(fill = guide_legend(
        label.position="left", label.hjust=1, reverse=TRUE, title.hjust = .5
    )) + theme(
      legend.position = c(.98, .98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill='white', color=line_color, size=.5),
      legend.title = element_text(size = 15, hjust = .5),
      legend.text = element_text(size = 11),
      legend.key.width = unit(.55, 'inches'),
      legend.key.height = unit(.3, 'inches'),
      legend.margin = margin(12, 12, 12, 12)
    )
  } else {
    fig <- fig + guides(fill = FALSE)
  }

  # Either save to file or return
  if(!is.null(save_fp)){
    png(save_fp, height = file_height, width = file_width, units='in', res=300)
    print(fig)
    dev.off()
    invisible()
  } else {
    return(fig)
  }
}

