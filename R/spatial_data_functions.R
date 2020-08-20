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
#'
#' @return sparse dsCMatrix representing adjacency between polygons in the
#'   `poly_sp` object
#'
#' @import spdep
#' @export
build_adjacency_matrix <- function(poly_sp, style='B', allow_zero_neighbors=TRUE){
  ## Generate adjacency matrix for the polygon
  adjmat <- spdep::nb2mat(
    neighbours = spdep::poly2nb(poly_sp),
    style = style,
    zero.policy = allow_zero_neighbors
  )
}
