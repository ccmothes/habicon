#' Rescale raster values to 0-1
#'
#' @param x The SpatRaster object to be rescaled
#'
#' @return A SpatRaster object
#'
#' @export
rescale_map <- function(x) {
  (x - min(terra::values(x), na.rm = TRUE)) / (max(terra::values(x), na.rm = TRUE) - min(terra::values(x), na.rm = TRUE))
}
