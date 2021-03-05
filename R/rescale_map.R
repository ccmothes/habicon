#' Rescale raster values to 0-1
#'
#' @param x A RasterLayer object to be rescaled
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
#'
rescale_map <- function(x){
  (x - min(values(x), na.rm = TRUE))/(max(values(x), na.rm = TRUE) - min(values(x), na.rm = TRUE))
}
