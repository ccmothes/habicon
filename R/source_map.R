
#' Create a source map, where all values other than focal value become NA
#'
#' @param x A RasterLayer object
#' @param value The cell value to create source surface from, default is 1
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
source_map <- function(x, value = 1){
  x[x != value] <- NA
  return(x)
}
