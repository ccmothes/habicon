#' Create a binary map
#'
#' This function creates a binary map of suitable/unsuitable habitat based on a user-specified threshold value
#'
#' @param x A SpatRaster object
#' @param threshold A value (or range of values) to use as suitability thresholds
#'
#' @return A SpatRaster object
#' @export
#'
#' @examples
bin_map <- function(x, threshold) {
  output <- vector("list", length = length(threshold))
  for (i in 1:length(threshold)) {

      output[[i]] <- classify(x,
                                matrix(
                                  c(-Inf, threshold[[i]], 0,
                                    threshold[[i]], Inf, 1),
                                  ncol = 3,
                                  byrow = TRUE,
                                ), right = FALSE)
      names(output[[i]]) <- paste("threshold", threshold[[i]])
  }
  if (length(output) > 1) {
    rast(output)
  } else {
    return(output[[i]])
  }
}
