

#' Create a binary map
#'
#' @param x A RasterLayer object
#' @param threshold
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
bin_map <- function(x, threshold) {
  output <- vector("list", length = length(threshold))
  for (i in 1:length(threshold)) {
    output[[i]] <- reclassify(x,
                              matrix(
                                c(-Inf, threshold[[i]], 0,
                                  threshold[[i]], Inf, 1),
                                ncol = 3,
                                byrow = TRUE,
                              ), right = FALSE)
    names(output[[i]]) <- paste("threshold", threshold[[i]])
  }
  if (length(output) > 1) {stack(output)} else {return(output[[i]])}
}
