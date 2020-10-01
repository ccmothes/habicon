#' Filter just the matrix of the corridor map. This removes corridor cells that are within identified habitat patches
#'
#' @param suit_bin A RasterLayer object of binary suitable (1) and unsuitable (0) habitat created using \code{\link[habicon]{bin_map}}
#' @param corr_bin A RasterLayer object of binary corridor (1) and non-corridor (0) values created using \code{\link[habicon]{bin_map}}
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
matrix_map <- function(suit_bin, corr_bin){
  suit_plus <- suit_bin + 1
  corr_source <- source_map(corr_bin)
  filter_patch <- corr_source + suit_plus
  matrix <- reclassify(filter_patch, matrix(c(-Inf, 2, 1, 2, Inf, NA), ncol = 3,
                                            byrow = TRUE), right = TRUE) %>%
    buffer(., width = res(corr_bin), doEdge = TRUE, dissolve = TRUE)
  return(matrix)


}
