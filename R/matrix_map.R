#' Filter just the matrix of the corridor map. This removes corridor cells that
#' are within identified habitat patches
#'
#' @param suit_bin A SpatRaster object of binary suitable (1) and unsuitable
#'   (0) habitat created using \code{\link{bin_map}}
#' @param corr_bin A SpatRaster object of binary corridor (1) and non-corridor
#'   (0) values created using \code{\link{bin_map}}
#'
#' @return A SpatRaster object
#' @export
matrix_map <- function(suit_bin, corr_bin) {
  suit_plus <- suit_bin + 1
  corr_source <- source_map(corr_bin)
  filter_patch <- corr_source + suit_plus
  matrix <- terra::classify(filter_patch, matrix(c(-Inf, 2, 1, 2, Inf, NA),
    ncol = 3,
    byrow = TRUE
  ), right = TRUE) %>%
    terra::buffer(., width = terra::res(corr_bin)[1])
  return(matrix)
}
