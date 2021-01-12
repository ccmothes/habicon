#' Calculate habitat patch priority
#'
#' @param suit A RasterLayer object of habitat suitability values, values range from 0-1
#' @param suit_bin A RasterLayer object of binary suitable (1) and unsuitable (0) habitat created using \code{\link[habicon]{bin_map}}
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
patch_priority <- function(suit, suit_bin){
  out <- raster(suit)
  #id patches
  landscape_suit <-
    landscapemetrics::get_patches(suit_bin)$'1' # ID discrete habitat patches (of class '1')

  # get patch area and associate with patch ID
  patch_area <-
    landscapemetrics::lsm_p_area(landscape_suit) %>%   #units are hectares
    dplyr::select(class, value) %>% dplyr::rename(patch = class, area_ha = value)

  # get patch quality and overall patch score
  patch_char <-
    zonal(suit, landscape_suit, fun = 'mean') %>% tabularaster::as_tibble() %>%
    dplyr::rename(quality = mean, patch = zone) %>%
    dplyr::left_join(patch_area, by = 'patch') %>% dplyr::mutate(area_sqm = area_ha * 10000) %>%
    dplyr::mutate(score = area_sqm*quality)

  # for each patch cell, assign value that corresponds with patch ID
  for (i in 1:ncell(suit_bin)) {
    if (is.na(suit_bin[i]) | suit_bin[i] == 0) {
      # keep non-patch cells as NA
      out[i] <- NA
    } else {
      patch_id <- landscape_suit[i]

      out[i] <- dplyr::filter(patch_char, patch == patch_id) %>% dplyr::pull(score)
    }
  }

  return(out)

}
