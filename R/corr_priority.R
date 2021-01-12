#' Corridor priority calculation
#'
#' @param suit A RasterLayer object of habitat suitability values, values range from 0-1
#' @param suit_bin A RasterLayer object of binary suitable (1) and unsuitable (0) habitat created using \code{\link[habicon]{bin_map}}
#' @param corr A RasterLayer object of connectivity values representing conductance to movement, values range from 0-1
#' @param corr_bin A RasterLayer object of binary corridor (1) and non-corridor (0) values created using \code{\link[habicon]{bin_map}}
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
corr_priority <-
  function(suit, suit_bin, corr, corr_bin) {
    #all rasters must be able to stack, need to add warning

    out <- raster(suit)

    #convert corr to resistance values (currently conductance, must invert for 'R' metric later)
    corr_resist <- calc(corr, function(x) {((x - max(x)) * -1) + min(x)})

    landscape_suit <-
      landscapemetrics::get_patches(suit_bin)$'1' # ID discrete habitat patches (of class '1')

    patch_edge <-
      landscapemetrics::get_boundaries(landscape_suit,
                     consider_boundary = TRUE,
                     return_raster = TRUE)[[1]] %>%
      reclassify(matrix(c(-Inf, 0, NA), ncol = 3)) %>%
      mask(landscape_suit, .) #associate edge cell values with their patch ID

    landscape_corr <-
      landscapemetrics::get_patches(corr_bin, class = 1, return_raster = TRUE)$'1'# ID individual corridors

    #use package tabularaster to ID all patches each corridor connects

    patch_corr <- stack(landscape_corr, landscape_suit)

    pc <-
      tabularaster::as_tibble(patch_corr) %>% dplyr::filter(cellvalue != 0 | !is.na(cellvalue)) %>%
      dplyr::filter(duplicated(cellindex) |
               duplicated(cellindex, fromLast = TRUE))  %>% #keep both duplicated values (i.e. patch and corridor)
      tidyr::spread(dimindex, cellvalue) %>% dplyr::rename(corridor = '1', patch = '2') %>%
      dplyr::distinct(corridor, patch)


    ## get patch area and associate with patch ID
    patch_area <-
      landscapemetrics::lsm_p_area(landscape_suit) %>%   #units are hectares
      dplyr::select(class, value) %>% dplyr::rename(patch = class, area_ha = value)

    ##quality
    patch_char <-
      zonal(suit, landscape_suit, fun = 'mean') %>% tabularaster::as_tibble() %>%
      dplyr::rename(quality = mean, patch = zone) %>%
      dplyr::left_join(patch_area, by = 'patch') %>%
      dplyr::mutate(area_sqm = area_ha * 10000)

    # corr edge cells for 'width' calculations
    corr_edge <-
      landscapemetrics::get_boundaries(landscape_corr,
                     consider_boundary = FALSE,
                     return_raster = TRUE)[[1]] %>%
      reclassify(., c(-Inf, 0, NA), ncol = 3) %>%
      mask(landscape_corr, .)

    # create transition layer for costDistance calc from point to patch

    tran <-
      transition(mask(corr, corr_bin), function(x)
        mean(x), 16)  # remember transition works on conductance values
    tran <- geoCorrection(tran)

    for (i in 1:ncell(corr_bin)) {
      if (is.na(corr_bin[i]) | corr_bin[i] == 0) {
        # keep NA cells as NA
        out[i] <- NA
      } else {
        #get corridor ID
        corr_id <- landscape_corr[i]

        # filter out corridor i
        corr_bin_i <-
          calc(landscape_corr, function(x) {
            x[x != corr_id] <-  NA
            return(x)
          })

        #corr i edge cells for distance calc
        corr_edge_i <-
          tabularaster::as_tibble(corr_edge) %>% dplyr::filter(cellvalue == corr_id) %>% dplyr::pull(cellindex) %>%
          xyFromCell(corr_edge, .)

        #get resistance value (make sure corr_resist is resistance values, not conductance!)
        R <- corr_resist[i]

        #turn pixel into point for distance calculations
        p <- xyFromCell(corr_bin, cell = i)

        #calculate nearest distance to corr edge
        W <-
          min(pointDistance(p, corr_edge_i, lonlat = FALSE)) #filter edge points to corr i


        # # get patches corridor connects
        pc_i <- dplyr::filter(pc, corridor == corr_id)

        if (nrow(pc_i) < 2 ) { #switch this back to == 0?
          out[i] <- NA
        } else {
          patch_eq <- vector("list", length = nrow(pc_i))
          for (k in 1:nrow(pc_i)) {
            #get patch id value
            pid <- as.numeric(paste(pc_i[k, 2]))

            #filter patch edge cells to corridor of interest and turn into points
            patch_edge_k <- mask(patch_edge, corr_bin_i) %>%
              tabularaster::as_tibble() %>% dplyr::filter(cellvalue == pid) %>% dplyr::pull(cellindex) %>%
              xyFromCell(patch_edge, .)


            #get area (Ai) and quality (Qi) from patch k and costDistance to patch edge
            patch_eq[[k]] <-
              (dplyr::filter(patch_char, patch == pid)$area_sqm) *
              (dplyr::filter(patch_char, patch == pid)$quality) *
              exp(-min(costDistance(tran, p[1, 1:2], patch_edge_k))) #what are these distance units?

          }
          #final product is sum of Ak*Qke^Dk for all K
          # use that final sum value in the final equation to get priority value for pixel i
          patch_eq <- do.call('sum', patch_eq)

          out[i] <- patch_eq - R * exp(-W)
        }
      }
    }
    return(out)

  }
